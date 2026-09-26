// Optional CUDA implementation for the RFCP algorithm.
//
// It lives in inst/cuda so R source builds distribute it without requiring
// CUDA. The local CUDA build script compiles it explicitly with nvcc.
//
// Optimized Algorithm Overview (100% GPU Execution, Zero PCIe stalls)
// ====================================================================
// For each tile of queries in [query_begin, query_end):
//
//   Phase A — QT recurrence & Distance calculation (GPU):
//     Double-buffered QT rows advance diagonally on the device.
//     For each query, Euclidean distances are computed into single-query
//     device buffers `d_dist_aa` and `d_dist_ab`.
//
//   Phase B — Warp/Block-Level Top-K Reduction with Exclusion (GPU):
//     Instead of sorting millions of elements via Radix Sort, a specialized
//     reduction kernel finds the K non-overlapping nearest neighbors iteratively.
//     Each reduction step scans the distance buffer in parallel, loads the K
//     already-accepted candidate indices into shared memory, skips candidates
//     within the exclusion radius [idx - r, idx + r], and reduces the global
//     minimum via warp shuffles (__shfl_down_sync) and block-level reduction.
//     Complexity per query is O(K * N / num_threads), avoiding all CUB memory
//     overhead and multi-pass radix sort passes.
//
//   Phase C — RMSC Contrast Computation (GPU):
//     A fast device kernel calculates the paired contrast score and RMS value
//     directly into `d_rms_profile[out_i]` and `d_paired_frequency[out_i]`.
//
//   Phase D — Single Bulk Transfer & Post-Loop Plato Retrieval:
//     Zero host-device synchronisations occur during the entire query loop.
//     Only the final `rms_profile`, `paired_frequency`, and optional profile
//     matrices are transferred back in bulk at the end.
//     The global Plato is found in a single fast host scan. If Plato distance
//     vectors are requested, they are computed once for the winning query and
//     transferred, eliminating redundant 70MB PCIe transfers inside the loop.

#include "rfcp_gpu_poc.h"

#include <cuda_runtime.h>
#include <math_constants.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace matrixprofiler
{
  namespace gpu
  {
    namespace
    {

      constexpr int kTPB = 256; // threads per block
      constexpr int kNaInteger = std::numeric_limits<int>::min(); // Matches R's NA_INTEGER
      constexpr std::uint32_t kMaxSupportedK = 64;
      constexpr std::uint32_t kNumReductionBlocks = 256;

      bool cuda_ok(cudaError_t s) { return s == cudaSuccess; }

      // -----------------------------------------------------------------------
      // Comparator for reduction: smaller distance first, then smaller index
      // -----------------------------------------------------------------------
      __device__ inline void reduce_candidate(
          double &best_d, int &best_idx, double val_d, int val_idx)
      {
        if (val_d < best_d || (val_d == best_d && val_idx < best_idx))
        {
          best_d   = val_d;
          best_idx = val_idx;
        }
      }

      // -----------------------------------------------------------------------
      // Kernel: initialise QT_aa row for anchor query.
      //   qt_aa[j] = sum_{k=0}^{m-1} pos[anchor+k] * pos[j+k]
      // -----------------------------------------------------------------------
      __global__ void rfcp_init_qt_aa_kernel(
          const double *pos,
          std::uint32_t window_size, std::uint32_t anchor,
          std::uint32_t pos_profiles, double *qt_aa)
      {
        std::uint32_t const j = blockIdx.x * blockDim.x + threadIdx.x;
        if (j >= pos_profiles)
          return;

        double sum = 0.0;
        const double *p_anc = pos + anchor;
        const double *p_ref = pos + j;
        for (std::uint32_t k = 0; k < window_size; ++k)
          sum += p_anc[k] * p_ref[k];
        qt_aa[j] = sum;
      }

      // -----------------------------------------------------------------------
      // Kernel: initialise QT_ab row for anchor query.
      //   qt_ab[j] = sum_{k=0}^{m-1} pos[anchor+k] * neg[j+k]
      // -----------------------------------------------------------------------
      __global__ void rfcp_init_qt_ab_kernel(
          const double *pos, const double *neg,
          std::uint32_t window_size, std::uint32_t anchor,
          std::uint32_t neg_profiles, double *qt_ab)
      {
        std::uint32_t const j = blockIdx.x * blockDim.x + threadIdx.x;
        if (j >= neg_profiles)
          return;

        double sum = 0.0;
        const double *p_anc = pos + anchor;
        const double *p_ref = neg + j;
        for (std::uint32_t k = 0; k < window_size; ++k)
          sum += p_anc[k] * p_ref[k];
        qt_ab[j] = sum;
      }

      // -----------------------------------------------------------------------
      // Kernel: precompute AA boundary values (j=0) for all queries in a tile.
      // -----------------------------------------------------------------------
      __global__ void rfcp_precompute_aa_boundary_kernel(
          const double *pos,
          std::uint32_t window_size, std::uint32_t tile_begin,
          std::uint32_t tile_count, double *aa_boundary)
      {
        std::uint32_t const t = blockIdx.x * blockDim.x + threadIdx.x;
        if (t >= tile_count)
          return;

        std::uint32_t const q = tile_begin + t;
        double sum = 0.0;
        const double *p_q = pos + q;
        for (std::uint32_t k = 0; k < window_size; ++k)
          sum += p_q[k] * pos[k];
        aa_boundary[t] = sum;
      }

      // -----------------------------------------------------------------------
      // Kernel: precompute AB boundary values (j=0) for all queries in a tile.
      // -----------------------------------------------------------------------
      __global__ void rfcp_precompute_ab_boundary_kernel(
          const double *pos, const double *neg,
          std::uint32_t window_size, std::uint32_t tile_begin,
          std::uint32_t tile_count, double *ab_boundary)
      {
        std::uint32_t const t = blockIdx.x * blockDim.x + threadIdx.x;
        if (t >= tile_count)
          return;

        std::uint32_t const q = tile_begin + t;
        double sum = 0.0;
        const double *p_q = pos + q;
        for (std::uint32_t k = 0; k < window_size; ++k)
          sum += p_q[k] * neg[k];
        ab_boundary[t] = sum;
      }

      // -----------------------------------------------------------------------
      // Kernel: advance QT_aa row by one query step (q -> q+1).
      // -----------------------------------------------------------------------
      __global__ void rfcp_advance_qt_aa_kernel(
          const double *pos,
          std::uint32_t window_size, std::uint32_t pos_profiles,
          std::uint32_t q, std::uint32_t tile_begin,
          const double *aa_boundary,
          const double *qt_prev, double *qt_curr)
      {
        std::uint32_t const j = blockIdx.x * blockDim.x + threadIdx.x;
        if (j >= pos_profiles)
          return;

        if (j == 0)
        {
          qt_curr[0] = aa_boundary[q - tile_begin];
        }
        else
        {
          double const term_old = pos[q - 1] * pos[j - 1];
          double const term_new = pos[q + window_size - 1] * pos[j + window_size - 1];
          qt_curr[j] = qt_prev[j - 1] - term_old + term_new;
        }
      }

      // -----------------------------------------------------------------------
      // Kernel: advance QT_ab row by one query step (q -> q+1).
      // -----------------------------------------------------------------------
      __global__ void rfcp_advance_qt_ab_kernel(
          const double *pos, const double *neg,
          std::uint32_t window_size, std::uint32_t neg_profiles,
          std::uint32_t q, std::uint32_t tile_begin,
          const double *ab_boundary,
          const double *qt_prev, double *qt_curr)
      {
        std::uint32_t const j = blockIdx.x * blockDim.x + threadIdx.x;
        if (j >= neg_profiles)
          return;

        if (j == 0)
        {
          qt_curr[0] = ab_boundary[q - tile_begin];
        }
        else
        {
          double const term_old = pos[q - 1] * neg[j - 1];
          double const term_new = pos[q + window_size - 1] * neg[j + window_size - 1];
          qt_curr[j] = qt_prev[j - 1] - term_old + term_new;
        }
      }

      // -----------------------------------------------------------------------
      // Kernel: compute Euclidean distances for AA join (self-join).
      // -----------------------------------------------------------------------
      __global__ void rfcp_distances_aa_kernel(
          const double *qt_aa, std::uint32_t pos_profiles,
          const double *mean_pos, const double *sig_pos, const int *valid_pos,
          std::uint32_t window_size, std::uint32_t q,
          std::uint32_t exclusion_radius, double *dist_aa,
          std::uint32_t *cand_count)
      {
        if (cand_count && blockIdx.x == 0 && threadIdx.x == 0)
          *cand_count = 0;

        std::uint32_t const j = blockIdx.x * blockDim.x + threadIdx.x;
        if (j >= pos_profiles)
          return;

        // Self-exclusion zone [q - exclusion_radius, q + exclusion_radius]
        std::uint32_t const left =
            q >= exclusion_radius ? q - exclusion_radius : 0;
        std::uint32_t const right =
            q + exclusion_radius < pos_profiles ? q + exclusion_radius : pos_profiles - 1;

        if (j >= left && j <= right)
        {
          dist_aa[j] = 1e300;
          return;
        }

        if (!valid_pos[j])
        {
          dist_aa[j] = 1e300;
          return;
        }

        double const mean_q = mean_pos[q];
        double const sig_q  = sig_pos[q];
        double const mean_j = mean_pos[j];
        double const sig_j  = sig_pos[j];

        double const cov  = qt_aa[j] - static_cast<double>(window_size) * mean_q * mean_j;
        double const corr = fmax(-1.0, fmin(1.0, cov * sig_q * sig_j));
        double const d2   = 2.0 * static_cast<double>(window_size) * (1.0 - corr);
        dist_aa[j]        = sqrt(fmax(0.0, d2));
      }

      // -----------------------------------------------------------------------
      // Kernel: compute Euclidean distances for AB join.
      // -----------------------------------------------------------------------
      __global__ void rfcp_distances_ab_kernel(
          const double *qt_ab, std::uint32_t neg_profiles,
          const double *mean_pos, const double *sig_pos,
          const double *mean_neg, const double *sig_neg, const int *valid_neg,
          std::uint32_t window_size, std::uint32_t q, double *dist_ab,
          std::uint32_t *cand_count)
      {
        if (cand_count && blockIdx.x == 0 && threadIdx.x == 0)
          *cand_count = 0;

        std::uint32_t const j = blockIdx.x * blockDim.x + threadIdx.x;
        if (j >= neg_profiles)
          return;

        if (!valid_neg[j])
        {
          dist_ab[j] = 1e300;
          return;
        }

        double const mean_q = mean_pos[q];
        double const sig_q  = sig_pos[q];
        double const mean_j = mean_neg[j];
        double const sig_j  = sig_neg[j];

        double const cov  = qt_ab[j] - static_cast<double>(window_size) * mean_q * mean_j;
        double const corr = fmax(-1.0, fmin(1.0, cov * sig_q * sig_j));
        double const d2   = 2.0 * static_cast<double>(window_size) * (1.0 - corr);
        dist_ab[j]        = sqrt(fmax(0.0, d2));
      }

      // -----------------------------------------------------------------------
      // Kernel: Find next nearest non-overlapping candidate via Block Reduction
      //
      // Scans `dist` array of size `n`, enforces non-overlapping exclusion zone
      // with all previously selected candidates in `accepted_idx[0..k-1]`, and
      // writes the per-block minimum to `block_min_dist` and `block_min_idx`.
      // -----------------------------------------------------------------------
      __global__ void rfcp_find_next_candidate_kernel(
          const double * __restrict__ dist,
          std::uint32_t n,
          const int * __restrict__ accepted_idx,
          std::uint32_t k,
          std::uint32_t exclusion_radius,
          double * __restrict__ block_min_dist,
          int * __restrict__ block_min_idx)
      {
        __shared__ int s_acc[kMaxSupportedK];
        if (threadIdx.x < k)
          s_acc[threadIdx.x] = accepted_idx[threadIdx.x];
        __syncthreads();

        double thread_min_d   = 1e300;
        int    thread_min_idx = -1;

        std::uint32_t const idx    = blockIdx.x * blockDim.x + threadIdx.x;
        std::uint32_t const stride = gridDim.x * blockDim.x;

        for (std::uint32_t j = idx; j < n; j += stride)
        {
          double const d = dist[j];
          if (d < thread_min_d)
          {
            bool excluded = false;
            #pragma unroll 4
            for (std::uint32_t a = 0; a < k; ++a)
            {
              int diff = static_cast<int>(j) - s_acc[a];
              if (diff < 0) diff = -diff;
              if (diff <= static_cast<int>(exclusion_radius))
              {
                excluded = true;
                break;
              }
            }
            if (!excluded)
            {
              thread_min_d   = d;
              thread_min_idx = static_cast<int>(j);
            }
          }
        }

        // Warp reduction within each 32-thread warp
        #pragma unroll
        for (int offset = 16; offset > 0; offset /= 2)
        {
          double const other_d   = __shfl_down_sync(0xffffffff, thread_min_d, offset);
          int    const other_idx = __shfl_down_sync(0xffffffff, thread_min_idx, offset);
          reduce_candidate(thread_min_d, thread_min_idx, other_d, other_idx);
        }

        __shared__ double s_min_d[8]; // 256 threads = 8 warps
        __shared__ int    s_min_idx[8];
        int const lane = threadIdx.x % 32;
        int const warp = threadIdx.x / 32;

        if (lane == 0)
        {
          s_min_d[warp]   = thread_min_d;
          s_min_idx[warp] = thread_min_idx;
        }
        __syncthreads();

        if (warp == 0)
        {
          double b_min_d   = (lane < 8) ? s_min_d[lane]   : 1e300;
          int    b_min_idx = (lane < 8) ? s_min_idx[lane] : -1;
          #pragma unroll
          for (int offset = 4; offset > 0; offset /= 2)
          {
            double const other_d   = __shfl_down_sync(0xffffffff, b_min_d, offset);
            int    const other_idx = __shfl_down_sync(0xffffffff, b_min_idx, offset);
            reduce_candidate(b_min_d, b_min_idx, other_d, other_idx);
          }
          if (lane == 0)
          {
            block_min_dist[blockIdx.x] = b_min_d;
            block_min_idx[blockIdx.x]  = b_min_idx;
          }
        }
      }

      // -----------------------------------------------------------------------
      // Kernel: Finalize Candidate Selection across Blocks
      //
      // Launched with 1 block of 256 threads. Reduces `num_blocks` block minima
      // to determine Candidate `k` and updates `cand_count`.
      // -----------------------------------------------------------------------
      __global__ void rfcp_finalize_candidate_kernel(
          const double * __restrict__ block_min_dist,
          const int * __restrict__ block_min_idx,
          std::uint32_t num_blocks,
          std::uint32_t k,
          double * __restrict__ cand_dist,
          int * __restrict__ cand_idx,
          std::uint32_t * __restrict__ cand_count)
      {
        int const tid = threadIdx.x;
        double d = (tid < num_blocks) ? block_min_dist[tid] : 1e300;
        int  idx = (tid < num_blocks) ? block_min_idx[tid]  : -1;

        #pragma unroll
        for (int offset = 16; offset > 0; offset /= 2)
        {
          double const other_d   = __shfl_down_sync(0xffffffff, d, offset);
          int    const other_idx = __shfl_down_sync(0xffffffff, idx, offset);
          reduce_candidate(d, idx, other_d, other_idx);
        }

        __shared__ double s_min_d[8];
        __shared__ int    s_min_idx[8];
        int const lane = tid % 32;
        int const warp = tid / 32;

        if (lane == 0)
        {
          s_min_d[warp]   = d;
          s_min_idx[warp] = idx;
        }
        __syncthreads();

        if (warp == 0)
        {
          double b_d   = (lane < 8) ? s_min_d[lane]   : 1e300;
          int    b_idx = (lane < 8) ? s_min_idx[lane] : -1;
          #pragma unroll
          for (int offset = 4; offset > 0; offset /= 2)
          {
            double const other_d   = __shfl_down_sync(0xffffffff, b_d, offset);
            int    const other_idx = __shfl_down_sync(0xffffffff, b_idx, offset);
            reduce_candidate(b_d, b_idx, other_d, other_idx);
          }
          if (lane == 0)
          {
            if (b_d < 1e200 && b_idx >= 0)
            {
              cand_dist[k]  = b_d;
              cand_idx[k]   = b_idx;
              cand_count[0] = k + 1;
            }
            else
            {
              cand_dist[k] = CUDART_NAN;
              cand_idx[k]  = -1;
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // Kernel: Compute Query RMSC and populate output profiles
      //
      // Computes paired contrast RMS and populates output buffers for query q.
      // -----------------------------------------------------------------------
      __global__ void rfcp_compute_query_rmsc_kernel(
          const double * __restrict__ cand_dist_aa,
          const int    * __restrict__ cand_idx_aa,
          const std::uint32_t * __restrict__ cand_count_aa,
          const double * __restrict__ cand_dist_ab,
          const int    * __restrict__ cand_idx_ab,
          const std::uint32_t * __restrict__ cand_count_ab,
          std::uint32_t out_i,
          std::uint32_t max_freq,
          double clip,
          int is_valid,
          double * __restrict__ d_rms_profile,
          int    * __restrict__ d_paired_frequency,
          double * __restrict__ d_profile_aa,
          double * __restrict__ d_profile_ab,
          double * __restrict__ d_profile_rfcp,
          int    * __restrict__ d_profile_index_aa,
          int    * __restrict__ d_profile_index_ab)
      {
        if (!is_valid)
        {
          if (threadIdx.x == 0)
          {
            d_rms_profile[out_i]      = CUDART_NAN;
            d_paired_frequency[out_i] = 0;
          }
          if (d_profile_aa)
          {
            std::uint64_t const col = static_cast<std::uint64_t>(out_i) * max_freq;
            for (std::uint32_t r = threadIdx.x; r < max_freq; r += blockDim.x)
            {
              d_profile_aa[col + r]       = CUDART_NAN;
              d_profile_ab[col + r]       = CUDART_NAN;
              d_profile_rfcp[col + r]     = CUDART_NAN;
              d_profile_index_aa[col + r] = kNaInteger;
              d_profile_index_ab[col + r] = kNaInteger;
            }
          }
          return;
        }

        __shared__ std::uint32_t s_aa_count;
        __shared__ std::uint32_t s_ab_count;
        __shared__ std::uint32_t s_paired;

        if (threadIdx.x == 0)
        {
          s_aa_count = cand_count_aa ? cand_count_aa[0] : 0;
          s_ab_count = cand_count_ab ? cand_count_ab[0] : 0;
          s_paired   = (s_aa_count < s_ab_count) ? s_aa_count : s_ab_count;
          d_paired_frequency[out_i] = static_cast<int>(s_paired);

          double sq_sum = 0.0;
          for (std::uint32_t r = 0; r < s_paired; ++r)
          {
            double const aa = (cand_dist_aa[r] < clip) ? cand_dist_aa[r] : clip;
            double const ab = (cand_dist_ab[r] < clip) ? cand_dist_ab[r] : clip;
            double const c  = (ab > aa) ? ((ab - aa) / clip) : 0.0;
            sq_sum += c * c;
          }
          double const rms = sqrt(sq_sum / static_cast<double>(max_freq));
          d_rms_profile[out_i] = rms;
        }
        __syncthreads();

        if (d_profile_aa)
        {
          std::uint64_t const col = static_cast<std::uint64_t>(out_i) * max_freq;
          for (std::uint32_t r = threadIdx.x; r < max_freq; r += blockDim.x)
          {
            double const aa_r = (r < s_aa_count) ? cand_dist_aa[r] : CUDART_NAN;
            double const ab_r = (r < s_ab_count) ? cand_dist_ab[r] : CUDART_NAN;
            d_profile_aa[col + r]       = aa_r;
            d_profile_ab[col + r]       = ab_r;
            d_profile_index_aa[col + r] = (r < s_aa_count) ? (cand_idx_aa[r] + 1) : kNaInteger;
            d_profile_index_ab[col + r] = (r < s_ab_count) ? (cand_idx_ab[r] + 1) : kNaInteger;

            if (!isnan(aa_r) && !isnan(ab_r))
            {
              double const aa_clip = (aa_r < clip) ? aa_r : clip;
              double const ab_clip = (ab_r < clip) ? ab_r : clip;
              d_profile_rfcp[col + r] = (ab_clip > aa_clip) ? ((ab_clip - aa_clip) / clip) : 0.0;
            }
            else
            {
              d_profile_rfcp[col + r] = CUDART_NAN;
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // Compute tile size matching the CPU rfcp.cpp heuristic.
      // -----------------------------------------------------------------------
      std::uint32_t rfcp_gpu_tile_size(std::uint32_t window_size,
                                       std::uint32_t query_count)
      {
        std::uint64_t const scaled = 16ULL * window_size;
        std::uint32_t const base =
            static_cast<std::uint32_t>(std::max<std::uint64_t>(
                1024ULL, std::min<std::uint64_t>(32768ULL, scaled)));
        return std::min(base, query_count);
      }

      // -----------------------------------------------------------------------
      // Core GPU implementation.
      // -----------------------------------------------------------------------
      RfcpGpuStatus rfcp_gpu_compute_impl(
          const double *positive,    const double *mean_positive,
          const double *sig_positive, const int *valid_positive,
          std::uint32_t positive_size,
          const double *negative,    const double *mean_negative,
          const double *sig_negative, const int *valid_negative,
          std::uint32_t negative_size,
          std::uint32_t window_size, std::uint32_t max_freq,
          std::uint32_t exclusion_radius,
          std::uint32_t query_begin, std::uint32_t query_end,
          double *rms_profile, int *paired_frequency,
          std::uint32_t *plato_index, double *plato_rms,
          double *plato_dist_aa, double *plato_dist_ab,
          double *profile_aa,       double *profile_ab,
          double *profile_rfcp,
          int    *profile_index_aa, int *profile_index_ab)
      {
        if (max_freq > kMaxSupportedK || max_freq < 1)
          return RfcpGpuStatus::invalid_argument;

        std::uint32_t const pos_profiles = positive_size - window_size + 1;
        std::uint32_t const neg_profiles = negative_size - window_size + 1;
        std::uint32_t const query_count  = query_end - query_begin;
        std::uint32_t const tile_size    = rfcp_gpu_tile_size(window_size, query_count);

        // -- Device allocations --
        double *d_pos         = nullptr;
        double *d_neg         = nullptr;
        double *d_mean_pos    = nullptr;
        double *d_sig_pos     = nullptr;
        double *d_mean_neg    = nullptr;
        double *d_sig_neg     = nullptr;
        int    *d_valid_pos   = nullptr;
        int    *d_valid_neg   = nullptr;
        double *d_qt_aa_0     = nullptr;
        double *d_qt_aa_1     = nullptr;
        double *d_qt_ab_0     = nullptr;
        double *d_qt_ab_1     = nullptr;
        double *d_aa_boundary = nullptr;
        double *d_ab_boundary = nullptr;

        // Current query distance buffers
        double *d_dist_aa     = nullptr;
        double *d_dist_ab     = nullptr;

        // Reduction auxiliary buffers
        double *d_block_min_dist = nullptr;
        int    *d_block_min_idx  = nullptr;
        double *d_cand_dist_aa   = nullptr;
        int    *d_cand_idx_aa    = nullptr;
        std::uint32_t *d_cand_count_aa = nullptr;
        double *d_cand_dist_ab   = nullptr;
        int    *d_cand_idx_ab    = nullptr;
        std::uint32_t *d_cand_count_ab = nullptr;

        // Device output arrays
        double *d_rms_profile       = nullptr;
        int    *d_paired_frequency  = nullptr;
        double *d_profile_aa        = nullptr;
        double *d_profile_ab        = nullptr;
        double *d_profile_rfcp      = nullptr;
        int    *d_profile_index_aa  = nullptr;
        int    *d_profile_index_ab  = nullptr;

        auto cleanup = [&]()
        {
          cudaFree(d_pos);             cudaFree(d_neg);
          cudaFree(d_mean_pos);        cudaFree(d_sig_pos);
          cudaFree(d_mean_neg);        cudaFree(d_sig_neg);
          cudaFree(d_valid_pos);       cudaFree(d_valid_neg);
          cudaFree(d_qt_aa_0);         cudaFree(d_qt_aa_1);
          cudaFree(d_qt_ab_0);         cudaFree(d_qt_ab_1);
          cudaFree(d_aa_boundary);     cudaFree(d_ab_boundary);
          cudaFree(d_dist_aa);         cudaFree(d_dist_ab);
          cudaFree(d_block_min_dist);  cudaFree(d_block_min_idx);
          cudaFree(d_cand_dist_aa);    cudaFree(d_cand_idx_aa);
          cudaFree(d_cand_count_aa);
          cudaFree(d_cand_dist_ab);    cudaFree(d_cand_idx_ab);
          cudaFree(d_cand_count_ab);
          cudaFree(d_rms_profile);     cudaFree(d_paired_frequency);
          cudaFree(d_profile_aa);      cudaFree(d_profile_ab);
          cudaFree(d_profile_rfcp);    cudaFree(d_profile_index_aa);
          cudaFree(d_profile_index_ab);
        };

        auto fail = [&]() -> RfcpGpuStatus
        {
          cleanup();
          return RfcpGpuStatus::cuda_error;
        };

        std::size_t const pos_bytes       = static_cast<std::size_t>(positive_size) * sizeof(double);
        std::size_t const neg_bytes       = static_cast<std::size_t>(negative_size) * sizeof(double);
        std::size_t const pos_prof_bytes  = static_cast<std::size_t>(pos_profiles)  * sizeof(double);
        std::size_t const neg_prof_bytes  = static_cast<std::size_t>(neg_profiles)  * sizeof(double);
        std::size_t const pos_prof_ibytes = static_cast<std::size_t>(pos_profiles)  * sizeof(int);
        std::size_t const neg_prof_ibytes = static_cast<std::size_t>(neg_profiles)  * sizeof(int);
        std::size_t const tile_bytes      = static_cast<std::size_t>(tile_size)     * sizeof(double);

        // Allocate base buffers
        if (!cuda_ok(cudaMalloc(&d_pos,         pos_bytes))       ||
            !cuda_ok(cudaMalloc(&d_neg,         neg_bytes))       ||
            !cuda_ok(cudaMalloc(&d_mean_pos,    pos_prof_bytes))  ||
            !cuda_ok(cudaMalloc(&d_sig_pos,     pos_prof_bytes))  ||
            !cuda_ok(cudaMalloc(&d_mean_neg,    neg_prof_bytes))  ||
            !cuda_ok(cudaMalloc(&d_sig_neg,     neg_prof_bytes))  ||
            !cuda_ok(cudaMalloc(&d_valid_pos,   pos_prof_ibytes)) ||
            !cuda_ok(cudaMalloc(&d_valid_neg,   neg_prof_ibytes)) ||
            !cuda_ok(cudaMalloc(&d_qt_aa_0,     pos_prof_bytes))  ||
            !cuda_ok(cudaMalloc(&d_qt_aa_1,     pos_prof_bytes))  ||
            !cuda_ok(cudaMalloc(&d_qt_ab_0,     neg_prof_bytes))  ||
            !cuda_ok(cudaMalloc(&d_qt_ab_1,     neg_prof_bytes))  ||
            !cuda_ok(cudaMalloc(&d_aa_boundary, tile_bytes))      ||
            !cuda_ok(cudaMalloc(&d_ab_boundary, tile_bytes))      ||
            !cuda_ok(cudaMalloc(&d_dist_aa,     pos_prof_bytes))  ||
            !cuda_ok(cudaMalloc(&d_dist_ab,     neg_prof_bytes)))
        {
          return fail();
        }

        // Allocate reduction buffers
        if (!cuda_ok(cudaMalloc(&d_block_min_dist, kNumReductionBlocks * sizeof(double))) ||
            !cuda_ok(cudaMalloc(&d_block_min_idx,  kNumReductionBlocks * sizeof(int)))    ||
            !cuda_ok(cudaMalloc(&d_cand_dist_aa,   max_freq * sizeof(double)))            ||
            !cuda_ok(cudaMalloc(&d_cand_idx_aa,    max_freq * sizeof(int)))               ||
            !cuda_ok(cudaMalloc(&d_cand_count_aa,  sizeof(std::uint32_t)))               ||
            !cuda_ok(cudaMalloc(&d_cand_dist_ab,   max_freq * sizeof(double)))            ||
            !cuda_ok(cudaMalloc(&d_cand_idx_ab,    max_freq * sizeof(int)))               ||
            !cuda_ok(cudaMalloc(&d_cand_count_ab,  sizeof(std::uint32_t))))
        {
          return fail();
        }

        // Allocate device output buffers
        if (!cuda_ok(cudaMalloc(&d_rms_profile,      query_count * sizeof(double))) ||
            !cuda_ok(cudaMalloc(&d_paired_frequency, query_count * sizeof(int))))
        {
          return fail();
        }

        std::uint64_t const profile_cells = static_cast<std::uint64_t>(query_count) * max_freq;
        if (profile_aa)
        {
          if (!cuda_ok(cudaMalloc(&d_profile_aa,       profile_cells * sizeof(double))) ||
              !cuda_ok(cudaMalloc(&d_profile_ab,       profile_cells * sizeof(double))) ||
              !cuda_ok(cudaMalloc(&d_profile_rfcp,     profile_cells * sizeof(double))) ||
              !cuda_ok(cudaMalloc(&d_profile_index_aa, profile_cells * sizeof(int)))    ||
              !cuda_ok(cudaMalloc(&d_profile_index_ab, profile_cells * sizeof(int))))
          {
            return fail();
          }
        }

        // Copy input arrays to device
        if (!cuda_ok(cudaMemcpy(d_pos,       positive,       pos_bytes,       cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(d_neg,       negative,       neg_bytes,       cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(d_mean_pos,  mean_positive,  pos_prof_bytes,  cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(d_sig_pos,   sig_positive,   pos_prof_bytes,  cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(d_mean_neg,  mean_negative,  neg_prof_bytes,  cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(d_sig_neg,   sig_negative,   neg_prof_bytes,  cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(d_valid_pos, valid_positive, pos_prof_ibytes, cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(d_valid_neg, valid_negative, neg_prof_ibytes, cudaMemcpyHostToDevice)))
        {
          return fail();
        }

        double const clip = std::sqrt(2.0 * static_cast<double>(window_size));
        dim3 const blk(kTPB);
        dim3 const grd_pos((pos_profiles + kTPB - 1) / kTPB);
        dim3 const grd_neg((neg_profiles + kTPB - 1) / kTPB);

        std::uint32_t const num_blocks_aa = std::max<std::uint32_t>(1U, std::min<std::uint32_t>(kNumReductionBlocks, (pos_profiles + kTPB - 1) / kTPB));
        std::uint32_t const num_blocks_ab = std::max<std::uint32_t>(1U, std::min<std::uint32_t>(kNumReductionBlocks, (neg_profiles + kTPB - 1) / kTPB));

        // --- Main tiled loop ---
        for (std::uint32_t tb = query_begin; tb < query_end; tb += tile_size)
        {
          std::uint32_t const te         = std::min(tb + tile_size, query_end);
          std::uint32_t const tile_count = te - tb;

          // Precompute boundary values (j=0) for queries in this tile
          dim3 const grd_tile((tile_count + kTPB - 1) / kTPB);
          rfcp_precompute_aa_boundary_kernel<<<grd_tile, blk>>>(
              d_pos, window_size, tb, tile_count, d_aa_boundary);
          rfcp_precompute_ab_boundary_kernel<<<grd_tile, blk>>>(
              d_pos, d_neg, window_size, tb, tile_count, d_ab_boundary);
          if (!cuda_ok(cudaGetLastError()))
            return fail();

          // Initialise QT for tile anchor (q = tb)
          rfcp_init_qt_aa_kernel<<<grd_pos, blk>>>(
              d_pos, window_size, tb, pos_profiles,
              (tb % 2 == 0) ? d_qt_aa_0 : d_qt_aa_1);
          rfcp_init_qt_ab_kernel<<<grd_neg, blk>>>(
              d_pos, d_neg, window_size, tb, neg_profiles,
              (tb % 2 == 0) ? d_qt_ab_0 : d_qt_ab_1);
          if (!cuda_ok(cudaGetLastError()))
            return fail();

          // Iterate through queries in tile
          for (std::uint32_t q = tb; q < te; ++q)
          {
            std::uint32_t const out_i = q - query_begin;
            double *qt_aa_curr = (q % 2 == 0) ? d_qt_aa_0 : d_qt_aa_1;
            double *qt_aa_prev = (q % 2 == 0) ? d_qt_aa_1 : d_qt_aa_0;
            double *qt_ab_curr = (q % 2 == 0) ? d_qt_ab_0 : d_qt_ab_1;
            double *qt_ab_prev = (q % 2 == 0) ? d_qt_ab_1 : d_qt_ab_0;

            // Advance QT (skip for tile anchor)
            if (q > tb)
            {
              rfcp_advance_qt_aa_kernel<<<grd_pos, blk>>>(
                  d_pos, window_size, pos_profiles, q, tb,
                  d_aa_boundary, qt_aa_prev, qt_aa_curr);
              rfcp_advance_qt_ab_kernel<<<grd_neg, blk>>>(
                  d_pos, d_neg, window_size, neg_profiles, q, tb,
                  d_ab_boundary, qt_ab_prev, qt_ab_curr);
            }

            if (!valid_positive[q])
            {
              rfcp_compute_query_rmsc_kernel<<<1, kTPB>>>(
                  nullptr, nullptr, nullptr,
                  nullptr, nullptr, nullptr,
                  out_i, max_freq, clip, 0,
                  d_rms_profile, d_paired_frequency,
                  d_profile_aa, d_profile_ab, d_profile_rfcp,
                  d_profile_index_aa, d_profile_index_ab);
              continue;
            }

            // Compute AA distances and reset AA candidate count
            rfcp_distances_aa_kernel<<<grd_pos, blk>>>(
                qt_aa_curr, pos_profiles,
                d_mean_pos, d_sig_pos, d_valid_pos,
                window_size, q, exclusion_radius, d_dist_aa, d_cand_count_aa);

            // Iterative greedy candidate reduction for AA
            for (std::uint32_t k = 0; k < max_freq; ++k)
            {
              rfcp_find_next_candidate_kernel<<<num_blocks_aa, kTPB>>>(
                  d_dist_aa, pos_profiles, d_cand_idx_aa, k, exclusion_radius,
                  d_block_min_dist, d_block_min_idx);
              rfcp_finalize_candidate_kernel<<<1, kTPB>>>(
                  d_block_min_dist, d_block_min_idx, num_blocks_aa, k,
                  d_cand_dist_aa, d_cand_idx_aa, d_cand_count_aa);
            }

            // Compute AB distances and reset AB candidate count
            rfcp_distances_ab_kernel<<<grd_neg, blk>>>(
                qt_ab_curr, neg_profiles,
                d_mean_pos, d_sig_pos,
                d_mean_neg, d_sig_neg, d_valid_neg,
                window_size, q, d_dist_ab, d_cand_count_ab);

            // Iterative greedy candidate reduction for AB
            for (std::uint32_t k = 0; k < max_freq; ++k)
            {
              rfcp_find_next_candidate_kernel<<<num_blocks_ab, kTPB>>>(
                  d_dist_ab, neg_profiles, d_cand_idx_ab, k, exclusion_radius,
                  d_block_min_dist, d_block_min_idx);
              rfcp_finalize_candidate_kernel<<<1, kTPB>>>(
                  d_block_min_dist, d_block_min_idx, num_blocks_ab, k,
                  d_cand_dist_ab, d_cand_idx_ab, d_cand_count_ab);
            }

            // Compute RMSC and record profiles for this query
            rfcp_compute_query_rmsc_kernel<<<1, kTPB>>>(
                d_cand_dist_aa, d_cand_idx_aa, d_cand_count_aa,
                d_cand_dist_ab, d_cand_idx_ab, d_cand_count_ab,
                out_i, max_freq, clip, 1,
                d_rms_profile, d_paired_frequency,
                d_profile_aa, d_profile_ab, d_profile_rfcp,
                d_profile_index_aa, d_profile_index_ab);
          }

          if (!cuda_ok(cudaGetLastError()))
            return fail();
        } // tile loop

        // Synchronize and copy final results to host
        if (!cuda_ok(cudaMemcpy(rms_profile, d_rms_profile, query_count * sizeof(double), cudaMemcpyDeviceToHost)) ||
            !cuda_ok(cudaMemcpy(paired_frequency, d_paired_frequency, query_count * sizeof(int), cudaMemcpyDeviceToHost)))
        {
          return fail();
        }

        if (profile_aa)
        {
          if (!cuda_ok(cudaMemcpy(profile_aa, d_profile_aa, profile_cells * sizeof(double), cudaMemcpyDeviceToHost)) ||
              !cuda_ok(cudaMemcpy(profile_ab, d_profile_ab, profile_cells * sizeof(double), cudaMemcpyDeviceToHost)) ||
              !cuda_ok(cudaMemcpy(profile_rfcp, d_profile_rfcp, profile_cells * sizeof(double), cudaMemcpyDeviceToHost)) ||
              !cuda_ok(cudaMemcpy(profile_index_aa, d_profile_index_aa, profile_cells * sizeof(int), cudaMemcpyDeviceToHost)) ||
              !cuda_ok(cudaMemcpy(profile_index_ab, d_profile_index_ab, profile_cells * sizeof(int), cudaMemcpyDeviceToHost)))
          {
            return fail();
          }
        }

        // Determine best Plato on host from rms_profile
        double        best_rms = -1.0;
        std::uint32_t best_q   = query_begin;
        bool          best_set = false;

        for (std::uint32_t i = 0; i < query_count; ++i)
        {
          double const r = rms_profile[i];
          if (!std::isnan(r) && r >= 0.0 && (!best_set || r > best_rms))
          {
            best_set = true;
            best_rms = r;
            best_q   = query_begin + i;
          }
        }

        if (plato_index)
          *plato_index = best_q;
        if (plato_rms)
          *plato_rms = best_set ? best_rms
                                : std::numeric_limits<double>::quiet_NaN();

        // If Plato distance profiles are requested, compute and copy for winning query best_q
        if (best_set && (plato_dist_aa || plato_dist_ab))
        {
          if (plato_dist_aa)
          {
            rfcp_init_qt_aa_kernel<<<grd_pos, blk>>>(
                d_pos, window_size, best_q, pos_profiles, d_qt_aa_0);
            rfcp_distances_aa_kernel<<<grd_pos, blk>>>(
                d_qt_aa_0, pos_profiles,
                d_mean_pos, d_sig_pos, d_valid_pos,
                window_size, best_q, exclusion_radius, d_dist_aa, nullptr);
            if (!cuda_ok(cudaMemcpy(plato_dist_aa, d_dist_aa, pos_prof_bytes, cudaMemcpyDeviceToHost)))
              return fail();
          }
          if (plato_dist_ab)
          {
            rfcp_init_qt_ab_kernel<<<grd_neg, blk>>>(
                d_pos, d_neg, window_size, best_q, neg_profiles, d_qt_ab_0);
            rfcp_distances_ab_kernel<<<grd_neg, blk>>>(
                d_qt_ab_0, neg_profiles,
                d_mean_pos, d_sig_pos,
                d_mean_neg, d_sig_neg, d_valid_neg,
                window_size, best_q, d_dist_ab, nullptr);
            if (!cuda_ok(cudaMemcpy(plato_dist_ab, d_dist_ab, neg_prof_bytes, cudaMemcpyDeviceToHost)))
              return fail();
          }
        }
        else if (!best_set)
        {
          if (plato_dist_aa)
            std::fill(plato_dist_aa, plato_dist_aa + pos_profiles, 1e300);
          if (plato_dist_ab)
            std::fill(plato_dist_ab, plato_dist_ab + neg_profiles, 1e300);
        }

        cleanup();
        return RfcpGpuStatus::ok;
      }

    } // anonymous namespace

    // -------------------------------------------------------------------------
    // Public API
    // -------------------------------------------------------------------------

    bool rfcp_gpu_available()
    {
      int count = 0;
      return cuda_ok(cudaGetDeviceCount(&count)) && count > 0;
    }

    RfcpGpuStatus rfcp_gpu_compute(
        const double *positive,    const double *mean_positive,
        const double *sig_positive, const int *valid_positive,
        std::uint32_t positive_size,
        const double *negative,    const double *mean_negative,
        const double *sig_negative, const int *valid_negative,
        std::uint32_t negative_size,
        std::uint32_t window_size, std::uint32_t max_freq,
        std::uint32_t exclusion_radius,
        std::uint32_t query_begin, std::uint32_t query_end,
        double *rms_profile, int *paired_frequency,
        std::uint32_t *plato_index, double *plato_rms,
        double *plato_dist_aa, double *plato_dist_ab,
        double *profile_aa,   double *profile_ab,
        double *profile_rfcp,
        int *profile_index_aa, int *profile_index_ab)
    {
      if (!positive || !mean_positive || !sig_positive || !valid_positive ||
          !negative || !mean_negative || !sig_negative || !valid_negative ||
          !rms_profile || !paired_frequency ||
          window_size < 2 ||
          positive_size <= static_cast<std::uint64_t>(window_size) ||
          negative_size <= static_cast<std::uint64_t>(window_size) ||
          query_begin >= query_end)
      {
        return RfcpGpuStatus::invalid_argument;
      }
      if (!rfcp_gpu_available())
        return RfcpGpuStatus::cuda_unavailable;

      return rfcp_gpu_compute_impl(
          positive, mean_positive, sig_positive, valid_positive, positive_size,
          negative, mean_negative, sig_negative, valid_negative, negative_size,
          window_size, max_freq, exclusion_radius, query_begin, query_end,
          rms_profile, paired_frequency, plato_index, plato_rms,
          plato_dist_aa, plato_dist_ab,
          profile_aa, profile_ab, profile_rfcp, profile_index_aa, profile_index_ab);
    }

  } // namespace gpu
} // namespace matrixprofiler
