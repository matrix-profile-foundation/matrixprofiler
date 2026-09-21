// Optional CUDA implementation for the RFCP algorithm using CUB.
//
// It lives in inst/cuda so R source builds distribute it without requiring
// CUDA.  The local CUDA build script compiles it explicitly with nvcc.
//
// Algorithm overview (100% GPU Execution via CUB)
// ===============================================
// For each tile of queries in [query_begin, query_end):
//
//   Phase A — QT recurrence & Distance calculation (GPU):
//     Double-buffered QT rows advance diagonally on the device.
//     For each query in the batch, AA and AB Euclidean distances are computed
//     into device batch buffers `d_batch_dist_aa` and `d_batch_dist_ab`.
//
//   Phase B — Segmented Radix Sort via CUB (GPU):
//     `cub::DeviceSegmentedRadixSort::SortPairs` sorts all queries in the
//     batch simultaneously in parallel on the GPU without any PCIe transfers.
//
//   Phase C — Greedy Neighbor Selection & RMSC Contrast (GPU):
//     A specialized CUDA kernel scans the sorted candidates per query,
//     enforces the non-overlapping exclusion zone [idx - r, idx + r],
//     and calculates the contrast score and RMS value directly on the device.
//
//   Phase D — Bulk Asynchronous Transfer (GPU -> Host):
//     Only the final `rms_profile`, `paired_frequency`, and optional profile
//     matrices are transferred back to host memory, eliminating >99.9% of
//     PCIe memory traffic compared to the CPU-sorting approach.

#include "rfcp_gpu_poc.h"

#include <cuda_runtime.h>
#include <math_constants.h>
#include <cub/cub.cuh>

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

      bool cuda_ok(cudaError_t s) { return s == cudaSuccess; }

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
          std::uint32_t exclusion_radius, double *dist_aa)
      {
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
          std::uint32_t window_size, std::uint32_t q, double *dist_ab)
      {
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
      // Kernel: initialise segment indices [0, 1, ..., n-1] across batch
      // -----------------------------------------------------------------------
      __global__ void rfcp_init_indices_kernel(int *indices, std::uint32_t n, std::uint32_t total)
      {
        std::uint32_t const idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < total)
          indices[idx] = static_cast<int>(idx % n);
      }

      // -----------------------------------------------------------------------
      // Kernel: fill buffer with a constant double (e.g. 1e300 for invalid queries)
      // -----------------------------------------------------------------------
      __global__ void rfcp_fill_kernel(double *buffer, double val, std::uint32_t total)
      {
        std::uint32_t const idx = blockIdx.x * blockDim.x + threadIdx.x;
        if (idx < total)
          buffer[idx] = val;
      }

      // -----------------------------------------------------------------------
      // Kernel: GPU Greedy Neighbor Selection & Contrast / RMS Computation
      //
      // One thread per query in the batch. Scans the CUB-sorted distance/index
      // arrays, enforces the exclusion radius, computes contrast & RMSC, and
      // writes directly to device outputs.
      // -----------------------------------------------------------------------
      __global__ void rfcp_gpu_greedy_selection_kernel(
          const double * __restrict__ sorted_dist_aa,
          const int    * __restrict__ sorted_idx_aa,
          std::uint32_t pos_profiles,
          const double * __restrict__ sorted_dist_ab,
          const int    * __restrict__ sorted_idx_ab,
          std::uint32_t neg_profiles,
          const int    * __restrict__ valid_positive,
          std::uint32_t query_begin,
          std::uint32_t batch_begin,
          std::uint32_t batch_count,
          std::uint32_t max_freq,
          std::uint32_t exclusion_radius,
          double clip,
          double * __restrict__ d_rms_profile,
          int    * __restrict__ d_paired_frequency,
          double * __restrict__ d_profile_aa,
          double * __restrict__ d_profile_ab,
          double * __restrict__ d_profile_rfcp,
          int    * __restrict__ d_profile_index_aa,
          int    * __restrict__ d_profile_index_ab,
          double * __restrict__ d_batch_rms,
          std::uint32_t * __restrict__ d_batch_q)
      {
        std::uint32_t const b = blockIdx.x * blockDim.x + threadIdx.x;
        if (b >= batch_count)
          return;

        std::uint32_t const q = batch_begin + b;
        std::uint32_t const out_i = q - query_begin;

        if (!valid_positive[q])
        {
          d_rms_profile[out_i]      = CUDART_NAN;
          d_paired_frequency[out_i] = 0;
          d_batch_rms[b] = -1.0;
          d_batch_q[b]   = q;

          if (d_profile_aa)
          {
            std::uint64_t const col = static_cast<std::uint64_t>(out_i) * max_freq;
            for (std::uint32_t r = 0; r < max_freq; ++r)
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

        // Scan AA sorted segment
        double aa_dists[kMaxSupportedK];
        int    aa_idxs[kMaxSupportedK];
        std::uint32_t aa_count = 0;
        const double *seg_dist_aa = sorted_dist_aa + static_cast<std::uint64_t>(b) * pos_profiles;
        const int    *seg_idx_aa  = sorted_idx_aa  + static_cast<std::uint64_t>(b) * pos_profiles;

        for (std::uint32_t j = 0; j < pos_profiles && aa_count < max_freq; ++j)
        {
          double const d = seg_dist_aa[j];
          if (d >= 1e200 || isnan(d))
            break;
          int const idx = seg_idx_aa[j];

          bool excluded = false;
          for (std::uint32_t a = 0; a < aa_count; ++a)
          {
            int diff = idx - aa_idxs[a];
            if (diff < 0) diff = -diff;
            if (diff <= static_cast<int>(exclusion_radius))
            {
              excluded = true;
              break;
            }
          }
          if (!excluded)
          {
            aa_dists[aa_count] = d;
            aa_idxs[aa_count]  = idx;
            ++aa_count;
          }
        }

        // Scan AB sorted segment
        double ab_dists[kMaxSupportedK];
        int    ab_idxs[kMaxSupportedK];
        std::uint32_t ab_count = 0;
        const double *seg_dist_ab = sorted_dist_ab + static_cast<std::uint64_t>(b) * neg_profiles;
        const int    *seg_idx_ab  = sorted_idx_ab  + static_cast<std::uint64_t>(b) * neg_profiles;

        for (std::uint32_t j = 0; j < neg_profiles && ab_count < max_freq; ++j)
        {
          double const d = seg_dist_ab[j];
          if (d >= 1e200 || isnan(d))
            break;
          int const idx = seg_idx_ab[j];

          bool excluded = false;
          for (std::uint32_t a = 0; a < ab_count; ++a)
          {
            int diff = idx - ab_idxs[a];
            if (diff < 0) diff = -diff;
            if (diff <= static_cast<int>(exclusion_radius))
            {
              excluded = true;
              break;
            }
          }
          if (!excluded)
          {
            ab_dists[ab_count] = d;
            ab_idxs[ab_count]  = idx;
            ++ab_count;
          }
        }

        std::uint32_t const paired = (aa_count < ab_count) ? aa_count : ab_count;
        d_paired_frequency[out_i] = static_cast<int>(paired);

        // Compute RMSC
        double sq_sum = 0.0;
        for (std::uint32_t r = 0; r < paired; ++r)
        {
          double const aa = (aa_dists[r] < clip) ? aa_dists[r] : clip;
          double const ab = (ab_dists[r] < clip) ? ab_dists[r] : clip;
          double const c  = (ab > aa) ? ((ab - aa) / clip) : 0.0;
          sq_sum += c * c;
        }
        double const rms = sqrt(sq_sum / static_cast<double>(max_freq));
        d_rms_profile[out_i] = rms;
        d_batch_rms[b] = rms;
        d_batch_q[b]   = q;

        // Profile matrices if requested
        if (d_profile_aa)
        {
          std::uint64_t const col = static_cast<std::uint64_t>(out_i) * max_freq;
          for (std::uint32_t r = 0; r < max_freq; ++r)
          {
            double const aa_r = (r < aa_count) ? aa_dists[r] : CUDART_NAN;
            double const ab_r = (r < ab_count) ? ab_dists[r] : CUDART_NAN;
            d_profile_aa[col + r]       = aa_r;
            d_profile_ab[col + r]       = ab_r;
            d_profile_index_aa[col + r] = (r < aa_count) ? (aa_idxs[r] + 1) : kNaInteger;
            d_profile_index_ab[col + r] = (r < ab_count) ? (ab_idxs[r] + 1) : kNaInteger;

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
      // Core GPU implementation with CUB.
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
        if (max_freq > kMaxSupportedK)
          return RfcpGpuStatus::invalid_argument;

        std::uint32_t const pos_profiles = positive_size - window_size + 1;
        std::uint32_t const neg_profiles = negative_size - window_size + 1;
        std::uint32_t const query_count  = query_end - query_begin;
        std::uint32_t const tile_size    = rfcp_gpu_tile_size(window_size, query_count);

        // Compute adaptive CUB batch size: bound batch memory to ~16M items (~128MB)
        std::uint32_t const max_ref = std::max(pos_profiles, neg_profiles);
        std::uint32_t const max_cub_batch = std::max<std::uint32_t>(1U, std::min<std::uint32_t>(64U, 16777216U / max_ref));
        std::uint32_t const B = std::min<std::uint32_t>(query_count, max_cub_batch);

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

        // CUB batch buffers
        double *d_batch_dist_aa  = nullptr;
        double *d_batch_dist_ab  = nullptr;
        double *d_sorted_dist_aa = nullptr;
        double *d_sorted_dist_ab = nullptr;
        int    *d_batch_idx_aa   = nullptr;
        int    *d_batch_idx_ab   = nullptr;
        int    *d_sorted_idx_aa  = nullptr;
        int    *d_sorted_idx_ab  = nullptr;
        int    *d_offsets_aa     = nullptr;
        int    *d_offsets_ab     = nullptr;
        double *d_batch_rms      = nullptr;
        std::uint32_t *d_batch_q = nullptr;

        void   *d_temp_storage_aa = nullptr;
        std::size_t temp_bytes_aa = 0;
        void   *d_temp_storage_ab = nullptr;
        std::size_t temp_bytes_ab = 0;

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
          cudaFree(d_pos);          cudaFree(d_neg);
          cudaFree(d_mean_pos);     cudaFree(d_sig_pos);
          cudaFree(d_mean_neg);     cudaFree(d_sig_neg);
          cudaFree(d_valid_pos);    cudaFree(d_valid_neg);
          cudaFree(d_qt_aa_0);      cudaFree(d_qt_aa_1);
          cudaFree(d_qt_ab_0);      cudaFree(d_qt_ab_1);
          cudaFree(d_aa_boundary);  cudaFree(d_ab_boundary);
          cudaFree(d_batch_dist_aa);  cudaFree(d_batch_dist_ab);
          cudaFree(d_sorted_dist_aa); cudaFree(d_sorted_dist_ab);
          cudaFree(d_batch_idx_aa);   cudaFree(d_batch_idx_ab);
          cudaFree(d_sorted_idx_aa);  cudaFree(d_sorted_idx_ab);
          cudaFree(d_offsets_aa);     cudaFree(d_offsets_ab);
          cudaFree(d_batch_rms);      cudaFree(d_batch_q);
          cudaFree(d_temp_storage_aa); cudaFree(d_temp_storage_ab);
          cudaFree(d_rms_profile);    cudaFree(d_paired_frequency);
          cudaFree(d_profile_aa);     cudaFree(d_profile_ab);
          cudaFree(d_profile_rfcp);   cudaFree(d_profile_index_aa);
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

        std::size_t const batch_pos_dbytes = static_cast<std::size_t>(B) * pos_prof_bytes;
        std::size_t const batch_neg_dbytes = static_cast<std::size_t>(B) * neg_prof_bytes;
        std::size_t const batch_pos_ibytes = static_cast<std::size_t>(B) * pos_prof_ibytes;
        std::size_t const batch_neg_ibytes = static_cast<std::size_t>(B) * neg_prof_ibytes;

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
            !cuda_ok(cudaMalloc(&d_ab_boundary, tile_bytes)))
        {
          return fail();
        }

        // Allocate CUB batch buffers
        if (!cuda_ok(cudaMalloc(&d_batch_dist_aa,  batch_pos_dbytes)) ||
            !cuda_ok(cudaMalloc(&d_batch_dist_ab,  batch_neg_dbytes)) ||
            !cuda_ok(cudaMalloc(&d_sorted_dist_aa, batch_pos_dbytes)) ||
            !cuda_ok(cudaMalloc(&d_sorted_dist_ab, batch_neg_dbytes)) ||
            !cuda_ok(cudaMalloc(&d_batch_idx_aa,   batch_pos_ibytes)) ||
            !cuda_ok(cudaMalloc(&d_batch_idx_ab,   batch_neg_ibytes)) ||
            !cuda_ok(cudaMalloc(&d_sorted_idx_aa,  batch_pos_ibytes)) ||
            !cuda_ok(cudaMalloc(&d_sorted_idx_ab,  batch_neg_ibytes)) ||
            !cuda_ok(cudaMalloc(&d_offsets_aa,     (B + 1) * sizeof(int))) ||
            !cuda_ok(cudaMalloc(&d_offsets_ab,     (B + 1) * sizeof(int))) ||
            !cuda_ok(cudaMalloc(&d_batch_rms,      B * sizeof(double))) ||
            !cuda_ok(cudaMalloc(&d_batch_q,        B * sizeof(std::uint32_t))))
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

        // Setup CUB segment offsets
        std::vector<int> h_offsets_aa(B + 1);
        std::vector<int> h_offsets_ab(B + 1);
        for (std::uint32_t i = 0; i <= B; ++i)
        {
          h_offsets_aa[i] = static_cast<int>(i * pos_profiles);
          h_offsets_ab[i] = static_cast<int>(i * neg_profiles);
        }
        if (!cuda_ok(cudaMemcpy(d_offsets_aa, h_offsets_aa.data(), (B + 1) * sizeof(int), cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(d_offsets_ab, h_offsets_ab.data(), (B + 1) * sizeof(int), cudaMemcpyHostToDevice)))
        {
          return fail();
        }

        // Initialize initial segment indices: [0, 1, ..., N-1] per segment
        rfcp_init_indices_kernel<<<(B * pos_profiles + kTPB - 1) / kTPB, kTPB>>>(
            d_batch_idx_aa, pos_profiles, B * pos_profiles);
        rfcp_init_indices_kernel<<<(B * neg_profiles + kTPB - 1) / kTPB, kTPB>>>(
            d_batch_idx_ab, neg_profiles, B * neg_profiles);
        if (!cuda_ok(cudaGetLastError()))
          return fail();

        // Query and allocate CUB temporary storage
        cub::DeviceSegmentedRadixSort::SortPairs(
            d_temp_storage_aa, temp_bytes_aa,
            d_batch_dist_aa, d_sorted_dist_aa,
            d_batch_idx_aa, d_sorted_idx_aa,
            B * pos_profiles, B,
            d_offsets_aa, d_offsets_aa + 1);
        if (!cuda_ok(cudaMalloc(&d_temp_storage_aa, temp_bytes_aa)))
          return fail();

        cub::DeviceSegmentedRadixSort::SortPairs(
            d_temp_storage_ab, temp_bytes_ab,
            d_batch_dist_ab, d_sorted_dist_ab,
            d_batch_idx_ab, d_sorted_idx_ab,
            B * neg_profiles, B,
            d_offsets_ab, d_offsets_ab + 1);
        if (!cuda_ok(cudaMalloc(&d_temp_storage_ab, temp_bytes_ab)))
          return fail();

        // Plato tracking on host
        std::vector<double>        h_batch_rms(B);
        std::vector<std::uint32_t> h_batch_q(B);
        double        best_rms = 0.0;
        std::uint32_t best_q   = 0;
        bool          best_set = false;
        double const  clip     = std::sqrt(2.0 * static_cast<double>(window_size));

        dim3 const blk(kTPB);
        dim3 const grd_pos((pos_profiles + kTPB - 1) / kTPB);
        dim3 const grd_neg((neg_profiles + kTPB - 1) / kTPB);

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

          // Iterate through tile in sub-batches of size B
          for (std::uint32_t bb = tb; bb < te; bb += B)
          {
            std::uint32_t const be          = std::min(bb + B, te);
            std::uint32_t const batch_count = be - bb;

            // Compute distances for each query in the batch
            for (std::uint32_t q = bb; q < be; ++q)
            {
              std::uint32_t const b = q - bb;
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
                if (!cuda_ok(cudaGetLastError()))
                  return fail();
              }

              double *dist_aa_ptr = d_batch_dist_aa + static_cast<std::uint64_t>(b) * pos_profiles;
              double *dist_ab_ptr = d_batch_dist_ab + static_cast<std::uint64_t>(b) * neg_profiles;

              if (valid_positive[q])
              {
                rfcp_distances_aa_kernel<<<grd_pos, blk>>>(
                    qt_aa_curr, pos_profiles,
                    d_mean_pos, d_sig_pos, d_valid_pos,
                    window_size, q, exclusion_radius, dist_aa_ptr);
                rfcp_distances_ab_kernel<<<grd_neg, blk>>>(
                    qt_ab_curr, neg_profiles,
                    d_mean_pos, d_sig_pos,
                    d_mean_neg, d_sig_neg, d_valid_neg,
                    window_size, q, dist_ab_ptr);
              }
              else
              {
                // Fill invalid query window with large values
                rfcp_fill_kernel<<<(pos_profiles + kTPB - 1) / kTPB, kTPB>>>(dist_aa_ptr, 1e300, pos_profiles);
                rfcp_fill_kernel<<<(neg_profiles + kTPB - 1) / kTPB, kTPB>>>(dist_ab_ptr, 1e300, neg_profiles);
              }
              if (!cuda_ok(cudaGetLastError()))
                return fail();
            }

            // Reset segment indices for CUB
            rfcp_init_indices_kernel<<<(batch_count * pos_profiles + kTPB - 1) / kTPB, kTPB>>>(
                d_batch_idx_aa, pos_profiles, batch_count * pos_profiles);
            rfcp_init_indices_kernel<<<(batch_count * neg_profiles + kTPB - 1) / kTPB, kTPB>>>(
                d_batch_idx_ab, neg_profiles, batch_count * neg_profiles);
            if (!cuda_ok(cudaGetLastError()))
              return fail();

            // Segmented Radix Sort on GPU via CUB for AA
            cub::DeviceSegmentedRadixSort::SortPairs(
                d_temp_storage_aa, temp_bytes_aa,
                d_batch_dist_aa, d_sorted_dist_aa,
                d_batch_idx_aa, d_sorted_idx_aa,
                batch_count * pos_profiles, batch_count,
                d_offsets_aa, d_offsets_aa + 1);
            if (!cuda_ok(cudaGetLastError()))
              return fail();

            // Segmented Radix Sort on GPU via CUB for AB
            cub::DeviceSegmentedRadixSort::SortPairs(
                d_temp_storage_ab, temp_bytes_ab,
                d_batch_dist_ab, d_sorted_dist_ab,
                d_batch_idx_ab, d_sorted_idx_ab,
                batch_count * neg_profiles, batch_count,
                d_offsets_ab, d_offsets_ab + 1);
            if (!cuda_ok(cudaGetLastError()))
              return fail();

            // Launch GPU Greedy Selection Kernel
            dim3 const grd_batch((batch_count + kTPB - 1) / kTPB);
            rfcp_gpu_greedy_selection_kernel<<<grd_batch, blk>>>(
                d_sorted_dist_aa, d_sorted_idx_aa, pos_profiles,
                d_sorted_dist_ab, d_sorted_idx_ab, neg_profiles,
                d_valid_pos, query_begin, bb, batch_count,
                max_freq, exclusion_radius, clip,
                d_rms_profile, d_paired_frequency,
                d_profile_aa, d_profile_ab, d_profile_rfcp,
                d_profile_index_aa, d_profile_index_ab,
                d_batch_rms, d_batch_q);
            if (!cuda_ok(cudaGetLastError()))
              return fail();

            // Copy batch RMS and Q values to host to track the global Plato
            if (!cuda_ok(cudaMemcpy(h_batch_rms.data(), d_batch_rms, batch_count * sizeof(double), cudaMemcpyDeviceToHost)) ||
                !cuda_ok(cudaMemcpy(h_batch_q.data(), d_batch_q, batch_count * sizeof(std::uint32_t), cudaMemcpyDeviceToHost)))
            {
              return fail();
            }

            for (std::uint32_t b = 0; b < batch_count; ++b)
            {
              double const r = h_batch_rms[b];
              if (r >= 0.0 && (!best_set || r > best_rms))
              {
                best_set = true;
                best_rms = r;
                best_q   = h_batch_q[b];

                if (plato_dist_aa)
                {
                  if (!cuda_ok(cudaMemcpy(plato_dist_aa,
                                          d_batch_dist_aa + static_cast<std::uint64_t>(b) * pos_profiles,
                                          pos_prof_bytes, cudaMemcpyDeviceToHost)))
                    return fail();
                }
                if (plato_dist_ab)
                {
                  if (!cuda_ok(cudaMemcpy(plato_dist_ab,
                                          d_batch_dist_ab + static_cast<std::uint64_t>(b) * neg_profiles,
                                          neg_prof_bytes, cudaMemcpyDeviceToHost)))
                    return fail();
                }
              }
            }
          } // sub-batch loop
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

        if (plato_index)
          *plato_index = best_q;
        if (plato_rms)
          *plato_rms = best_set ? best_rms
                                : std::numeric_limits<double>::quiet_NaN();

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
