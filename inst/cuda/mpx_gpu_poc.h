#ifndef MATRIXPROFILER_MPX_GPU_POC_H
#define MATRIXPROFILER_MPX_GPU_POC_H

#include <cstdint>

namespace matrixprofiler
{
  namespace gpu
  {

    enum class MpxGpuStatus : int
    {
      ok = 0,
      invalid_argument = 1,
      cuda_unavailable = 2,
      cuda_error = 3,
    };

    // Computes an exact finite, nonconstant-window MPX self-join in correlation
    // space. The caller owns output buffers of length data_size - window_size +
    // 1. Indices are one-based, matching matrixprofiler's R-facing convention.
    MpxGpuStatus mpx_gpu_self_join(const double *data, std::uint32_t data_size, std::uint32_t window_size,
                                   std::uint32_t exclusion_zone, double *matrix_profile, int *profile_index);

    // Computes an exact self-join from package-precomputed, stable statistics.
    // `data` must be recurrence-safe (non-finite samples replaced with zero),
    // and `valid_window` controls which profile updates are eligible.
    MpxGpuStatus mpx_gpu_self_join_precomputed(const double *data, const double *mean, const double *inverse_norm,
                                               const int *valid_window, std::uint32_t data_size,
                                               std::uint32_t window_size, std::uint32_t exclusion_zone,
                                               bool direct_dot_products, double *matrix_profile, int *profile_index);

    // Computes both sides of an AB join from package-precomputed statistics.
    // The two directional passes avoid device-wide reductions when updating
    // profile indices; outputs retain matrixprofiler's one-based convention.
    MpxGpuStatus mpx_gpu_ab_join_precomputed(
        const double *data_a, const double *mean_a, const double *inverse_norm_a, const int *valid_a,
        std::uint32_t data_size_a, const double *data_b, const double *mean_b, const double *inverse_norm_b,
        const int *valid_b, std::uint32_t data_size_b, std::uint32_t window_size, bool direct_dot_products,
        double *matrix_profile_a, int *profile_index_a, double *matrix_profile_b, int *profile_index_b);

    bool mpx_gpu_available();

  } // namespace gpu
} // namespace matrixprofiler

#endif // MATRIXPROFILER_MPX_GPU_POC_H
