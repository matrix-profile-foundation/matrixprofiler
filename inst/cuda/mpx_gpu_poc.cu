// Optional CUDA implementation for the MPX proof of concept.
//
// It lives in inst/cuda so R source builds distribute it without requiring
// CUDA. The local CUDA build script compiles it explicitly with nvcc.

#include "mpx_gpu_poc.h"

#include <cuda_runtime.h>

#include <cmath>
#include <limits>
#include <vector>

namespace matrixprofiler
{
  namespace gpu
  {
    namespace
    {

      constexpr int threads_per_block = 256;

      bool cuda_ok(cudaError_t status) { return status == cudaSuccess; }

      __global__ void initialize_qt_kernel(const double *data, std::uint32_t window_size, std::uint32_t profile_length,
                                           double *qt_initial)
      {
        std::uint32_t const index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index >= profile_length)
          return;

        double dot_product = 0.0;
        for (std::uint32_t offset = 0; offset < window_size; offset++)
        {
          dot_product += data[offset] * data[index + offset];
        }
        qt_initial[index] = dot_product;
      }

      __global__ void update_profile_kernel(const double *data, const double *mean, const double *inverse_norm,
                                            const int *valid_window, const double *qt_previous,
                                            const double *qt_initial, double *qt_current, std::uint32_t profile_length,
                                            std::uint32_t window_size, std::uint32_t query,
                                            std::uint32_t exclusion_zone, bool direct_dot_products, double *profile,
                                            int *profile_index)
      {
        std::uint32_t const index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index >= profile_length)
          return;

        double dot_product;
        if (direct_dot_products)
        {
          dot_product = 0.0;
          for (std::uint32_t offset = 0; offset < window_size; offset++)
          {
            dot_product += (data[index + offset] - mean[index]) * (data[query + offset] - mean[query]);
          }
        }
        else if (query == 0)
        {
          dot_product = qt_initial[index];
        }
        else if (index == 0)
        {
          dot_product = qt_initial[query];
        }
        else
        {
          dot_product = qt_previous[index - 1] - data[index - 1] * data[query - 1] +
                        data[index + window_size - 1] * data[query + window_size - 1];
        }
        qt_current[index] = dot_product;

        if (!valid_window[index] || !valid_window[query] || inverse_norm[index] == 0.0 || inverse_norm[query] == 0.0 ||
            (index > query ? index - query : query - index) <= exclusion_zone)
        {
          return;
        }

        double const covariance = direct_dot_products
                                      ? dot_product
                                      : dot_product - static_cast<double>(window_size) * mean[index] * mean[query];
        double correlation = covariance * inverse_norm[index] * inverse_norm[query];
        correlation = fmax(-1.0, fmin(1.0, correlation));
        if (correlation > profile[index])
        {
          profile[index] = correlation;
          profile_index[index] = static_cast<int>(query + 1);
        }
      }

      MpxGpuStatus mpx_gpu_self_join_impl(const double *data, const double *mean, const double *inverse_norm,
                                          const int *valid_window, std::uint32_t data_size,
                                          std::uint32_t window_size, std::uint32_t exclusion_zone,
                                          bool direct_dot_products, double *matrix_profile, int *profile_index)
      {
        std::uint32_t const profile_length = data_size - window_size + 1;
        double *device_data = nullptr;
        double *device_mean = nullptr;
        double *device_inverse_norm = nullptr;
        int *device_valid_window = nullptr;
        double *device_qt_initial = nullptr;
        double *device_qt_even = nullptr;
        double *device_qt_odd = nullptr;
        double *device_profile = nullptr;
        int *device_profile_index = nullptr;
        auto cleanup = [&]()
        {
          cudaFree(device_data);
          cudaFree(device_mean);
          cudaFree(device_inverse_norm);
          cudaFree(device_valid_window);
          cudaFree(device_qt_initial);
          cudaFree(device_qt_even);
          cudaFree(device_qt_odd);
          cudaFree(device_profile);
          cudaFree(device_profile_index);
        };

        std::size_t const data_bytes = static_cast<std::size_t>(data_size) * sizeof(double);
        std::size_t const profile_bytes = static_cast<std::size_t>(profile_length) * sizeof(double);
        std::size_t const index_bytes = static_cast<std::size_t>(profile_length) * sizeof(int);
        if (!cuda_ok(cudaMalloc(&device_data, data_bytes)) || !cuda_ok(cudaMalloc(&device_mean, profile_bytes)) ||
            !cuda_ok(cudaMalloc(&device_inverse_norm, profile_bytes)) ||
            !cuda_ok(cudaMalloc(&device_valid_window, index_bytes)) ||
            !cuda_ok(cudaMalloc(&device_qt_initial, profile_bytes)) ||
            !cuda_ok(cudaMalloc(&device_qt_even, profile_bytes)) || !cuda_ok(cudaMalloc(&device_qt_odd, profile_bytes)) ||
            !cuda_ok(cudaMalloc(&device_profile, profile_bytes)) ||
            !cuda_ok(cudaMalloc(&device_profile_index, index_bytes)))
        {
          cleanup();
          return MpxGpuStatus::cuda_error;
        }

        if (!cuda_ok(cudaMemcpy(device_data, data, data_bytes, cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(device_mean, mean, profile_bytes, cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(device_inverse_norm, inverse_norm, profile_bytes, cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(device_valid_window, valid_window, index_bytes, cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemset(device_profile_index, 0xff, index_bytes)))
        {
          cleanup();
          return MpxGpuStatus::cuda_error;
        }

        std::vector<double> initial_profile(profile_length, -std::numeric_limits<double>::infinity());
        if (!cuda_ok(cudaMemcpy(device_profile, initial_profile.data(), profile_bytes, cudaMemcpyHostToDevice)))
        {
          cleanup();
          return MpxGpuStatus::cuda_error;
        }

        dim3 const block(threads_per_block);
        dim3 const grid((profile_length + threads_per_block - 1) / threads_per_block);
        initialize_qt_kernel<<<grid, block>>>(device_data, window_size, profile_length, device_qt_initial);
        if (!cuda_ok(cudaGetLastError()))
        {
          cleanup();
          return MpxGpuStatus::cuda_error;
        }

        for (std::uint32_t query = 0; query < profile_length; query++)
        {
          double const *previous = query % 2 == 0 ? device_qt_odd : device_qt_even;
          double *current = query % 2 == 0 ? device_qt_even : device_qt_odd;
          update_profile_kernel<<<grid, block>>>(device_data, device_mean, device_inverse_norm, device_valid_window,
                                                 previous, device_qt_initial, current, profile_length, window_size,
                                                 query, exclusion_zone, direct_dot_products, device_profile,
                                                 device_profile_index);
          if (!cuda_ok(cudaGetLastError()))
          {
            cleanup();
            return MpxGpuStatus::cuda_error;
          }
        }

        if (!cuda_ok(cudaMemcpy(matrix_profile, device_profile, profile_bytes, cudaMemcpyDeviceToHost)) ||
            !cuda_ok(cudaMemcpy(profile_index, device_profile_index, index_bytes, cudaMemcpyDeviceToHost)) ||
            !cuda_ok(cudaDeviceSynchronize()))
        {
          cleanup();
          return MpxGpuStatus::cuda_error;
        }

        cleanup();
        return MpxGpuStatus::ok;
      }

      __global__ void initialize_qt_ab_kernel(const double *data_a, const double *data_b, std::uint32_t window_size,
                                              std::uint32_t profile_length_a, double *qt_initial)
      {
        std::uint32_t const index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index >= profile_length_a)
          return;

        double dot_product = 0.0;
        for (std::uint32_t offset = 0; offset < window_size; offset++)
        {
          dot_product += data_a[index + offset] * data_b[offset];
        }
        qt_initial[index] = dot_product;
      }

      __global__ void initialize_qt_first_ab_kernel(const double *data_a, const double *data_b,
                                                    std::uint32_t window_size, std::uint32_t profile_length_b,
                                                    double *qt_first)
      {
        std::uint32_t const query = blockIdx.x * blockDim.x + threadIdx.x;
        if (query >= profile_length_b)
          return;

        double dot_product = 0.0;
        for (std::uint32_t offset = 0; offset < window_size; offset++)
        {
          dot_product += data_a[offset] * data_b[query + offset];
        }
        qt_first[query] = dot_product;
      }

      __global__ void update_profile_ab_kernel(
          const double *data_a, const double *mean_a, const double *inverse_norm_a, const int *valid_a,
          const double *data_b, const double *mean_b, const double *inverse_norm_b, const int *valid_b,
          const double *qt_previous, const double *qt_initial, const double *qt_first, double *qt_current,
          std::uint32_t profile_length_a, std::uint32_t window_size, std::uint32_t query, bool direct_dot_products, double *profile,
          int *profile_index)
      {
        std::uint32_t const index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index >= profile_length_a)
          return;

        double dot_product;
        if (direct_dot_products)
        {
          dot_product = 0.0;
          for (std::uint32_t offset = 0; offset < window_size; offset++)
          {
            dot_product += (data_a[index + offset] - mean_a[index]) * (data_b[query + offset] - mean_b[query]);
          }
        }
        else if (query == 0)
        {
          dot_product = qt_initial[index];
        }
        else if (index == 0)
        {
          dot_product = qt_first[query];
        }
        else
        {
          dot_product = qt_previous[index - 1] - data_a[index - 1] * data_b[query - 1] +
                        data_a[index + window_size - 1] * data_b[query + window_size - 1];
        }
        qt_current[index] = dot_product;

        if (!valid_a[index] || !valid_b[query] || inverse_norm_a[index] == 0.0 || inverse_norm_b[query] == 0.0)
          return;

        double const covariance = direct_dot_products
                                      ? dot_product
                                      : dot_product - static_cast<double>(window_size) * mean_a[index] * mean_b[query];
        double correlation = covariance * inverse_norm_a[index] * inverse_norm_b[query];
        correlation = fmax(-1.0, fmin(1.0, correlation));
        if (correlation > profile[index])
        {
          profile[index] = correlation;
          profile_index[index] = static_cast<int>(query + 1);
        }
      }

      MpxGpuStatus mpx_gpu_ab_direction_impl(
          const double *data_a, const double *mean_a, const double *inverse_norm_a, const int *valid_a,
          std::uint32_t data_size_a, const double *data_b, const double *mean_b, const double *inverse_norm_b,
          const int *valid_b, std::uint32_t data_size_b, std::uint32_t window_size, bool direct_dot_products,
          double *matrix_profile, int *profile_index)
      {
        std::uint32_t const profile_length_a = data_size_a - window_size + 1;
        std::uint32_t const profile_length_b = data_size_b - window_size + 1;
        double *device_data_a = nullptr;
        double *device_mean_a = nullptr;
        double *device_inverse_norm_a = nullptr;
        int *device_valid_a = nullptr;
        double *device_data_b = nullptr;
        double *device_mean_b = nullptr;
        double *device_inverse_norm_b = nullptr;
        int *device_valid_b = nullptr;
        double *device_qt_initial = nullptr;
        double *device_qt_first = nullptr;
        double *device_qt_even = nullptr;
        double *device_qt_odd = nullptr;
        double *device_profile = nullptr;
        int *device_profile_index = nullptr;
        auto cleanup = [&]()
        {
          cudaFree(device_data_a);
          cudaFree(device_mean_a);
          cudaFree(device_inverse_norm_a);
          cudaFree(device_valid_a);
          cudaFree(device_data_b);
          cudaFree(device_mean_b);
          cudaFree(device_inverse_norm_b);
          cudaFree(device_valid_b);
          cudaFree(device_qt_initial);
          cudaFree(device_qt_first);
          cudaFree(device_qt_even);
          cudaFree(device_qt_odd);
          cudaFree(device_profile);
          cudaFree(device_profile_index);
        };

        std::size_t const data_a_bytes = static_cast<std::size_t>(data_size_a) * sizeof(double);
        std::size_t const data_b_bytes = static_cast<std::size_t>(data_size_b) * sizeof(double);
        std::size_t const profile_a_bytes = static_cast<std::size_t>(profile_length_a) * sizeof(double);
        std::size_t const profile_a_index_bytes = static_cast<std::size_t>(profile_length_a) * sizeof(int);
        if (!cuda_ok(cudaMalloc(&device_data_a, data_a_bytes)) ||
            !cuda_ok(cudaMalloc(&device_mean_a, profile_a_bytes)) ||
            !cuda_ok(cudaMalloc(&device_inverse_norm_a, profile_a_bytes)) ||
            !cuda_ok(cudaMalloc(&device_valid_a, profile_a_index_bytes)) ||
            !cuda_ok(cudaMalloc(&device_data_b, data_b_bytes)) ||
            !cuda_ok(cudaMalloc(&device_mean_b, static_cast<std::size_t>(profile_length_b) * sizeof(double))) ||
            !cuda_ok(cudaMalloc(&device_inverse_norm_b, static_cast<std::size_t>(profile_length_b) * sizeof(double))) ||
            !cuda_ok(cudaMalloc(&device_valid_b, static_cast<std::size_t>(profile_length_b) * sizeof(int))) ||
            !cuda_ok(cudaMalloc(&device_qt_initial, profile_a_bytes)) ||
            !cuda_ok(cudaMalloc(&device_qt_first, static_cast<std::size_t>(profile_length_b) * sizeof(double))) ||
            !cuda_ok(cudaMalloc(&device_qt_even, profile_a_bytes)) || !cuda_ok(cudaMalloc(&device_qt_odd, profile_a_bytes)) ||
            !cuda_ok(cudaMalloc(&device_profile, profile_a_bytes)) ||
            !cuda_ok(cudaMalloc(&device_profile_index, profile_a_index_bytes)))
        {
          cleanup();
          return MpxGpuStatus::cuda_error;
        }

        if (!cuda_ok(cudaMemcpy(device_data_a, data_a, data_a_bytes, cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(device_mean_a, mean_a, profile_a_bytes, cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(device_inverse_norm_a, inverse_norm_a, profile_a_bytes, cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(device_valid_a, valid_a, profile_a_index_bytes, cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(device_data_b, data_b, data_b_bytes, cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(device_mean_b, mean_b, static_cast<std::size_t>(profile_length_b) * sizeof(double),
                                cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(device_inverse_norm_b, inverse_norm_b,
                                static_cast<std::size_t>(profile_length_b) * sizeof(double), cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemcpy(device_valid_b, valid_b, static_cast<std::size_t>(profile_length_b) * sizeof(int),
                                cudaMemcpyHostToDevice)) ||
            !cuda_ok(cudaMemset(device_profile_index, 0xff, profile_a_index_bytes)))
        {
          cleanup();
          return MpxGpuStatus::cuda_error;
        }

        std::vector<double> initial_profile(profile_length_a, -std::numeric_limits<double>::infinity());
        if (!cuda_ok(cudaMemcpy(device_profile, initial_profile.data(), profile_a_bytes, cudaMemcpyHostToDevice)))
        {
          cleanup();
          return MpxGpuStatus::cuda_error;
        }

        dim3 const block(threads_per_block);
        dim3 const grid((profile_length_a + threads_per_block - 1) / threads_per_block);
        initialize_qt_ab_kernel<<<grid, block>>>(device_data_a, device_data_b, window_size, profile_length_a,
                                                 device_qt_initial);
        dim3 const first_grid((profile_length_b + threads_per_block - 1) / threads_per_block);
        initialize_qt_first_ab_kernel<<<first_grid, block>>>(device_data_a, device_data_b, window_size, profile_length_b,
                                                             device_qt_first);
        if (!cuda_ok(cudaGetLastError()))
        {
          cleanup();
          return MpxGpuStatus::cuda_error;
        }

        for (std::uint32_t query = 0; query < profile_length_b; query++)
        {
          double const *previous = query % 2 == 0 ? device_qt_odd : device_qt_even;
          double *current = query % 2 == 0 ? device_qt_even : device_qt_odd;
          update_profile_ab_kernel<<<grid, block>>>(
              device_data_a, device_mean_a, device_inverse_norm_a, device_valid_a, device_data_b, device_mean_b,
              device_inverse_norm_b, device_valid_b, previous, device_qt_initial, device_qt_first, current, profile_length_a,
              window_size, query, direct_dot_products, device_profile, device_profile_index);
          if (!cuda_ok(cudaGetLastError()))
          {
            cleanup();
            return MpxGpuStatus::cuda_error;
          }
        }

        if (!cuda_ok(cudaMemcpy(matrix_profile, device_profile, profile_a_bytes, cudaMemcpyDeviceToHost)) ||
            !cuda_ok(cudaMemcpy(profile_index, device_profile_index, profile_a_index_bytes, cudaMemcpyDeviceToHost)) ||
            !cuda_ok(cudaDeviceSynchronize()))
        {
          cleanup();
          return MpxGpuStatus::cuda_error;
        }

        cleanup();
        return MpxGpuStatus::ok;
      }

    } // namespace

    bool mpx_gpu_available()
    {
      int device_count = 0;
      return cuda_ok(cudaGetDeviceCount(&device_count)) && device_count > 0;
    }

    MpxGpuStatus mpx_gpu_self_join(const double *data, std::uint32_t data_size, std::uint32_t window_size,
                                   std::uint32_t exclusion_zone, double *matrix_profile, int *profile_index)
    {
      if (data == nullptr || matrix_profile == nullptr || profile_index == nullptr || window_size < 2 ||
          window_size >= data_size)
      {
        return MpxGpuStatus::invalid_argument;
      }
      if (!mpx_gpu_available())
        return MpxGpuStatus::cuda_unavailable;

      for (std::uint32_t index = 0; index < data_size; index++)
      {
        if (!std::isfinite(data[index]))
          return MpxGpuStatus::invalid_argument;
      }

      double const origin = data[0];
      std::vector<double> centered_data(data_size);
      for (std::uint32_t index = 0; index < data_size; index++)
      {
        centered_data[index] = data[index] - origin;
      }

      std::uint32_t const profile_length = data_size - window_size + 1;
      std::vector<double> means(profile_length);
      std::vector<double> inverse_norms(profile_length);
      std::vector<int> valid_window(profile_length, 1);
      for (std::uint32_t index = 0; index < profile_length; index++)
      {
        double mean = 0.0;
        for (std::uint32_t offset = 0; offset < window_size; offset++)
          mean += centered_data[index + offset];
        mean /= window_size;
        means[index] = mean;

        double squared_norm = 0.0;
        for (std::uint32_t offset = 0; offset < window_size; offset++)
        {
          double const centered = centered_data[index + offset] - mean;
          squared_norm += centered * centered;
        }
        inverse_norms[index] = squared_norm > 0.0 ? 1.0 / std::sqrt(squared_norm) : 0.0;
      }

      return mpx_gpu_self_join_impl(centered_data.data(), means.data(), inverse_norms.data(), valid_window.data(),
                                    data_size, window_size, exclusion_zone, false, matrix_profile, profile_index);
    }

    MpxGpuStatus mpx_gpu_self_join_precomputed(const double *data, const double *mean, const double *inverse_norm,
                                               const int *valid_window, std::uint32_t data_size,
                                               std::uint32_t window_size, std::uint32_t exclusion_zone,
                                               bool direct_dot_products, double *matrix_profile, int *profile_index)
    {
      if (data == nullptr || mean == nullptr || inverse_norm == nullptr || valid_window == nullptr ||
          matrix_profile == nullptr || profile_index == nullptr || window_size < 2 || window_size >= data_size)
      {
        return MpxGpuStatus::invalid_argument;
      }
      if (!mpx_gpu_available())
        return MpxGpuStatus::cuda_unavailable;
      return mpx_gpu_self_join_impl(data, mean, inverse_norm, valid_window, data_size, window_size, exclusion_zone,
                                    direct_dot_products, matrix_profile, profile_index);
    }

    MpxGpuStatus mpx_gpu_ab_join_precomputed(
        const double *data_a, const double *mean_a, const double *inverse_norm_a, const int *valid_a,
        std::uint32_t data_size_a, const double *data_b, const double *mean_b, const double *inverse_norm_b,
        const int *valid_b, std::uint32_t data_size_b, std::uint32_t window_size, bool direct_dot_products,
        double *matrix_profile_a, int *profile_index_a, double *matrix_profile_b, int *profile_index_b)
    {
      if (data_a == nullptr || mean_a == nullptr || inverse_norm_a == nullptr || valid_a == nullptr || data_b == nullptr ||
          mean_b == nullptr || inverse_norm_b == nullptr || valid_b == nullptr || matrix_profile_a == nullptr ||
          profile_index_a == nullptr || matrix_profile_b == nullptr || profile_index_b == nullptr || window_size < 2 ||
          window_size >= data_size_a || window_size >= data_size_b)
      {
        return MpxGpuStatus::invalid_argument;
      }
      if (!mpx_gpu_available())
        return MpxGpuStatus::cuda_unavailable;

      MpxGpuStatus const ab = mpx_gpu_ab_direction_impl(
          data_a, mean_a, inverse_norm_a, valid_a, data_size_a, data_b, mean_b, inverse_norm_b, valid_b, data_size_b,
          window_size, direct_dot_products, matrix_profile_a, profile_index_a);
      if (ab != MpxGpuStatus::ok)
        return ab;
      return mpx_gpu_ab_direction_impl(data_b, mean_b, inverse_norm_b, valid_b, data_size_b, data_a, mean_a,
                                       inverse_norm_a, valid_a, data_size_a, window_size, direct_dot_products,
                                       matrix_profile_b, profile_index_b);
    }

  } // namespace gpu
} // namespace matrixprofiler
