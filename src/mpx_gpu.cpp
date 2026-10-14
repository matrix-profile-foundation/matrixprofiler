#include "mpx_gpu.h"
#include "windowfunc.h"

#include <cfloat>
#include <cmath>
#include <limits>

#ifdef MATRIXPROFILER_ENABLE_CUDA
#include "../inst/cuda/mpx_gpu_poc.h"
#endif

// Optional CUDA proof-of-concept entry point. It is intentionally separate
// from mpx_rcpp: default package builds remain CUDA-free and CRAN-compatible.
// [[Rcpp::export]]
List mpx_gpu_rcpp(NumericVector data_ref, uint64_t window_size, double ez, bool idxs, bool euclidean) {
#ifndef MATRIXPROFILER_ENABLE_CUDA
  Rcpp::stop("CUDA support is not enabled; rebuild with MATRIXPROFILER_ENABLE_CUDA=1");
#else
  uint64_t const data_size = data_ref.length();
  if (window_size < 2 || window_size >= data_size) {
    Rcpp::stop("window_size must leave at least two subsequences");
  }
  if (data_size > std::numeric_limits<uint32_t>::max()) {
    Rcpp::stop("CUDA MPX currently supports at most 2^32 - 1 samples");
  }
  if (!std::isfinite(ez) || ez < 0.0) {
    Rcpp::stop("ez must be finite and non-negative");
  }
  for (R_xlen_t index = 0; index < data_ref.size(); index++) {
    if (!std::isfinite(data_ref[index])) {
      Rcpp::stop("CUDA MPX proof of concept supports finite inputs only");
    }
  }

  uint32_t const profile_length = static_cast<uint32_t>(data_size - window_size + 1);
  uint32_t const exclusion_radius = static_cast<uint32_t>(round(window_size * ez + DBL_EPSILON));
  NumericVector matrix_profile(profile_length);
  IntegerVector profile_index(profile_length, -1);
  matrixprofiler::gpu::MpxGpuStatus const status = matrixprofiler::gpu::mpx_gpu_self_join(
      data_ref.begin(), static_cast<uint32_t>(data_size), static_cast<uint32_t>(window_size), exclusion_radius,
      matrix_profile.begin(), profile_index.begin());
  if (status == matrixprofiler::gpu::MpxGpuStatus::cuda_unavailable) {
    Rcpp::stop("no CUDA-capable GPU is available");
  }
  if (status != matrixprofiler::gpu::MpxGpuStatus::ok) {
    Rcpp::stop("CUDA MPX computation failed");
  }

  for (uint32_t index = 0; index < profile_length; index++) {
    matrix_profile[index] = std::max(-1.0, std::min(1.0, matrix_profile[index]));
    if (euclidean) {
      matrix_profile[index] = sqrt(std::max(0.0, 2.0 * window_size * (1.0 - matrix_profile[index])));
    }
  }

  if (idxs) {
    return List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("profile_index") = profile_index,
                        Rcpp::Named("partial") = false);
  }
  return List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("partial") = false);
#endif
}

// Optional CUDA NA-aware self-join. Stable statistics and validity semantics
// are shared with mpx_na_rcpp; the quadratic recurrence runs on the GPU.
// [[Rcpp::export]]
List mpx_na_gpu_rcpp(NumericVector data_ref, uint64_t window_size, double ez, bool idxs, bool euclidean) {
#ifndef MATRIXPROFILER_ENABLE_CUDA
  Rcpp::stop("CUDA support is not enabled; rebuild with MATRIXPROFILER_ENABLE_CUDA=1");
#else
  uint64_t const data_size = data_ref.length();
  if (window_size < 2 || window_size >= data_size) {
    Rcpp::stop("window_size must leave at least two subsequences");
  }
  if (data_size > std::numeric_limits<uint32_t>::max()) {
    Rcpp::stop("CUDA MPX currently supports at most 2^32 - 1 samples");
  }
  if (!std::isfinite(ez) || ez < 0.0) {
    Rcpp::stop("ez must be finite and non-negative");
  }

  uint32_t const profile_length = static_cast<uint32_t>(data_size - window_size + 1);
  uint32_t const exclusion_radius = static_cast<uint32_t>(round(window_size * ez + DBL_EPSILON));
  double min_value = R_PosInf;
  double max_value = R_NegInf;
  double max_abs = 0.0;
  for (R_xlen_t index = 0; index < data_ref.size(); index++) {
    if (std::isfinite(data_ref[index])) {
      min_value = std::min(min_value, data_ref[index]);
      max_value = std::max(max_value, data_ref[index]);
      max_abs = std::max(max_abs, std::abs(data_ref[index]));
    }
  }
  double const dynamic_range = max_value - min_value;
  bool const direct_dot_products = max_abs / std::max(1.0, dynamic_range) > 1000.0;
  List const stats = muinvn_na(data_ref, static_cast<uint32_t>(window_size));
  NumericVector data = stats["data"];
  NumericVector mean = stats["avg"];
  NumericVector inverse_norm = stats["sig"];
  LogicalVector valid_window = stats["valid_window"];
  NumericVector matrix_profile(profile_length);
  IntegerVector profile_index(profile_length, NA_INTEGER);

  matrixprofiler::gpu::MpxGpuStatus const status = matrixprofiler::gpu::mpx_gpu_self_join_precomputed(
      data.begin(), mean.begin(), inverse_norm.begin(), valid_window.begin(), static_cast<uint32_t>(data_size),
      static_cast<uint32_t>(window_size), exclusion_radius, direct_dot_products, matrix_profile.begin(),
      profile_index.begin());
  if (status == matrixprofiler::gpu::MpxGpuStatus::cuda_unavailable) {
    Rcpp::stop("no CUDA-capable GPU is available");
  }
  if (status != matrixprofiler::gpu::MpxGpuStatus::ok) {
    Rcpp::stop("CUDA MPX computation failed");
  }

  for (uint32_t index = 0; index < profile_length; index++) {
    if (!valid_window[index] || !std::isfinite(matrix_profile[index])) {
      matrix_profile[index] = NA_REAL;
      profile_index[index] = NA_INTEGER;
      continue;
    }
    matrix_profile[index] = std::max(-1.0, std::min(1.0, matrix_profile[index]));
    if (euclidean) {
      matrix_profile[index] = sqrt(std::max(0.0, 2.0 * window_size * (1.0 - matrix_profile[index])));
    }
  }

  if (idxs) {
    return List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("profile_index") = profile_index,
                        Rcpp::Named("valid_window") = valid_window, Rcpp::Named("partial") = false);
  }
  return List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("valid_window") = valid_window,
                      Rcpp::Named("partial") = false);
#endif
}

// Optional CUDA NA-aware AB join. Both directions are evaluated on the GPU so
// the returned profiles have the same shape as mpxab_na_rcpp_parallel().
// [[Rcpp::export]]
List mpxab_na_gpu_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, bool idxs,
                       bool euclidean) {
#ifndef MATRIXPROFILER_ENABLE_CUDA
  Rcpp::stop("CUDA support is not enabled; rebuild with MATRIXPROFILER_ENABLE_CUDA=1");
#else
  uint64_t const data_size = data_ref.length();
  uint64_t const query_size = query_ref.length();
  if (window_size < 2 || window_size >= data_size || window_size >= query_size) {
    Rcpp::stop("window_size must leave at least two subsequences in each series");
  }
  if (data_size > std::numeric_limits<uint32_t>::max() || query_size > std::numeric_limits<uint32_t>::max()) {
    Rcpp::stop("CUDA MPX currently supports at most 2^32 - 1 samples per series");
  }

  auto requires_direct_dot_products = [](const NumericVector &series) {
    double min_value = R_PosInf;
    double max_value = R_NegInf;
    double max_abs = 0.0;
    for (R_xlen_t index = 0; index < series.size(); index++) {
      if (std::isfinite(series[index])) {
        min_value = std::min(min_value, series[index]);
        max_value = std::max(max_value, series[index]);
        max_abs = std::max(max_abs, std::abs(series[index]));
      }
    }
    return max_abs / std::max(1.0, max_value - min_value) > 1000.0;
  };

  List const stats_a = muinvn_na(data_ref, static_cast<uint32_t>(window_size));
  List const stats_b = muinvn_na(query_ref, static_cast<uint32_t>(window_size));
  NumericVector data_a = stats_a["data"];
  NumericVector mean_a = stats_a["avg"];
  NumericVector inverse_norm_a = stats_a["sig"];
  LogicalVector valid_a = stats_a["valid_window"];
  NumericVector data_b = stats_b["data"];
  NumericVector mean_b = stats_b["avg"];
  NumericVector inverse_norm_b = stats_b["sig"];
  LogicalVector valid_b = stats_b["valid_window"];
  uint32_t const profile_length_a = static_cast<uint32_t>(data_size - window_size + 1);
  uint32_t const profile_length_b = static_cast<uint32_t>(query_size - window_size + 1);
  NumericVector matrix_profile(profile_length_a);
  NumericVector mpb(profile_length_b);
  IntegerVector profile_index(profile_length_a, NA_INTEGER);
  IntegerVector pib(profile_length_b, NA_INTEGER);

  matrixprofiler::gpu::MpxGpuStatus const status = matrixprofiler::gpu::mpx_gpu_ab_join_precomputed(
      data_a.begin(), mean_a.begin(), inverse_norm_a.begin(), valid_a.begin(), static_cast<uint32_t>(data_size),
      data_b.begin(), mean_b.begin(), inverse_norm_b.begin(), valid_b.begin(), static_cast<uint32_t>(query_size),
      static_cast<uint32_t>(window_size),
      requires_direct_dot_products(data_ref) || requires_direct_dot_products(query_ref), matrix_profile.begin(),
      profile_index.begin(), mpb.begin(), pib.begin());
  if (status == matrixprofiler::gpu::MpxGpuStatus::cuda_unavailable) {
    Rcpp::stop("no CUDA-capable GPU is available");
  }
  if (status != matrixprofiler::gpu::MpxGpuStatus::ok) {
    Rcpp::stop("CUDA MPX AB computation failed");
  }

  auto finalize = [window_size, euclidean](NumericVector &profile, IntegerVector &index, const LogicalVector &valid) {
    for (R_xlen_t position = 0; position < profile.size(); position++) {
      if (!valid[position] || !std::isfinite(profile[position])) {
        profile[position] = NA_REAL;
        index[position] = NA_INTEGER;
        continue;
      }
      profile[position] = std::max(-1.0, std::min(1.0, profile[position]));
      if (euclidean) {
        profile[position] = sqrt(std::max(0.0, 2.0 * window_size * (1.0 - profile[position])));
      }
    }
  };
  finalize(matrix_profile, profile_index, valid_a);
  finalize(mpb, pib, valid_b);

  if (idxs) {
    return List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("profile_index") = profile_index,
                        Rcpp::Named("mpb") = mpb, Rcpp::Named("pib") = pib, Rcpp::Named("valid_window_a") = valid_a,
                        Rcpp::Named("valid_window_b") = valid_b, Rcpp::Named("partial") = false);
  }
  return List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("mpb") = mpb,
                      Rcpp::Named("valid_window_a") = valid_a, Rcpp::Named("valid_window_b") = valid_b,
                      Rcpp::Named("partial") = false);
#endif
}

// Optional CUDA finite AB join. Constant windows remain unsupported by this
// entry point; callers needing them should use mpxab_na_gpu_rcpp().
// [[Rcpp::export]]
List mpxab_gpu_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, bool idxs, bool euclidean) {
  for (R_xlen_t index = 0; index < data_ref.size(); index++) {
    if (!std::isfinite(data_ref[index]))
      Rcpp::stop("CUDA MPX proof of concept supports finite inputs only");
  }
  for (R_xlen_t index = 0; index < query_ref.size(); index++) {
    if (!std::isfinite(query_ref[index]))
      Rcpp::stop("CUDA MPX proof of concept supports finite inputs only");
  }
  List result = mpxab_na_gpu_rcpp(data_ref, query_ref, window_size, idxs, euclidean);
  if (idxs) {
    return List::create(Rcpp::Named("matrix_profile") = result["matrix_profile"],
                        Rcpp::Named("profile_index") = result["profile_index"], Rcpp::Named("mpb") = result["mpb"],
                        Rcpp::Named("pib") = result["pib"], Rcpp::Named("partial") = result["partial"]);
  }
  return List::create(Rcpp::Named("matrix_profile") = result["matrix_profile"], Rcpp::Named("mpb") = result["mpb"],
                      Rcpp::Named("partial") = result["partial"]);
}

//' Check if package was compiled with CUDA support
//'
//' @description
//' Queries whether `matrixprofiler` was compiled with NVIDIA CUDA GPU acceleration
//' (`MATRIXPROFILER_ENABLE_CUDA=1`).
//'
//' @details
//' By default, `matrixprofiler` is compiled without CUDA to maintain portability and
//' compliance with CRAN policies, running all algorithms on multi-threaded CPU engines
//' powered by `RcppThread`.
//'
//' If an NVIDIA GPU and the CUDA Toolkit are available on the host system, the package
//' can be recompiled with CUDA acceleration for supported algorithms (such as \code{\link{rfcp}}
//' and \code{\link{mpx}}) using the helper script \code{inst/cuda/build_cuda_rcpp.sh}
//' or by setting the environment variables:
//' \itemize{
//'   \item \code{MATRIXPROFILER_ENABLE_CUDA=1}
//'   \item \code{MATRIXPROFILER_CUDA_OBJECT} pointing to precompiled CUDA object files (\code{.o})
//'   \item \code{MATRIXPROFILER_CUDA_LIBS} pointing to the CUDA runtime library (e.g. \code{-lcudart})
//' }
//'
//' When `has_cuda()` returns `FALSE`, attempting to execute an algorithm with
//' `backend = "cuda"` will stop with an informative error instructing how to enable CUDA.
//'
//' @return A logical value: `TRUE` if CUDA support was enabled at compile time, `FALSE` otherwise.
//' @export
//' @examples
//' has_cuda()
// [[Rcpp::export]]
bool has_cuda() {
#ifdef MATRIXPROFILER_ENABLE_CUDA
  return true;
#else
  return false;
#endif
}
