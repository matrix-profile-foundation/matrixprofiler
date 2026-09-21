#include "rfcp_gpu.h"
#include "windowfunc.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifdef MATRIXPROFILER_ENABLE_CUDA
#include "../inst/cuda/rfcp_gpu_poc.h"
#endif

// Optional CUDA entry point for the RFCP algorithm.  It is intentionally
// separate from rfcp_na_segmented_native_rcpp_parallel: default package
// builds remain CUDA-free and CRAN-compatible.
//
// query_begin is zero-based and query_end is exclusive.  The public R wrapper
// (rfcp.R) applies the one-based → zero-based normalisation before calling
// this function, matching the convention used by the CPU parallel path.
//
// The returned list has the same structure as rfcp_na_segmented_native_rcpp_parallel
// so that the R wrapper can be backend-transparent.
// [[Rcpp::export]]
List rfcp_gpu_rcpp(NumericVector positive_ref, NumericVector negative_ref,
                   uint64_t window_size, uint32_t max_freq,
                   uint32_t exclusion_radius,
                   uint64_t query_begin, uint64_t query_end,
                   bool return_profiles) {
#ifndef MATRIXPROFILER_ENABLE_CUDA
  Rcpp::stop("CUDA support is not enabled; rebuild with MATRIXPROFILER_ENABLE_CUDA=1");
#else
  uint64_t const positive_size = positive_ref.length();
  uint64_t const negative_size = negative_ref.length();
  if (window_size < 2 || window_size >= positive_size || window_size >= negative_size) {
    Rcpp::stop("window_size must fit both series and be at least 2");
  }
  if (positive_size > std::numeric_limits<uint32_t>::max() ||
      negative_size > std::numeric_limits<uint32_t>::max()) {
    Rcpp::stop("RFCP GPU currently supports series shorter than 2^32 samples");
  }
  if (max_freq < 1) {
    Rcpp::stop("max_freq must be at least 1");
  }

  uint32_t const pos_profiles = static_cast<uint32_t>(positive_size - window_size + 1);
  uint32_t const neg_profiles = static_cast<uint32_t>(negative_size - window_size + 1);
  if (query_begin >= query_end || query_end > pos_profiles) {
    Rcpp::stop("query range must be a non-empty subset of positive window starts");
  }

  uint32_t const query_count = static_cast<uint32_t>(query_end - query_begin);
  uint64_t const profile_cells = static_cast<uint64_t>(query_count) * max_freq;
  constexpr uint64_t rfcp_max_profile_cells = 5000000;
  if (return_profiles && profile_cells > rfcp_max_profile_cells) {
    Rcpp::stop("return_profiles exceeds the limit of 5,000,000 rank-query cells");
  }

  // Precompute stable, recurrence-safe statistics — identical to the CPU NA path.
  List const pos_stats = muinvn_na(positive_ref, static_cast<uint32_t>(window_size));
  List const neg_stats = muinvn_na(negative_ref, static_cast<uint32_t>(window_size));
  NumericVector positive       = pos_stats["data"];
  NumericVector mean_positive  = pos_stats["avg"];
  NumericVector sig_positive   = pos_stats["sig"];
  LogicalVector valid_positive = pos_stats["valid_window"];
  NumericVector negative       = neg_stats["data"];
  NumericVector mean_negative  = neg_stats["avg"];
  NumericVector sig_negative   = neg_stats["sig"];
  LogicalVector valid_negative = neg_stats["valid_window"];

  NumericVector rms_profile(query_count, NA_REAL);
  IntegerVector paired_frequency(query_count, 0);

  NumericVector profile_aa;
  NumericVector profile_ab;
  NumericVector profile_rfcp;
  IntegerVector profile_index_aa;
  IntegerVector profile_index_ab;
  if (return_profiles) {
    profile_aa = NumericVector(profile_cells, NA_REAL);
    profile_ab = NumericVector(profile_cells, NA_REAL);
    profile_rfcp = NumericVector(profile_cells, NA_REAL);
    profile_index_aa = IntegerVector(profile_cells, NA_INTEGER);
    profile_index_ab = IntegerVector(profile_cells, NA_INTEGER);
  }

  // Plato distance buffers: filled by the GPU when a new best query is found.
  std::vector<double> plato_dist_aa(pos_profiles, 1e300);
  std::vector<double> plato_dist_ab(neg_profiles, 1e300);
  std::uint32_t plato_idx = 0;
  double        plato_rms_val = std::numeric_limits<double>::quiet_NaN();

  matrixprofiler::gpu::RfcpGpuStatus const status = matrixprofiler::gpu::rfcp_gpu_compute(
      positive.begin(), mean_positive.begin(), sig_positive.begin(), valid_positive.begin(),
      static_cast<uint32_t>(positive_size),
      negative.begin(), mean_negative.begin(), sig_negative.begin(), valid_negative.begin(),
      static_cast<uint32_t>(negative_size),
      static_cast<uint32_t>(window_size), max_freq, exclusion_radius,
      static_cast<uint32_t>(query_begin), static_cast<uint32_t>(query_end),
      rms_profile.begin(), paired_frequency.begin(),
      &plato_idx, &plato_rms_val,
      plato_dist_aa.data(), plato_dist_ab.data(),
      return_profiles ? profile_aa.begin() : nullptr,
      return_profiles ? profile_ab.begin() : nullptr,
      return_profiles ? profile_rfcp.begin() : nullptr,
      return_profiles ? profile_index_aa.begin() : nullptr,
      return_profiles ? profile_index_ab.begin() : nullptr);

  if (status == matrixprofiler::gpu::RfcpGpuStatus::cuda_unavailable) {
    Rcpp::stop("no CUDA-capable GPU is available");
  }
  if (status != matrixprofiler::gpu::RfcpGpuStatus::ok) {
    Rcpp::stop("CUDA RFCP computation failed");
  }

  bool const plato_available = std::isfinite(plato_rms_val);

  // Effective max_freq: largest paired-frequency observed among completed queries.
  uint32_t effective_max_freq = 0;
  for (uint32_t i = 0; i < query_count; ++i) {
    effective_max_freq =
        std::max(effective_max_freq, static_cast<uint32_t>(paired_frequency[i]));
  }

  // Plato rank vectors: re-run greedy neighbour selection from saved distance arrays.
  NumericVector rfmp_aa_at_plato(max_freq, NA_REAL);
  NumericVector rfmp_ab_at_plato(max_freq, NA_REAL);
  NumericVector rfcp_at_plato(max_freq, NA_REAL);
  IntegerVector rfmpi_aa_at_plato(max_freq, NA_INTEGER);
  IntegerVector rfmpi_ab_at_plato(max_freq, NA_INTEGER);

  if (plato_available) {
    constexpr double kInfThreshold = 1e200;
    double const clip = std::sqrt(2.0 * static_cast<double>(window_size));

    // Generic greedy selection lambda (mirrors select_neighbors_host in the .cu)
    auto select = [&](const std::vector<double> &dists, uint32_t n_ref,
                      NumericVector &out_d, IntegerVector &out_i) {
      std::vector<std::pair<double, uint32_t>> cands;
      cands.reserve(n_ref);
      for (uint32_t j = 0; j < n_ref; ++j)
        if (dists[j] < kInfThreshold)
          cands.push_back({dists[j], j});
      std::sort(cands.begin(), cands.end());

      std::vector<bool> excl(n_ref, false);
      uint32_t acc = 0;
      for (auto const &cand : cands) {
        if (acc >= max_freq) break;
        uint32_t const idx = cand.second;
        if (excl[idx]) continue;
        out_d[acc] = cand.first;
        out_i[acc] = static_cast<int>(idx + 1); // 1-based
        ++acc;
        uint32_t const l = idx >= exclusion_radius ? idx - exclusion_radius : 0;
        uint32_t const r = idx + exclusion_radius < n_ref ? idx + exclusion_radius : n_ref - 1;
        for (uint32_t k = l; k <= r; ++k) excl[k] = true;
      }
    };

    select(plato_dist_aa, pos_profiles, rfmp_aa_at_plato, rfmpi_aa_at_plato);
    select(plato_dist_ab, neg_profiles, rfmp_ab_at_plato, rfmpi_ab_at_plato);

    // RFCP contrasts at plato
    for (uint32_t r = 0; r < max_freq; ++r) {
      if (!ISNA(rfmp_aa_at_plato[r]) && !ISNA(rfmp_ab_at_plato[r])) {
        double const aa = std::min((double)rfmp_aa_at_plato[r], clip);
        double const ab = std::min((double)rfmp_ab_at_plato[r], clip);
        rfcp_at_plato[r] = std::max(0.0, (ab - aa) / clip);
      }
    }
  }

  List result = List::create(
      Rcpp::Named("rms_profile")          = rms_profile,
      Rcpp::Named("valid_window_positive") = valid_positive,
      Rcpp::Named("valid_window_negative") = valid_negative,
      Rcpp::Named("plato_index")           = plato_available
                                               ? Rcpp::wrap(static_cast<double>(plato_idx + 1))
                                               : Rcpp::wrap(NA_REAL),
      Rcpp::Named("plato_rms")             = plato_available ? plato_rms_val : NA_REAL,
      Rcpp::Named("plato_twin_index")      = plato_available && rfmpi_aa_at_plato[0] != NA_INTEGER
                                               ? Rcpp::wrap(rfmpi_aa_at_plato[0])
                                               : Rcpp::wrap(NA_INTEGER),
      Rcpp::Named("rfcp_at_plato")         = rfcp_at_plato,
      Rcpp::Named("rfmp_aa_at_plato")      = rfmp_aa_at_plato,
      Rcpp::Named("rfmp_ab_at_plato")      = rfmp_ab_at_plato,
      Rcpp::Named("rfmpi_aa_at_plato")     = rfmpi_aa_at_plato,
      Rcpp::Named("rfmpi_ab_at_plato")     = rfmpi_ab_at_plato,
      Rcpp::Named("requested_max_freq")    = static_cast<int>(max_freq),
      Rcpp::Named("effective_max_freq")    = static_cast<int>(effective_max_freq),
      Rcpp::Named("query_begin")           = static_cast<double>(query_begin + 1),
      Rcpp::Named("query_end")             = static_cast<double>(query_end),
      Rcpp::Named("completed_query_end")   = static_cast<double>(query_end),
      Rcpp::Named("partial")               = false,
      Rcpp::Named("query_rows_per_task")   = 1,
      Rcpp::Named("n_tasks")               = static_cast<int>(query_count),
      Rcpp::Named("effective_workers")     = 1);

  if (return_profiles) {
    IntegerVector const dimensions = IntegerVector::create(max_freq, query_count);
    profile_aa.attr("dim") = dimensions;
    profile_ab.attr("dim") = dimensions;
    profile_rfcp.attr("dim") = dimensions;
    profile_index_aa.attr("dim") = dimensions;
    profile_index_ab.attr("dim") = dimensions;
    result["rfmp_aa"] = profile_aa;
    result["rfmp_ab"] = profile_ab;
    result["rfcp"] = profile_rfcp;
    result["rfmpi_aa"] = profile_index_aa;
    result["rfmpi_ab"] = profile_index_ab;
  }

  return result;
#endif
}
