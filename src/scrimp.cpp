// This is a personal academic project. Dear PVS-Studio, please check it.

#include "mathtools.h" // math first to fix OSX error
#include "scrimp.h"
#include "mpx.h"
#include "mpx_parallel.h"
#include "windowfunc.h"
#include <cfloat>
#include <cmath>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

namespace {

void validate_scrimp_inputs(const NumericVector &data, const NumericVector &query, uint32_t window_size, double ez,
                            double s_size) {
  if (window_size < 2 || window_size >= static_cast<uint32_t>(data.size()) ||
      window_size >= static_cast<uint32_t>(query.size())) {
    Rcpp::stop("window_size must leave at least two subsequences in each series");
  }
  if (!std::isfinite(ez) || ez < 0.0) Rcpp::stop("ez must be finite and non-negative");
  if (!std::isfinite(s_size) || s_size < 0.0 || s_size > 1.0) Rcpp::stop("s_size must be between 0 and 1");
}

List scrimp_result(const List &profile, double ez) {
  return List::create(_["matrix_profile"] = profile["matrix_profile"],
                      _["profile_index"] = profile["profile_index"], _["partial"] = profile["partial"],
                      _["ez"] = ez);
}

List scrimp_ab_result(const List &profile) {
  return List::create(_["matrix_profile"] = profile["matrix_profile"],
                      _["profile_index"] = profile["profile_index"], _["mpb"] = profile["mpb"],
                      _["pib"] = profile["pib"], _["partial"] = profile["partial"]);
}

void apply_pre_scrimp(List &profile, const NumericVector &data_ref, uint32_t window_size, double ez,
                      double pre_scrimp, bool progress) {
  if (pre_scrimp <= 0.0) return;

  List const stats = muinvn_na(data_ref, window_size);
  NumericVector const data = stats["data"];
  NumericVector const mean = stats["avg"];
  NumericVector const inverse_norm = stats["sig"];
  LogicalVector const valid = stats["valid_window"];
  NumericVector matrix_profile = profile["matrix_profile"];
  IntegerVector profile_index = profile["profile_index"];
  uint32_t const profile_size = matrix_profile.size();
  uint32_t const exclusion_zone = std::round(window_size * ez + DBL_EPSILON);
  uint32_t const step = std::max<uint32_t>(1, std::floor(window_size * pre_scrimp + DBL_EPSILON));
  uint32_t const samples = (profile_size + step - 1) / step;
  Progress p(samples, progress);

  for (uint32_t i = 0; i < profile_size; i += step) {
    RcppThread::checkUserInterrupt();
    p.increment();
    if (!valid[i]) continue;
    for (uint32_t j = 0; j < profile_size; j++) {
      if (!valid[j] || std::abs(static_cast<int64_t>(i) - static_cast<int64_t>(j)) <=
                           static_cast<int64_t>(exclusion_zone))
        continue;
      double covariance = 0.0;
      for (uint32_t k = 0; k < window_size; k++) {
        covariance += (data[i + k] - mean[i]) * (data[j + k] - mean[j]);
      }
      double correlation = covariance * inverse_norm[i] * inverse_norm[j];
      correlation = std::max(-1.0, std::min(1.0, correlation));
      double const distance = std::sqrt(std::max(0.0, 2.0 * window_size * (1.0 - correlation)));
      if (!std::isfinite(matrix_profile[j]) || distance < matrix_profile[j]) {
        matrix_profile[j] = distance;
        profile_index[j] = i + 1;
      }
      if (!std::isfinite(matrix_profile[i]) || distance < matrix_profile[i]) {
        matrix_profile[i] = distance;
        profile_index[i] = j + 1;
      }
    }
  }
}

} // namespace

// [[Rcpp::export]]
List scrimp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                 double s_size, double pre_scrimp, bool progress) {
  validate_scrimp_inputs(data_ref, query_ref, window_size, ez, s_size);
  if (!std::isfinite(pre_scrimp) || pre_scrimp < 0.0 || pre_scrimp > 1.0) {
    Rcpp::stop("pre_scrimp must be between 0 and 1");
  }

  // SCRIMP and MPX traverse the same diagonal correlation recurrence. Sharing
  // the NA-aware implementation prevents the former SCRIMP recurrence from
  // dividing by zero for constant windows and from indexing invalid windows.
  List profile = mpx_na_rcpp(data_ref, window_size, ez, s_size, true, true, progress);
  apply_pre_scrimp(profile, data_ref, window_size, ez, pre_scrimp, progress);
  return scrimp_result(profile, ez);
}

// [[Rcpp::export]]
List scrimp_rcpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                          double s_size, bool progress) {
  validate_scrimp_inputs(data_ref, query_ref, window_size, ez, s_size);
  List const profile = mpx_na_rcpp_parallel(data_ref, window_size, ez, s_size, true, true, progress);
  return scrimp_result(profile, ez);
}

// [[Rcpp::export]]
List scrimpab_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double s_size,
                   bool progress) {
  validate_scrimp_inputs(data_ref, query_ref, window_size, 0.0, s_size);
  List const profile = mpxab_na_rcpp(data_ref, query_ref, window_size, s_size, true, true, progress);
  return scrimp_ab_result(profile);
}
