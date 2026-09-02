// This is a personal academic project. Dear PVS-Studio, please check it.

#include "mathtools.h" // math first to fix OSX error
#include "stamp.h"
#include "mass.h"
#include <cfloat>
#include <cmath>
#include <vector>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

using namespace RcppParallel;

#if RCPP_PARALLEL_USE_TBB
#include "tbb/mutex.h"
#else
#include "rcpp_parallel_fix.h"
#include "tthread/tinythread.h"
#endif

namespace {

void validate_stamp_inputs(const NumericVector &data, const NumericVector &query, uint32_t window_size, double ez,
                           double s_size) {
  if (window_size < 2 || window_size > static_cast<uint32_t>(data.size()) ||
      window_size > static_cast<uint32_t>(query.size())) {
    Rcpp::stop("window_size is incompatible with data or query");
  }
  if (!std::isfinite(ez) || ez < 0.0) Rcpp::stop("ez must be finite and non-negative");
  if (!std::isfinite(s_size) || s_size < 0.0 || s_size > 1.0) Rcpp::stop("s_size must be between 0 and 1");
}

LogicalVector valid_windows(const NumericVector &x, uint32_t window_size) {
  uint32_t const profile_size = x.size() - window_size + 1;
  LogicalVector valid(profile_size, true);
  uint32_t invalid = 0;
  for (uint32_t i = 0; i < static_cast<uint32_t>(x.size()); i++) {
    if (!std::isfinite(x[i])) invalid++;
    if (i >= window_size && !std::isfinite(x[i - window_size])) invalid--;
    if (i >= window_size - 1) valid[i - window_size + 1] = invalid == 0;
  }
  return valid;
}

uint64_t requested_work(uint64_t total, double s_size, bool &partial) {
  if (s_size > 0.0 && s_size < 1.0) {
    uint64_t const work = static_cast<uint64_t>(std::round(total * s_size + DBL_EPSILON));
    partial = work < total;
    return work;
  }
  return total;
}

struct StampWorker : public Worker {
private:
  const RVector<double> data;
  const RVector<double> query;
  const RVector<int> order;
  const uint32_t window_size;
  const uint64_t data_size;
  const RVector<double> data_mean;
  const RVector<double> data_sd;
  const RVector<double> query_mean;
  const RVector<double> query_sd;
  const RVector<int> valid_data;
  const RVector<int> valid_query;
  const uint64_t exclusion_zone;
  const uint32_t grain;
  RVector<double> matrix_profile;
  RVector<int> profile_index;

#if RCPP_PARALLEL_USE_TBB
  tbb::mutex mutex;
#else
  tthread::mutex mutex;
#endif

public:
  StampWorker(const NumericVector &data, const NumericVector &query, const IntegerVector &order, uint32_t window_size,
              const NumericVector &data_mean, const NumericVector &data_sd, const NumericVector &query_mean,
              const NumericVector &query_sd, const LogicalVector &valid_data, const LogicalVector &valid_query,
              uint64_t exclusion_zone, uint32_t grain, const NumericVector &matrix_profile,
              const IntegerVector &profile_index)
      : data(data), query(query), order(order), window_size(window_size), data_size(data.size()), data_mean(data_mean),
        data_sd(data_sd), query_mean(query_mean), query_sd(query_sd), valid_data(valid_data), valid_query(valid_query),
        exclusion_zone(exclusion_zone), grain(grain), matrix_profile(matrix_profile), profile_index(profile_index) {}

  void operator()(std::size_t begin, std::size_t end) override {
    std::vector<double> local_profile(matrix_profile.size(), R_PosInf);
    std::vector<int> local_index(profile_index.size(), -1);

    for (std::size_t position = begin; position < end; position++) {
      uint32_t const query_index = order[position];
      if (!valid_query[query_index] || !std::isfinite(query_sd[query_index]) || query_sd[query_index] < DBL_EPSILON)
        continue;

      Mass3Result result = mass3_cpp_result(query.begin() + query_index, data.begin(), data_size, window_size,
                                            data_mean.begin(), data_sd.begin(), query_mean[query_index],
                                            query_sd[query_index], grain);
      for (uint32_t data_index = 0; data_index < matrix_profile.size(); data_index++) {
        if (!valid_data[data_index] || !std::isfinite(data_sd[data_index]) || data_sd[data_index] < DBL_EPSILON)
          continue;
        if (exclusion_zone > 0 &&
            std::abs(static_cast<int64_t>(data_index) - static_cast<int64_t>(query_index)) <=
                static_cast<int64_t>(exclusion_zone))
          continue;
        double distance = result.distance_profile[data_index];
        if (!std::isfinite(distance)) continue;
        distance = std::max(0.0, distance);
        if (distance < local_profile[data_index]) {
          local_profile[data_index] = distance;
          local_index[data_index] = query_index + 1;
        }
      }
    }

    std::lock_guard<decltype(mutex)> lock(mutex);
    for (uint32_t i = 0; i < matrix_profile.size(); i++) {
      if (local_profile[i] < matrix_profile[i]) {
        matrix_profile[i] = local_profile[i];
        profile_index[i] = local_index[i];
      }
    }
  }
};

} // namespace

// [[Rcpp::export]]
List stamp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                double s_size, bool progress) {
  validate_stamp_inputs(data_ref, query_ref, window_size, ez, s_size);
  uint64_t const profile_size = data_ref.size() - window_size + 1;
  uint64_t const query_profile_size = query_ref.size() - window_size + 1;
  uint64_t const exclusion_zone = std::round(window_size * ez + DBL_EPSILON);
  bool partial = false;

  LogicalVector const valid_data = valid_windows(data_ref, window_size);
  LogicalVector const valid_query = valid_windows(query_ref, window_size);
  NumericVector data = clone(data_ref);
  NumericVector query = clone(query_ref);
  for (R_xlen_t i = 0; i < data.size(); i++) if (!std::isfinite(data[i])) data[i] = 0.0;
  for (R_xlen_t i = 0; i < query.size(); i++) if (!std::isfinite(query[i])) query[i] = 0.0;

  List const pre = mass_pre_rcpp(data, query, window_size);
  NumericVector data_mean = pre["data_mean"];
  NumericVector data_sd = pre["data_sd"];
  NumericVector query_mean = pre["query_mean"];
  NumericVector query_sd = pre["query_sd"];
  IntegerVector order = Rcpp::Range(0, query_profile_size - 1);
  order = sample(order, order.size());
  uint64_t const work_size = requested_work(query_profile_size, s_size, partial);
  uint32_t const grain = set_k_rcpp(window_size, data.size(), window_size);
  NumericVector matrix_profile(profile_size, R_PosInf);
  IntegerVector profile_index(profile_size, -1);
  Progress p(work_size, progress);

  try {
    for (uint64_t position = 0; position < work_size; position++) {
      RcppThread::checkUserInterrupt();
      p.increment();
      uint32_t const query_index = order[position];
      if (!valid_query[query_index] || !std::isfinite(query_sd[query_index]) || query_sd[query_index] < DBL_EPSILON)
        continue;
      Mass3Result result = mass3_cpp_result(query.begin() + query_index, data.begin(), data.size(), window_size,
                                            data_mean.begin(), data_sd.begin(), query_mean[query_index],
                                            query_sd[query_index], grain);
      for (uint32_t data_index = 0; data_index < profile_size; data_index++) {
        if (!valid_data[data_index] || !std::isfinite(data_sd[data_index]) || data_sd[data_index] < DBL_EPSILON)
          continue;
        if (exclusion_zone > 0 &&
            std::abs(static_cast<int64_t>(data_index) - static_cast<int64_t>(query_index)) <=
                static_cast<int64_t>(exclusion_zone))
          continue;
        double distance = result.distance_profile[data_index];
        if (std::isfinite(distance) && distance < matrix_profile[data_index]) {
          matrix_profile[data_index] = std::max(0.0, distance);
          profile_index[data_index] = query_index + 1;
        }
      }
    }
  } catch (RcppThread::UserInterruptException &) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
  }

  return List::create(_["matrix_profile"] = sqrt(matrix_profile), _["profile_index"] = profile_index,
                      _["partial"] = partial, _["ez"] = ez);
}

// [[Rcpp::export]]
List stamp_rcpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                         double s_size, bool progress) {
  validate_stamp_inputs(data_ref, query_ref, window_size, ez, s_size);
  uint64_t const profile_size = data_ref.size() - window_size + 1;
  uint64_t const query_profile_size = query_ref.size() - window_size + 1;
  uint64_t const exclusion_zone = std::round(window_size * ez + DBL_EPSILON);
  bool partial = false;

  LogicalVector const valid_data = valid_windows(data_ref, window_size);
  LogicalVector const valid_query = valid_windows(query_ref, window_size);
  NumericVector data = clone(data_ref);
  NumericVector query = clone(query_ref);
  for (R_xlen_t i = 0; i < data.size(); i++) if (!std::isfinite(data[i])) data[i] = 0.0;
  for (R_xlen_t i = 0; i < query.size(); i++) if (!std::isfinite(query[i])) query[i] = 0.0;

  List const pre = mass_pre_rcpp(data, query, window_size);
  IntegerVector order = Rcpp::Range(0, query_profile_size - 1);
  order = sample(order, order.size());
  uint64_t const work_size = requested_work(query_profile_size, s_size, partial);
  if (work_size == 0) order = IntegerVector(0);
  else if (work_size < query_profile_size) order = order[Rcpp::Range(0, work_size - 1)];
  NumericVector matrix_profile(profile_size, R_PosInf);
  IntegerVector profile_index(profile_size, -1);
  uint32_t const grain = set_k_rcpp(window_size, data.size(), window_size);

  StampWorker worker(data, query, order, window_size, pre["data_mean"], pre["data_sd"], pre["query_mean"],
                     pre["query_sd"], valid_data, valid_query, exclusion_zone, grain, matrix_profile, profile_index);
  if (work_size > 0) {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, work_size, worker, 1);
#else
    RcppParallel2::ttParallelFor(0, work_size, worker, 1);
#endif
  }
  Progress p(work_size, progress);
  for (uint64_t i = 0; i < work_size; i++) p.increment();

  return List::create(_["matrix_profile"] = sqrt(matrix_profile), _["profile_index"] = profile_index,
                      _["partial"] = partial, _["ez"] = ez);
}
