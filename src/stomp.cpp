// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#include "mathtools.h" // math first to fix OSX error
#include "stomp.h"
#include "mass.h"
#include <numeric>
#include <mutex>
#include <cfloat> // DBL_EPSILON when STRICT_R_HEADERS
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

// [[Rcpp::export]]
List stomp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                bool progress, bool left_right_profile) {
  bool partial = false;
  double const exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint32_t const data_size = data_ref.length();
  uint32_t const query_size = query_ref.length();
  uint32_t const matrix_profile_size = data_size - window_size + 1;
  uint32_t const num_queries = query_size - window_size + 1;

  // check skip position
  LogicalVector skip_location(matrix_profile_size, FALSE);

  for (uint64_t i = 0; i < matrix_profile_size; i++) {
    NumericVector const range = data_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_location[i] = TRUE;
    }
  }

  NumericVector data = data_ref;
  NumericVector query = query_ref;

  data[is_na(data)] = 0;
  data[is_infinite(data)] = 0;
  query[is_na(query)] = 0;
  query[is_infinite(query)] = 0;

  NumericVector matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector profile_index(matrix_profile_size, -1);
  NumericVector right_matrix_profile;
  IntegerVector right_profile_index;
  NumericVector left_matrix_profile;
  IntegerVector left_profile_index;

  if (left_right_profile) {
    right_matrix_profile = NumericVector(matrix_profile_size, R_PosInf);
    right_profile_index = IntegerVector(matrix_profile_size, -1);
    left_matrix_profile = NumericVector(matrix_profile_size, R_PosInf);
    left_profile_index = IntegerVector(matrix_profile_size, -1);
  }

  uint32_t const grain = find_best_k_rcpp(data, query, window_size);

  List pre = mass_pre_rcpp(data, query, window_size);
  // NOLINTNEXTLINE(readability-suspicious-call-argument): it is really swapped
  List rpre = mass_pre_rcpp(query, data, window_size); //-V764
  List nn =
      mass3_rcpp(query[Range(0, window_size - 1)], data, pre["data_size"], rpre["window_size"], pre["data_mean"],
                 pre["data_sd"], as<NumericVector>(pre["query_mean"])[0], as<NumericVector>(pre["query_sd"])[0], grain);

  ///// This is needed for JOIN similarity
  List rnn = mass3_rcpp(data[Range(0, window_size - 1)], query, query_size, rpre["window_size"], rpre["data_mean"],
                        rpre["data_sd"], as<NumericVector>(rpre["query_mean"])[0],
                        as<NumericVector>(rpre["query_sd"])[0], grain);

  NumericVector first_product = rnn["last_product"];

  ///// This is needed for JOIN similarity
  IntegerVector const lp_range = Range(1, (data_size - window_size));
  IntegerVector const lp2_range = Range(0, (data_size - window_size - 1));
  IntegerVector const dt_range = Range(window_size, data_size - 1);
  NumericVector distance_profile;
  NumericVector last_product;
  double drop_value = query[0];

  IntegerVector const order = Range(0, num_queries - 1);

  Progress p(100, progress);

  uint32_t const num_progress =
      ceil((double)order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

  try {
    for (int32_t const &i : order) {

      if ((i % num_progress) == 0) {
        RcppThread::checkUserInterrupt();
        p.increment();
      }
      // compute the distance profile
      NumericVector query_window = query[Range(i, (i + window_size - 1))];
      if (i == 0) {
        distance_profile = nn["distance_profile"];
        last_product = nn["last_product"];
      } else {
        last_product[lp_range] = (NumericVector)(as<NumericVector>(last_product[lp2_range]) -
                                                 as<NumericVector>(data[lp2_range]) * drop_value +
                                                 as<NumericVector>(data[dt_range]) * query_window[window_size - 1]);
        last_product[0] = first_product[i];
        distance_profile =
            2 * (window_size - (last_product - window_size * as<NumericVector>(pre["data_mean"]) *
                                                   as<NumericVector>(pre["query_mean"])[i]) /
                                   (as<NumericVector>(pre["data_sd"]) * as<NumericVector>(pre["query_sd"])[i]));
        distance_profile[distance_profile < 0] = 0;
      }

      // distance_profile = sqrt(distance_profile);
      drop_value = query_window[0];

      // apply exclusion zone
      if (exclusion_zone > 0) {
        uint32_t const exc_st = MAX(0, i - exclusion_zone);
        uint32_t const exc_ed = MIN(matrix_profile_size - 1, i + exclusion_zone);
        IntegerVector const dp_range = Range(exc_st, exc_ed);
        distance_profile[dp_range] = R_PosInf;
      }

      distance_profile[as<NumericVector>(pre["data_sd"]) < DBL_EPSILON] = R_PosInf;
      if (skip_location[i] == TRUE || as<NumericVector>(pre["query_sd"])[i] < DBL_EPSILON) {
        distance_profile.fill(R_PosInf);
      }

      distance_profile[skip_location] = R_PosInf;

      if (left_right_profile) {
        // left matrix_profile
        LogicalVector ldx(matrix_profile_size, FALSE);
        ldx[Range(i, matrix_profile_size - 1)] = (distance_profile[Range(i, matrix_profile_size - 1)] <
                                                  left_matrix_profile[Range(i, matrix_profile_size - 1)]);
        left_matrix_profile[ldx] = distance_profile[ldx];
        left_profile_index[which_cpp(ldx)] = i;

        // right matrix_profile
        LogicalVector rdx(matrix_profile_size, FALSE);
        rdx[Range(0, i)] = (distance_profile[Range(0, i)] < right_matrix_profile[Range(0, i)]);
        right_matrix_profile[rdx] = distance_profile[rdx];
        right_profile_index[which_cpp(rdx)] = i;
      }

      LogicalVector const idx = (distance_profile < matrix_profile);
      matrix_profile[idx] = distance_profile[idx];
      profile_index[which_cpp(idx)] = i + 1;
    }

    if (left_right_profile) {
      left_matrix_profile = sqrt(left_matrix_profile);
      right_matrix_profile = sqrt(right_matrix_profile);
    }

    matrix_profile = sqrt(matrix_profile);
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  if (left_right_profile) {
    return (List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("profile_index") = profile_index,
                         Rcpp::Named("right_matrix_profile") = right_matrix_profile,
                         Rcpp::Named("right_profile_index") = right_profile_index,
                         Rcpp::Named("left_matrix_profile") = left_matrix_profile,
                         Rcpp::Named("left_profile_index") = left_profile_index, Rcpp::Named("partial") = partial,
                         Rcpp::Named("ez") = ez));
  } else {
    return (List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("profile_index") = profile_index,
                         Rcpp::Named("partial") = partial, Rcpp::Named("ez") = ez));
  }
}

struct StompWorker : public Worker {
private:
  // input
  const RVector<double> data_ref;
  const RVector<double> window_ref;
  const uint64_t w_size;
  const uint64_t d_size;
  const RVector<double> d_mean;
  const RVector<double> d_std;
  const RVector<double> q_mean;
  const RVector<double> q_std;
  const RVector<int> skip_location;
  const RVector<double> first_product;
  const uint64_t ez;

  RVector<double> mp;
  RVector<int> pi;
  const bool left_right_profile;
  RVector<double> rmp;
  RVector<int> rpi;
  RVector<double> lmp;
  RVector<int> lpi;

#if RCPP_PARALLEL_USE_TBB
  tbb::spin_mutex m;
#else
  tthread::mutex m;
#endif

public:
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  StompWorker(const NumericVector &data_ref, const NumericVector &window_ref, const uint64_t w_size,
              const uint64_t d_size, const NumericVector &d_mean, const NumericVector &d_std,
              const NumericVector &q_mean, const NumericVector &q_std, const IntegerVector &skip_location,
              const NumericVector &first_product, const uint64_t ez, const NumericVector &mp, const IntegerVector &pi,
              const bool left_right_profile, const NumericVector &rmp, const IntegerVector &rpi,
              const NumericVector &lmp, const IntegerVector &lpi)
      : data_ref(data_ref), window_ref(window_ref), w_size(w_size), d_size(d_size), d_mean(d_mean), d_std(d_std),
        q_mean(q_mean), q_std(q_std), skip_location(skip_location), first_product(first_product), ez(ez), mp(mp), pi(pi),
        left_right_profile(left_right_profile), rmp(rmp), rpi(rpi), lmp(lmp), lpi(lpi) {}

  ~StompWorker() override = default;

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) override {
    // begin and end are the query window

    uint64_t const grain = set_k_cpp(w_size * 2, d_size, w_size);
    Mass3Result nn = mass3_cpp_result(window_ref.begin() + begin, data_ref.begin(), d_size, w_size, d_mean.begin(),
                                      d_std.begin(), q_mean[begin], q_std[begin], grain);
    std::vector<double> distance_profile = std::move(nn.distance_profile);
    std::vector<double> last_product = std::move(nn.last_product);

    std::vector<double> matrix_profile(d_mean.size(), R_PosInf);
    std::vector<int> profile_index(d_mean.size(), -1);
    std::vector<double> right_matrix_profile(left_right_profile ? d_mean.size() : 0, R_PosInf);
    std::vector<int> right_profile_index(left_right_profile ? d_mean.size() : 0, -1);
    std::vector<double> left_matrix_profile(left_right_profile ? d_mean.size() : 0, R_PosInf);
    std::vector<int> left_profile_index(left_right_profile ? d_mean.size() : 0, -1);
    double drop_value = 0;
    uint32_t exc_st = 0;
    uint32_t exc_ed = 0;

    for (uint64_t i = begin; i < end; i++) {

        // compute the distance profile
        if (i > begin) {

          for (uint64_t j = d_mean.size() - 1; j > 0; j--) {
            last_product[j] = last_product[j - 1] - data_ref[j - 1] * drop_value +
                              data_ref[w_size + j - 1] * window_ref[i + w_size - 1];
          }
          last_product[0] = first_product[i];
        }

        if (ez > 0) {
          exc_st = MAX(0, i > ez ? (i - ez) : 0);
          exc_ed = MIN(d_std.size() - 1, i + ez);
        }

        for (uint64_t j = 0; j < d_mean.size(); j++) {
          if (skip_location[j] == 1 || d_std[j] < DBL_EPSILON || q_std[i] < DBL_EPSILON) {
            distance_profile[j] = R_PosInf;
          } else if (ez == 0 || j < exc_st || j > exc_ed) {
            double const dp = 2 * (w_size - (last_product[j] - w_size * d_mean[j] * q_mean[i]) / (d_std[j] * q_std[i]));
            distance_profile[j] = (dp > 0) ? dp : 0;
          } else {
            distance_profile[j] = R_PosInf;
          }
        }

        drop_value = window_ref[i];

      for (uint64_t j = 0; j < d_mean.size(); j++) {
        if (distance_profile[j] < matrix_profile[j]) {
          matrix_profile[j] = distance_profile[j];
          profile_index[j] = i + 1;
        }
        if (left_right_profile && j >= i && distance_profile[j] < left_matrix_profile[j]) {
          left_matrix_profile[j] = distance_profile[j];
          left_profile_index[j] = i;
        }
        if (left_right_profile && j <= i && distance_profile[j] < right_matrix_profile[j]) {
          right_matrix_profile[j] = distance_profile[j];
          right_profile_index[j] = i;
        }
      }
    }

    std::lock_guard<decltype(m)> lock(m);
    for (uint64_t j = 0; j < d_mean.size(); j++) {
      if (matrix_profile[j] < mp[j]) {
        mp[j] = matrix_profile[j];
        pi[j] = profile_index[j];
      }
      if (left_right_profile && right_matrix_profile[j] < rmp[j]) {
        rmp[j] = right_matrix_profile[j];
        rpi[j] = right_profile_index[j];
      }
      if (left_right_profile && left_matrix_profile[j] < lmp[j]) {
        lmp[j] = left_matrix_profile[j];
        lpi[j] = left_profile_index[j];
      }
    }
  }
};

// [[Rcpp::export]]
List stomp_rcpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                         bool progress, bool left_right_profile) {
  uint64_t const exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint64_t const data_size = data_ref.length();
  uint64_t const query_size = query_ref.length();
  uint64_t const matrix_profile_size = data_size - window_size + 1;
  uint64_t const num_queries = query_size - window_size + 1;
  bool partial = false;

  // check skip position
  IntegerVector skip_location(matrix_profile_size, 0);

  for (uint64_t i = 0; i < matrix_profile_size; i++) {
    NumericVector const range = data_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_location[i] = 1;
    }
  }

  NumericVector data = data_ref;
  NumericVector query = query_ref;

  data[is_na(data)] = 0;
  data[is_infinite(data)] = 0;
  query[is_na(query)] = 0;
  query[is_infinite(query)] = 0;

  NumericVector const matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector const profile_index(matrix_profile_size, -1);
  NumericVector const right_matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector const right_profile_index(matrix_profile_size, -1);
  NumericVector const left_matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector const left_profile_index(matrix_profile_size, -1);

  uint64_t grain = set_k_rcpp(256, data_size, window_size);

  ///// This is needed for JOIN similarity
  // NOLINTNEXTLINE(readability-suspicious-call-argument): it is really swapped
  List rpre = mass_pre_rcpp(query, data, window_size); //-V764
  List rnn = mass3_rcpp(data[Range(0, window_size - 1)], query, query_size, rpre["window_size"], rpre["data_mean"],
                        rpre["data_sd"], as<NumericVector>(rpre["query_mean"])[0],
                        as<NumericVector>(rpre["query_sd"])[0], grain);

  NumericVector const first_product = rnn["last_product"];

  List pre = mass_pre_rcpp(data, query, window_size);

  Progress p(100, progress);

  StompWorker stomp_worker(data, query, pre["window_size"], data_size, pre["data_mean"], pre["data_sd"],
                           pre["query_mean"], pre["query_sd"], skip_location, first_product, exclusion_zone,
                           matrix_profile, profile_index, left_right_profile, right_matrix_profile, right_profile_index,
                           left_matrix_profile, left_profile_index);

  grain = set_k_rcpp(1024, num_queries, window_size);

  // call parallelFor to do the work
  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, num_queries, stomp_worker, 2 * grain);
#else
    RcppParallel2::ttParallelFor(0, num_queries, stomp_worker, 2 * grain);
#endif

    p.increment(100);

  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results "
             "were returned.";
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  if (left_right_profile) {
    return List::create(Rcpp::Named("matrix_profile") = sqrt(matrix_profile),
                        Rcpp::Named("profile_index") = profile_index,
                        Rcpp::Named("right_matrix_profile") = sqrt(right_matrix_profile),
                        Rcpp::Named("right_profile_index") = right_profile_index,
                        Rcpp::Named("left_matrix_profile") = sqrt(left_matrix_profile),
                        Rcpp::Named("left_profile_index") = left_profile_index, Rcpp::Named("partial") = partial,
                        Rcpp::Named("ez") = ez);
  }

  return List::create(Rcpp::Named("matrix_profile") = sqrt(matrix_profile),
                      Rcpp::Named("profile_index") = profile_index, Rcpp::Named("partial") = partial,
                      Rcpp::Named("ez") = ez);
}
