#include "stomp.h"

#include "fft.h"
#include "mass.h"
#include "math.h"
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
#include "tthread/fast_mutex.h"
#endif

// [[Rcpp::export]]
List stomp_rcpp(const NumericVector data_ref, const NumericVector query_ref,
                uint32_t window_size, double ez = 0.5, bool progress = false) {
  bool partial = false;
  double exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint32_t data_size = data_ref.length();
  uint32_t query_size = query_ref.length();
  uint32_t matrix_profile_size = data_size - window_size + 1;
  uint32_t num_queries = query_size - window_size + 1;

  // check skip position
  LogicalVector skip_location(matrix_profile_size);

  for (uint32_t i = 0; i < matrix_profile_size; i++) {
    NumericVector range = data_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_location[i] = true;
    }
  }

  NumericVector data = data_ref;
  NumericVector query = query_ref;

  data[is_na(data)] = 0;
  data[is_infinite(data)] = 0;
  query[is_na(query)] = 0;
  query[is_infinite(query)] = 0;

  NumericVector matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector profile_index(matrix_profile_size, R_NegInf);

  uint32_t k = find_best_k_rcpp(data, query, window_size);

  List pre = mass_pre_rcpp(data, query, window_size);
  List rpre = mass_pre_rcpp(query, data, window_size);
  List nn = mass3_rcpp(query[Range(0, window_size - 1)], data, pre["data_size"],
                       rpre["window_size"], pre["data_mean"], pre["data_sd"],
                       as<NumericVector>(pre["query_mean"])[0],
                       as<NumericVector>(pre["query_sd"])[0], k);

  ///// This is needed for JOIN similarity
  List rnn = mass3_rcpp(data[Range(0, window_size - 1)], query, query_size,
                        rpre["window_size"], rpre["data_mean"], rpre["data_sd"],
                        as<NumericVector>(rpre["query_mean"])[0],
                        as<NumericVector>(rpre["query_sd"])[0], k);

  NumericVector first_product = rnn["last_product"];

  ///// This is needed for JOIN similarity
  IntegerVector lp_range = Range(1, (data_size - window_size));
  IntegerVector lp2_range = Range(0, (data_size - window_size - 1));
  IntegerVector dt_range = Range(window_size, data_size - 1);
  NumericVector distance_profile;
  NumericVector last_product;
  double drop_value = query[0];

  IntegerVector order = Range(0, num_queries - 1);

  Progress p(num_queries, progress);

  try {
    for (int32_t i : order) {
      p.increment();

      if (i % 100 == 0) {
        RcppThread::checkUserInterrupt();
      }

      // compute the distance profile
      NumericVector query_window = query[Range(i, (i + window_size - 1))];
      if (i == 0) {
        distance_profile = nn["distance_profile"];
        last_product = nn["last_product"];
      } else {
        last_product[lp_range] = (NumericVector)(
            as<NumericVector>(last_product[lp2_range]) -
            as<NumericVector>(data[lp2_range]) * drop_value +
            as<NumericVector>(data[dt_range]) * query_window[window_size - 1]);
        last_product[0] = first_product[i];
        distance_profile =
            2 *
            (window_size -
             (last_product - window_size * as<NumericVector>(pre["data_mean"]) *
                                 as<NumericVector>(pre["query_mean"])[i]) /
                 (as<NumericVector>(pre["data_sd"]) *
                  as<NumericVector>(pre["query_sd"])[i]));
      }
      distance_profile[distance_profile < 0] = 0;
      // distance_profile = sqrt(distance_profile);
      drop_value = query_window[0];

      // apply exclusion zone
      if (exclusion_zone > 0) {
        uint32_t exc_st = MAX(0, i - exclusion_zone);
        uint32_t exc_ed = MIN(matrix_profile_size - 1, i + exclusion_zone);
        IntegerVector dp_range = Range(exc_st, exc_ed);
        distance_profile[dp_range] = R_PosInf;
      }

      distance_profile[as<NumericVector>(pre["data_sd"]) < DBL_EPSILON] =
          R_PosInf;
      if (skip_location[i] ||
          as<NumericVector>(pre["query_sd"])[i] < DBL_EPSILON) {
        distance_profile.fill(R_PosInf);
      }
      distance_profile[skip_location] = R_PosInf;

      LogicalVector idx = (distance_profile < matrix_profile);
      matrix_profile[idx] = distance_profile[idx];
      profile_index[which(idx)] = i + 1;
    }

    matrix_profile = sqrt(matrix_profile);
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results "
             "were returned.";
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("matrix_profile") = matrix_profile,
                       Rcpp::Named("profile_index") = profile_index,
                       Rcpp::Named("partial") = partial,
                       Rcpp::Named("ez") = ez));
}

struct StompWorker : public Worker {
  // input
  const RVector<double> data_ref;
  const RVector<double> window_ref;
  const uint64_t w_size;
  const uint64_t d_size;
  const RVector<double> d_mean;
  const RVector<double> d_std;
  const RVector<double> q_mean;
  const RVector<double> q_std;
  const IntegerVector skip_location;
  const RVector<double> first_product;
  const IntegerVector lp_range;
  const IntegerVector lp2_range;
  const IntegerVector dt_range;
  const uint64_t ez;

  Progress *p;

  RVector<double> mp;
  RVector<int> pi;

  tbb::mutex m;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  StompWorker(const NumericVector data_ref, const NumericVector window_ref,
              const uint64_t w_size, const uint64_t d_size,
              const NumericVector d_mean, const NumericVector d_std,
              const NumericVector q_mean, const NumericVector q_std,
              const LogicalVector skip_location,
              const NumericVector first_product, const IntegerVector lp_range,
              const IntegerVector lp2_range, const IntegerVector dt_range,
              const uint64_t ez, Progress *p, NumericVector mp,
              IntegerVector pi)
      : data_ref(data_ref), window_ref(window_ref), w_size(w_size),
        d_size(d_size), d_mean(d_mean), d_std(d_std), q_mean(q_mean),
        q_std(q_std), skip_location(skip_location),
        first_product(first_product), lp_range(lp_range), lp2_range(lp2_range),
        dt_range(dt_range), ez(ez), p(p), mp(mp), pi(pi) {}

  ~StompWorker() {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    // begin and end are the query window

    uint64_t chunk = (end - begin);

    if (<= w_size) {
      ::Rf_error("Chunk size is too small.");
    }

    Rcout << "Chunk size: " << (end - begin) << std::endl;

    uint64_t k = set_k_rcpp(w_size * 2, (end - begin), w_size);

    NumericVector query(window_ref.begin() + begin,
                        window_ref.begin() + begin + w_size - 1);
    NumericVector data(data_ref.begin(), data_ref.end());
    NumericVector data_mean(d_mean.begin(), d_mean.end());
    NumericVector data_std(d_std.begin(), d_std.end());
    float qm = q_mean[begin];
    float qs = q_std[begin];

    List nn = mass3_rcpp_parallel(query, data, d_size, w_size, data_mean,
                                  data_std, qm, qs, k);

    // NumericVector distance_profile;
    // NumericVector last_product;
    // double drop_value = 0;

    // // index of sliding window
    // try {
    // for (uint64_t i = begin; i < end; i++) {
    //   // compute the distance profile
    //   NumericVector query(window_ref.begin() + i,
    //                       window_ref.begin() + (i + w_size - 1));
    //
    //   if (i == begin) {
    //     distance_profile = nn["distance_profile"];
    //     last_product = nn["last_product"];
    //   } else {
    //     last_product[lp_range] = (NumericVector)(
    //         as<NumericVector>(last_product[lp2_range]) -
    //         as<NumericVector>(data[lp2_range]) * drop_value +
    //         as<NumericVector>(data[dt_range]) * query[w_size - 1]);
    //     last_product[0] = first_product[i];
    //     distance_profile =
    //         2 *
    //         (w_size - (last_product -
    //                    w_size * as<NumericVector>(data_mean) * q_mean[i]) /
    //                       (as<NumericVector>(data_std) * q_std[i]));
    //   }
    //
    //   distance_profile[distance_profile < 0] = 0;
    //   //   // distance_profile = sqrt(distance_profile);
    //   drop_value = query[0];
    //
    //   // apply exclusion zone
    //   if (ez > 0) {
    //     uint32_t exc_st = MAX(0, ((int64_t)i - (int64_t)ez));
    //     uint32_t exc_ed = MIN((int64_t)(data_std.size() - 1), (int64_t)(i +
    //     ez)); IntegerVector dp_range = Range(exc_st, exc_ed);
    //     distance_profile[dp_range] = R_PosInf;
    //   }
    //
    //   distance_profile[as<NumericVector>(data_std) < DBL_EPSILON] = R_PosInf;
    //   if (skip_location[i] || q_std[i] < DBL_EPSILON) {
    //     distance_profile.fill(R_PosInf);
    //   }
    //   distance_profile[skip_location] = R_PosInf;
    //
    //   LogicalVector idx = (distance_profile < data_mean);
    //   data_mean[idx] = distance_profile[idx];
    //   // pi[which(idx)] = i + 1;
    //   // }
    // }
    // } catch (RcppThread::UserInterruptException &e) {
    //   Rcout << "Computation interrupted by the user." << std::endl;
    //   Rcout << "Please wait for other threads to stop." << std::endl;
    //   throw;
    // } catch (...) {
    //   ::Rf_error("c++ exception (1unknown reason)");
    // }
  }
};

// [[Rcpp::export]]
List stomp_rcpp_parallel(const NumericVector data_ref,
                         const NumericVector query_ref, uint32_t window_size,
                         double ez = 0.5, bool progress = false) {
  uint64_t exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();
  uint64_t matrix_profile_size = data_size - window_size + 1;
  uint64_t num_queries = query_size - window_size + 1;
  bool partial = false;

  // TODO: check skip position
  LogicalVector skip_location(matrix_profile_size);

  for (uint64_t i = 0; i < matrix_profile_size; i++) {
    NumericVector range = data_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_location[i] = true;
    }
  }

  NumericVector data = data_ref;
  NumericVector query = query_ref;

  data[is_na(data)] = 0;
  data[is_infinite(data)] = 0;
  query[is_na(query)] = 0;
  query[is_infinite(query)] = 0;

  NumericVector matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector profile_index(matrix_profile_size, R_NegInf);

  uint64_t k = set_k_rcpp(256, data_size, window_size);

  ///// This is needed for JOIN similarity
  List rpre = mass_pre_rcpp(query, data, window_size);
  List rnn = mass3_rcpp(data[Range(0, window_size - 1)], query, query_size,
                        rpre["window_size"], rpre["data_mean"], rpre["data_sd"],
                        as<NumericVector>(rpre["query_mean"])[0],
                        as<NumericVector>(rpre["query_sd"])[0], k);

  NumericVector first_product = rnn["last_product"];

  List pre = mass_pre_rcpp(data, query, window_size);

  ///// This is needed for JOIN similarity
  IntegerVector lp_range = Range(1, (data_size - window_size));
  IntegerVector lp2_range = Range(0, (data_size - window_size - 1));
  IntegerVector dt_range = Range(window_size, data_size - 1);

  IntegerVector order = Range(0, num_queries - 1);

  uint64_t parts =
      pow(2, ceil(log2(data_size / k))) * matrix_profile_size / window_size;

  Progress p(parts, progress);

  StompWorker stomp_worker(data, query, pre["window_size"], data_size,
                           pre["data_mean"], pre["data_sd"], pre["query_mean"],
                           pre["query_sd"], skip_location, first_product,
                           lp_range, lp2_range, dt_range, exclusion_zone, &p,
                           matrix_profile, profile_index);

  k = set_k_rcpp(256, num_queries, window_size);

  // call parallelFor to do the work
  try {
    RcppParallel::parallelFor(0, num_queries, stomp_worker, k);
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results "
             "were returned.";
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("matrix_profile") = sqrt(matrix_profile),
                       Rcpp::Named("profile_index") = profile_index,
                       Rcpp::Named("partial") = partial,
                       Rcpp::Named("ez") = ez));
}
