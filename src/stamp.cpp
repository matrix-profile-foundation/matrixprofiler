#include "math.h" // math first to fix OSX error
#include "stamp.h"
#include "fft.h"
#include "mass.h"
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
List stamp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                double s_size, bool progress) {

  bool partial = false;
  int64_t exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();
  uint64_t matrix_profile_size = data_size - window_size + 1;
  uint64_t num_queries = query_size - window_size + 1;

  // check skip position
  LogicalVector skip_location(matrix_profile_size);

  // TODO: data or query?
  for (uint64_t i = 0; i < matrix_profile_size; i++) {
    NumericVector range = data_ref[Range(i, (i + window_size - 1))];
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
  List pre = mass_pre_rcpp(data, query, window_size);

  IntegerVector order = Range(0, num_queries - 1);

  uint32_t k = find_best_k_rcpp(data, query, window_size);

  Progress p(num_queries, progress);

  order = sample(order, order.size());

  uint64_t stop = 0;

  if (s_size < 1.0) {
    stop = round(order.size() * s_size + DBL_EPSILON);
  }

  try {
    uint64_t j = 1;
    for (int32_t i : order) {
      RcppThread::checkUserInterrupt();
      p.increment();

      List nn =
          mass3_rcpp(query[Range(i, i + window_size - 1)], data, pre["data_size"], pre["window_size"], pre["data_mean"],
                     pre["data_sd"], as<NumericVector>(pre["query_mean"])[i], as<NumericVector>(pre["query_sd"])[i], k);

      NumericVector distance_profile = as<NumericVector>(nn["distance_profile"]);

      // apply exclusion zone
      if (exclusion_zone > 0) {
        int64_t exc_st = MAX(0, (i - exclusion_zone));
        int64_t exc_ed = MIN(matrix_profile_size - 1, (uint64_t)(i + exclusion_zone));
        IntegerVector dp_range = Range(exc_st, exc_ed);
        distance_profile[dp_range] = R_PosInf;
      }

      distance_profile[as<NumericVector>(pre["data_sd"]) < DBL_EPSILON] = R_PosInf;
      if (skip_location[i] || as<NumericVector>(pre["query_sd"])[i] < DBL_EPSILON) {
        distance_profile.fill(R_PosInf);
      }
      distance_profile[skip_location] = R_PosInf;

      // normal matrix_profile
      LogicalVector idx = (distance_profile < matrix_profile);
      matrix_profile[idx] = distance_profile[idx];
      profile_index[which_cpp(idx)] = i + 1;

      if (stop > 0 && j >= stop) {
        partial = true;
        break;
      }

      j++;
    }

  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("matrix_profile") = sqrt(matrix_profile),
                       Rcpp::Named("profile_index") = profile_index, Rcpp::Named("partial") = partial,
                       Rcpp::Named("ez") = ez));
}

struct StampWorker : public Worker {
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
  const uint64_t ez;

  Progress *p;

  RVector<double> mp;
  RVector<int> pi;

#if RCPP_PARALLEL_USE_TBB
  tbb::mutex m;
#else
  tthread::mutex m;
#endif

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  StampWorker(const NumericVector data_ref, const NumericVector window_ref, const uint64_t w_size,
              const uint64_t d_size, const NumericVector d_mean, const NumericVector d_std, const NumericVector q_mean,
              const NumericVector q_std, const IntegerVector skip_location, const uint64_t ez, Progress *p,
              NumericVector mp, IntegerVector pi)
      : data_ref(data_ref), window_ref(window_ref), w_size(w_size), d_size(d_size), d_mean(d_mean), d_std(d_std),
        q_mean(q_mean), q_std(q_std), skip_location(skip_location), ez(ez), p(p), mp(mp), pi(pi) {}

  ~StampWorker() {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    // begin and end are the indexes of data.
    uint64_t i, j;
    uint64_t start_ez, end_ez;
    double dp;

    uint64_t chunk = (end - begin);

    if (chunk <= w_size) {
      Rcout << "Chunk size is too small (" << chunk << ") for a window size of " << w_size << std::endl;
      return;
    }

    FFT::fftw *fft = new FFT::fftw();

    // index of sliding window
    try {
      for (uint64_t w = 0; w < q_mean.size(); w++) {

        if (w % w_size == 0) {
          RcppThread::checkUserInterrupt();
          m.lock();
          p->increment();
          m.unlock();
        }

        // exclusion zone
        if (w > ez) {
          start_ez = w - ez;
        } else {
          start_ez = 0;
        }

        if (w < (d_size - ez)) {
          end_ez = w + ez;
        } else {
          end_ez = d_size - 1;
        }

        uint64_t jump = end - begin;
        uint64_t k = jump + w_size - 1;
        uint64_t s = pow(2, (ceil(log2(k))));

        std::vector<std::complex<double>> rev_query(s);

        for (i = w_size, j = 0; i > 0; i--, j++) {
          rev_query[i - 1] = std::complex<double>(window_ref[w + j], 0.0);
        }
        std::vector<std::complex<double>> Y = fft->fft(rev_query, false);

        if (end > d_size - w_size) { // Last
          jump = (d_size - w_size + 1) /*end*/ - begin;
          if (jump > 1000000) {
            Rcout << "Error on jump" << std::endl;
            Rcout << "begin: " << begin << " end: " << end << std::endl;
            return;
          }

          k = end - begin;
        }

        std::vector<std::complex<double>> data(s);

        for (uint64_t i = 0; i < k; i++) {
          data[i] = std::complex<double>(data_ref[begin + i], 0.0);
        }

        std::vector<std::complex<double>> X = fft->fft(data, false);
        std::vector<std::complex<double>> Z(X.size());
        std::transform(X.begin(), X.end(), Y.begin(), Z.begin(), std::multiplies<std::complex<double>>());
        std::vector<std::complex<double>> z = fft->fft(Z, true);

        for (uint64_t i = 0; i < jump; i++) {
          if (skip_location[begin + i] == 1 || d_std[begin + i] < DBL_EPSILON || q_std[w] < DBL_EPSILON) {
            dp = R_PosInf;
          } else if (ez == 0 || (begin + i) < start_ez || end_ez < (begin + i)) {
            dp = 2 * (w_size - (z[k - jump + i].real() - w_size * d_mean[begin + i] * q_mean[w]) /
                                   (d_std[begin + i] * q_std[w]));
            if (dp < 0) {
              dp = 0;
            }
            if (dp < mp[begin + i]) {
              mp[begin + i] = dp;
              pi[begin + i] = w + 1;
            }
          }
        }
      }
      delete (fft);
    } catch (RcppThread::UserInterruptException &e) {
      delete (fft);
      Rcout << "Computation interrupted by the user." << std::endl;
      Rcout << "Please wait for other threads to stop." << std::endl;
      throw;
    }
  }
};

// [[Rcpp::export]]
List stamp_rcpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                         double s_size, bool progress) {

  uint64_t data_size = data_ref.length();
  uint64_t matrix_profile_size = data_size - window_size + 1;
  uint64_t exclusion_zone = round(window_size * ez + DBL_EPSILON);
  bool partial = false;

  IntegerVector skip_location(matrix_profile_size, 0);

  for (uint64_t i = 0; i < matrix_profile_size; i++) {
    NumericVector range = data_ref[Range(i, (i + window_size - 1))];
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

  NumericVector matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector profile_index(matrix_profile_size, -1);
  List pre = mass_pre_rcpp(data, query, window_size);

  uint64_t k = set_k_rcpp(window_size, data_size, window_size);

  uint64_t jobs = pow(2, (ceil(log2((double)data_size / (double)k)) - 1));
  uint64_t steps = ceil((double)matrix_profile_size / (double)window_size);

  Progress p(jobs * steps, progress);

  StampWorker stamp_worker(data, query, pre["window_size"], data_size, pre["data_mean"], pre["data_sd"],
                           pre["query_mean"], pre["query_sd"], skip_location, exclusion_zone, &p, matrix_profile,
                           profile_index);

  // call parallelFor to do the work
  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, data.size(), stamp_worker, 2 * k);
#else
    RcppParallel2::ttParallelFor(0, data.size(), stamp_worker, 2 * k);
#endif
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results "
             "were returned.";
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("matrix_profile") = sqrt(matrix_profile),
                       Rcpp::Named("profile_index") = profile_index, Rcpp::Named("partial") = partial,
                       Rcpp::Named("ez") = ez));
}
