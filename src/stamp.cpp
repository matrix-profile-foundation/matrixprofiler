#include "stamp.h"
#include "mass.h"
#include "math.h"
#include "fft.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

#if RCPP_PARALLEL_USE_TBB
#include "tbb/mutex.h"
#else
#include "tthread/fast_mutex.h"
#endif

// [[Rcpp::export]]
List stamp_cpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
               double ez) {

  double exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();
  uint64_t matrix_profile_size = data_size - window_size + 1;
  uint64_t num_queries = query_size - window_size + 1;


// check skip position
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
  List pre = mass_pre_rcpp(data, query, window_size);

  IntegerVector order = Range(0, num_queries - 1);
  order = sample(order, num_queries);

  uint32_t k = find_best_k_rcpp(data, query, window_size);

  try {
    for (int32_t i : order) {
      List nn = mass3_rcpp(query[Range(i, i + window_size - 1)], data, pre["data_size"], pre["window_size"],
                           pre["data_mean"], pre["data_sd"], as<NumericVector>(pre["query_mean"])[i],
                           as<NumericVector>(pre["query_sd"])[i], k);

      NumericVector distance_profile = sqrt(as<NumericVector>(nn["distance_profile"]));

      // apply exclusion zone
      if (exclusion_zone > 0) {
        uint64_t exc_st = MAX(0, i - exclusion_zone);
        uint64_t exc_ed = MIN(matrix_profile_size - 1, i + exclusion_zone);
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
      profile_index[which(idx)] = i + 1;
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(
            Rcpp::Named("matrix_profile") = matrix_profile,
            Rcpp::Named("profile_index") = profile_index
          ));
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

#if RCPP_PARALLEL_USE_TBB
  tbb::mutex m;
#else
  tthread::fast_mutex m;
#endif
  // output
  RVector<double> mp;
  RVector<int> pi;

  bool exit = false;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  StampWorker(const NumericVector data_ref, const NumericVector window_ref, const uint64_t w_size, const uint64_t d_size,
              const NumericVector d_mean, const NumericVector d_std, const NumericVector q_mean,
              const NumericVector q_std, NumericVector mp, IntegerVector pi) :
    data_ref(data_ref), window_ref(window_ref), w_size(w_size), d_size(d_size), d_mean(d_mean), d_std(d_std),
    q_mean(q_mean), q_std(q_std), mp(mp), pi(pi) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    // begin and end are the indexes of sliding window.

    uint64_t i, j;
    double dp;

    FFT::fftw *fft = new FFT::fftw();

    try {
      for (uint64_t w = 0; w < d_mean.size() && !exit; w++) {

        if (w % 100 == 0) {
          Rcpp::checkUserInterrupt();
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
          dp = 2 * (w_size - (z[k - jump + i].real() - w_size * d_mean[begin + i] * q_mean[w])
                    / (d_std[begin + i] * q_std[w]));

          mp[begin + i] = MIN(mp[begin + i], dp);
        }
      }

      delete(fft);
    } catch (Rcpp::internal::InterruptedException &ex) {
      delete(fft);
      Rcout << "Thread terminated.\n";
      m.lock();
      exit = true;
      m.unlock();
    }
  }
};

// [[Rcpp::export]]
NumericVector stamp_cpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
                                 double ez) {

  // double exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();
  uint64_t matrix_profile_size = data_size - window_size + 1;
  uint64_t num_queries = query_size - window_size + 1;


  // check skip position
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
  List pre = mass_pre_rcpp(data, query, window_size);

  IntegerVector order = Range(0, num_queries - 1);
  order = sample(order, num_queries);

  StampWorker stamp_worker(data, query, pre["window_size"], data.size(),
                           pre["data_mean"], pre["data_sd"], pre["query_mean"], pre["query_sd"],
                           matrix_profile, profile_index);


  // call parallelFor to do the work
  try {
    RcppParallel::parallelFor(0, data.size(), stamp_worker);
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated.\n";
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (sqrt(matrix_profile));
}
