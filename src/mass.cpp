#include "math.h" // math first to fix OSX error
#include "mass.h"
#include "fft.h"
#include "windowfunc.h"
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
using namespace RcppParallel;
#include <Rcpp/Benchmark/Timer.h>
#if RCPP_PARALLEL_USE_TBB
#include "tbb/mutex.h"
#else
#include "rcpp_parallel_fix.h"
#include "tthread/tinythread.h"
#endif

// #### MASS functions ####
//[[Rcpp::export]]
List mass_weighted_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint32_t data_size,
                        uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd,
                        double query_mean, double query_sd, const NumericVector data_pre, const NumericVector weight,
                        const bool normalized) {
  NumericVector distance_profile;
  NumericVector last_product;

  NumericVector query;
  NumericVector rev_query(data_fft.length());
  NumericVector rev_weight(data_fft.length());

  // normalized query
  if (normalized) {
    query = (query_window - query_mean) / query_sd;
  } else {
    query = query_window;
  }

  double sumwy = sum(query * weight);
  double sumwy2 = sum(weight * query * query);

  std::reverse_copy(query.begin(), query.end(), rev_query.begin());
  std::reverse_copy(weight.begin(), weight.end(), rev_weight.begin());

  ComplexVector prod = data_fft * fft_rcpp(rev_weight * rev_query);
  NumericVector data_queryw = Re(fft_rcpp(prod, true));

  IntegerVector range_d = Range(window_size - 1, data_size - 1);

  last_product = data_queryw[range_d];

  distance_profile = data_pre - 2 * (last_product - sumwy * data_mean) / data_sd + sumwy2;
  distance_profile[distance_profile < 0] = 0;

  return (List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product));
}

//[[Rcpp::export]]
List mass_absolute_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint32_t data_size,
                        uint32_t window_size, const NumericVector sumx2, double sumy2) {
  NumericVector distance_profile;
  NumericVector last_product;

  // pre-process query for fft
  try {
    NumericVector rev_query(data_fft.length()); // use data_fft because is padded to a power of 2
    std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());

    // compute the product
    ComplexVector prod = data_fft * fft_rcpp(rev_query);
    NumericVector z = Re(fft_rcpp(prod, true));
    // compute the distance profile
    IntegerVector range_z = Range(window_size - 1, data_size - 1);
    IntegerVector range_d = Range(0, data_size - window_size);
    last_product = z[range_z];
    distance_profile = as<NumericVector>(sumx2[range_d]) - 2 * last_product + sumy2;
    distance_profile[distance_profile < 0] = 0;
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product));
}

//[[Rcpp::export]]
List mass2_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint64_t data_size,
                uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd, double query_mean,
                double query_sd) {
  NumericVector distance_profile;
  NumericVector last_product;
  double q_mean, q_std = 0;

  // compute query_window stats -- O(d_size)
  q_mean = query_mean;
  q_std = query_sd;

  // pre-process query for fft
  try {
    NumericVector rev_query(data_fft.length()); // use data_fft because is padded to a power of 2
    std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());

    // compute the product
    ComplexVector prod = data_fft * fft_rcpp(rev_query);
    NumericVector z = Re(fft_rcpp(prod, true));
    // compute the distance profile
    last_product = z[Range(window_size - 1, data_size - 1)];

    distance_profile = 2 * (window_size - (last_product - window_size * data_mean * q_mean) / (data_sd * q_std));
    distance_profile[distance_profile < 0] = 0;
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product));
}

//[[Rcpp::export]]
List mass3_rcpp(const NumericVector query_window, const NumericVector data_ref, uint64_t data_size,
                uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd, double query_mean,
                double query_sd, uint32_t k) {
  // data_ref is the long time series
  // query_window is the query
  // k is the size of pieces, preferably a power of two
  // data_ref is the data, query_window is the query
  //

  uint32_t w_size = window_size;
  uint64_t d_size = data_size;
  NumericVector dist(data_mean.length());
  NumericVector::iterator dist_it = dist.begin();
  NumericVector last(data_mean.length());
  NumericVector::iterator last_it = last.begin();

  k = set_k_rcpp(k, d_size, w_size);
  double q_mean, q_std = 0;

  // compute query_window stats -- O(d_size)
  q_mean = query_mean;
  q_std = query_sd;

  // compute data stats -- O(d_size)
  NumericVector d_mean = data_mean;
  NumericVector d_std = data_sd;

  NumericVector rev_query(k);
  std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());
  ComplexVector Y = fft_rcpp(rev_query);

  uint64_t j = 0;
  uint64_t jump = k - w_size + 1;
  uint64_t seq_end = d_size - k;
  ComplexVector Z;
  NumericVector z;
  NumericVector d;

  try {
    for (j = 0; j <= seq_end; j = j + jump) {
      // The main trick of getting dot products in O(d_size log d_size) time
      uint64_t idx_begin = j;
      uint64_t idx_end = j + k - w_size;                         // idx_begin + jump - 1
      ComplexVector X = fft_rcpp(data_ref[Range(j, j + k - 1)]); // idx_begin:(idx_begin + jump + w_size - 2)

      Z = X * Y;

      z = Re(fft_rcpp(Z, true));

      d = 2 * (w_size - (z[Range(w_size - 1, k - 1)] - w_size * d_mean[Range(idx_begin, idx_end)] * q_mean) /
                            (d_std[Range(idx_begin, idx_end)] * q_std));
      std::copy(d.begin(), d.end(), dist_it + j);
      std::copy(z.begin() + w_size - 1, z.begin() + k, last_it + j);
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  jump = d_size - j;

  try {
    if (jump >= w_size) {
      uint64_t idx_begin = j;
      uint64_t idx_end = d_size - w_size;

      if ((jump - (w_size - 1) + j) > (uint64_t)data_mean.length()) {
        Rcout << "DEBUG: error." << std::endl;
      } else {
        ComplexVector X = fft_rcpp(data_ref[Range(j, d_size - 1)]);
        Y = fft_rcpp(rev_query[Range(0, jump - 1)]);

        Z = X * Y;
        z = Re(fft_rcpp(Z, true));

        d = 2 * (w_size - (z[Range(w_size - 1, jump - 1)] - w_size * d_mean[Range(idx_begin, idx_end)] * q_mean) /
                              (d_std[Range(idx_begin, idx_end)] * q_std));

        std::copy(d.begin(), d.end(), dist_it + j);
        std::copy(z.begin() + w_size - 1, z.begin() + jump, last_it + j);
      }
    }

    dist[dist < 0] = 0;
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("distance_profile") = dist, Rcpp::Named("last_product") = last));
}

struct MassWorker : public Worker {
  // input
  const RVector<double> data_ref;
  const RVector<double> window_ref;
  const uint64_t w_size;
  const uint64_t d_size;
  const RVector<double> d_mean;
  const RVector<double> d_std;
  const double q_mean;
  const double q_std;

#if RCPP_PARALLEL_USE_TBB
  tbb::mutex m;
#else
  tthread::mutex m;
#endif

  std::vector<std::complex<double>> Y;
  // output
  RVector<double> dp;
  RVector<double> lp;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  MassWorker(const NumericVector data_ref, const NumericVector window_ref, const uint64_t w_size, const uint64_t d_size,
             const NumericVector d_mean, const NumericVector d_std, const double q_mean, const double q_std,
             NumericVector dp, NumericVector lp)
      : data_ref(data_ref), window_ref(window_ref), w_size(w_size), d_size(d_size), d_mean(d_mean), d_std(d_std),
        q_mean(q_mean), q_std(q_std), dp(dp), lp(lp) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    uint64_t jump = end - begin;
    uint64_t k = jump + w_size - 1;
    uint64_t s = pow(2, (ceil(log2(k))));

    FFT::fftw *fft = new FFT::fftw();

    if (end > d_size - w_size) { // Last
      jump = (d_size - w_size + 1) /*end*/ - begin;
      k = end - begin;
    }

    std::vector<std::complex<double>> data(s);

    for (uint64_t i = 0; i < k; i++) {
      data[i] = std::complex<double>(data_ref[begin + i], 0.0);
    }

    m.lock();
    if (Y.size() == 0) {
      std::vector<std::complex<double>> rev_query(s);

      uint64_t j = 0;
      for (uint64_t i = w_size; i > 0; i--, j++) {
        rev_query[i - 1] = std::complex<double>(window_ref[j], 0.0);
      }
      Y = fft->fft(rev_query, false);
    }
    m.unlock();

    std::vector<std::complex<double>> X = fft->fft(data, false);
    std::vector<std::complex<double>> Z(X.size());
    std::transform(X.begin(), X.end(), Y.begin(), Z.begin(), std::multiplies<std::complex<double>>());
    std::vector<std::complex<double>> z = fft->fft(Z, true);

    for (uint64_t i = 0; i < jump; i++) {
      dp[begin + i] =
          2 * (w_size - (z[k - jump + i].real() - w_size * d_mean[begin + i] * q_mean) / (d_std[begin + i] * q_std));
      lp[begin + i] = z[k - jump + i].real();
    }
    delete (fft);
  }
};

//[[Rcpp::export]]
List mass3_rcpp_parallel(const NumericVector query_window, const NumericVector data_ref, uint64_t data_size,
                         uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd,
                         double query_mean, double query_sd, uint16_t k) {
  try {
    k = set_k_rcpp(k, data_size, window_size);

    // allocate the output matrix
    NumericVector distance_profile(data_mean.length());
    NumericVector last_product(data_mean.length());

    // SquareRoot functor (pass input and output matrixes)
    MassWorker mass_worker(data_ref, query_window, window_size, data_size, data_mean, data_sd, query_mean, query_sd,
                           distance_profile, last_product);

    // call parallelFor to do the work
    try {
#if RCPP_PARALLEL_USE_TBB
      RcppParallel::parallelFor(0, data_size, mass_worker, k);
#else
      RcppParallel2::ttParallelFor(0, data_size, mass_worker, k);
#endif

      distance_profile[distance_profile < 0] = 0;

    } catch (RcppThread::UserInterruptException &ex) {
      Rcout << "Process terminated.\n";
    } catch (...) {
      ::Rf_error("c++ exception (unknown reason)");
    }

    return (
        List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product));
  } catch (std::out_of_range &ex) {
    ::Rf_error("out_of_range\\n");
  } catch (...) {
    ::Rf_error("Interrupted\\n");
  }
}

//[[Rcpp::export]]
uint32_t set_k_rcpp(uint32_t k, uint64_t data_size, uint64_t window_size) {
  try {
    if (k > data_size) {
      k = pow(2, ceil(log2(sqrt((double)data_size)))); // added double inside sqrt to avoid ambiguity on Solaris
    }

    if (k <= window_size) {
      k = pow(2, (ceil(log2(window_size)) + 1));
      if (k > data_size) {
        k = data_size;
      }
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (k);
}

//[[Rcpp::export]]
uint32_t find_best_k_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size) {
  uint64_t data_size = data_ref.length();
  uint32_t k = set_k_rcpp(window_size, data_size, window_size); // Set baseline
  uint64_t best_time = (int)pow((double)2, (double)50); // added double inside sqrt to avoid ambiguity on Solaris
  uint32_t best_k = k;
  List pre = mass_pre_rcpp(data_ref, query_ref, window_size);
  Timer timer;

  try {
    for (uint16_t j = 0; j < 10; j++) {
      uint64_t tictoc = timer.now();
      for (uint16_t i = 0; i < 10; i++) {
        List nn = mass3_rcpp(query_ref[Range(i, i + window_size - 1)], data_ref, pre["data_size"], pre["window_size"],
                             pre["data_mean"], pre["data_sd"], as<NumericVector>(pre["query_mean"])[i],
                             as<NumericVector>(pre["query_sd"])[i], k);
      }
      uint64_t time_res = timer.now() - tictoc;
      if (time_res < best_time) {
        best_time = time_res;
        best_k = k;
        k = k * 2;
        if (k > data_size) {
          break;
        }
      } else {
        break;
      }
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return best_k;
}

// #### MASS PRE functions ####

//[[Rcpp::export]]
List mass_pre_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size) {
  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();

  List data_avgsd = movmean_std_rcpp(data_ref, window_size); // precompute moving average and SD

  uint32_t pad_size = pow(2, ceil(log2(data_size))); // padded to a power of 2
  NumericVector data_padded(pad_size);

  std::copy(data_ref.begin(), data_ref.end(), data_padded.begin());

  ComplexVector data_fft = fft_rcpp(data_padded); // precompute fft of data

  NumericVector query_mean;
  NumericVector query_sd;

  try {
    if (query_size > 0) {
      List query_avgsd = movmean_std_rcpp(query_ref, window_size); // precompute moving average and SD
      query_mean = query_avgsd["avg"];
      query_sd = query_avgsd["sd"];
    } else {
      query_mean = data_avgsd["avg"];
      query_sd = data_avgsd["sd"];
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("data_fft") = data_fft, Rcpp::Named("data_size") = data_size,
                       Rcpp::Named("window_size") = window_size, Rcpp::Named("data_mean") = data_avgsd["avg"],
                       Rcpp::Named("data_sd") = data_avgsd["sd"], Rcpp::Named("query_mean") = query_mean,
                       Rcpp::Named("query_sd") = query_sd));
}

//[[Rcpp::export]]
List mass_pre_abs_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size) {
  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();

  uint32_t pad_size = pow(2, ceil(log2(data_size))); // padded to a power of 2
  NumericVector data_padded(pad_size);

  std::copy(data_ref.begin(), data_ref.end(), data_padded.begin());

  ComplexVector data_fft = fft_rcpp(data_padded); // precompute fft of data
  NumericVector sumx2 = movsum_ogita_rcpp(data_ref * data_ref, window_size);
  NumericVector sumy2;
  try {
    if (query_size > 0) {
      sumy2 = movsum_ogita_rcpp(query_ref * query_ref, window_size);
    } else {
      sumy2 = sumx2;
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("data_fft") = data_fft, Rcpp::Named("window_size") = window_size,
                       Rcpp::Named("data_size") = data_size, Rcpp::Named("sumx2") = sumx2,
                       Rcpp::Named("sumy2") = sumy2));
}

//[[Rcpp::export]]
List mass_pre_weighted_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
                            const NumericVector weight) {
  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();

  List data_avgsd = movmean_std_rcpp(data_ref, window_size); // precompute moving average and SD

  uint32_t pad_size = pow(2, ceil(log2(data_size))); // padded to a power of 2
  NumericVector data_padded(pad_size);
  NumericVector rev_weight(pad_size);
  double sumw = sum(weight);

  std::reverse_copy(weight.begin(), weight.end(), rev_weight.begin());
  std::copy(data_ref.begin(), data_ref.end(), data_padded.begin());

  ComplexVector data_fft = fft_rcpp(data_padded); // precompute fft of data
  ComplexVector w_fft = fft_rcpp(rev_weight);
  ComplexVector data_fft_w_fft = data_fft * w_fft;

  NumericVector query_mean;
  NumericVector query_sd;

  try {
    if (query_size > 0) {
      List query_avgsd = movmean_std_rcpp(query_ref, window_size); // precompute moving average and SD
      query_mean = query_avgsd["avg"];
      query_sd = query_avgsd["sd"];
    } else {
      query_mean = data_avgsd["avg"];
      query_sd = data_avgsd["sd"];
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  IntegerVector range_s = Range(window_size - 1, data_size - 1);

  NumericVector data_w = Re(fft_rcpp(data_fft_w_fft, true));
  ComplexVector data2_fft = fft_rcpp(data_padded * data_padded);
  ComplexVector data2_fft_w_fft = data2_fft * w_fft;
  NumericVector data2_w = Re(fft_rcpp(data2_fft_w_fft, true));
  NumericVector sumxw2 = data2_w[range_s];
  NumericVector sumxw = data_w[range_s];

  NumericVector data_mean = data_avgsd["avg"];
  NumericVector data_sd = data_avgsd["sd"];

  NumericVector data_pre = (sumxw2 - 2 * sumxw * data_mean + sumw * (data_mean * data_mean)) / (data_sd * data_sd);

  return (List::create(
      Rcpp::Named("data_fft") = data_fft, Rcpp::Named("data_pre") = data_pre, Rcpp::Named("data_size") = data_size,
      Rcpp::Named("window_size") = window_size, Rcpp::Named("data_mean") = data_mean, Rcpp::Named("data_sd") = data_sd,
      Rcpp::Named("query_mean") = query_mean, Rcpp::Named("query_sd") = query_sd, Rcpp::Named("weight") = weight));
}
