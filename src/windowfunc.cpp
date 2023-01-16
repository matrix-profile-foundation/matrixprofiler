// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

/*===========================================================================*/
/* Original source from caTools library:                                     */
/* runfunc - running window functions                                        */
/* Copyright (C) 2005 Jarek Tuszynski                                        */
/* Distributed under GNU General Public License version 3                    */
/* Converted to Rcpp by Francisco Bischoff                                   */
/*===========================================================================*/

// Supports:
//               NA/NaN  -Inf/Inf  Edge
// movmin          Unk     Unk      No
// movmax          Unk     Unk      No

#include "mathtools.h" // math first to fix OSX error
#include "windowfunc.h"

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

#if RCPP_PARALLEL_USE_TBB
#include "tbb/mutex.h"
#else
#include "rcpp_parallel_fix.h"
#include "tthread/tinythread.h"
#endif

//[[Rcpp::export]]
NumericVector movmean_rcpp(const NumericVector data, const uint32_t window_size) {
  uint32_t const data_size = data.length();

  NumericVector out(data_size - window_size + 1);

  double n = 0.0;
  double sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    sum = sum + data[i];
    n = n + 1;

    if (i >= window_size) {
      sum = sum - data[i - window_size];
      n = n - 1;
    }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = sum / n;
    }
  }
  return (out);
}

//[[Rcpp::export]]
NumericVector movstd_rcpp(const NumericVector data, const uint32_t window_size) {

  NumericVector const mu = movsum_ogita_rcpp(data, window_size) / window_size;
  NumericVector const data2_sum = movsum_ogita_rcpp(data * data, window_size);
  NumericVector const data_var = (data2_sum / window_size) - (mu * mu);
  NumericVector const data_sd = sqrt(data_var);

  return (data_sd);
}

//[[Rcpp::export]]
List movmean_std_rcpp(const NumericVector data, const uint32_t window_size) {

  NumericVector const data_sum = movsum_ogita_rcpp(data, window_size);
  NumericVector const data_mean = data_sum / window_size;
  NumericVector const data2 = (data * data);
  NumericVector const data2_sum = movsum_ogita_rcpp(data2, window_size);

  NumericVector const data_var = (data2_sum / window_size) - (data_mean * data_mean); // variance
  NumericVector const data_sd = sqrt(data_var);
  NumericVector const data_sig = sqrt(1 / (data_var * window_size));

  return (List::create(Rcpp::Named("avg") = data_mean, Rcpp::Named("sd") = data_sd, Rcpp::Named("sig") = data_sig,
                       Rcpp::Named("sum") = data_sum, Rcpp::Named("sqrsum") = data2_sum));
}

//[[Rcpp::export]]
NumericVector movvar_rcpp(const NumericVector data, const uint32_t window_size) {

  NumericVector const mu = movsum_ogita_rcpp(data, window_size) / window_size;
  NumericVector const data2_sum = movsum_ogita_rcpp(data * data, window_size);
  NumericVector const data_var = (data2_sum / window_size) - (mu * mu);

  return (data_var);
}

//[[Rcpp::export]]
NumericVector movvar2_rcpp(const NumericVector data, uint32_t window_size) {

  uint32_t const data_size = data.length();

  NumericVector out(data_size - window_size + 1);

  double n = 0.0;
  double data_sum = 0.0;
  double data2_sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    data_sum = data_sum + data[i];
    data2_sum = data2_sum + (data[i] * data[i]);
    n = n + 1;
    if (i >= window_size) {
      data_sum = data_sum - data[i - window_size];
      data2_sum = data2_sum - (data[i - window_size] * data[i - window_size]);
      n = n - 1;
    }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = data2_sum / n - ((data_sum * data_sum) / (n * n));
    }
  }
  return (out);
}

//[[Rcpp::export]]
NumericVector movsum_rcpp(NumericVector data, uint32_t window_size) {
  uint32_t const data_size = data.length();

  NumericVector out(data_size - window_size + 1);

  double sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    sum = sum + data[i];

    if (i >= window_size) {
      sum = sum - data[i - window_size];
    }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = sum;
    }
  }
  return (out);
}

//[[Rcpp::export]]
NumericVector movsum_ogita_rcpp(const NumericVector data, uint32_t window_size) {
  NumericVector res(data.length() - window_size + 1, 0);
  double accum = data[0];
  double resid = 0.0;

  for (uint32_t i = 1; i < window_size; i++) {
    double const m = data[i];
    double const p = accum;
    accum = accum + m;
    double const q = accum - p;
    resid = resid + ((p - (accum - q)) + (m - q));
  }

  if (resid > 0.001) {
    Function const warning("warning");
    warning("Residual value is large. Some precision may be lost. res = %f\n", resid);
  }

  res[0] = accum + resid;

  for (int64_t i = window_size; i < data.length(); i++) {
    double const m = data[i - window_size];
    double const n = data[i];
    double const p = accum - m;
    double const q = p - accum;
    double const r = resid + ((accum - (p - q)) - (m + q));
    accum = p + n;
    double const t = accum - p;
    resid = r + ((p - (accum - t)) + (n - t));
    res[i - window_size + 1] = accum + resid;
  }

  return (res);
}

//[[Rcpp::export]]
double precision_test_rcpp(std::vector<double> dd) {

  std::vector<double> data(dd);

  double const sum = std::accumulate(data.begin(), data.end(), 0.0);
  double const mean = sum / static_cast<double>(data.size());

  for (double &i : data) {
    i = i - mean;
  }

  double const out = std::accumulate(data.begin(), data.end(), 0.0);

  return (out);
}

//[[Rcpp::export]]
NumericVector movmin_rcpp(const NumericVector data, uint32_t window_size) {
  uint32_t const data_size = data.length();

  if (window_size <= 1) {
    return data;
  }

  if (window_size > data_size) {
    window_size = data_size;
  }

  uint32_t i = 0UL, d = 0UL, k = 0UL;
  double min_res = 0.0, res_out = 0.0;
  NumericVector out(data_size - window_size + 1);

  k = 0;
  d = 0;

  min_res = res_out = R_PosInf;
  for (i = window_size - 1; i < data_size; i++) {
    // if point coming out of the window was window's min than we need to
    // recalculate 'min_res'
    if (res_out == min_res) {
      // find minimum over a window of length 'window_size'
      min_res = *std::min_element(data.begin() + d, data.begin() + d + window_size); // cpp11, faster
    } else {
      // if point coming out of the window was NOT window min than min of
      // window's first
      //  'window_size - 1' points is still 'min_res', so we have to add a
      //  single point
      min_res = MIN(min_res, data[d + window_size - 1]);
    }

    res_out = data[d++];                                   // store point comming out of the window for future use
                                                           // and move window
    out[k++] = (min_res == R_PosInf ? R_NaReal : min_res); // save 'min_res' and move window
  }

  return out;
}

//[[Rcpp::export]]
NumericVector movmax_rcpp(const NumericVector data, uint32_t window_size) {
  uint32_t const data_size = data.length();

  if (window_size <= 1) {
    return data;
  }

  if (window_size > data_size) {
    window_size = data_size;
  }

  uint32_t i = 0UL, d = 0UL, k = 0UL;
  double max_res = 0.0, res_out = 0.0;
  NumericVector out(data_size - window_size + 1);

  k = 0;
  d = 0;

  max_res = res_out = R_NegInf;
  for (i = window_size - 1; i < data_size; i++) {
    // if point coming out of the window was window's max than we need to
    // recalculate 'max_res'
    if (res_out == max_res) {
      // find maximum over a window of length 'window_size'
      max_res = *std::max_element(data.begin() + d, data.begin() + d + window_size); // cpp11, faster
    } else {
      // if point coming out of the window was NOT window max than max of
      // window's first
      //  'window_size - 1' points is still 'max_res', so we have to add a
      //  single point
      max_res = MAX(max_res, data[d + window_size - 1]);
    }

    res_out = data[d++];                                   // store point comming out of the window for future use
                                                           // and move window
    out[k++] = (max_res == R_NegInf ? R_NaReal : max_res); // save 'max_res' and move window
  }

  return out;
}

// #### ONLINE ALGORITHMS --------

//[[Rcpp::export]]
NumericVector movmean_weighted_rcpp(const NumericVector data, uint32_t window_size, double eps) {

  uint32_t const data_size = data.length();

  double const w = window_size;

  double const alpha = pow(eps, 1 / w);
  NumericVector out(data_size - window_size + 1);

  double n = 0.0;
  double sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    sum = sum * alpha + data[i];
    n = n * alpha + 1;

    if (i >= window_size) {
      sum = sum - data[i - window_size] * pow(alpha, window_size - 1);
      n = n - 1 * pow(alpha, window_size - 1);
    }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = sum / n;
    }
  }
  return (out);
}

//[[Rcpp::export]]
NumericVector movmean_fading_rcpp(const NumericVector data, uint32_t window_size, double eps) {

  uint32_t const data_size = data.length();

  double const w = window_size;

  double const alpha = pow(eps, 1 / w);
  NumericVector out(data_size - window_size + 1);

  double n = 0.0;
  double sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    sum = sum * alpha + data[i];
    n = n * alpha + 1;

    // if (i >= window_size) {
    //   sum = sum - data[i - window_size] * pow(alpha, window_size - 1);
    //   n = n - 1 * pow(alpha, window_size - 1);
    // }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = sum / n;
    }
  }
  return (out);
}

//[[Rcpp::export]]
NumericVector movsum_weighted_rcpp(NumericVector data, uint32_t window_size, double eps) {
  uint32_t const data_size = data.length();

  double const w = window_size;

  double const alpha = pow(eps, 1 / w);
  NumericVector out(data_size - window_size + 1);

  double sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    sum = sum * alpha + data[i];

    if (i >= window_size) {
      sum = sum - data[i - window_size] * pow(alpha, window_size - 1);
    }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = sum;
    }
  }
  return (out);
}

//[[Rcpp::export]]
NumericVector movsum_fading_rcpp(NumericVector data, uint32_t window_size, double eps) {

  uint32_t const data_size = data.length();

  double const w = window_size;

  double const alpha = pow(eps, 1 / w);
  NumericVector out(data_size - window_size + 1);

  double sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    sum = sum * alpha + data[i];

    // if (i >= window_size) {
    //   sum = sum - data[i - window_size] * pow(alpha, window_size - 1);
    // }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = sum;
    }
  }
  return (out);
}

//[[Rcpp::export]]
NumericVector movvar_weighted_rcpp(const NumericVector data, uint32_t window_size, double eps) {

  uint32_t const data_size = data.length();

  double const w = window_size;

  double const alpha = pow(eps, 1 / w);
  NumericVector out(data_size - window_size + 1);

  double n = 0.0;
  double data_sum = 0.0;
  double data2_sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    data_sum = data_sum * alpha + data[i];
    data2_sum = data2_sum * alpha + (data[i] * data[i]);
    n = n * alpha + 1;

    if (i >= window_size) {
      data_sum = data_sum - data[i - window_size] * pow(alpha, window_size - 1);
      data2_sum = data2_sum - (data[i - window_size] * data[i - window_size]) * pow(alpha, window_size - 1);
      n = n - 1 * pow(alpha, window_size - 1);
    }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = data2_sum / n - ((data_sum * data_sum) / (n * n));
    }
  }
  return (out);
}

//[[Rcpp::export]]
NumericVector movvar_fading_rcpp(const NumericVector data, uint32_t window_size, double eps) {

  uint32_t const data_size = data.length();

  double const w = window_size;

  double const alpha = pow(eps, 1 / w);
  NumericVector out(data_size - window_size + 1);

  double n = 0.0;
  double data_sum = 0.0;
  double data2_sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    data_sum = data_sum * alpha + data[i];
    data2_sum = data2_sum * alpha + (data[i] * data[i]);
    n = n * alpha + 1;

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = data2_sum / n - ((data_sum * data_sum) / (n * n));
    }
  }
  return (out);
}

//[[Rcpp::export]]
List muinvn_rcpp(const NumericVector data, uint32_t window_size) {
  // Functions here are based on the work in
  // Ogita et al, Accurate Sum and Dot Product
  // results here are a moving average and stable inverse centered norm based
  // on Accurate Sum and Dot Product, Ogita et al

  // NumericVector sig(data.length() - window_size + 1, 0);
  NumericVector const mu = movsum_ogita_rcpp(data, window_size) / window_size;
  NumericVector const data2_sum = movsum_ogita_rcpp(data * data, window_size);
  NumericVector const sig = 1 / sqrt(data2_sum - mu * mu * window_size);

  // std is equals to 1 / (sig * sqrt(w))
  // sig is equals to 1 / (std * sqrt(w))

  return (List::create(Rcpp::Named("avg") = mu, Rcpp::Named("sig") = sig));
}

struct MuinWorker : public Worker {

private:
  // input
  const RVector<double> data2_sum;
  const RVector<double> mu;
  const uint32_t w_size;
  // output
  RVector<double> sig;

public:
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  MuinWorker(const NumericVector &data2_sum, const NumericVector &mu, const uint32_t w_size, const NumericVector &sig)
      : data2_sum(data2_sum), mu(mu), w_size(w_size), sig(sig) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) override {
    for (uint32_t i = begin; i < end; i++) {
      sig[i] = 1 / sqrt(data2_sum[i] - mu[i] * mu[i] * w_size);
    }
  }
};

//[[Rcpp::export]]
List muinvn_rcpp_parallel(const NumericVector data, uint32_t window_size) {
  // Functions here are based on the work in
  // Ogita et al, Accurate Sum and Dot Product
  // results here are a moving average and stable inverse centered norm based
  // on Accurate Sum and Dot Product, Ogita et al

  NumericVector const sig(data.length() - window_size + 1);
  NumericVector const mu = movsum_ogita_rcpp(data, window_size) / window_size;
  NumericVector const data2_sum = movsum_ogita_rcpp(data * data, window_size);

  MuinWorker muin_worker(data2_sum, mu, window_size, sig);

  // call parallelFor to do the work
  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, mu.length(), muin_worker);
#else
    RcppParallel2::ttParallelFor(0, mu.length(), muin_worker);
#endif
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated.\n";
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("avg") = mu, Rcpp::Named("sig") = sig));
}

// Counts number of zero-crossings
//
// Count the number of zero-crossings from the supplied time-domain input vector. A simple method is
// applied here that can be easily ported to a real-time system that would minimize the number of
// if-else conditionals.
//
// @param data a `vector` of `numeric`.
//
// @return Returns the amount of zero-crossings in the input signal.
// @author sparafucile17 06/27/04
// @references <https://www.dsprelated.com/showcode/179.php>
// @keywords internal
// @noRd

//[[Rcpp::export]]
IntegerVector zero_crossing_rcpp(const NumericVector data, const uint32_t window_size) {

  uint32_t const profile_size = data.size() - window_size + 1;
  NumericVector norm_data = znorm_rcpp(data);
  IntegerVector crossings(profile_size);

  for (uint64_t j = 0; j < profile_size; j++) {
    uint32_t count = 0;

    for (uint32_t i = j + 1; i < (j + window_size - 1); i++) {
      // Any time you multiply to adjacent values that have a sign difference
      // the result will always be negative.  When the signs are identical,
      // the product will always be positive.
      if ((norm_data[i] * norm_data[i - 1]) < 0) {
        count++;
      }
    }

    crossings[j] = count;
  }

  return (crossings);
}
