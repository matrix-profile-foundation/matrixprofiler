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

#include "windowfunc.h"
#include "math.h" // math first to fix OSX error

//[[Rcpp::export]]
NumericVector movmean_rcpp(const NumericVector data, const uint32_t window_size) {
  uint32_t data_size = data.length();

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
NumericVector movvar_rcpp(const NumericVector data, const uint32_t window_size) {

  // Improve the numerical analysis by subtracting off the series mean
  // this has no effect on the standard deviation.
  NumericVector data_zeromean = data - mean(data);

  NumericVector data_sum =
      cumsum(diff2_lag(data_zeromean, window_size, sum(data_zeromean[Range(0, (window_size - 1))])));
  NumericVector data_mean = data_sum / window_size;

  NumericVector data2 = pow(data_zeromean, 2);
  NumericVector data2_sum = cumsum(diff2_lag(data2, window_size, sum(data2[Range(0, (window_size - 1))])));
  NumericVector data_var = (data2_sum / window_size) - pow(data_mean, 2); // variance

  return (data_var);
}

//[[Rcpp::export]]
NumericVector movvar2_rcpp(const NumericVector data, uint32_t window_size) {

  uint32_t data_size = data.length();

  NumericVector out(data_size - window_size + 1);

  double n = 0.0;
  double data_sum = 0.0;
  double data2_sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    data_sum = data_sum + data[i];
    data2_sum = data2_sum + pow(data[i], 2);
    n = n + 1;
    if (i >= window_size) {
      data_sum = data_sum - data[i - window_size];
      data2_sum = data2_sum - pow(data[i - window_size], 2);
      n = n - 1;
    }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = data2_sum / n - (pow(data_sum, 2) / pow(n, 2));
    }
  }
  return (out);
}

//[[Rcpp::export]]
NumericVector movsum_rcpp(NumericVector data, uint32_t window_size) {
  uint32_t data_size = data.length();

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
NumericVector movsum_ogita_rcpp(const NumericVector a, uint32_t w) {
  NumericVector res(a.length() - w + 1, 0);
  double accum = a[0];
  double resid = 0.0;

  for (uint32_t i = 1; i < w; i++) {
    double m = a[i];
    double p = accum;
    accum = accum + m;
    double q = accum - p;
    resid = resid + ((p - (accum - q)) + (m - q));
  }

  if (resid > 0.001) {
    Rf_warning("Residual value is large. Some precision may be lost. res = %f\n", resid);
  }

  res[0] = accum + resid;

  for (uint32_t i = w; i < a.length(); i++) {
    double m = a[i - w];
    double n = a[i];
    double p = accum - m;
    double q = p - accum;
    double r = resid + ((accum - (p - q)) - (m + q));
    accum = p + n;
    double t = accum - p;
    resid = r + ((p - (accum - t)) + (n - t));
    res[i - w + 1] = accum + resid;
  }

  return (res);
}

//[[Rcpp::export]]
double precision_test_rcpp(std::vector<double> d) {

  std::vector<double> data(d);

  double sum = std::accumulate(data.begin(), data.end(), 0.0);
  double mean = sum / data.size();

  for (uint32_t i = 0; i < data.size(); i++) {
    data[i] = data[i] - mean;
  }

  double out = std::accumulate(data.begin(), data.end(), 0.0);

  return (out);
}

//[[Rcpp::export]]
NumericVector movmin_rcpp(const NumericVector data, uint32_t window_size) {
  uint32_t data_size = data.length();

  if (window_size <= 1) {
    return data;
  }

  if (window_size > data_size) {
    window_size = data_size;
  }

  uint32_t i, d, k;
  double min_res, res_out;
  NumericVector out(data_size - window_size + 1);

  k = 0;
  d = 0;

  min_res = res_out = R_PosInf;
  for (i = window_size - 1; i < data_size; i++) {
    // if point comining out of the window was window's min than we need to
    // recalculate 'min_res'
    if (res_out == min_res) {
      // find minimum over a window of length 'window_size'
      min_res = *std::min_element(data.begin() + d, data.begin() + d + window_size); // cpp11, faster
    } else {
      // if point comining out of the window was NOT window min than min of
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
  uint32_t data_size = data.length();

  if (window_size <= 1) {
    return data;
  }

  if (window_size > data_size) {
    window_size = data_size;
  }

  uint32_t i, d, k;
  double max_res, res_out;
  NumericVector out(data_size - window_size + 1);

  k = 0;
  d = 0;

  max_res = res_out = R_NegInf;
  for (i = window_size - 1; i < data_size; i++) {
    // if point comining out of the window was window's max than we need to
    // recalculate 'max_res'
    if (res_out == max_res) {
      // find maximum over a window of length 'window_size'
      max_res = *std::max_element(data.begin() + d, data.begin() + d + window_size); // cpp11, faster
    } else {
      // if point comining out of the window was NOT window max than max of
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
// ONLINE ALGORITHMS

//[[Rcpp::export]]
NumericVector movmean_weighted_rcpp(const NumericVector data, uint32_t window_size, double eps) {

  uint32_t data_size = data.length();

  double w = window_size;

  double alpha = pow(eps, 1 / w);
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

  uint32_t data_size = data.length();

  double w = window_size;

  double alpha = pow(eps, 1 / w);
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
  uint32_t data_size = data.length();

  double w = window_size;

  double alpha = pow(eps, 1 / w);
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

  uint32_t data_size = data.length();

  double w = window_size;

  double alpha = pow(eps, 1 / w);
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

  uint32_t data_size = data.length();

  double w = window_size;

  double alpha = pow(eps, 1 / w);
  NumericVector out(data_size - window_size + 1);

  double n = 0.0;
  double data_sum = 0.0;
  double data2_sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    data_sum = data_sum * alpha + data[i];
    data2_sum = data2_sum * alpha + pow(data[i], 2);
    n = n * alpha + 1;

    if (i >= window_size) {
      data_sum = data_sum - data[i - window_size] * pow(alpha, window_size - 1);
      data2_sum = data2_sum - pow(data[i - window_size], 2) * pow(alpha, window_size - 1);
      n = n - 1 * pow(alpha, window_size - 1);
    }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = data2_sum / n - (pow(data_sum, 2) / pow(n, 2));
    }
  }
  return (out);
}

//[[Rcpp::export]]
NumericVector movvar_fading_rcpp(const NumericVector data, uint32_t window_size, double eps) {

  uint32_t data_size = data.length();

  double w = window_size;

  double alpha = pow(eps, 1 / w);
  NumericVector out(data_size - window_size + 1);

  double n = 0.0;
  double data_sum = 0.0;
  double data2_sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    data_sum = data_sum * alpha + data[i];
    data2_sum = data2_sum * alpha + pow(data[i], 2);
    n = n * alpha + 1;

    // if (i >= window_size) {
    //   data_sum = data_sum - data[i - window_size] * pow(alpha, window_size -
    //   1); data2_sum = data2_sum - pow(data[i - window_size], 2) * pow(alpha,
    //   window_size - 1); n = n - 1 * pow(alpha, window_size - 1);
    // }

    if (i >= (window_size - 1)) {
      out[i - (window_size - 1)] = data2_sum / n - (pow(data_sum, 2) / pow(n, 2));
    }
  }
  return (out);
}
