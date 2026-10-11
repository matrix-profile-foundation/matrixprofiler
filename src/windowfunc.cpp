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
#include <limits>

// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

class SlidingNonFiniteTracker {
public:
  SlidingNonFiniteTracker(const double *data, uint32_t n, uint32_t w)
      : data_(data), n_(n), w_(w), count_(0) {}

  void update(uint32_t i) {
    if (!std::isfinite(data_[i])) {
      count_++;
    }
    if (i >= w_ && !std::isfinite(data_[i - w_])) {
      count_--;
    }
  }

  uint32_t count() const { return count_; }
  bool has_na() const { return count_ > 0; }
  uint32_t valid_count() const { return w_ - count_; }

private:
  const double *data_;
  uint32_t n_;
  uint32_t w_;
  uint32_t count_;
};

//[[Rcpp::export]]
NumericVector movmean_rcpp(const NumericVector data, const uint32_t window_size, int na_mode) {
  uint32_t const data_size = data.length();
  NumericVector out(data_size - window_size + 1);

  bool has_non_finite = false;
  for (uint32_t i = 0; i < data_size; ++i) {
    if (!std::isfinite(data[i])) {
      has_non_finite = true;
      break;
    }
  }

  if (!has_non_finite) {
    double n = 0.0;
    double sum = 0.0;
    for (uint32_t i = 0; i < data_size; i++) {
      sum += data[i];
      n += 1.0;
      if (i >= window_size) {
        sum -= data[i - window_size];
        n -= 1.0;
      }
      if (i >= (window_size - 1)) {
        out[i - (window_size - 1)] = sum / n;
      }
    }
    return out;
  }

  if (na_mode == 2) {
    Rcpp::stop("Data contains missing/non-finite values.");
  }

  const double *data_ptr = data.begin();
  SlidingNonFiniteTracker tracker(data_ptr, data_size, window_size);
  double sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    if (std::isfinite(data_ptr[i])) {
      sum += data_ptr[i];
    }
    if (i >= window_size && std::isfinite(data_ptr[i - window_size])) {
      sum -= data_ptr[i - window_size];
    }
    tracker.update(i);

    if (i >= (window_size - 1)) {
      uint32_t out_idx = i - (window_size - 1);
      if (na_mode == 0) { // window_na
        out[out_idx] = tracker.has_na() ? NA_REAL : (sum / window_size);
      } else { // remove
        uint32_t v_n = tracker.valid_count();
        out[out_idx] = (v_n == 0) ? NA_REAL : (sum / v_n);
      }
    }
  }
  return out;
}

//[[Rcpp::export]]
NumericVector movsum_rcpp(NumericVector data, uint32_t window_size, int na_mode) {
  uint32_t const data_size = data.length();
  NumericVector out(data_size - window_size + 1);

  bool has_non_finite = false;
  for (uint32_t i = 0; i < data_size; ++i) {
    if (!std::isfinite(data[i])) {
      has_non_finite = true;
      break;
    }
  }

  if (!has_non_finite) {
    double sum = 0.0;
    for (uint32_t i = 0; i < data_size; i++) {
      sum += data[i];
      if (i >= window_size) {
        sum -= data[i - window_size];
      }
      if (i >= (window_size - 1)) {
        out[i - (window_size - 1)] = sum;
      }
    }
    return out;
  }

  if (na_mode == 2) {
    Rcpp::stop("Data contains missing/non-finite values.");
  }

  const double *data_ptr = data.begin();
  SlidingNonFiniteTracker tracker(data_ptr, data_size, window_size);
  double sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    if (std::isfinite(data_ptr[i])) {
      sum += data_ptr[i];
    }
    if (i >= window_size && std::isfinite(data_ptr[i - window_size])) {
      sum -= data_ptr[i - window_size];
    }
    tracker.update(i);

    if (i >= (window_size - 1)) {
      uint32_t out_idx = i - (window_size - 1);
      if (na_mode == 0) { // window_na
        out[out_idx] = tracker.has_na() ? NA_REAL : sum;
      } else { // remove
        out[out_idx] = (tracker.valid_count() == 0) ? NA_REAL : sum;
      }
    }
  }
  return out;
}

//[[Rcpp::export]]
NumericVector movsum_ogita_rcpp(const NumericVector data, uint32_t window_size, int na_mode) {
  uint32_t const data_size = data.length();
  NumericVector res(data_size - window_size + 1, 0.0);

  bool has_non_finite = false;
  for (uint32_t i = 0; i < data_size; ++i) {
    if (!std::isfinite(data[i])) {
      has_non_finite = true;
      break;
    }
  }

  if (!has_non_finite) {
    double accum = data[0];
    double resid = 0.0;

    for (uint32_t i = 1; i < window_size; i++) {
      double m = data[i];
      double p = accum;
      accum = accum + m;
      double q = accum - p;
      resid = resid + ((p - (accum - q)) + (m - q));
    }

    if (resid > 0.001) {
      Rf_warning("Residual value is large. Some precision may be lost. res = %f\n", resid);
    }

    res[0] = accum + resid;

    for (int64_t i = window_size; i < data_size; i++) {
      double m = data[i - window_size];
      double n = data[i];
      double p = accum - m;
      double q = p - accum;
      double r = resid + ((accum - (p - q)) - (m + q));
      accum = p + n;
      double t = accum - p;
      resid = r + ((p - (accum - t)) + (n - t));
      res[i - window_size + 1] = accum + resid;
    }
    return res;
  }

  if (na_mode == 2) {
    Rcpp::stop("Data contains missing/non-finite values.");
  }

  NumericVector clean_data(data_size);
  const double *data_ptr = data.begin();
  double *clean_ptr = clean_data.begin();
  for (uint32_t i = 0; i < data_size; ++i) {
    clean_ptr[i] = std::isfinite(data_ptr[i]) ? data_ptr[i] : 0.0;
  }

  double accum = clean_ptr[0];
  double resid = 0.0;

  for (uint32_t i = 1; i < window_size; i++) {
    double m = clean_ptr[i];
    double p = accum;
    accum = accum + m;
    double q = accum - p;
    resid = resid + ((p - (accum - q)) + (m - q));
  }

  res[0] = accum + resid;

  for (int64_t i = window_size; i < data_size; i++) {
    double m = clean_ptr[i - window_size];
    double n = clean_ptr[i];
    double p = accum - m;
    double q = p - accum;
    double r = resid + ((accum - (p - q)) - (m + q));
    accum = p + n;
    double t = accum - p;
    resid = r + ((p - (accum - t)) + (n - t));
    res[i - window_size + 1] = accum + resid;
  }

  SlidingNonFiniteTracker tracker(data_ptr, data_size, window_size);
  for (uint32_t i = 0; i < data_size; ++i) {
    tracker.update(i);
    if (i >= window_size - 1) {
      uint32_t out_idx = i - window_size + 1;
      if (na_mode == 0) { // window_na
        if (tracker.has_na()) {
          res[out_idx] = NA_REAL;
        }
      } else { // remove
        if (tracker.valid_count() == 0) {
          res[out_idx] = NA_REAL;
        }
      }
    }
  }

  return res;
}

//[[Rcpp::export]]
NumericVector movvar2_rcpp(const NumericVector data, uint32_t window_size, int na_mode) {
  uint32_t const data_size = data.length();
  NumericVector out(data_size - window_size + 1);

  bool has_non_finite = false;
  for (uint32_t i = 0; i < data_size; ++i) {
    if (!std::isfinite(data[i])) {
      has_non_finite = true;
      break;
    }
  }

  if (!has_non_finite) {
    double n = 0.0;
    double data_sum = 0.0;
    double data2_sum = 0.0;

    for (uint32_t i = 0; i < data_size; i++) {
      data_sum = data_sum + data[i];
      data2_sum = data2_sum + (data[i] * data[i]);
      n = n + 1.0;
      if (i >= window_size) {
        data_sum = data_sum - data[i - window_size];
        data2_sum = data2_sum - (data[i - window_size] * data[i - window_size]);
        n = n - 1.0;
      }

      if (i >= (window_size - 1)) {
        double var = data2_sum / n - ((data_sum * data_sum) / (n * n));
        out[i - (window_size - 1)] = (var > 0.0) ? var : 0.0;
      }
    }
    return out;
  }

  if (na_mode == 2) {
    Rcpp::stop("Data contains missing/non-finite values.");
  }

  const double *data_ptr = data.begin();
  SlidingNonFiniteTracker tracker(data_ptr, data_size, window_size);
  double data_sum = 0.0;
  double data2_sum = 0.0;

  for (uint32_t i = 0; i < data_size; i++) {
    if (std::isfinite(data_ptr[i])) {
      data_sum += data_ptr[i];
      data2_sum += data_ptr[i] * data_ptr[i];
    }
    if (i >= window_size && std::isfinite(data_ptr[i - window_size])) {
      data_sum -= data_ptr[i - window_size];
      data2_sum -= data_ptr[i - window_size] * data_ptr[i - window_size];
    }
    tracker.update(i);

    if (i >= (window_size - 1)) {
      uint32_t out_idx = i - (window_size - 1);
      if (na_mode == 0) { // window_na
        if (tracker.has_na()) {
          out[out_idx] = NA_REAL;
        } else {
          double var = data2_sum / window_size - ((data_sum * data_sum) / (static_cast<double>(window_size) * window_size));
          out[out_idx] = (var > 0.0) ? var : 0.0;
        }
      } else { // remove
        uint32_t v_n = tracker.valid_count();
        if (v_n <= 1) {
          out[out_idx] = NA_REAL;
        } else {
          double var = data2_sum / v_n - ((data_sum * data_sum) / (static_cast<double>(v_n) * v_n));
          out[out_idx] = (var > 0.0) ? var : 0.0;
        }
      }
    }
  }
  return out;
}

//[[Rcpp::export]]
NumericVector movvar_rcpp(const NumericVector data, const uint32_t window_size, int na_mode) {
  uint32_t const data_size = data.length();

  bool has_non_finite = false;
  for (uint32_t i = 0; i < data_size; ++i) {
    if (!std::isfinite(data[i])) {
      has_non_finite = true;
      break;
    }
  }

  if (!has_non_finite) {
    NumericVector mu = movsum_ogita_rcpp(data, window_size, 0) / window_size;
    NumericVector data2_sum = movsum_ogita_rcpp(data * data, window_size, 0);
    NumericVector data_var = (data2_sum / window_size) - (mu * mu);
    for (int i = 0; i < data_var.length(); ++i) {
      if (data_var[i] < 0.0) data_var[i] = 0.0;
    }
    return data_var;
  }

  if (na_mode == 2) {
    Rcpp::stop("Data contains missing/non-finite values.");
  }

  return movvar2_rcpp(data, window_size, na_mode);
}

//[[Rcpp::export]]
NumericVector movstd_rcpp(const NumericVector data, const uint32_t window_size, int na_mode) {
  NumericVector data_var = movvar_rcpp(data, window_size, na_mode);
  NumericVector data_sd(data_var.length());
  for (int i = 0; i < data_var.length(); ++i) {
    if (ISNA(data_var[i]) || ISNAN(data_var[i])) {
      data_sd[i] = NA_REAL;
    } else {
      data_sd[i] = std::sqrt(data_var[i]);
    }
  }
  return data_sd;
}

//[[Rcpp::export]]
List movmean_std_rcpp(const NumericVector data, const uint32_t window_size, int na_mode) {
  uint32_t const data_size = data.length();
  uint32_t const profile_len = data_size - window_size + 1;

  bool has_non_finite = false;
  for (uint32_t i = 0; i < data_size; ++i) {
    if (!std::isfinite(data[i])) {
      has_non_finite = true;
      break;
    }
  }

  if (!has_non_finite) {
    NumericVector data_sum = movsum_ogita_rcpp(data, window_size, 0);
    NumericVector data_mean = data_sum / window_size;
    NumericVector data2 = (data * data);
    NumericVector data2_sum = movsum_ogita_rcpp(data2, window_size, 0);

    NumericVector data_var = (data2_sum / window_size) - (data_mean * data_mean);
    NumericVector data_sd(profile_len);
    NumericVector data_sig(profile_len);
    for (uint32_t i = 0; i < profile_len; ++i) {
      double v = (data_var[i] > 0.0) ? data_var[i] : 0.0;
      data_sd[i] = std::sqrt(v);
      data_sig[i] = (v > DBL_EPSILON) ? (1.0 / std::sqrt(v * window_size)) : R_PosInf;
    }

    return List::create(
      Rcpp::Named("avg") = data_mean, Rcpp::Named("sd") = data_sd, Rcpp::Named("sig") = data_sig,
      Rcpp::Named("sum") = data_sum, Rcpp::Named("sqrsum") = data2_sum
    );
  }

  if (na_mode == 2) {
    Rcpp::stop("Data contains missing/non-finite values.");
  }

  NumericVector data_mean(profile_len);
  NumericVector data_sd(profile_len);
  NumericVector data_sig(profile_len);
  NumericVector data_sum(profile_len);
  NumericVector data2_sum(profile_len);

  const double *data_ptr = data.begin();
  SlidingNonFiniteTracker tracker(data_ptr, data_size, window_size);
  double sum1 = 0.0;
  double sum2 = 0.0;

  for (uint32_t i = 0; i < data_size; ++i) {
    if (std::isfinite(data_ptr[i])) {
      sum1 += data_ptr[i];
      sum2 += data_ptr[i] * data_ptr[i];
    }
    if (i >= window_size && std::isfinite(data_ptr[i - window_size])) {
      sum1 -= data_ptr[i - window_size];
      sum2 -= data_ptr[i - window_size] * data_ptr[i - window_size];
    }
    tracker.update(i);

    if (i >= (window_size - 1)) {
      uint32_t out_idx = i - (window_size - 1);
      if (na_mode == 0) { // window_na
        if (tracker.has_na()) {
          data_mean[out_idx] = NA_REAL;
          data_sd[out_idx] = NA_REAL;
          data_sig[out_idx] = NA_REAL;
          data_sum[out_idx] = NA_REAL;
          data2_sum[out_idx] = NA_REAL;
        } else {
          double m = sum1 / window_size;
          double v = sum2 / window_size - m * m;
          if (v < 0.0) v = 0.0;
          double s = std::sqrt(v);
          data_mean[out_idx] = m;
          data_sd[out_idx] = s;
          data_sig[out_idx] = (v > DBL_EPSILON) ? (1.0 / std::sqrt(v * window_size)) : R_PosInf;
          data_sum[out_idx] = sum1;
          data2_sum[out_idx] = sum2;
        }
      } else { // remove
        uint32_t v_n = tracker.valid_count();
        if (v_n <= 1) {
          data_mean[out_idx] = (v_n == 1) ? sum1 : NA_REAL;
          data_sd[out_idx] = NA_REAL;
          data_sig[out_idx] = NA_REAL;
          data_sum[out_idx] = (v_n == 1) ? sum1 : NA_REAL;
          data2_sum[out_idx] = (v_n == 1) ? sum2 : NA_REAL;
        } else {
          double m = sum1 / v_n;
          double v = sum2 / v_n - m * m;
          if (v < 0.0) v = 0.0;
          double s = std::sqrt(v);
          data_mean[out_idx] = m;
          data_sd[out_idx] = s;
          data_sig[out_idx] = (v > DBL_EPSILON) ? (1.0 / std::sqrt(v * v_n)) : R_PosInf;
          data_sum[out_idx] = sum1;
          data2_sum[out_idx] = sum2;
        }
      }
    }
  }

  return List::create(
    Rcpp::Named("avg") = data_mean, Rcpp::Named("sd") = data_sd, Rcpp::Named("sig") = data_sig,
    Rcpp::Named("sum") = data_sum, Rcpp::Named("sqrsum") = data2_sum
  );
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

inline double window_min_finite(const double *ptr, uint32_t w) {
  double m = R_PosInf;
  for (uint32_t j = 0; j < w; ++j) {
    if (std::isfinite(ptr[j]) && ptr[j] < m) {
      m = ptr[j];
    }
  }
  return (m == R_PosInf) ? NA_REAL : m;
}

inline double window_max_finite(const double *ptr, uint32_t w) {
  double m = R_NegInf;
  for (uint32_t j = 0; j < w; ++j) {
    if (std::isfinite(ptr[j]) && ptr[j] > m) {
      m = ptr[j];
    }
  }
  return (m == R_NegInf) ? NA_REAL : m;
}

//[[Rcpp::export]]
NumericVector movmin_rcpp(const NumericVector data, uint32_t window_size, int na_mode) {
  uint32_t const data_size = data.length();

  if (window_size <= 1) {
    return data;
  }

  if (window_size > data_size) {
    window_size = data_size;
  }

  bool has_non_finite = false;
  for (uint32_t i = 0; i < data_size; ++i) {
    if (!std::isfinite(data[i])) {
      has_non_finite = true;
      break;
    }
  }

  if (!has_non_finite) {
    uint32_t i = 0UL, d = 0UL, k = 0UL;
    double min_res = 0.0, res_out = 0.0;
    NumericVector out(data_size - window_size + 1);

    k = 0;
    d = 0;

    min_res = res_out = R_PosInf;
    for (i = window_size - 1; i < data_size; i++) {
      if (res_out == min_res) {
        min_res = *std::min_element(data.begin() + d, data.begin() + d + window_size);
      } else {
        min_res = MIN(min_res, data[d + window_size - 1]);
      }

      res_out = data[d++];
      out[k++] = (min_res == R_PosInf ? R_NaReal : min_res);
    }

    return out;
  }

  if (na_mode == 2) {
    Rcpp::stop("Data contains missing/non-finite values.");
  }

  const double *data_ptr = data.begin();
  SlidingNonFiniteTracker tracker(data_ptr, data_size, window_size);
  NumericVector out(data_size - window_size + 1);
  double min_res = R_PosInf, res_out = R_PosInf;
  uint32_t d = 0;

  for (uint32_t i = 0; i < data_size; ++i) {
    tracker.update(i);
    if (i >= window_size - 1) {
      uint32_t k = i - window_size + 1;
      if (tracker.has_na()) {
        if (na_mode == 0) { // window_na
          out[k] = NA_REAL;
        } else { // remove
          out[k] = (tracker.valid_count() == 0) ? NA_REAL : window_min_finite(data_ptr + d, window_size);
        }
        min_res = res_out = R_PosInf;
        d++;
      } else {
        if (res_out == min_res) {
          min_res = *std::min_element(data_ptr + d, data_ptr + d + window_size);
        } else {
          min_res = MIN(min_res, data_ptr[d + window_size - 1]);
        }
        res_out = data_ptr[d++];
        out[k] = (min_res == R_PosInf ? NA_REAL : min_res);
      }
    }
  }

  return out;
}

//[[Rcpp::export]]
NumericVector movmax_rcpp(const NumericVector data, uint32_t window_size, int na_mode) {
  uint32_t const data_size = data.length();

  if (window_size <= 1) {
    return data;
  }

  if (window_size > data_size) {
    window_size = data_size;
  }

  bool has_non_finite = false;
  for (uint32_t i = 0; i < data_size; ++i) {
    if (!std::isfinite(data[i])) {
      has_non_finite = true;
      break;
    }
  }

  if (!has_non_finite) {
    uint32_t i = 0UL, d = 0UL, k = 0UL;
    double max_res = 0.0, res_out = 0.0;
    NumericVector out(data_size - window_size + 1);

    k = 0;
    d = 0;

    max_res = res_out = R_NegInf;
    for (i = window_size - 1; i < data_size; i++) {
      if (res_out == max_res) {
        max_res = *std::max_element(data.begin() + d, data.begin() + d + window_size);
      } else {
        max_res = MAX(max_res, data[d + window_size - 1]);
      }

      res_out = data[d++];
      out[k++] = (max_res == R_NegInf ? R_NaReal : max_res);
    }

    return out;
  }

  if (na_mode == 2) {
    Rcpp::stop("Data contains missing/non-finite values.");
  }

  const double *data_ptr = data.begin();
  SlidingNonFiniteTracker tracker(data_ptr, data_size, window_size);
  NumericVector out(data_size - window_size + 1);
  double max_res = R_NegInf, res_out = R_NegInf;
  uint32_t d = 0;

  for (uint32_t i = 0; i < data_size; ++i) {
    tracker.update(i);
    if (i >= window_size - 1) {
      uint32_t k = i - window_size + 1;
      if (tracker.has_na()) {
        if (na_mode == 0) { // window_na
          out[k] = NA_REAL;
        } else { // remove
          out[k] = (tracker.valid_count() == 0) ? NA_REAL : window_max_finite(data_ptr + d, window_size);
        }
        max_res = res_out = R_NegInf;
        d++;
      } else {
        if (res_out == max_res) {
          max_res = *std::max_element(data_ptr + d, data_ptr + d + window_size);
        } else {
          max_res = MAX(max_res, data_ptr[d + window_size - 1]);
        }
        res_out = data_ptr[d++];
        out[k] = (max_res == R_NegInf ? NA_REAL : max_res);
      }
    }
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


inline bool inverse_centered_norm(const double *data, uint32_t start, uint32_t window_size, double mean,
                                  double &inverse_norm) {
  if (!std::isfinite(mean)) {
    return false;
  }

  double scale = 0.0;
  double sum_squares = 1.0;
  for (uint32_t j = 0; j < window_size; j++) {
    double const centered = std::abs(data[start + j] - mean);
    if (centered == 0.0) {
      continue;
    }
    if (scale < centered) {
      double const ratio = scale / centered;
      sum_squares = 1.0 + sum_squares * ratio * ratio;
      scale = centered;
    } else {
      double const ratio = centered / scale;
      sum_squares += ratio * ratio;
    }
  }

  double const centered_norm = scale * sqrt(sum_squares);
  inverse_norm = 1.0 / centered_norm;
  return std::isfinite(centered_norm) && std::isfinite(inverse_norm);
}

// Precondition data for MPX joins containing non-finite samples. This follows
// the MATLAB implementation: non-finite samples are zeroed only for the
// recurrence, the rolling mean is compensated, and each centered norm is
// computed directly instead of being derived from rolling sums of squares.
List muinvn_na_impl(NumericVector data_ref, uint32_t window_size, bool parallel, uint32_t n_workers) {
  uint32_t const data_size = data_ref.length();
  uint32_t const profile_len = data_size - window_size + 1;
  NumericVector data = clone(data_ref);
  LogicalVector valid_window(profile_len, true);

  double const *const data_ref_ptr = data_ref.begin();
  double *const data_ptr = data.begin();
  int *const valid_window_ptr = valid_window.begin();

  uint32_t non_finite_count = 0;
  for (uint32_t i = 0; i < data_size; i++) {
    if (!std::isfinite(data_ptr[i])) {
      non_finite_count++;
      data_ptr[i] = 0.0;
    }
    if (i >= window_size && !std::isfinite(data_ref_ptr[i - window_size])) {
      non_finite_count--;
    }
    if (i >= window_size - 1) {
      valid_window_ptr[i - window_size + 1] = non_finite_count == 0;
    }
  }

  NumericVector mu = movsum_ogita_rcpp(data, window_size, 0) / window_size;
  NumericVector sig(profile_len, R_NaN);
  double const *const mu_ptr = mu.begin();
  double *const sig_ptr = sig.begin();

  uint64_t const norm_work = static_cast<uint64_t>(profile_len) * window_size;
  if (parallel && norm_work >= 1000000ULL) {
    uint32_t actual_workers = std::max<uint32_t>(1, n_workers);
    RcppThread::parallelFor(0, profile_len, [&](std::size_t i) {
      if (!valid_window_ptr[i]) {
        return;
      }

      double inverse_norm = std::numeric_limits<double>::quiet_NaN();
      if (!inverse_centered_norm(data_ptr, i, window_size, mu_ptr[i], inverse_norm)) {
        valid_window_ptr[i] = false;
        return;
      }
      sig_ptr[i] = inverse_norm;
    }, actual_workers);
  } else {
    for (uint32_t i = 0; i < profile_len; i++) {
      if (!valid_window_ptr[i]) {
        continue;
      }

      double inverse_norm = std::numeric_limits<double>::quiet_NaN();
      if (!inverse_centered_norm(data_ptr, i, window_size, mu_ptr[i], inverse_norm)) {
        valid_window_ptr[i] = false;
        continue;
      }
      sig_ptr[i] = inverse_norm;
    }
  }

  return List::create(Rcpp::Named("data") = data, Rcpp::Named("avg") = mu, Rcpp::Named("sig") = sig,
                      Rcpp::Named("valid_window") = valid_window);
}

List muinvn_na(NumericVector data_ref, uint32_t window_size) { return muinvn_na_impl(data_ref, window_size, false, 1); }

List muinvn_na_parallel(NumericVector data_ref, uint32_t window_size, uint32_t n_workers) {
  return muinvn_na_impl(data_ref, window_size, true, n_workers);
}

// A finite series can still contain constant windows, which the NA-aware
// implementation deliberately masks. Also avoid the ordinary MPX statistics
// path when the series has a large offset relative to its dynamic range: its
// rolling sum-of-squares calculation is then susceptible to cancellation.
// This O(n) guard keeps the finite fast path both fast and numerically safe.
bool mpx_finite_fast_path_safe(const NumericVector &data, uint32_t window_size) {
  uint32_t const data_size = data.length();
  if (window_size < 2 || window_size > data_size) {
    return false;
  }

  double const *const data_ptr = data.begin();
  double min_value = R_PosInf;
  double max_value = R_NegInf;
  double max_abs = 0.0;
  uint32_t equal_run = 1;
  for (uint32_t i = 0; i < data_size; i++) {
    if (!std::isfinite(data_ptr[i])) {
      return false;
    }
    min_value = std::min(min_value, data_ptr[i]);
    max_value = std::max(max_value, data_ptr[i]);
    max_abs = std::max(max_abs, std::abs(data_ptr[i]));
    if (i > 0 && data_ptr[i] == data_ptr[i - 1]) {
      equal_run++;
      if (equal_run >= window_size) {
        return false;
      }
    } else {
      equal_run = 1;
    }
  }

  double const dynamic_range = max_value - min_value;
  if (!std::isfinite(dynamic_range)) {
    return false;
  }
  double const scale = std::max(1.0, dynamic_range);
  constexpr double max_offset_to_range = 1000.0;
  return (max_abs / scale) <= max_offset_to_range;
}

//[[Rcpp::export]]
List muinvn_rcpp(const NumericVector data, uint32_t window_size, int na_mode) {
  uint32_t const data_size = data.length();
  if (window_size > data_size) {
    Rcpp::stop("'window_size' cannot be larger than data length.");
  }

  bool has_non_finite = false;
  for (uint32_t i = 0; i < data_size; ++i) {
    if (!std::isfinite(data[i])) {
      has_non_finite = true;
      break;
    }
  }

  if (!has_non_finite) {
    NumericVector mu = movsum_ogita_rcpp(data, window_size, 0) / window_size;
    NumericVector data2_sum = movsum_ogita_rcpp(data * data, window_size, 0);
    NumericVector sig = 1.0 / sqrt(data2_sum - mu * mu * window_size);

    return List::create(Rcpp::Named("avg") = mu, Rcpp::Named("sig") = sig);
  }

  if (na_mode == 2) {
    Rcpp::stop("Data contains missing/non-finite values.");
  }

  if (na_mode == 0) {
    List na_res = muinvn_na_impl(data, window_size, false, 1);
    NumericVector mu = na_res["avg"];
    NumericVector sig = na_res["sig"];
    LogicalVector valid = na_res["valid_window"];
    for (int i = 0; i < mu.length(); ++i) {
      if (!valid[i]) {
        mu[i] = NA_REAL;
        sig[i] = NA_REAL;
      }
    }
    return List::create(Rcpp::Named("avg") = mu, Rcpp::Named("sig") = sig);
  }

  // na_mode == 1: remove
  List stats = movmean_std_rcpp(data, window_size, 1);
  return List::create(Rcpp::Named("avg") = stats["avg"], Rcpp::Named("sig") = stats["sig"]);
}

//[[Rcpp::export]]
List muinvn_rcpp_parallel(const NumericVector data, uint32_t window_size, uint32_t n_workers, int na_mode) {
  uint32_t const data_size = data.length();
  if (window_size > data_size) {
    Rcpp::stop("'window_size' cannot be larger than data length.");
  }

  bool has_non_finite = false;
  for (uint32_t i = 0; i < data_size; ++i) {
    if (!std::isfinite(data[i])) {
      has_non_finite = true;
      break;
    }
  }

  if (!has_non_finite) {
    NumericVector sig(data_size - window_size + 1);
    NumericVector const mu = movsum_ogita_rcpp(data, window_size, 0) / window_size;
    NumericVector const data2_sum = movsum_ogita_rcpp(data * data, window_size, 0);

    double *const sig_ptr = sig.begin();
    double const *const mu_ptr = mu.begin();
    double const *const d2_ptr = data2_sum.begin();
    uint32_t const len = mu.length();

    uint32_t actual_workers = std::max<uint32_t>(1, n_workers);

    try {
      RcppThread::parallelFor(0, len, [&](std::size_t i) {
        sig_ptr[i] = 1.0 / std::sqrt(d2_ptr[i] - mu_ptr[i] * mu_ptr[i] * window_size);
      }, actual_workers);
    } catch (RcppThread::UserInterruptException &ex) {
      Rcout << "Process terminated.\n";
    } catch (...) {
      Rcpp::stop("c++ exception (unknown reason)");
    }

    return List::create(Rcpp::Named("avg") = mu, Rcpp::Named("sig") = sig);
  }

  if (na_mode == 2) {
    Rcpp::stop("Data contains missing/non-finite values.");
  }

  if (na_mode == 0) {
    List na_res = muinvn_na_impl(data, window_size, true, n_workers);
    NumericVector mu = na_res["avg"];
    NumericVector sig = na_res["sig"];
    LogicalVector valid = na_res["valid_window"];
    for (int i = 0; i < mu.length(); ++i) {
      if (!valid[i]) {
        mu[i] = NA_REAL;
        sig[i] = NA_REAL;
      }
    }
    return List::create(Rcpp::Named("avg") = mu, Rcpp::Named("sig") = sig);
  }

  // na_mode == 1: remove
  List stats = movmean_std_rcpp(data, window_size, 1);
  return List::create(Rcpp::Named("avg") = stats["avg"], Rcpp::Named("sig") = stats["sig"]);
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

  uint32_t profile_size = data.size() - window_size + 1;
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
