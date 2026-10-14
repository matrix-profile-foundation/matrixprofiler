// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#include "mathtools.h"
#include "mass.h"
#include "fft.h"
#include "windowfunc.h"
#include "mpx.h"
#include "mpx_parallel.h"
#include "mpdist.h"
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>
#include <complex>

using namespace Rcpp;

inline void sliding_min(const double* in, size_t n, size_t w, double* out) {
  if (w <= 1) {
    std::copy(in, in + n, out);
    return;
  }
  std::deque<size_t> dq;
  for (size_t i = 0; i < n; ++i) {
    while (!dq.empty() && dq.front() + w <= i) {
      dq.pop_front();
    }
    while (!dq.empty() && in[dq.back()] >= in[i]) {
      dq.pop_back();
    }
    dq.push_back(i);
    if (i >= w - 1) {
      out[i - (w - 1)] = in[dq.front()];
    }
  }
}

inline double cal_mp_dist_cpp(std::vector<double>& mp, double thr, size_t data_size) {
  if (mp.empty()) {
    return 0.0;
  }
  size_t k = static_cast<size_t>(std::ceil(thr * static_cast<double>(data_size)));
  if (k == 0) {
    k = 1;
  }
  if (k > mp.size()) {
    return *std::max_element(mp.begin(), mp.end());
  }
  std::nth_element(mp.begin(), mp.begin() + (k - 1), mp.end());
  return mp[k - 1];
}

// [[Rcpp::export]]
double mpdist_simple_rcpp(NumericVector ref_data, NumericVector query_data, uint32_t window_size, double thr = 0.05,
                          uint32_t n_workers = 1) {
  uint64_t const ref_size = ref_data.length();
  uint64_t const query_size = query_data.length();

  if (query_size == window_size) {
    List pre = mass_pre_rcpp(ref_data, query_data, window_size);
    NumericVector query_mean = pre["query_mean"];
    NumericVector query_sd = pre["query_sd"];
    List res = mass2_rcpp(pre["data_fft"], query_data, ref_size, window_size, pre["data_mean"], pre["data_sd"],
                          query_mean[0], query_sd[0], n_workers);
    NumericVector dp = res["distance_profile"];
    std::vector<double> dp_vec(dp.begin(), dp.end());
    double dist = cal_mp_dist_cpp(dp_vec, thr, ref_size);
    return std::sqrt(dist);
  }

  List ab_mp;
  if (n_workers > 1) {
    ab_mp = mpxab_rcpp_parallel(ref_data, query_data, window_size, 1.0, false, true, false, n_workers);
  } else {
    ab_mp = mpxab_rcpp(ref_data, query_data, window_size, 1.0, false, true, false);
  }

  NumericVector mpb = ab_mp["mpb"];
  std::vector<double> mp_vec(mpb.begin(), mpb.end());
  double dist = cal_mp_dist_cpp(mp_vec, thr, ref_size + query_size);
  return dist;
}

// [[Rcpp::export]]
NumericVector mpdist_vect_rcpp(NumericVector data, NumericVector query, uint32_t window_size, double thr = 0.05,
                              uint32_t n_workers = 1) {
  uint64_t const data_size = data.length();
  uint64_t const query_size = query.length();
  uint32_t const num_subseqs = query_size - window_size + 1;
  uint32_t const dist_profile_size = data_size - window_size + 1;
  uint32_t const mp_dist_length = data_size - query_size + 1;

  List data_avgsd = movmean_std_rcpp(data, window_size);
  NumericVector data_mean = data_avgsd["avg"];
  NumericVector data_sd = data_avgsd["sd"];

  List query_avgsd = movmean_std_rcpp(query, window_size);
  NumericVector query_mean = query_avgsd["avg"];
  NumericVector query_sd = query_avgsd["sd"];

  uint32_t const pad_size = static_cast<uint32_t>(pow(2, ceil(log2(static_cast<double>(data_size)))));
  std::vector<double> data_padded(pad_size, 0.0);
  std::copy(data.begin(), data.end(), data_padded.begin());

  FFT::fftw fft_engine;
  std::vector<std::complex<double>> data_fft = fft_engine.fft(data_padded, false);

  std::vector<double> mat(static_cast<size_t>(num_subseqs) * dist_profile_size);

  RcppThread::parallelFor(0, num_subseqs, [&](size_t i) {
    FFT::fftw local_fft;
    std::vector<double> rev_query(pad_size, 0.0);
    for (uint32_t j = 0; j < window_size; ++j) {
      rev_query[j] = query[i + window_size - 1 - j];
    }
    std::vector<std::complex<double>> q_fft = local_fft.fft(rev_query, false);
    std::vector<std::complex<double>> prod(pad_size);
    for (size_t p = 0; p < pad_size; ++p) {
      prod[p] = data_fft[p] * q_fft[p];
    }
    std::vector<std::complex<double>> z = local_fft.fft(prod, true);

    double q_m = query_mean[i];
    double q_s = query_sd[i];

    for (size_t j = 0; j < dist_profile_size; ++j) {
      double d_m = data_mean[j];
      double d_s = data_sd[j];
      double val;
      if (d_s <= 0.0 || q_s <= 0.0 || std::isnan(d_s) || std::isnan(q_s)) {
        val = R_PosInf;
      } else {
        double dot = z[window_size - 1 + j].real();
        val = 2.0 * (static_cast<double>(window_size) - (dot - static_cast<double>(window_size) * d_m * q_m) / (d_s * q_s));
        if (std::isnan(val)) {
          val = R_PosInf;
        } else if (val < 0.0) {
          val = 0.0;
        }
      }
      mat[i * dist_profile_size + j] = val;
    }
  }, n_workers);

  // Column minimums: all_right_histogram
  std::vector<double> all_right_histogram(dist_profile_size, R_PosInf);
  for (size_t i = 0; i < num_subseqs; ++i) {
    for (size_t j = 0; j < dist_profile_size; ++j) {
      double v = mat[i * dist_profile_size + j];
      if (v < all_right_histogram[j]) {
        all_right_histogram[j] = v;
      }
    }
  }

  // Sliding min for each row
  size_t slid_min_len = dist_profile_size - num_subseqs + 1;
  std::vector<double> mass_dist_slid_min(static_cast<size_t>(num_subseqs) * slid_min_len);

  RcppThread::parallelFor(0, num_subseqs, [&](size_t i) {
    sliding_min(&mat[i * dist_profile_size], dist_profile_size, num_subseqs, &mass_dist_slid_min[i * slid_min_len]);
  }, n_workers);

  // Final MPdist calculation
  NumericVector mpdist_array(mp_dist_length);
  size_t rec_len = 2 * num_subseqs;

  RcppThread::parallelFor(0, mp_dist_length, [&](size_t k) {
    std::vector<double> recreated_mp(rec_len);
    for (size_t i = 0; i < num_subseqs; ++i) {
      recreated_mp[i] = mass_dist_slid_min[i * slid_min_len + k];
    }
    for (size_t i = 0; i < num_subseqs; ++i) {
      recreated_mp[num_subseqs + i] = all_right_histogram[k + i];
    }

    double dist_val = cal_mp_dist_cpp(recreated_mp, thr, 2 * query_size);
    if (dist_val < 1e-14) {
      dist_val = 0.0;
    }
    mpdist_array[k] = std::sqrt(dist_val);
  }, n_workers);

  return mpdist_array;
}
