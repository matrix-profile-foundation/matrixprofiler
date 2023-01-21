// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#include "mathtools.h"
#include "fft.h"
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

IntegerVector seq(uint64_t start, uint64_t end) {

  if (start > end) {
    return rev(Range(end, start));
  }

  return Range(start, end);
}

IntegerVector seq_by(uint64_t start, uint64_t end, uint32_t by) {

  uint64_t const total_length = ceil((double)end / by);

  IntegerVector const result = Range(1, total_length) * by + start - by;

  return (result);
}

//[[Rcpp::export]]
double std_rcpp(const NumericVector data, const bool na_rm = false) {

  NumericVector the_data = data;

  // if there is NaN in vector the result will be NaN
  if (any(is_na(data))) {
    if (na_rm) {
      the_data = na_omit(data);
    } else {
      return NA_REAL;
    }
  }

  double const result = sqrt(sum(((the_data - mean(the_data)) * (the_data - mean(the_data)))) / the_data.length());

  return (result);
}

//[[Rcpp::export]]
NumericMatrix list_to_matrix(const List x) {
  int32_t const nlines = x.size();
  uint32_t colmax = 0;

  for (int32_t i = 0; i < nlines; i++) {
    uint32_t const currsize = as<NumericVector>(x[i]).size();
    if (colmax < currsize) {
      colmax = currsize;
    }
  }

  NumericMatrix m(nlines, colmax);

  for (int32_t i = 0; i < nlines; i++) {
    // int32_t line = nlines - i - 1;
    uint32_t const currsize = as<NumericVector>(x[i]).size();
    NumericMatrix::Row row = m(i, _);
    row = as<NumericVector>(x[i]); //-V519

    for (uint32_t j = currsize; j < colmax; j++) {
      row[j] = 0;
    }
  }

  return (m);
}

IntegerVector which_cpp(const LogicalVector x) {

  int const nx = x.size();
  std::vector<int> y;
  y.reserve(nx);

  for (int i = 0; i < nx; i++) {
    if (x[i] == true) {
      y.push_back(i);
    }
  }

  return wrap(y);
}

//[[Rcpp::export]]
int32_t mode_rcpp(const IntegerVector x) {

  // is slower than R implementation...
  IntegerVector ux = unique(x);
  int32_t const y = ux[which_max(table(match(x, ux)))];
  return y;
}

// Normalizes data for mean Zero and Standard Deviation One
//
// @inheritParams std
//
// @return Returns the normalized data
// @keywords internal
// @noRd
//

//[[Rcpp::export]]
NumericVector znorm_rcpp(const NumericVector data) {
  double const data_mean = mean(data);
  double const data_dev = sqrt(sum(((data - data_mean) * (data - data_mean))) / data.length());

  if (data_dev == NA_REAL || data_dev <= 0.01) {
    return (data - data_mean);
  } else {
    return ((data - data_mean) / (data_dev));
  }
}

//[[Rcpp::export]]
NumericVector normalize_rcpp(const NumericVector data, double min, double max) {
  double const min_val = ::min(data);
  double const max_val = ::max(data);

  double const a = (max - min) / (max_val - min_val);
  double const b = max - a * max_val;
  NumericVector norm_data = a * data + b;

  norm_data[norm_data < min] = min;
  norm_data[norm_data > max] = max; //-V519

  return (norm_data);
}

//[[Rcpp::export]]
IntegerVector binary_split_rcpp(const uint32_t n) {

  IntegerVector idxs(n);

  idxs[0] = 1; // We always begin by explore the first integer
  // After exploring the first integer, we begin splitting the interval 2:n

  std::deque<uint32_t> lb_list;
  std::deque<uint32_t> ub_list;

  lb_list.push_back(2);
  ub_list.push_back(n);

  uint32_t lb = 0;
  uint32_t ub = 0;
  uint32_t mid = 0;

  for (uint32_t i = 1; i < n; i++) {
    lb = lb_list.front();
    ub = ub_list.front();
    mid = (lb + ub) / 2; // integer division is automatically floor()
    lb_list.pop_front();
    ub_list.pop_front();

    idxs[i] = mid;

    if (lb == ub) {
      continue;
    } else {
      if (lb < mid) {
        lb_list.push_back(lb);
        ub_list.push_back(mid - 1);
      }

      if (ub > mid) {
        lb_list.push_back(mid + 1);
        ub_list.push_back(ub);
      }
    }
  }

  return (idxs);
}

//[[Rcpp::export]]
NumericVector ed_corr_rcpp(const NumericVector data, uint32_t window_size) {
  // this format is less prone to rounding problems
  NumericVector const res = (2 * window_size - data * data) / (2 * window_size);

  return (res);
}

//[[Rcpp::export]]
NumericVector corr_ed_rcpp(const NumericVector data, uint32_t window_size) {
  // this format is less prone to rounding problems
  NumericVector const res = sqrt(2 * window_size * (1 - ifelse(data > 1, 1, data)));

  return (res);
}

//[[Rcpp::export]]
double inner_product(const NumericVector a, const NumericVector b) {
  double const res = std::inner_product(a.begin(), a.end(), b.begin(), 0.0);

  return (res);
}

//[[Rcpp::export]]
double sum_of_squares(const NumericVector a) {
  double const res = std::inner_product(a.begin(), a.end(), a.begin(), 0.0);

  return (res);
}

// Version compatible with RCPP code
//[[Rcpp::export]]
ComplexVector fft_rcpp(const ComplexVector z, bool invert) {

  ComplexVector result;
  int const n = z.length();
  std::vector<std::complex<double>> zz(n);
  FFT::fftw *fft = new FFT::fftw();

  for (int i = 0; i < n; i++) {
    zz[i].real(z[i].r);
    zz[i].imag(z[i].i);
  }

  result = wrap(fft->fft(zz, invert));
  delete (fft);

  return result;
}

ComplexVector fft_rcpp(const NumericVector z, bool invert) {

  ComplexVector result;
  int const n = z.length();
  std::vector<std::complex<double>> zz(n);
  FFT::fftw *fft = new FFT::fftw();

  for (int i = 0; i < n; i++) {
    zz[i].real(z[i]);
    zz[i].imag(0.0);
  }

  result = wrap(fft->fft(zz, invert));
  delete (fft);

  return result;
}

std::vector<std::complex<double>> fft_rcpp(const std::vector<double> z, bool invert) {

  std::vector<std::complex<double>> result;
  int const n = z.size();
  std::vector<std::complex<double>> zz(n);
  FFT::fftw *fft = new FFT::fftw();

  for (int i = 0; i < n; i++) {
    zz[i].real(z[i]);
    zz[i].imag(0);
  }

  result = fft->fft(zz, invert);
  delete (fft);

  return result;
}

std::vector<std::complex<double>> fft_rcpp(const std::vector<std::complex<double>> z, bool invert) {

  std::vector<std::complex<double>> result;
  int const n = z.size();
  std::vector<std::complex<double>> zz(n);
  FFT::fftw *fft = new FFT::fftw();

  for (int i = 0; i < n; i++) {
    zz[i].real(z[i].real());
    zz[i].imag(z[i].imag());
  }

  result = fft->fft(zz, invert);
  delete (fft);

  return result;
}

std::vector<double> fft_rcpp_real(const std::vector<std::complex<double>> z, bool invert) {

  int const n = z.size();
  std::vector<double> result(n);
  std::vector<std::complex<double>> zz(n);
  std::vector<std::complex<double>> result_cplx;
  FFT::fftw *fft = new FFT::fftw();

  for (int i = 0; i < n; i++) {
    zz[i].real(z[i].real());
    zz[i].imag(z[i].imag());
  }

  result_cplx = fft->fft(zz, invert);
  delete (fft);

  for (int i = 0; i < n; i++) {
    result[i] = result_cplx[i].real();
  }

  return result;
}
