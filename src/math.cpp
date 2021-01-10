#include "math.h"
#include "fft.h"
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

IntegerVector seq(uint64_t start, uint64_t end) {

  if (start > end) {
    return rev(Range(end, start));
  }

  return Range(start, end);
}

//[[Rcpp::export]]
IntegerVector seq_by(uint64_t start, uint64_t end, uint32_t by) {

  uint64_t total_length = ceil((double)end / by);

  IntegerVector result = Range(1, total_length) * by + start - by;

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

  double result = sqrt(sum(pow((the_data - mean(the_data)), 2)) / the_data.length());

  return (result);
}

//[[Rcpp::export]]
NumericMatrix list_to_matrix(const List x) {
  int32_t nlines = x.size();
  uint32_t colmax = 0;

  for (int32_t i = 0; i < nlines; i++) {
    uint32_t currsize = as<NumericVector>(x[i]).size();
    if (colmax < currsize) {
      colmax = currsize;
    }
  }

  NumericMatrix m(nlines, colmax);

  for (int32_t i = 0; i < nlines; i++) {
    // int32_t line = nlines - i - 1;
    uint32_t currsize = as<NumericVector>(x[i]).size();
    NumericMatrix::Row row = m(i, _);
    row = as<NumericVector>(x[i]);

    for (uint32_t j = currsize; j < colmax; j++) {
      row[j] = 0;
    }
  }

  return (m);
}

// [[Rcpp::export]]
IntegerVector which(const LogicalVector x) {

  int nx = x.size();
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
NumericVector diff_lag(const NumericVector x, const uint32_t lag = 1) {
  uint32_t n = x.size();
  NumericVector out(n - lag);

  for (uint32_t i = 0; i < (n - lag); i++) {
    out[i] = x[i + lag] - x[i];
  }
  return out;
}

//[[Rcpp::export]]
NumericVector diff2_lag(const NumericVector x, const uint32_t lag = 1, const double v = 0.0) {
  uint32_t n = x.size();
  NumericVector out(n - lag + 1);

  out[0] = v;

  for (uint32_t i = 0; i < (n - lag); i++) {
    out[i + 1] = x[i + lag] - x[i];
  }
  return out;
}

//[[Rcpp::export]]
NumericVector fast_movsd_rcpp(const NumericVector data, const uint32_t window_size) {

  // Improve the numerical analysis by subtracting off the series mean
  // this has no effect on the standard deviation.
  NumericVector data_zeromean = data - mean(data);

  NumericVector data_sum =
      cumsum(diff2_lag(data_zeromean, window_size, sum(data_zeromean[Range(0, (window_size - 1))])));
  NumericVector data_mean = data_sum / window_size;

  NumericVector data2 = pow(data_zeromean, 2);
  NumericVector data2_sum = cumsum(diff2_lag(data2, window_size, sum(data2[Range(0, (window_size - 1))])));
  NumericVector data_sd2 = (data2_sum / window_size) - pow(data_mean, 2); // variance
  NumericVector data_sd = sqrt(data_sd2);

  return (data_sd);
}

//[[Rcpp::export]]
List fast_avg_sd_rcpp(const NumericVector data, const uint32_t window_size) {

  NumericVector mov_sum =
      cumsum(diff2_lag(data, window_size, sum(as<NumericVector>(data[Range(0, (window_size - 1))]))));
  NumericVector mov_mean = mov_sum / window_size;
  NumericVector data2 = pow(data, 2);
  NumericVector mov2_sum = cumsum(diff2_lag(data2, window_size, sum(data2[Range(0, (window_size - 1))])));

  // Improve the numerical analysis by subtracting off the series mean
  // this has no effect on the standard deviation.
  NumericVector data_zeromean = data - mean(data);

  NumericVector data_sum =
      cumsum(diff2_lag(data_zeromean, window_size, sum(data_zeromean[Range(0, (window_size - 1))])));
  NumericVector data_mean = data_sum / window_size;

  data2 = pow(data_zeromean, 2);
  NumericVector data2_sum = cumsum(diff2_lag(data2, window_size, sum(data2[Range(0, (window_size - 1))])));
  NumericVector data_sd2 = (data2_sum / window_size) - pow(data_mean, 2); // variance
  NumericVector data_sd = sqrt(data_sd2);
  NumericVector data_sig = sqrt(1 / (data_sd2 * window_size));

  return (List::create(Rcpp::Named("avg") = mov_mean, Rcpp::Named("sd") = data_sd, Rcpp::Named("sig") = data_sig,
                       Rcpp::Named("sum") = mov_sum, Rcpp::Named("sqrsum") = mov2_sum));
}

//[[Rcpp::export]]
int32_t mode_rcpp(const NumericVector x) {

  // is slower than R implementation...
  NumericVector ux = unique(x);
  int32_t y = ux[which_max(table(match(x, ux)))];
  return y;
}

//[[Rcpp::export]]
NumericVector znorm_rcpp(const NumericVector data) {
  double data_mean = mean(data);
  double data_dev = sqrt(sum(pow((data - data_mean), 2)) / data.length());

  if (data_dev == NA_REAL || data_dev <= 0.01) {
    return (data - data_mean);
  } else {
    return ((data - data_mean) / (data_dev));
  }
}

//[[Rcpp::export]]
NumericVector binary_split_rcpp(const uint32_t n) {

  NumericVector idxs(n);

  idxs[0] = 1; // We always begin by explore the first integer
  // After exploring the first integer, we begin splitting the interval 2:n

  std::deque<uint32_t> lb_list;
  std::deque<uint32_t> ub_list;

  lb_list.push_back(2);
  ub_list.push_back(n);

  uint32_t lb;
  uint32_t ub;
  uint32_t mid;

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
double inner_product(const NumericVector a, const NumericVector b) {
  double res = std::inner_product(a.begin(), a.end(), b.begin(), 0.0);

  return (res);
}

//[[Rcpp::export]]
double sum_of_squares(const NumericVector a) {
  double res = std::inner_product(a.begin(), a.end(), a.begin(), 0.0);

  return (res);
}

struct MuinWorker : public Worker {
  // input
  const RVector<double> a;
  const RVector<double> mu;
  const uint32_t w_size;
  // output
  RVector<double> sig;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  MuinWorker(const NumericVector a, const NumericVector mu, const uint32_t w_size, NumericVector sig)
      : a(a), mu(mu), w_size(w_size), sig(sig) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    uint32_t j, k;
    std::vector<double> b(w_size);

    for (uint32_t i = begin; i < end; i++) {

      for (j = i, k = 0; j < (i + w_size); j++, k++) {
        b[k] = a[j] - mu[i];
      }

      sig[i] = 1 / sqrt(std::inner_product(b.begin(), b.end(), b.begin(), 0.0));
    }
  }
};

//[[Rcpp::export]]
List muinvn_parallel_rcpp(const NumericVector a, uint32_t w) {
  // Functions here are based on the work in
  // Ogita et al, Accurate Sum and Dot Product
  // results here are a moving average and stable inverse centered norm based
  // on Accurate Sum and Dot Product, Ogita et al

  NumericVector sig(a.length() - w + 1);
  NumericVector mu = movsum_ogita_rcpp(a, w) / w;

  MuinWorker muin_worker(a, mu, w, sig);

  // call parallelFor to do the work
  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, mu.length(), muin_worker, 2 * w);
#else
    RcppParallel2::ttParallelFor(0, mu.length(), muin_worker, 2 * w);
#endif
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated.\n";
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("avg") = mu, Rcpp::Named("sig") = sig));
}

//[[Rcpp::export]]
List muinvn_rcpp(const NumericVector a, uint32_t w) {
  // Functions here are based on the work in
  // Ogita et al, Accurate Sum and Dot Product
  // results here are a moving average and stable inverse centered norm based
  // on Accurate Sum and Dot Product, Ogita et al

  NumericVector sig(a.length() - w + 1, 0);
  NumericVector mu = movsum_ogita_rcpp(a, w) / w;

  for (uint32_t i = 0; i < mu.length(); i++) {
    IntegerVector a_range = Range(i, i + w - 1);
    sig[i] = sum_of_squares(as<NumericVector>(a[a_range]) - mu[i]);
  }

  sig = 1 / sqrt(sig);

  return (List::create(Rcpp::Named("avg") = mu, Rcpp::Named("sig") = sig));
}

// Version compatible with RCPP code
//[[Rcpp::export]]
ComplexVector fft_rcpp(const ComplexVector z, bool invert) {

  ComplexVector result;
  int n = z.length();
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

std::vector<std::complex<double>> fft_rcpp(const std::vector<double> z, bool invert) {

  std::vector<std::complex<double>> result;
  int n = z.size();
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
  int n = z.size();
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

  int n = z.size();
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

ComplexVector fft_rcpp(const NumericVector z, bool invert) {

  ComplexVector result;
  int n = z.length();
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
