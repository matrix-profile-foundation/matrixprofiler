#include "mass.h"

//[[Rcpp::export]]
List mass_rcpp(ComplexVector data_fft, NumericVector query_window, uint64_t data_size,
               uint32_t window_size, NumericVector data_mean, NumericVector data_sd,
               double query_mean, double query_sd, Function fft) {

  // pre-process query for fft

  NumericVector rev_query(data_fft.length());
  std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());

  ComplexVector prod = data_fft * as<ComplexVector>(fft(rev_query));
  NumericVector z = Re(as<ComplexVector>(fft(prod, true))) / prod.length();
  // compute the distance profile
  // IntegerVector range(data_fft.length());
  // std::iota(range.begin(), range.end(), window_size - 1);
  // NumericVector last_product = z[range];
  NumericVector last_product = z[Range(window_size - 1, data_size - 1)];
  NumericVector distance_profile = 2 * (window_size - (last_product - window_size * data_mean * query_mean) / (data_sd * query_sd));
  distance_profile[distance_profile < 0] = 0;

  return (List::create(
            Rcpp::Named("distance_profile") = distance_profile,
            Rcpp::Named("last_product") = last_product
          ));
}

//[[Rcpp::export]]
List mass2_rcpp(ComplexVector data_fft, NumericVector query_window, uint64_t data_size,
                uint32_t window_size, NumericVector data_mean, NumericVector data_sd,
                double query_mean, double query_sd, Function fft) {

  // pre-process query for fft

  NumericVector rev_query(data_fft.length());
  std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());

  // compute the product
  ComplexVector prod = data_fft * as<ComplexVector>(fft(rev_query));
  NumericVector z = Re(as<ComplexVector>(fft(prod, true))) / prod.length();
  // compute the distance profile
  NumericVector last_product = z[Range(window_size - 1, data_size - 1)];
  NumericVector distance_profile = 2 * (window_size - (last_product - window_size * data_mean * query_mean) / (data_sd * query_sd));
  distance_profile[distance_profile < 0] = 0;

  return (List::create(
            Rcpp::Named("distance_profile") = distance_profile,
            Rcpp::Named("last_product") = last_product
          ));
}

//[[Rcpp::export]]
List mass3_rcpp(NumericVector query_window, NumericVector data,
                uint32_t window_size, uint64_t data_size, NumericVector data_mean,
                NumericVector data_sd, double query_mean, double query_sd, Function fft,
                uint32_t k = 1024) {

  // data is the long time series
  // query_window is the query
  // k is the size of pieces, preferably a power of two
  // data is the data, query_window is the query
  //

  if (k > data_size) {
    k = pow(2, ceil(log2(sqrt(data_size))));
  }

  if (k <= window_size) {
    k = pow(2,(ceil(log2(window_size)) + 1));
    if (k > data_size) {
      k = data_size;
    }
  }

  uint32_t w_size = window_size;
  uint64_t d_size = data_size;
  NumericVector dist(data_mean.length());
  NumericVector::iterator dist_it = dist.begin();
  NumericVector last(data_mean.length());
  NumericVector::iterator last_it = last.begin();


  // compute query_window stats -- O(d_size)
  double q_mean = query_mean;
  double q_std = query_sd;

  // compute data stats -- O(d_size)
  NumericVector d_mean = data_mean;
  NumericVector d_std = data_sd;

  // k = 4096; // assume k > w_size
  // k = pow2(nextpow2(sqrt(d_size)));

  NumericVector rev_query(k);
  std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());
  ComplexVector Y = fft(rev_query);

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
      uint64_t idx_end = j + k - w_size;
      ComplexVector X = fft(data[Range(j, j + k - 1)]);
      Z = X * Y;
      z = Re(as<ComplexVector>(fft(Z, true))) / Z.length();
      d = 2 * (w_size - (z[Range(w_size - 1, k - 1)] - w_size * d_mean[Range(idx_begin, idx_end)]
                         * q_mean) / (d_std[Range(idx_begin, idx_end)] * q_std));
      std::copy(d.begin(), d.end(), dist_it + j);
      std::copy(z.begin() + w_size - 1, z.begin() + k, last_it + j);
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  // j = j + jump;
  jump = d_size - j;

  try {
    if (jump >= w_size) {
      uint64_t idx_begin = j;
      uint64_t idx_end = d_size - w_size;

      if((jump - (w_size - 1) + j) > data_mean.length()) {
        Rcout << "DEBUG: error." << std::endl;
      } else {
        ComplexVector X = fft(data[Range(j, d_size - 1)]);
        Y = fft(rev_query[Range(0, jump - 1)]);
        Z = X * Y;
        z = Re(as<ComplexVector>(fft(Z, true))) / Z.length();
        d = 2 * (w_size - (z[Range(w_size - 1, jump - 1)] - w_size * d_mean[Range(idx_begin, idx_end)] * q_mean) / (d_std[Range(idx_begin, idx_end)] * q_std));
        std::copy(d.begin(), d.end(), dist_it + j);
        std::copy(z.begin() + w_size - 1, z.begin() + jump, last_it + j);
        Rcout << "DEBUG: Finished. k=" << k <<std::endl;
      }
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(
            Rcpp::Named("distance_profile") = dist,
            Rcpp::Named("last_product") = last
          ));
}
