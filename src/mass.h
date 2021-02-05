#ifndef __MASS__
#define __MASS__

#include <Rcpp.h>
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

List mass_weighted_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint32_t data_size,
                        uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd,
                        double query_mean, double query_sd, const NumericVector data_pre, const NumericVector weight,
                        const bool normalized = true);
List mass_absolute_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint32_t data_size,
                        uint32_t window_size, const NumericVector sumx2, double sumy2);
List mass2_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint64_t data_size,
                uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd, double query_mean,
                double query_sd);
List mass3_rcpp(const NumericVector query_window, const NumericVector data_ref, uint64_t data_size,
                uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd, double query_mean,
                double query_sd, uint32_t k = 4096);
List mass3_rcpp_parallel(const NumericVector query_window, const NumericVector data_ref, uint64_t data_size,
                         uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd,
                         double query_mean, double query_sd, uint16_t k = 8192);
uint32_t set_k_rcpp(uint32_t k, uint64_t data_size, uint64_t window_size);
uint32_t find_best_k_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size);
List mass_pre_rcpp(const NumericVector data, const NumericVector query, uint32_t window_size);
List mass_pre_abs_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size);
List mass_pre_weighted_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
                            const NumericVector weight);

// redeclaration for MacOS compilation
std::vector<std::complex<double>> fft_rcpp(const std::vector<double> z, bool invert);
std::vector<double> fft_rcpp_real(const std::vector<std::complex<double>> z, bool invert);

template <typename Iterator>
List mass3_cpp(const Iterator query_it, const Iterator data_it, const uint64_t data_size, const uint32_t window_size,
               const Iterator data_mean_it, const Iterator data_sd_it, const double query_mean, const double query_sd,
               uint32_t k) {

  uint32_t w_size = window_size;
  uint64_t d_size = data_size;
  uint64_t p_size = data_size - window_size + 1;

  std::vector<double> dist(p_size);
  std::vector<double>::iterator dist_it = dist.begin();
  std::vector<double> last(p_size);
  std::vector<double>::iterator last_it = last.begin();

  k = set_k_rcpp(k, d_size, w_size);

  // compute query_it stats -- O(d_size)
  double q_mean = query_mean;
  double q_std = query_sd;

  // compute data stats -- O(d_size)
  Iterator d_mean_it = data_mean_it;
  Iterator d_std_it = data_sd_it;

  std::vector<double> rev_query(k);
  std::reverse_copy(query_it, query_it + w_size, rev_query.begin());

  std::vector<std::complex<double>> Y = fft_rcpp(rev_query, false);

  uint64_t j = 0;
  uint64_t jump = k - w_size + 1;
  uint64_t seq_end = d_size - k;
  std::vector<std::complex<double>> Z(k);
  std::vector<double> z;
  std::vector<double> d(k - w_size + 1);

  try {
    for (j = 0; j <= seq_end; j = j + jump) {
      // The main trick of getting dot products in O(d_size log d_size) time
      uint64_t idx_begin = j;

      std::vector<double> data_chunk(data_it + j, data_it + j + k);
      std::vector<std::complex<double>> X = fft_rcpp(data_chunk, false);

      std::transform(X.begin(), X.end(), Y.begin(), Z.begin(), std::multiplies<std::complex<double>>());

      z = fft_rcpp_real(Z, true);

      for (uint64_t i = 0; i < (k - w_size + 1); i++) {
        d[i] = 2 * (w_size - (z[w_size - 1 + i] - w_size * *(d_mean_it + idx_begin + i) * q_mean) /
                                 (*(d_std_it + idx_begin + i) * q_std));
      }

      std::copy(d.begin(), d.begin() + k - w_size + 1, dist_it + j);
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

      if ((jump - (w_size - 1) + j) > (uint64_t)p_size) {
        Rcout << "DEBUG: error." << std::endl;
      } else {
        std::vector<double> data_chunk(data_it + j, data_it + d_size);
        std::vector<std::complex<double>> X = fft_rcpp(data_chunk, false);
        std::vector<double> rev_query_chunk(rev_query.begin(), rev_query.begin() + jump);
        Y = fft_rcpp(rev_query_chunk, false);

        Z = std::vector<std::complex<double>>(Y.size());

        std::transform(X.begin(), X.end(), Y.begin(), Z.begin(), std::multiplies<std::complex<double>>());

        z = fft_rcpp_real(Z, true);

        for (uint64_t i = 0; i < (jump - w_size + 1); i++) {
          d[i] = 2 * (w_size - (z[i + w_size - 1] - w_size * *(d_mean_it + idx_begin + i) * q_mean) /
                                   (*(d_std_it + idx_begin + i) * q_std));
        }

        std::copy(d.begin(), d.begin() + (jump - w_size + 1), dist_it + j);
        std::copy(z.begin() + w_size - 1, z.begin() + jump, last_it + j);
      }
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("distance_profile") = dist, Rcpp::Named("last_product") = last));
}

#endif // __MASS__
