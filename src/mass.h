#ifndef __MASS__
#define __MASS__

#include <Rcpp.h>
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

List mass_weighted_rcpp(ComplexVector data_fft, NumericVector query_window, uint32_t data_size, uint32_t window_size,
                        NumericVector data_mean, NumericVector data_sd, double query_mean, double query_sd,
                        NumericVector data_pre, NumericVector weight, bool normalized = true);
List mass_absolute_rcpp(ComplexVector data_fft, NumericVector query_window, uint32_t data_size, uint32_t window_size,
                        NumericVector sumx2, double sumy2);
List mass2_rcpp(ComplexVector data_fft, NumericVector query_window, uint64_t data_size, uint32_t window_size,
                NumericVector data_mean, NumericVector data_sd, double query_mean, double query_sd);
List mass3_rcpp(const NumericVector &query_window, const NumericVector &data_ref, uint64_t data_size,
                uint32_t window_size, const NumericVector &data_mean, const NumericVector &data_sd, double query_mean,
                double query_sd, uint32_t grain = 4096);
List mass3_rcpp_parallel(NumericVector query_window, NumericVector data_ref, uint64_t data_size, uint32_t window_size,
                         NumericVector data_mean, NumericVector data_sd, double query_mean, double query_sd,
                         uint16_t grain = 8192);
uint32_t set_k_rcpp(uint32_t grain, uint64_t data_size, uint64_t window_size);
uint32_t find_best_k_rcpp(NumericVector data_ref, NumericVector query_ref, uint32_t window_size);
List mass_pre_rcpp(NumericVector data, NumericVector query, uint32_t window_size);
List mass_pre_abs_rcpp(NumericVector data_ref, NumericVector query_ref, uint32_t window_size);
List mass_pre_weighted_rcpp(NumericVector data_ref, NumericVector query_ref, uint32_t window_size,
                            NumericVector weight);

// redeclaration for MacOS compilation
std::vector<std::complex<double>> fft_rcpp(std::vector<double> z, bool invert);
std::vector<double> fft_rcpp_real(std::vector<std::complex<double>> z, bool invert);

template <typename Iterator>
List mass3_cpp(const Iterator query_it, const Iterator data_it, const uint64_t data_size, const uint32_t window_size,
               const Iterator data_mean_it, const Iterator data_sd_it, const double query_mean, const double query_sd,
               uint32_t grain) {

  uint32_t w_size = window_size;
  uint64_t d_size = data_size;
  uint64_t const p_size = data_size - window_size + 1;

  std::vector<double> dist(p_size);
  std::vector<double>::iterator const dist_it = dist.begin();
  std::vector<double> last(p_size);
  std::vector<double>::iterator const last_it = last.begin();

  grain = set_k_rcpp(grain, d_size, w_size);

  // compute query_it stats -- O(d_size)
  double q_mean = query_mean;
  double q_std = query_sd;

  // compute data stats -- O(d_size)
  Iterator d_mean_it = data_mean_it;
  Iterator d_std_it = data_sd_it;

  std::vector<double> rev_query(grain);
  std::reverse_copy(query_it, query_it + w_size, rev_query.begin());

  std::vector<std::complex<double>> Y = fft_rcpp(rev_query, false);

  uint64_t j = 0;
  uint64_t jump = grain - w_size + 1;
  uint64_t const seq_end = d_size - grain;
  std::vector<std::complex<double>> Z(grain);
  std::vector<double> z;
  std::vector<double> d(grain - w_size + 1);

  try {
    for (j = 0; j <= seq_end; j = j + jump) {
      // The main trick of getting dot products in O(d_size log d_size) time
      uint64_t idx_begin = j;

      std::vector<double> data_chunk(data_it + j, data_it + j + grain);
      std::vector<std::complex<double>> X = fft_rcpp(data_chunk, false);

      std::transform(X.begin(), X.end(), Y.begin(), Z.begin(), std::multiplies<std::complex<double>>());

      z = fft_rcpp_real(Z, true);

      for (uint64_t i = 0; i < (grain - w_size + 1); i++) {
        d[i] = 2 * (w_size - (z[w_size - 1 + i] - w_size * *(d_mean_it + idx_begin + i) * q_mean) /
                                 (*(d_std_it + idx_begin + i) * q_std));
      }

      std::copy(d.begin(), d.begin() + grain - w_size + 1, dist_it + static_cast<int64_t>(j));
      std::copy(z.begin() + w_size - 1, z.begin() + grain, last_it + static_cast<int64_t>(j));
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    stop("c++ exception (unknown reason)");
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
        std::vector<double> const rev_query_chunk(rev_query.begin(), rev_query.begin() + static_cast<int64_t>(jump));
        Y = fft_rcpp(rev_query_chunk, false);

        Z = std::vector<std::complex<double>>(Y.size());

        std::transform(X.begin(), X.end(), Y.begin(), Z.begin(), std::multiplies<std::complex<double>>());

        z = fft_rcpp_real(Z, true);

        for (uint64_t i = 0; i < (jump - w_size + 1); i++) {
          d[i] = 2 * (w_size - (z[i + w_size - 1] - w_size * *(d_mean_it + idx_begin + i) * q_mean) /
                                   (*(d_std_it + idx_begin + i) * q_std));
        }

        std::copy(d.begin(), d.begin() + static_cast<int64_t>(jump - w_size + 1), dist_it + static_cast<int64_t>(j));
        std::copy(z.begin() + w_size - 1, z.begin() + static_cast<int64_t>(jump), last_it + static_cast<int64_t>(j));
      }
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    stop("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("distance_profile") = dist, Rcpp::Named("last_product") = last));
}

#endif // __MASS__
