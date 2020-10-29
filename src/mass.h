#ifndef __MASS__
#define __MASS__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
// using namespace RcppParallel;
List      mass_weighted_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint32_t data_size,
                             uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd,
                             double query_mean, double query_sd, const NumericVector data_pre,
                             const NumericVector weight);
List      mass_absolute_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint32_t data_size,
                             uint32_t window_size, const NumericVector sumx2, double sumy2);
List          mass2_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint64_t data_size,
                         uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd,
                         double query_mean, double query_sd);
List          mass3_rcpp(const NumericVector query_window, const NumericVector data_ref,
                         uint64_t data_size, uint32_t window_size, const NumericVector data_mean,
                         const NumericVector data_sd, double query_mean, double query_sd,
                         uint32_t k = 4096);
List          mass3_rcpp_parallel(const NumericVector query_window, const NumericVector data_ref,
                                  uint64_t data_size, uint32_t window_size, const NumericVector data_mean,
                                  const NumericVector data_sd, double query_mean, double query_sd, uint16_t k = 8192);
uint32_t      set_k_rcpp(uint32_t k, uint64_t data_size, uint64_t window_size);
uint32_t      find_best_k_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size);
List          mass_pre_rcpp(const NumericVector data, const NumericVector query, uint32_t window_size);
List          mass_pre_abs_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size);
List          mass_pre_weighted_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
                                     const NumericVector weight);

#endif // __MASS__
