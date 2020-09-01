#ifndef __MASS__
#define __MASS__

#include <Rcpp.h>
// #include <RcppParallel.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
// using namespace RcppParallel;

List      mass_weighted_rcpp(ComplexVector data_fft, NumericVector query_window, uint32_t data_size,
                        uint32_t window_size, NumericVector data_mean, NumericVector data_sd,
                        double query_mean, double query_sd, NumericVector data_pre,
                        NumericVector weight, Function fft);
List      mass_absolute_rcpp(ComplexVector data_fft, NumericVector query_window, uint32_t data_size,
                             uint32_t window_size, NumericVector sumx2, double sumy2, Function fft);
List      mass2_rcpp(ComplexVector data_fft, NumericVector query_window, uint64_t data_size,
                     uint32_t window_size, NumericVector data_mean, NumericVector data_sd,
                     double query_mean, double query_sd, Function fft);
List      mass3_rcpp(NumericVector query_window, NumericVector data_ref,
                     uint64_t data_size, uint32_t window_size, NumericVector data_mean,
                     NumericVector data_sd, double query_mean, double query_sd, Function fft,
                     uint32_t k = 1024);
uint32_t  set_k(uint32_t k, uint64_t data_size, uint64_t window_size);
uint32_t  find_best_k(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size);
List      mass_pre_rcpp(const NumericVector data, const NumericVector query, uint32_t window_size);
List      mass_pre_abs_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size);
List      mass_pre_weighted_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
                            const NumericVector weight);
#endif // __MASS__
