#ifndef __MASS__
#define __MASS__

#include <Rcpp.h>
// #include <RcppParallel.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
// using namespace RcppParallel;

List mass_rcpp(ComplexVector data_fft, NumericVector query_window, uint64_t data_size,
               uint32_t window_size, NumericVector data_mean, NumericVector data_sd,
               double query_mean, double query_sd, Function fft);
List mass2_rcpp(ComplexVector data_fft, NumericVector query_window, uint64_t data_size,
                uint32_t window_size, NumericVector data_mean, NumericVector data_sd,
                double query_mean, double query_sd, Function fft);
List mass3_rcpp(NumericVector query_window, NumericVector data,
                uint32_t window_size, uint64_t data_size, NumericVector data_mean,
                NumericVector data_sd, double query_mean, double query_sd, Function fft,
                uint32_t k = 1024);
List mass_pre_rcpp(const NumericVector data, const NumericVector query, uint32_t window_size, Function fft);

#endif // __MASS__
