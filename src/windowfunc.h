#ifndef __WINDOWFUNC__
#define __WINDOWFUNC__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

NumericVector movmean_rcpp(NumericVector data, uint32_t window_size);
NumericVector movstd_rcpp(NumericVector data, uint32_t window_size);
List movmean_std_rcpp(NumericVector data, uint32_t window_size);
NumericVector movvar_rcpp(NumericVector data, uint32_t window_size);
NumericVector movvar2_rcpp(NumericVector data, uint32_t window_size);
NumericVector movsum_rcpp(NumericVector data, uint32_t window_size);
NumericVector movsum_ogita_rcpp(NumericVector data, uint32_t window_size);
double precision_test_rcpp(std::vector<double> d);
NumericVector movmin_rcpp(NumericVector data, uint32_t window_size);
NumericVector movmax_rcpp(NumericVector data, uint32_t window_size);
NumericVector movmean_weighted_rcpp(NumericVector data, uint32_t window_size, double eps = 1);
NumericVector movmean_fading_rcpp(NumericVector data, uint32_t window_size, double eps = 0.24);
NumericVector movsum_weighted_rcpp(NumericVector data, uint32_t window_size, double eps = 1);
NumericVector movsum_fading_rcpp(NumericVector data, uint32_t window_size, double eps = 0.4);
NumericVector movvar_weighted_rcpp(NumericVector data, uint32_t window_size, double eps = 1);
NumericVector movvar_fading_rcpp(NumericVector data, uint32_t window_size, double eps = 0.19);
List muinvn_rcpp(NumericVector data, uint32_t window_size);
List muinvn_rcpp_parallel(NumericVector data, uint32_t window_size);
IntegerVector zero_crossing_rcpp(NumericVector data, uint32_t window_size);

#endif // __WINDOWFUNC__
