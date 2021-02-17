#ifndef __WINDOWFUNC__
#define __WINDOWFUNC__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

NumericVector movmean_rcpp(const NumericVector data, const uint32_t window_size);
NumericVector movstd_rcpp(const NumericVector data, const uint32_t window_size);
List movmean_std_rcpp(const NumericVector data, const uint32_t window_size);
NumericVector movvar_rcpp(const NumericVector data, const uint32_t window_size);
NumericVector movvar2_rcpp(const NumericVector data, uint32_t window_size);
NumericVector movsum_rcpp(NumericVector data, uint32_t window_size);
NumericVector movsum_ogita_rcpp(const NumericVector data, uint32_t window_size);
double precision_test_rcpp(std::vector<double> d);
NumericVector movmin_rcpp(const NumericVector data, uint32_t window_size);
NumericVector movmax_rcpp(const NumericVector data, uint32_t window_size);
NumericVector movmean_weighted_rcpp(const NumericVector data, uint32_t window_size, double eps = 1);
NumericVector movmean_fading_rcpp(const NumericVector data, uint32_t window_size, double eps = 0.24);
NumericVector movsum_weighted_rcpp(NumericVector data, uint32_t window_size, double eps = 1);
NumericVector movsum_fading_rcpp(NumericVector data, uint32_t window_size, double eps = 0.4);
NumericVector movvar_weighted_rcpp(const NumericVector data, uint32_t window_size, double eps = 1);
NumericVector movvar_fading_rcpp(const NumericVector data, uint32_t window_size, double eps = 0.19);
List muinvn_rcpp(const NumericVector data, uint32_t window_size);
List muinvn_rcpp_parallel(const NumericVector data, uint32_t window_size);
IntegerVector zero_crossing_rcpp(const NumericVector data, uint32_t window_size);

#endif // __WINDOWFUNC__
