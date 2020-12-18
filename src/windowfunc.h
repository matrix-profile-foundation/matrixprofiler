#ifndef __WINDOWFUNC__
#define __WINDOWFUNC__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector movmean_rcpp(const NumericVector data,
                           const uint32_t window_size);
NumericVector movvar_rcpp(const NumericVector data, const uint32_t window_size);
NumericVector movvar2_rcpp(const NumericVector data, uint32_t window_size);
NumericVector movsum_rcpp(NumericVector data, uint32_t window_size);
NumericVector movsum_ogita_rcpp(const NumericVector a, uint32_t w);
double precision_test_rcpp(std::vector<double> d);
NumericVector movmin_rcpp(const NumericVector data, uint32_t window_size);
NumericVector movmax_rcpp(const NumericVector data, uint32_t window_size);
NumericVector movmean_weighted_rcpp(const NumericVector data,
                                    uint32_t window_size, double eps = 1);
NumericVector movmean_fading_rcpp(const NumericVector data,
                                  uint32_t window_size, double eps = 0.24);
NumericVector movsum_weighted_rcpp(NumericVector data, uint32_t window_size,
                                   double eps = 1);
NumericVector movsum_fading_rcpp(NumericVector data, uint32_t window_size,
                                 double eps = 0.4);
NumericVector movvar_weighted_rcpp(const NumericVector data,
                                   uint32_t window_size, double eps = 1);
NumericVector movvar_fading_rcpp(const NumericVector data, uint32_t window_size,
                                 double eps = 0.19);

#endif // __WINDOWFUNC__
