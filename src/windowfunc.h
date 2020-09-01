#ifndef __WINDOWFUNC__
#define __WINDOWFUNC__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector movsum(NumericVector data, uint32_t window_size);
NumericVector movmin(const NumericVector data, uint32_t window_size);
NumericVector movmax(const NumericVector data, uint32_t window_size);
NumericVector mov_mean_online_rcpp(const NumericVector data, uint32_t window_size, double eps);

#endif // __WINDOWFUNC__
