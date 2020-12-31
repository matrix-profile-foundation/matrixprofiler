#ifndef __STAMP__
#define __STAMP__

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List stamp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint64_t window_size, double ez,
                bool progress);
List stamp_rcpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint64_t window_size, double ez,
                         bool progress);

#endif // __STAMP__
