#ifndef __STAMP__
#define __STAMP__

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List stamp_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double ez = 0.5,
                double s_size = 1.0, bool progress = false);
List stamp_rcpp_parallel(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double ez = 0.5,
                         double s_size = 1.0, bool progress = false);

#endif // __STAMP__
