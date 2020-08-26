#ifndef __STAMP__
#define __STAMP__

#include <Rcpp.h>
// #include <RcppParallel.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
// using namespace RcppParallel;

List stamp_cpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
               double ez = 0.5);


#endif // __STAMP__
