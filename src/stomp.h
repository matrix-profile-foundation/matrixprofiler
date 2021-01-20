#ifndef __STOMP__
#define __STOMP__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List stomp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint64_t window_size, double ez = 0.5,
                bool progress = false);

List stomp_cpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint64_t window_size,
                        double ez = 0.5, bool progress = false);

#endif // __STOMP__
