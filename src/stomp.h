#ifndef __STOMP__
#define __STOMP__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List stomp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint64_t window_size, double ez,
                bool progress);

List stomp_cpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint64_t window_size, double ez,
                        bool progress);

#endif // __STOMP__
