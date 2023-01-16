#ifndef __STOMP__
#define __STOMP__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List stomp_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double ez = 0.5,
                bool progress = false, bool left_right_profile = false);

List stomp_cpp_parallel(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double ez = 0.5,
                        bool progress = false, bool left_right_profile = false);

#endif // __STOMP__
