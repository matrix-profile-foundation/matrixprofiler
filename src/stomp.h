#ifndef __STOMP__
#define __STOMP__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List stomp_cpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
               double ez = 0.5);

#endif // __STOMP__
