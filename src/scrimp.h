#ifndef __SCRIMP__
#define __SCRIMP__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List scrimp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                 double pre_scrimp, bool progress);

#endif // __SCRIMP__
