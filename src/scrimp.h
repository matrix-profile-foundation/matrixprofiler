#ifndef __SCRIMP__
#define __SCRIMP__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List scrimp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez = 0.5,
                 double s_size = 1.0, double pre_scrimp = 0.25, bool progress = false);
List scrimp_rcpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
                          double ez = 0.5, double s_size = 1.0, bool progress = false);
List scrimpab_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
                   double s_size = 1.0, bool progress = false);
// List scrimpab_rcpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double
// s_size, bool progress = false);

#endif // __SCRIMP__
