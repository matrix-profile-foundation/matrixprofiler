#ifndef __MSTOMP__
#define __MSTOMP__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List mstomp_rcpp(const NumericMatrix data_ref, uint32_t window_size, double ez = 0.5,
                 Nullable<IntegerVector> must_dim = R_NilValue,
                 Nullable<IntegerVector> exc_dim = R_NilValue,
                 uint32_t n_workers = 1, bool progress = false);

#endif // __MSTOMP__
