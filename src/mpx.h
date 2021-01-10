#ifndef __MPX__
#define __MPX__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List mpx_rcpp(NumericVector data_ref, uint64_t window_size, double ez, bool idxs, bool euclidean, bool progress);
List mpx_rcpp_parallel(NumericVector data_ref, uint64_t window_size, double ez, bool idxs, bool euclidean,
                       bool progress);
List mpxab_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, bool idxs, bool euclidean,
                bool progress);
List mpxab_rcpp_parallel(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, bool idxs,
                         bool euclidean, bool progress);

#endif // __MPX__
