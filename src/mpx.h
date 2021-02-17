#ifndef __MPX__
#define __MPX__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List mpx_rcpp(NumericVector data_ref, uint64_t window_size, double ez = 0.5, double s_size = 1.0, bool idxs = true,
              bool euclidean = true, bool progress = false);
List mpxab_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double s_size = 1.0,
                bool idxs = true, bool euclidean = true, bool progress = false);
List mpx_rcpp_parallel(NumericVector data_ref, uint64_t window_size, double ez = 0.5, double s_size = 1.0,
                       bool idxs = true, bool euclidean = true, bool progress = false);
List mpxab_rcpp_parallel(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double s_size = 1.0,
                         bool idxs = true, bool euclidean = true, bool progress = false);

#endif // __MPX__
