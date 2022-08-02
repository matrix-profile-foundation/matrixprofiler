#ifndef __CONTRAST__
#define __CONTRAST__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List contrast_profile_rcpp(NumericVector negative_data, List positive_data, uint64_t window_size, List positive_matrix,
                           double ez = 0.5, double s_size = 1.0, uint8_t n_workers = 1, bool idxs = true,
                           bool euclidean = true, bool progress = false);

#endif // __CONTRAST__
