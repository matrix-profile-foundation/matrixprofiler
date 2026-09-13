#ifndef __RFCP__
#define __RFCP__

#include <Rcpp.h>

using namespace Rcpp;

List rfcp_na_segmented_native_rcpp_parallel(NumericVector positive_ref, NumericVector negative_ref,
                                             uint64_t window_size, uint32_t max_freq,
                                             uint32_t exclusion_radius, uint64_t query_begin,
                                             uint64_t query_end, bool return_profiles = false,
                                             bool progress = false);

#endif // __RFCP__
