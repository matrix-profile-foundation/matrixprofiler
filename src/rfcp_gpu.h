#ifndef MATRIXPROFILER_RFCP_GPU_H
#define MATRIXPROFILER_RFCP_GPU_H

#include <Rcpp.h>

using namespace Rcpp;

List rfcp_gpu_rcpp(NumericVector positive_ref, NumericVector negative_ref,
                   uint64_t window_size, uint32_t max_freq,
                   uint32_t exclusion_radius,
                   uint64_t query_begin, uint64_t query_end,
                   bool return_profiles = false);

#endif // MATRIXPROFILER_RFCP_GPU_H
