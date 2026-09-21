#ifndef MATRIXPROFILER_MPX_GPU_H
#define MATRIXPROFILER_MPX_GPU_H

#include <Rcpp.h>

using namespace Rcpp;

List mpx_gpu_rcpp(NumericVector data_ref, uint64_t window_size, double ez = 0.5, bool idxs = true,
                  bool euclidean = true);
List mpx_na_gpu_rcpp(NumericVector data_ref, uint64_t window_size, double ez = 0.5, bool idxs = true,
                     bool euclidean = true);
List mpxab_gpu_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, bool idxs = true,
                    bool euclidean = true);
List mpxab_na_gpu_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, bool idxs = true,
                       bool euclidean = true);

#endif // MATRIXPROFILER_MPX_GPU_H
