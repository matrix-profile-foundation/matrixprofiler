#ifndef __MPX__
#define __MPX__

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List mpxis_rcpp(NumericVector data_ref, uint64_t batch_size, List object, List stats, uint64_t history = 0,
                uint64_t mp_time_constraint = 0, bool progress = false, float threshold = -1.0);
List mpxi_rcpp(NumericVector new_data, List object, uint64_t history = 0, uint64_t mp_time_constraint = 0,
               bool progress = false);
List mpxiright_rcpp(NumericVector data_ref, uint64_t window_size, double ez = 0.5, uint64_t mp_time_constraint = 0,
                    double s_size = 1.0, bool idxs = true, bool euclidean = true, bool progress = false,
                    uint64_t start = 1, List old = R_NilValue);
List mpxileft_rcpp(NumericVector data_ref, uint64_t window_size, double ez = 0.5, double s_size = 1.0, bool idxs = true,
                   bool euclidean = true, bool progress = false, uint64_t start = 1, List old = R_NilValue);
List mpx_rcpp_new(NumericVector data_ref, uint64_t window_size, double ez = 0.5, uint64_t mp_time_constraint = 0,
                  double s_size = 1.0, bool idxs = true, bool euclidean = true, bool progress = false);
List mpx_rcpp(NumericVector data_ref, uint64_t window_size, double ez = 0.5, double s_size = 1.0, bool idxs = true,
              bool euclidean = true, bool progress = false);
List mpxab_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double s_size = 1.0,
                bool idxs = true, bool euclidean = true, bool progress = false);
#endif // __MPX__
