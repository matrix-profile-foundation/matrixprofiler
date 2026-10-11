#ifndef __MPDIST__
#define __MPDIST__

#include <Rcpp.h>
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

using namespace Rcpp;

NumericVector mpdist_vect_rcpp(NumericVector data, NumericVector query, uint32_t window_size, double thr,
                              uint32_t n_workers);

double mpdist_simple_rcpp(NumericVector ref_data, NumericVector query_data, uint32_t window_size, double thr,
                          uint32_t n_workers);

#endif // __MPDIST__
