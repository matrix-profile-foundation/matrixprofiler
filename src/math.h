#ifndef __MATH__
#define __MATH__

#include "fft.h"
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

#ifndef MIN
#define MIN(y, x) ((x) < (y) && (x) == (x) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(y, x) ((x) > (y) && (x) == (x) ? (x) : (y))
#endif

IntegerVector seq(uint64_t start, uint64_t end);
IntegerVector seq_by(uint64_t start, uint64_t end, uint32_t by);
double std_rcpp(const NumericVector data, const bool na_rm);
NumericMatrix list_to_matrix(const List x); // unused?
IntegerVector which_cpp(const LogicalVector x);
int32_t mode_rcpp(const IntegerVector x);
NumericVector znorm_rcpp(const NumericVector data);
NumericVector normalize_rcpp(const NumericVector data, double min = 0, double max = 1);
IntegerVector binary_split_rcpp(const uint32_t n);
NumericVector ed_corr_rcpp(const NumericVector data, uint32_t window_size);
NumericVector corr_ed_rcpp(const NumericVector data, uint32_t window_size);
double inner_product(const NumericVector a, const NumericVector b);
double sum_of_squares(const NumericVector a);
ComplexVector fft_rcpp(const ComplexVector z, bool invert = false);
ComplexVector fft_rcpp(const NumericVector z, bool invert = false);
std::vector<std::complex<double>> fft_rcpp(const std::vector<double> z, bool invert = false);
std::vector<std::complex<double>> fft_rcpp(const std::vector<std::complex<double>> z, bool invert = false);
std::vector<double> fft_rcpp_real(const std::vector<std::complex<double>> z, bool invert = false);

#endif // __MATH__
