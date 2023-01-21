#ifndef __MATHTOOLS__
#define __MATHTOOLS__

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
double std_rcpp(NumericVector data, bool na_rm);
NumericMatrix list_to_matrix(List x); // unused?
IntegerVector which_cpp(LogicalVector x);
int32_t mode_rcpp(IntegerVector x);
NumericVector znorm_rcpp(NumericVector data);
NumericVector normalize_rcpp(NumericVector data, double min = 0, double max = 1);
IntegerVector binary_split_rcpp(uint32_t n);
NumericVector ed_corr_rcpp(NumericVector data, uint32_t window_size);
NumericVector corr_ed_rcpp(NumericVector data, uint32_t window_size);
double inner_product(NumericVector a, NumericVector b);
double sum_of_squares(NumericVector a);
ComplexVector fft_rcpp(ComplexVector z, bool invert = false);
ComplexVector fft_rcpp(NumericVector z, bool invert = false);
std::vector<std::complex<double>> fft_rcpp(std::vector<double> z, bool invert = false);
std::vector<std::complex<double>> fft_rcpp(std::vector<std::complex<double>> z, bool invert = false);
std::vector<double> fft_rcpp_real(std::vector<std::complex<double>> z, bool invert = false);

#endif // __MATHTOOLS__
