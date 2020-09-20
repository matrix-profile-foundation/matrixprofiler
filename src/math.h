#ifndef __MATH__
#define __MATH__

#include <Rcpp.h>

using namespace Rcpp;

#ifndef MIN
#define MIN(y, x) ((x) < (y) && (x) == (x) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(y, x) ((x) > (y) && (x) == (x) ? (x) : (y))
#endif

double        std_rcpp(const NumericVector data, const bool na_rm);
NumericMatrix list_to_matrix(const List x); // unnused?
IntegerVector which(const LogicalVector x);
NumericVector diff_lag(const NumericVector x, const uint32_t lag); // unnused?
NumericVector diff2_lag(const NumericVector x, const uint32_t lag, const double v);
NumericVector fast_movsd_rcpp(const NumericVector data, const uint32_t window_size);
List          fast_avg_sd_rcpp(const NumericVector data, const uint32_t window_size);
int32_t       mode_rcpp(const NumericVector x);
NumericVector znorm_rcpp(const NumericVector data);
NumericVector binary_split_rcpp(const uint32_t n);
double        inner_product(const NumericVector a, const NumericVector b);
double        sum_of_squares(const NumericVector a);
NumericVector sum2s_rcpp(const NumericVector a, uint32_t w);
List          muinvn_rcpp(const NumericVector a, uint32_t w);
ComplexVector fft_rcpp(const ComplexVector z, bool invert = false);
ComplexVector fft_rcpp(const NumericVector z, bool invert = false);
ComplexVector fft_rcpp_debug(const ComplexVector z);
ComplexVector fft_rcpp_debug2(const ComplexVector z, Function stats_fft);
ComplexVector fft_rcpp_debug3(const ComplexVector z, bool invert);
std::vector<std::complex<double>> fftw_rcpp_debug(const std::vector<std::complex<double>> z);

#endif // __MATH__
