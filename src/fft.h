#ifndef __FFT__
#define __FFT__

#include <vector>
#include <complex>

std::vector<std::complex<double>> fftw(std::vector<std::complex<double>> z, bool inverse = false);
std::vector<std::complex<double>> fftw(std::vector<double> z, bool inverse = false);

#endif // __FFT__
