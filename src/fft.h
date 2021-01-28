#ifndef __FFT__
#define __FFT__
#include <complex>
#include <vector>

namespace FFT {

#ifndef MIN
#define MIN(y, x) ((x) < (y) && (x) == (x) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(y, x) ((x) > (y) && (x) == (x) ? (x) : (y))
#endif

#ifndef INT_MAX
#define INT_MAX 2147483647
#endif

#ifndef M_SQRT_3
#define M_SQRT_3 1.732050807568877293527446341506 /* sqrt(3) */
#endif
#ifndef M_PI_4
#define M_PI_4 0.785398163397448309615660845820 /* pi/4 */
#endif

typedef struct complex {
  double r;
  double i;
} complex_t;

class fftw {

public:
  // Constructor
  fftw();

  // Destructor
  ~fftw();

  std::vector<std::complex<double>> fft(std::vector<double> z, bool inverse);
  std::vector<std::complex<double>> fft(std::vector<std::complex<double>> z, bool inverse);

private:
  int old_n = 0;
  int nfac[20];
  int m_fac;
  int kt;
  int maxf;
  int maxp;
  double *work = nullptr;
  int *iwork = nullptr;
  complex_t *cplx = nullptr;

  /* non-API, but used by package RandomFields */
  void fft_factor(int n, int *pmaxf, int *pmaxp);
  int fft_work(double *a, double *b, int nseg, int n, int nspn, int isn, double *work, int *iwork);
  void fftmx(double *a, double *b, int ntot, int n, int nspan, int isn, int m, int kt, double *at, double *ck,
             double *bt, double *sk, int *np, int *nfac);
};
} // namespace FFT

#endif // __FFT__
