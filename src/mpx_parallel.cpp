// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#include "mathtools.h" // math first to fix OSX error
#include "mpx_parallel.h"
#include "windowfunc.h"
#include <cfloat> // DBL_EPSILON when STRICT_R_HEADERS
#include <mutex>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

#if RCPP_PARALLEL_USE_TBB
#include "tbb/mutex.h"
#else
#include "rcpp_parallel_fix.h"
#include "tthread/tinythread.h"
#endif

struct MatrixProfileP : public Worker {
  // input
private:
  const RVector<double> data_ref;
  const uint64_t window_size;
  const RVector<int> compute_order;
  const RVector<double> df;
  const RVector<double> dg;
  const RVector<double> mu;
  const RVector<double> sig;
  const RVector<double> ww;

  // output
  RVector<double> mp;
  RVector<int> mpi;

#if RCPP_PARALLEL_USE_TBB
  tbb::spin_mutex m;
#else
  tthread::mutex m;
#endif

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
public:
  MatrixProfileP(const NumericVector &data_ref, const uint64_t window_size, const IntegerVector &compute_order,
                 const NumericVector &df, const NumericVector &dg, const NumericVector &mmu, const NumericVector &sig,
                 const NumericVector &ww, const NumericVector &mp, const IntegerVector &mpi)
      : data_ref(data_ref), window_size(window_size), compute_order(compute_order), df(df), dg(dg), mu(mmu), sig(sig),
        ww(ww), mp(mp), mpi(mpi) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) override { // exclusion_zone:profile_len
    double c = 0.0, c_cmp = 0.0;
    uint32_t off_max = 0UL, off_diag = 0UL, offset = 0UL;
    uint32_t const n = data_ref.length();
    std::vector<double> aa(window_size);

    std::vector<double> mpp(mp.size(), -1.0);
    std::vector<int> mpip(mp.size(), -1);

    for (std::size_t dd = begin; dd < end; dd++) {
        uint32_t const diag = compute_order[dd];

        for (uint64_t i = 0; i < window_size; i++) {
          aa[i] = data_ref[diag + i] - mu[diag];
        }

        c = std::inner_product(aa.begin(), aa.end(), ww.begin(), 0.0);

        off_max = (n - window_size - diag + 1);

        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          c = c + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = c * sig[offset] * sig[off_diag];
          if (c_cmp > mpp[offset]) {
            mpp[offset] = c_cmp;
            mpip[offset] = off_diag + 1;
          }
          if (c_cmp > mpp[off_diag]) {
            mpp[off_diag] = c_cmp;
            mpip[off_diag] = offset + 1;
          }
        }
      }

      std::lock_guard<decltype(m)> lock(m);
      for (uint64_t i = 0; i < mp.size(); i++) {
        if (mpp[i] > mp[i]) {
          mp[i] = mpp[i];
          mpi[i] = mpip[i];
        }
    }
  }
};

struct MatrixProfilePNA : public Worker {
private:
  const RVector<double> data;
  const uint64_t window_size;
  const RVector<int> compute_order;
  const RVector<double> df;
  const RVector<double> dg;
  const RVector<double> mu;
  const RVector<double> sig;
  const RVector<double> first_window;
  const RVector<int> valid_window;
  RVector<double> mp;
  RVector<int> mpi;

#if RCPP_PARALLEL_USE_TBB
  tbb::spin_mutex mutex;
#else
  tthread::mutex mutex;
#endif

public:
  MatrixProfilePNA(const NumericVector &data, uint64_t window_size, const IntegerVector &compute_order,
                   const NumericVector &df, const NumericVector &dg, const NumericVector &mu,
                   const NumericVector &sig, const NumericVector &first_window, const LogicalVector &valid_window,
                   const NumericVector &mp, const IntegerVector &mpi)
      : data(data), window_size(window_size), compute_order(compute_order), df(df), dg(dg), mu(mu), sig(sig),
        first_window(first_window), valid_window(valid_window), mp(mp), mpi(mpi) {}

  void operator()(std::size_t begin, std::size_t end) override {
    uint32_t const profile_len = mp.size();
    std::vector<double> centered_window(window_size);
    std::vector<double> local_mp(profile_len, R_NegInf);
    std::vector<int> local_mpi(profile_len, NA_INTEGER);

    double const *const data_ptr = data.begin();
    int const *const compute_order_ptr = compute_order.begin();
    double const *const df_ptr = df.begin();
    double const *const dg_ptr = dg.begin();
    double const *const mu_ptr = mu.begin();
    double const *const sig_ptr = sig.begin();
    double const *const first_window_ptr = first_window.begin();
    int const *const valid_window_ptr = valid_window.begin();
    double *const mp_ptr = mp.begin();
    int *const mpi_ptr = mpi.begin();

    for (std::size_t order_index = begin; order_index < end; order_index++) {
      uint32_t const diag = compute_order_ptr[order_index];
      for (uint64_t i = 0; i < window_size; i++) {
        centered_window[i] = data_ptr[diag + i] - mu_ptr[diag];
      }

      double covariance =
          std::inner_product(centered_window.begin(), centered_window.end(), first_window_ptr, 0.0);
      uint32_t const diagonal_length = profile_len - diag;

      for (uint32_t offset = 0; offset < diagonal_length; offset++) {
        uint32_t const off_diag = offset + diag;
        covariance = covariance + df_ptr[offset] * dg_ptr[off_diag] + df_ptr[off_diag] * dg_ptr[offset];

        if (!valid_window_ptr[offset] || !valid_window_ptr[off_diag]) {
          continue;
        }

        double const correlation = covariance * sig_ptr[offset] * sig_ptr[off_diag];
        if (!std::isfinite(correlation)) {
          continue;
        }

        if (correlation > local_mp[offset]) {
          local_mp[offset] = correlation;
          local_mpi[offset] = off_diag + 1;
        }
        if (correlation > local_mp[off_diag]) {
          local_mp[off_diag] = correlation;
          local_mpi[off_diag] = offset + 1;
        }
      }
    }

    std::lock_guard<decltype(mutex)> lock(mutex);
    for (uint32_t i = 0; i < profile_len; i++) {
      if (local_mp[i] > mp_ptr[i]) {
        mp_ptr[i] = local_mp[i];
        mpi_ptr[i] = local_mpi[i];
      }
    }
  }
};

// [[Rcpp::export]]
List mpx_rcpp_parallel(NumericVector data_ref, uint64_t window_size, double ez, double s_size, bool idxs,
                       bool euclidean, bool progress) {

  uint64_t const exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

  try {
    // matrix profile using cross correlation,
    bool partial = false;
    uint32_t const n = data_ref.length();

    List msd = muinvn_rcpp(data_ref, window_size);

    NumericVector mu = msd["avg"];
    NumericVector const sig = msd["sig"];

    const uint32_t profile_len = n - window_size + 1;
    NumericVector mp(profile_len, -1.0);
    IntegerVector const mpi(profile_len, -1);

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector df = 0.5 * (data_ref[Range(window_size, n - 1)] - data_ref[Range(0, n - window_size - 1)]);
    df.push_front(0);
    NumericVector dg = (data_ref[Range(window_size, n - 1)] - mu[Range(1, profile_len - 1)]) +
                       (data_ref[Range(0, n - window_size - 1)] - mu[Range(0, n - window_size - 1)]);
    dg.push_front(0);

    NumericVector const ww = (data_ref[Range(0, window_size - 1)] - mu[0]);

    IntegerVector compute_order = Range(exclusion_zone, profile_len - 1);
    compute_order = sample(compute_order, compute_order.size());

    Progress p(100, progress);

    uint64_t work_size = compute_order.size();
    if (s_size > 0.0 && s_size < 1.0) {
      work_size = static_cast<uint64_t>(round(work_size * s_size + DBL_EPSILON));
      partial = work_size < static_cast<uint64_t>(compute_order.size());
    }

    MatrixProfileP matrix_profile(data_ref, window_size, compute_order, df, dg, mu, sig, ww, mp, mpi);

    try {
#if RCPP_PARALLEL_USE_TBB
      RcppParallel::parallelFor(0, work_size, matrix_profile, 4 * window_size);
#else
      RcppParallel2::ttParallelFor(0, work_size, matrix_profile, 4 * window_size);
#endif
      p.increment(100);
    } catch (RcppThread::UserInterruptException &e) {
      partial = true;
      Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
    } catch (...) {
      Rcpp::stop("c++ exception (unknown reason)");
    }

    // to do ed
    mp[mp > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      mp = sqrt(2 * window_size * (1 - mp));
    }

    if (idxs) {
      return (List::create(Rcpp::Named("matrix_profile") = mp, Rcpp::Named("profile_index") = mpi,
                           Rcpp::Named("partial") = partial));
    } else {
      return (List::create(Rcpp::Named("matrix_profile") = mp, Rcpp::Named("partial") = partial));
    }
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }
}

// Parallel MPX self-join with explicit support for non-finite samples.
// [[Rcpp::export]]
List mpx_na_rcpp_parallel(NumericVector data_ref, uint64_t window_size, double ez, double s_size, bool idxs,
                          bool euclidean, bool progress) {
  uint64_t const data_size = data_ref.length();
  if (window_size < 2 || window_size >= data_size) {
    Rcpp::stop("window_size must leave at least two subsequences");
  }
  if (!std::isfinite(ez) || ez < 0.0) {
    Rcpp::stop("ez must be finite and non-negative");
  }
  if (!std::isfinite(s_size) || s_size < 0.0 || s_size > 1.0) {
    Rcpp::stop("s_size must be between 0 and 1");
  }

  uint64_t const exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;
  uint32_t const n = data_size;
  uint32_t const profile_len = data_size - window_size + 1;
  bool partial = false;

  List const stats = muinvn_na_parallel(data_ref, window_size);
  NumericVector data = stats["data"];
  NumericVector mu = stats["avg"];
  NumericVector sig = stats["sig"];
  LogicalVector valid_window = stats["valid_window"];

  NumericVector mp(profile_len, R_NegInf);
  IntegerVector mpi(profile_len, NA_INTEGER);

  NumericVector df = 0.5 * (data[Range(window_size, n - 1)] - data[Range(0, n - window_size - 1)]);
  df.push_front(0);
  NumericVector dg = (data[Range(window_size, n - 1)] - mu[Range(1, profile_len - 1)]) +
                     (data[Range(0, n - window_size - 1)] - mu[Range(0, profile_len - 2)]);
  dg.push_front(0);
  NumericVector const first_window = data[Range(0, window_size - 1)] - mu[0];

  IntegerVector compute_order;
  if (exclusion_zone < profile_len) {
    compute_order = Range(exclusion_zone, profile_len - 1);
    compute_order = sample(compute_order, compute_order.size());
  }

  uint64_t work_size = compute_order.size();
  if (s_size > 0.0 && s_size < 1.0) {
    work_size = static_cast<uint64_t>(round(work_size * s_size + DBL_EPSILON));
    partial = work_size < static_cast<uint64_t>(compute_order.size());
  }

  Progress progress_bar(100, progress);
  MatrixProfilePNA matrix_profile(data, window_size, compute_order, df, dg, mu, sig, first_window, valid_window, mp,
                                  mpi);
  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, work_size, matrix_profile, 4 * window_size);
#else
    RcppParallel2::ttParallelFor(0, work_size, matrix_profile, 4 * window_size);
#endif
    progress_bar.increment(100);
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  double *const mp_ptr = mp.begin();
  int *const mpi_ptr = mpi.begin();
  int const *const valid_window_ptr = valid_window.begin();
  for (uint32_t i = 0; i < profile_len; i++) {
    if (!valid_window_ptr[i] || !std::isfinite(mp_ptr[i])) {
      mp_ptr[i] = NA_REAL;
      mpi_ptr[i] = NA_INTEGER;
      continue;
    }
    mp_ptr[i] = std::max(-1.0, std::min(1.0, mp_ptr[i]));
    if (euclidean) {
      mp_ptr[i] = sqrt(std::max(0.0, 2.0 * window_size * (1.0 - mp_ptr[i])));
    }
  }

  if (idxs) {
    return List::create(Rcpp::Named("matrix_profile") = mp, Rcpp::Named("profile_index") = mpi,
                        Rcpp::Named("valid_window") = valid_window, Rcpp::Named("partial") = partial);
  }
  return List::create(Rcpp::Named("matrix_profile") = mp, Rcpp::Named("valid_window") = valid_window,
                      Rcpp::Named("partial") = partial);
}

// ##### Parallel version #####

struct MatrixProfilePAB : public Worker {

private:
  // input
  const RVector<double> data_ref;
  const RVector<double> query_ref;
  const uint64_t window_size;
  const RVector<double> df_a;
  const RVector<double> df_b;
  const RVector<double> dg_a;
  const RVector<double> dg_b;
  const RVector<double> mu_a;
  const RVector<double> mu_b;
  const RVector<double> sig_a;
  const RVector<double> sig_b;
  const RVector<double> ww_a;
  const RVector<double> ww_b;
  const RVector<int> compute_ordera;
  const RVector<int> compute_orderb;
  // output
  RVector<double> mp_a;
  RVector<double> mp_b;
  RVector<int> mpi_a;
  RVector<int> mpi_b;

  // AB == 0, BA == 1
  uint8_t ab_ba = 0;

#if RCPP_PARALLEL_USE_TBB
  tbb::spin_mutex m;
  tbb::mutex m1;
#else
  tthread::mutex m;
  tthread::mutex m1;
#endif

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
public:
  MatrixProfilePAB(const NumericVector &data_ref, const NumericVector &query_ref, const uint64_t window_size,
                   const NumericVector &df_a, const NumericVector &df_b, const NumericVector &dg_a,
                   const NumericVector &dg_b, const NumericVector &mu_a, const NumericVector &mu_b,
                   const NumericVector &sig_a, const NumericVector &sig_b, const NumericVector &ww_a,
                   const NumericVector &ww_b, const IntegerVector &compute_ordera, const IntegerVector &compute_orderb,
                   const NumericVector &mp_a, const NumericVector &mp_b, const IntegerVector &mpi_a,
                   const IntegerVector &mpi_b)
      : data_ref(data_ref), query_ref(query_ref), window_size(window_size), df_a(df_a), df_b(df_b), dg_a(dg_a),
        dg_b(dg_b), mu_a(mu_a), mu_b(mu_b), sig_a(sig_a), sig_b(sig_b), ww_a(ww_a), ww_b(ww_b),
        compute_ordera(compute_ordera), compute_orderb(compute_orderb), mp_a(mp_a), mp_b(mp_b), mpi_a(mpi_a),
        mpi_b(mpi_b) {}

  void set_ab() { this->ab_ba = 0; }
  void set_ba() { this->ab_ba = 1; }

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) override { // exclusion_zone:profile_len
    double c = 0.0, c_cmp = 0.0;
    uint32_t off_max = 0UL, off_diag = 0UL, offset = 0UL;
    uint32_t const a_len = data_ref.length();
    uint32_t const b_len = query_ref.length();
    std::vector<double> inn(window_size);
    std::vector<double> mpp_a(mp_a.size(), -1.0);
    std::vector<int> mpip_a(mp_a.size(), -1);
    std::vector<double> mpp_b(mp_b.size(), -1.0);
    std::vector<int> mpip_b(mp_b.size(), -1);

    if (ab_ba == 0) {

        for (std::size_t diag = begin; diag < end; diag++) {

          for (uint64_t i = 0; i < window_size; i++) {
            inn[i] = data_ref[diag + i] - mu_a[diag];
          }

          off_max = MIN(a_len - window_size - diag + 1, b_len - window_size + 1);
          c = std::inner_product(inn.begin(), inn.end(), ww_b.begin(), 0.0);

          // off_max = (a_len - window_size - diag + 1);

          for (offset = 0; offset < off_max; offset++) {
            off_diag = offset + diag;
            c = c + df_a[off_diag] * dg_b[offset] + dg_a[off_diag] * df_b[offset];
            c_cmp = c * sig_b[offset] * sig_a[off_diag];

            if (c_cmp > mpp_b[offset]) {
              mpp_b[offset] = c_cmp;
              mpip_b[offset] = off_diag + 1;
            }
            if (c_cmp > mpp_a[off_diag]) {
              mpp_a[off_diag] = c_cmp;
              mpip_a[off_diag] = offset + 1;
            }
          }
        }
      } else {
        for (std::size_t diag = begin; diag < end; diag++) {

          for (uint64_t i = 0; i < window_size; i++) {
            inn[i] = query_ref[diag + i] - mu_b[diag];
          }

          off_max = MIN(b_len - window_size - diag + 1, a_len - window_size + 1);

          c = std::inner_product(inn.begin(), inn.end(), ww_a.begin(), 0.0);

          // off_max = (b_len - window_size - diag + 1);

          for (offset = 0; offset < off_max; offset++) {
            off_diag = offset + diag;
            c = c + df_b[off_diag] * dg_a[offset] + dg_b[off_diag] * df_a[offset];
            c_cmp = c * sig_a[offset] * sig_b[off_diag];

            if (c_cmp > mpp_a[offset]) {
              mpp_a[offset] = c_cmp;
              mpip_a[offset] = off_diag + 1;
            }
            if (c_cmp > mpp_b[off_diag]) {
              mpp_b[off_diag] = c_cmp;
              mpip_b[off_diag] = offset + 1;
            }
          }
        }
      }

      std::lock_guard<decltype(m1)> lock(m1);
      for (uint32_t i = 0; i < mp_a.size(); i++) {
        if (mpp_a[i] > mp_a[i]) {
          mp_a[i] = mpp_a[i];
          mpi_a[i] = mpip_a[i];
        }
      }
      for (uint32_t i = 0; i < mp_b.size(); i++) {
        if (mpp_b[i] > mp_b[i]) {
          mp_b[i] = mpp_b[i];
          mpi_b[i] = mpip_b[i];
        }
    }
  }
};

struct MatrixProfilePABNA : public Worker {
private:
  const RVector<double> data;
  const RVector<double> query;
  const uint64_t window_size;
  const RVector<double> df_a;
  const RVector<double> df_b;
  const RVector<double> dg_a;
  const RVector<double> dg_b;
  const RVector<double> mu_a;
  const RVector<double> mu_b;
  const RVector<double> sig_a;
  const RVector<double> sig_b;
  const RVector<double> first_a;
  const RVector<double> first_b;
  const RVector<int> valid_a;
  const RVector<int> valid_b;
  const RVector<int> order_a;
  const RVector<int> order_b;
  RVector<double> mp_a;
  RVector<double> mp_b;
  RVector<int> mpi_a;
  RVector<int> mpi_b;
  bool ba = false;

#if RCPP_PARALLEL_USE_TBB
  tbb::spin_mutex mutex;
#else
  tthread::mutex mutex;
#endif

public:
  MatrixProfilePABNA(const NumericVector &data, const NumericVector &query, uint64_t window_size,
                     const NumericVector &df_a, const NumericVector &df_b, const NumericVector &dg_a,
                     const NumericVector &dg_b, const NumericVector &mu_a, const NumericVector &mu_b,
                     const NumericVector &sig_a, const NumericVector &sig_b, const NumericVector &first_a,
                     const NumericVector &first_b, const LogicalVector &valid_a, const LogicalVector &valid_b,
                     const IntegerVector &order_a, const IntegerVector &order_b, const NumericVector &mp_a,
                     const NumericVector &mp_b, const IntegerVector &mpi_a, const IntegerVector &mpi_b)
      : data(data), query(query), window_size(window_size), df_a(df_a), df_b(df_b), dg_a(dg_a), dg_b(dg_b),
        mu_a(mu_a), mu_b(mu_b), sig_a(sig_a), sig_b(sig_b), first_a(first_a), first_b(first_b), valid_a(valid_a),
        valid_b(valid_b), order_a(order_a), order_b(order_b), mp_a(mp_a), mp_b(mp_b), mpi_a(mpi_a), mpi_b(mpi_b) {}

  void set_ba() { ba = true; }

  void operator()(std::size_t begin, std::size_t end) override {
    uint32_t const profile_len_a = mp_a.size();
    uint32_t const profile_len_b = mp_b.size();
    std::vector<double> centered_window(window_size);
    std::vector<double> local_mp_a(profile_len_a, R_NegInf);
    std::vector<double> local_mp_b(profile_len_b, R_NegInf);
    std::vector<int> local_mpi_a(profile_len_a, NA_INTEGER);
    std::vector<int> local_mpi_b(profile_len_b, NA_INTEGER);

    double const *const data_ptr = data.begin();
    double const *const query_ptr = query.begin();
    double const *const df_a_ptr = df_a.begin();
    double const *const df_b_ptr = df_b.begin();
    double const *const dg_a_ptr = dg_a.begin();
    double const *const dg_b_ptr = dg_b.begin();
    double const *const mu_a_ptr = mu_a.begin();
    double const *const mu_b_ptr = mu_b.begin();
    double const *const sig_a_ptr = sig_a.begin();
    double const *const sig_b_ptr = sig_b.begin();
    double const *const first_a_ptr = first_a.begin();
    double const *const first_b_ptr = first_b.begin();
    int const *const valid_a_ptr = valid_a.begin();
    int const *const valid_b_ptr = valid_b.begin();
    int const *const order_a_ptr = order_a.begin();
    int const *const order_b_ptr = order_b.begin();
    double *const mp_a_ptr = mp_a.begin();
    double *const mp_b_ptr = mp_b.begin();
    int *const mpi_a_ptr = mpi_a.begin();
    int *const mpi_b_ptr = mpi_b.begin();

    if (!ba) {
      for (std::size_t order_index = begin; order_index < end; order_index++) {
        uint32_t const diag = order_a_ptr[order_index];
        for (uint64_t i = 0; i < window_size; i++) {
          centered_window[i] = data_ptr[diag + i] - mu_a_ptr[diag];
        }
        double covariance = std::inner_product(centered_window.begin(), centered_window.end(), first_b_ptr, 0.0);
        uint32_t const diagonal_length = MIN(profile_len_a - diag, profile_len_b);

        for (uint32_t offset = 0; offset < diagonal_length; offset++) {
          uint32_t const off_diag = offset + diag;
          covariance =
              covariance + df_a_ptr[off_diag] * dg_b_ptr[offset] + dg_a_ptr[off_diag] * df_b_ptr[offset];
          if (!valid_a_ptr[off_diag] || !valid_b_ptr[offset]) {
            continue;
          }
          double const correlation = covariance * sig_a_ptr[off_diag] * sig_b_ptr[offset];
          if (!std::isfinite(correlation)) {
            continue;
          }
          if (correlation > local_mp_a[off_diag]) {
            local_mp_a[off_diag] = correlation;
            local_mpi_a[off_diag] = offset + 1;
          }
          if (correlation > local_mp_b[offset]) {
            local_mp_b[offset] = correlation;
            local_mpi_b[offset] = off_diag + 1;
          }
        }
      }
    } else {
      for (std::size_t order_index = begin; order_index < end; order_index++) {
        uint32_t const diag = order_b_ptr[order_index];
        for (uint64_t i = 0; i < window_size; i++) {
          centered_window[i] = query_ptr[diag + i] - mu_b_ptr[diag];
        }
        double covariance = std::inner_product(centered_window.begin(), centered_window.end(), first_a_ptr, 0.0);
        uint32_t const diagonal_length = MIN(profile_len_b - diag, profile_len_a);

        for (uint32_t offset = 0; offset < diagonal_length; offset++) {
          uint32_t const off_diag = offset + diag;
          covariance =
              covariance + df_b_ptr[off_diag] * dg_a_ptr[offset] + dg_b_ptr[off_diag] * df_a_ptr[offset];
          if (!valid_b_ptr[off_diag] || !valid_a_ptr[offset]) {
            continue;
          }
          double const correlation = covariance * sig_b_ptr[off_diag] * sig_a_ptr[offset];
          if (!std::isfinite(correlation)) {
            continue;
          }
          if (correlation > local_mp_b[off_diag]) {
            local_mp_b[off_diag] = correlation;
            local_mpi_b[off_diag] = offset + 1;
          }
          if (correlation > local_mp_a[offset]) {
            local_mp_a[offset] = correlation;
            local_mpi_a[offset] = off_diag + 1;
          }
        }
      }
    }

    std::lock_guard<decltype(mutex)> lock(mutex);
    for (uint32_t i = 0; i < profile_len_a; i++) {
      if (local_mp_a[i] > mp_a_ptr[i]) {
        mp_a_ptr[i] = local_mp_a[i];
        mpi_a_ptr[i] = local_mpi_a[i];
      }
    }
    for (uint32_t i = 0; i < profile_len_b; i++) {
      if (local_mp_b[i] > mp_b_ptr[i]) {
        mp_b_ptr[i] = local_mp_b[i];
        mpi_b_ptr[i] = local_mpi_b[i];
      }
    }
  }
};

// [[Rcpp::export]]
List mpxab_rcpp_parallel(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double s_size,
                         bool idxs, bool euclidean, bool progress) {

  try {
    // matrix profile using cross correlation,
    bool partial = false;
    uint32_t a_len = data_ref.length();
    uint32_t b_len = query_ref.length();

    List msd_a = muinvn_rcpp_parallel(data_ref, window_size);
    List msd_b = muinvn_rcpp_parallel(query_ref, window_size);

    NumericVector mu_a = msd_a["avg"];
    NumericVector sig_a = msd_a["sig"];
    NumericVector mu_b = msd_b["avg"];
    NumericVector sig_b = msd_b["sig"];

    uint32_t const profile_len_a = a_len - window_size + 1;
    uint32_t const profile_len_b = b_len - window_size + 1;

    NumericVector mp_a(profile_len_a, -1.0);
    NumericVector mp_b(profile_len_b, -1.0);

    IntegerVector mpi_a(profile_len_a, -1);
    IntegerVector mpi_b(profile_len_b, -1);

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.
    NumericVector df_a = 0.5 * (data_ref[Range(window_size, a_len - 1)] - data_ref[Range(0, a_len - window_size - 1)]);
    df_a.push_front(0);
    NumericVector dg_a = (data_ref[Range(window_size, a_len - 1)] - mu_a[Range(1, profile_len_a - 1)]) +
                         (data_ref[Range(0, a_len - window_size - 1)] - mu_a[Range(0, a_len - window_size - 1)]);
    dg_a.push_front(0);
    NumericVector df_b =
        0.5 * (query_ref[Range(window_size, b_len - 1)] - query_ref[Range(0, b_len - window_size - 1)]);
    df_b.push_front(0);
    NumericVector dg_b = (query_ref[Range(window_size, b_len - 1)] - mu_b[Range(1, profile_len_b - 1)]) +
                         (query_ref[Range(0, b_len - window_size - 1)] - mu_b[Range(0, b_len - window_size - 1)]);
    dg_b.push_front(0);

    NumericVector ww_a = (data_ref[Range(0, window_size - 1)] - mu_a[0]);
    NumericVector ww_b = (query_ref[Range(0, window_size - 1)] - mu_b[0]);

    Progress p(100, progress);

    IntegerVector compute_ordera = Range(0, profile_len_a - 1);
    compute_ordera = sample(compute_ordera, compute_ordera.size());
    IntegerVector compute_orderb = Range(0, profile_len_b - 1);
    compute_orderb = sample(compute_orderb, compute_orderb.size());

    MatrixProfilePAB matrix_profile(data_ref, query_ref, window_size, df_a, df_b, dg_a, dg_b, mu_a, mu_b, sig_a, sig_b,
                                    ww_a, ww_b, compute_ordera, compute_orderb, mp_a, mp_b, mpi_a, mpi_b);

    try {
#if RCPP_PARALLEL_USE_TBB
      RcppParallel::parallelFor(0, profile_len_a, matrix_profile, 4 * window_size);
#else
      RcppParallel2::ttParallelFor(0, profile_len_a, matrix_profile, 4 * window_size);
#endif
      p.increment(50);
    } catch (RcppThread::UserInterruptException &e) {
      partial = true;
      Rcout << "Process AB terminated by the user successfully, partial results were returned." << std::endl;
    } catch (...) {
      Rcpp::stop("c++ exception (unknown reason)");
    }

    // switch worker to BA
    matrix_profile.set_ba();

    try {
#if RCPP_PARALLEL_USE_TBB
      RcppParallel::parallelFor(1, profile_len_b, matrix_profile, 4 * window_size);
#else
      RcppParallel2::ttParallelFor(1, profile_len_b, matrix_profile, 4 * window_size);
#endif
      p.increment(50);
    } catch (RcppThread::UserInterruptException &e) {
      partial = true;
      Rcout << "Process BA terminated by the user successfully, partial results were returned." << std::endl;
    } catch (...) {
      Rcpp::stop("c++ exception (unknown reason)");
    }

    // to do ed
    mp_a[mp_a > 1.0] = 1.0;
    mp_b[mp_b > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      mp_a = sqrt(2 * window_size * (1 - mp_a));
      mp_b = sqrt(2 * window_size * (1 - mp_b));
    }

    if (idxs) {
      return (List::create(Rcpp::Named("matrix_profile") = mp_a, Rcpp::Named("profile_index") = mpi_a,
                           Rcpp::Named("mpb") = mp_b, Rcpp::Named("pib") = mpi_b, Rcpp::Named("partial") = partial));
    } else {
      return (List::create(Rcpp::Named("matrix_profile") = mp_a, Rcpp::Named("mpb") = mp_b,
                           Rcpp::Named("partial") = partial));
    }
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }
}

// Parallel MPX AB-join with explicit support for non-finite samples.
// [[Rcpp::export]]
List mpxab_na_rcpp_parallel(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double s_size,
                            bool idxs, bool euclidean, bool progress) {
  uint64_t const data_size = data_ref.length();
  uint64_t const query_size = query_ref.length();
  if (window_size < 2 || window_size >= data_size || window_size >= query_size) {
    Rcpp::stop("window_size must leave at least two subsequences in each series");
  }
  if (!std::isfinite(s_size) || s_size < 0.0 || s_size > 1.0) {
    Rcpp::stop("s_size must be between 0 and 1");
  }

  uint32_t const profile_len_a = data_size - window_size + 1;
  uint32_t const profile_len_b = query_size - window_size + 1;
  bool partial = false;

  List const stats_a = muinvn_na_parallel(data_ref, window_size);
  List const stats_b = muinvn_na_parallel(query_ref, window_size);
  NumericVector data = stats_a["data"];
  NumericVector query = stats_b["data"];
  NumericVector mu_a = stats_a["avg"];
  NumericVector mu_b = stats_b["avg"];
  NumericVector sig_a = stats_a["sig"];
  NumericVector sig_b = stats_b["sig"];
  LogicalVector valid_a = stats_a["valid_window"];
  LogicalVector valid_b = stats_b["valid_window"];

  NumericVector mp_a(profile_len_a, R_NegInf);
  NumericVector mp_b(profile_len_b, R_NegInf);
  IntegerVector mpi_a(profile_len_a, NA_INTEGER);
  IntegerVector mpi_b(profile_len_b, NA_INTEGER);

  NumericVector df_a = 0.5 * (data[Range(window_size, data_size - 1)] -
                              data[Range(0, data_size - window_size - 1)]);
  df_a.push_front(0);
  NumericVector dg_a = (data[Range(window_size, data_size - 1)] - mu_a[Range(1, profile_len_a - 1)]) +
                       (data[Range(0, data_size - window_size - 1)] - mu_a[Range(0, profile_len_a - 2)]);
  dg_a.push_front(0);

  NumericVector df_b = 0.5 * (query[Range(window_size, query_size - 1)] -
                              query[Range(0, query_size - window_size - 1)]);
  df_b.push_front(0);
  NumericVector dg_b = (query[Range(window_size, query_size - 1)] - mu_b[Range(1, profile_len_b - 1)]) +
                       (query[Range(0, query_size - window_size - 1)] - mu_b[Range(0, profile_len_b - 2)]);
  dg_b.push_front(0);

  NumericVector const first_a = data[Range(0, window_size - 1)] - mu_a[0];
  NumericVector const first_b = query[Range(0, window_size - 1)] - mu_b[0];
  IntegerVector order_a = Range(0, profile_len_a - 1);
  IntegerVector order_b = Range(0, profile_len_b - 1);
  order_a = sample(order_a, order_a.size());
  order_b = sample(order_b, order_b.size());

  uint64_t work_a = order_a.size();
  uint64_t work_b = order_b.size();
  if (s_size > 0.0 && s_size < 1.0) {
    work_a = static_cast<uint64_t>(round(work_a * s_size + DBL_EPSILON));
    work_b = static_cast<uint64_t>(round(work_b * s_size + DBL_EPSILON));
    partial = work_a < static_cast<uint64_t>(order_a.size()) || work_b < static_cast<uint64_t>(order_b.size());
  }

  Progress progress_bar(100, progress);
  MatrixProfilePABNA matrix_profile(data, query, window_size, df_a, df_b, dg_a, dg_b, mu_a, mu_b, sig_a, sig_b,
                                    first_a, first_b, valid_a, valid_b, order_a, order_b, mp_a, mp_b, mpi_a, mpi_b);

  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, work_a, matrix_profile, 4 * window_size);
#else
    RcppParallel2::ttParallelFor(0, work_a, matrix_profile, 4 * window_size);
#endif
    progress_bar.increment(50);
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process AB terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  matrix_profile.set_ba();
  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, work_b, matrix_profile, 4 * window_size);
#else
    RcppParallel2::ttParallelFor(0, work_b, matrix_profile, 4 * window_size);
#endif
    progress_bar.increment(50);
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process BA terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  double *const mp_a_ptr = mp_a.begin();
  double *const mp_b_ptr = mp_b.begin();
  int *const mpi_a_ptr = mpi_a.begin();
  int *const mpi_b_ptr = mpi_b.begin();
  int const *const valid_a_ptr = valid_a.begin();
  int const *const valid_b_ptr = valid_b.begin();
  for (uint32_t i = 0; i < profile_len_a; i++) {
    if (!valid_a_ptr[i] || !std::isfinite(mp_a_ptr[i])) {
      mp_a_ptr[i] = NA_REAL;
      mpi_a_ptr[i] = NA_INTEGER;
      continue;
    }
    mp_a_ptr[i] = std::max(-1.0, std::min(1.0, mp_a_ptr[i]));
    if (euclidean) {
      mp_a_ptr[i] = sqrt(std::max(0.0, 2.0 * window_size * (1.0 - mp_a_ptr[i])));
    }
  }
  for (uint32_t i = 0; i < profile_len_b; i++) {
    if (!valid_b_ptr[i] || !std::isfinite(mp_b_ptr[i])) {
      mp_b_ptr[i] = NA_REAL;
      mpi_b_ptr[i] = NA_INTEGER;
      continue;
    }
    mp_b_ptr[i] = std::max(-1.0, std::min(1.0, mp_b_ptr[i]));
    if (euclidean) {
      mp_b_ptr[i] = sqrt(std::max(0.0, 2.0 * window_size * (1.0 - mp_b_ptr[i])));
    }
  }

  if (idxs) {
    return List::create(Rcpp::Named("matrix_profile") = mp_a, Rcpp::Named("profile_index") = mpi_a,
                        Rcpp::Named("mpb") = mp_b, Rcpp::Named("pib") = mpi_b,
                        Rcpp::Named("valid_window_a") = valid_a, Rcpp::Named("valid_window_b") = valid_b,
                        Rcpp::Named("partial") = partial);
  }
  return List::create(Rcpp::Named("matrix_profile") = mp_a, Rcpp::Named("mpb") = mp_b,
                      Rcpp::Named("valid_window_a") = valid_a, Rcpp::Named("valid_window_b") = valid_b,
                      Rcpp::Named("partial") = partial);
}
