// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#include "mathtools.h" // math first to fix OSX error
#include "mpx_parallel.h"
#include "windowfunc.h"
#include <cfloat> // DBL_EPSILON when STRICT_R_HEADERS
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <thread>
#include <unordered_set>
#include <vector>

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

// Worker instrumentation is deliberately opt-in.  The hot loops only update
// local counters when MATRIXPROFILER_PROFILE_WORKERS is enabled; aggregation
// happens once per scheduler task, outside the matrix-profile arithmetic.
struct MatrixProfileWorkerTaskStats {
  uint64_t diagonals = 0;
  uint64_t pairs_visited = 0;
  uint64_t invalid_pairs = 0;
  uint64_t nonfinite_pairs = 0;
  uint64_t scratch_initialized = 0;
  uint64_t scratch_reset = 0;
  uint64_t merge_values = 0;
  uint64_t scratch_ns = 0;
  uint64_t compute_ns = 0;
  uint64_t mutex_wait_ns = 0;
  uint64_t merge_ns = 0;
  uint64_t task_ns = 0;
};

struct MatrixProfileWorkerPhaseStats {
  uint64_t tasks = 0;
  uint64_t diagonals = 0;
  uint64_t pairs_visited = 0;
  uint64_t invalid_pairs = 0;
  uint64_t nonfinite_pairs = 0;
  uint64_t scratch_initialized = 0;
  uint64_t scratch_reset = 0;
  uint64_t merge_values = 0;
  uint64_t scratch_ns = 0;
  uint64_t compute_ns = 0;
  uint64_t mutex_wait_ns = 0;
  uint64_t merge_ns = 0;
  uint64_t task_ns = 0;
  std::unordered_set<std::size_t> thread_ids;
};

struct MatrixProfileWorkerDiagnostics {
  const bool enabled;
  mutable std::mutex mutex;
  MatrixProfileWorkerPhaseStats ab;
  MatrixProfileWorkerPhaseStats ba;

  explicit MatrixProfileWorkerDiagnostics(bool enabled_) : enabled(enabled_) {}

  void record(bool is_ba, const MatrixProfileWorkerTaskStats &task) {
    if (!enabled) {
      return;
    }
    std::lock_guard<std::mutex> lock(mutex);
    MatrixProfileWorkerPhaseStats &phase = is_ba ? ba : ab;
    phase.tasks++;
    phase.diagonals += task.diagonals;
    phase.pairs_visited += task.pairs_visited;
    phase.invalid_pairs += task.invalid_pairs;
    phase.nonfinite_pairs += task.nonfinite_pairs;
    phase.scratch_initialized += task.scratch_initialized;
    phase.scratch_reset += task.scratch_reset;
    phase.merge_values += task.merge_values;
    phase.scratch_ns += task.scratch_ns;
    phase.compute_ns += task.compute_ns;
    phase.mutex_wait_ns += task.mutex_wait_ns;
    phase.merge_ns += task.merge_ns;
    phase.task_ns += task.task_ns;
    phase.thread_ids.insert(std::hash<std::thread::id>{}(std::this_thread::get_id()));
  }

  static Rcpp::List phase_list(const MatrixProfileWorkerPhaseStats &phase) {
    const double billion = 1000000000.0;
    return Rcpp::List::create(
        Rcpp::Named("tasks") = static_cast<double>(phase.tasks),
        Rcpp::Named("diagonals") = static_cast<double>(phase.diagonals),
        Rcpp::Named("pairs_visited") = static_cast<double>(phase.pairs_visited),
        Rcpp::Named("invalid_pairs") = static_cast<double>(phase.invalid_pairs),
        Rcpp::Named("nonfinite_pairs") = static_cast<double>(phase.nonfinite_pairs),
        Rcpp::Named("scratch_initialized") = static_cast<double>(phase.scratch_initialized),
        Rcpp::Named("scratch_reset") = static_cast<double>(phase.scratch_reset),
        Rcpp::Named("merge_values") = static_cast<double>(phase.merge_values),
        Rcpp::Named("scratch_seconds") = static_cast<double>(phase.scratch_ns) / billion,
        Rcpp::Named("compute_seconds") = static_cast<double>(phase.compute_ns) / billion,
        Rcpp::Named("mutex_wait_seconds") = static_cast<double>(phase.mutex_wait_ns) / billion,
        Rcpp::Named("merge_seconds") = static_cast<double>(phase.merge_ns) / billion,
        Rcpp::Named("task_seconds") = static_cast<double>(phase.task_ns) / billion,
        Rcpp::Named("threads_observed") = static_cast<int>(phase.thread_ids.size()));
  }

  Rcpp::List as_list() const {
    if (!enabled) {
      return Rcpp::List::create(Rcpp::Named("enabled") = false);
    }
    std::lock_guard<std::mutex> lock(mutex);
    return Rcpp::List::create(Rcpp::Named("enabled") = true, Rcpp::Named("ab") = phase_list(ab),
                              Rcpp::Named("ba") = phase_list(ba));
  }
};

static bool matrix_profile_worker_diagnostics_enabled() {
  const char *value = std::getenv("MATRIXPROFILER_PROFILE_WORKERS");
  if (value == nullptr || value[0] == '\0') {
    return false;
  }
  return value[0] == '1' || value[0] == 't' || value[0] == 'T' || value[0] == 'y' || value[0] == 'Y';
}

// Benchmark-only override for the invalid-block heuristic. The default is
// deliberately conservative; "legacy" disables block restarts and "force"
// (or numeric zero) selects them for every NA-aware call. This is read on the
// R thread before a worker is launched and is never consulted by a worker.
static double matrix_profile_invalid_block_threshold() {
  const double default_threshold = 0.01;
  const char *value = std::getenv("MATRIXPROFILER_INVALID_BLOCK_THRESHOLD");
  if (value == nullptr || value[0] == '\0' || std::strcmp(value, "auto") == 0) {
    return default_threshold;
  }
  if (std::strcmp(value, "legacy") == 0 || std::strcmp(value, "off") == 0 || std::strcmp(value, "none") == 0) {
    return -1.0;
  }
  if (std::strcmp(value, "force") == 0) {
    return 0.0;
  }
  char *end = nullptr;
  double threshold = std::strtod(value, &end);
  if (end == value || *end != '\0' || !std::isfinite(threshold) || threshold < 0.0 || threshold > 1.0) {
    return default_threshold;
  }
  return threshold;
}

using MatrixProfileClock = std::chrono::steady_clock;

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
  const bool keep_indices;

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
                 const NumericVector &ww, const NumericVector &mp, const IntegerVector &mpi, bool idxs)
      : data_ref(data_ref), window_size(window_size), compute_order(compute_order), df(df), dg(dg), mu(mmu), sig(sig),
        ww(ww), mp(mp), mpi(mpi), keep_indices(idxs) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) override { // exclusion_zone:profile_len
    double c = 0.0, c_cmp = 0.0;
    uint32_t off_max = 0UL, off_diag = 0UL, offset = 0UL;
    uint32_t const n = data_ref.length();
    std::vector<double> aa(window_size);

    std::vector<double> mpp(mp.size(), -1.0);
    std::vector<int> mpip;
    if (keep_indices) {
      mpip.assign(mp.size(), -1);
    }

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
            if (keep_indices) {
              mpip[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mpp[off_diag]) {
            mpp[off_diag] = c_cmp;
            if (keep_indices) {
              mpip[off_diag] = offset + 1;
            }
          }
        }
      }

      std::lock_guard<decltype(m)> lock(m);
      for (uint64_t i = 0; i < mp.size(); i++) {
        if (mpp[i] > mp[i]) {
          mp[i] = mpp[i];
          if (keep_indices) {
            mpi[i] = mpip[i];
          }
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
  const bool keep_indices;

#if RCPP_PARALLEL_USE_TBB
  tbb::spin_mutex mutex;
#else
  tthread::mutex mutex;
#endif

public:
  MatrixProfilePNA(const NumericVector &data, uint64_t window_size, const IntegerVector &compute_order,
                   const NumericVector &df, const NumericVector &dg, const NumericVector &mu,
                   const NumericVector &sig, const NumericVector &first_window, const LogicalVector &valid_window,
                   const NumericVector &mp, const IntegerVector &mpi, bool idxs)
      : data(data), window_size(window_size), compute_order(compute_order), df(df), dg(dg), mu(mu), sig(sig),
        first_window(first_window), valid_window(valid_window), mp(mp), mpi(mpi), keep_indices(idxs) {}

  void operator()(std::size_t begin, std::size_t end) override {
    uint32_t const profile_len = mp.size();
    std::vector<double> centered_window(window_size);
    std::vector<double> local_mp(profile_len, R_NegInf);
    std::vector<int> local_mpi;
    if (keep_indices) {
      local_mpi.assign(profile_len, NA_INTEGER);
    }

    double const *const data_ptr = data.begin();
    int const *const compute_order_ptr = compute_order.begin();
    double const *const df_ptr = df.begin();
    double const *const dg_ptr = dg.begin();
    double const *const mu_ptr = mu.begin();
    double const *const sig_ptr = sig.begin();
    double const *const first_window_ptr = first_window.begin();
    int const *const valid_window_ptr = valid_window.begin();
    double *const mp_ptr = mp.begin();
    int *const mpi_ptr = keep_indices ? mpi.begin() : nullptr;

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
          if (keep_indices) {
            local_mpi[offset] = off_diag + 1;
          }
        }
        if (correlation > local_mp[off_diag]) {
          local_mp[off_diag] = correlation;
          if (keep_indices) {
            local_mpi[off_diag] = offset + 1;
          }
        }
      }
    }

    std::lock_guard<decltype(mutex)> lock(mutex);
    for (uint32_t i = 0; i < profile_len; i++) {
      if (local_mp[i] > mp_ptr[i]) {
        mp_ptr[i] = local_mp[i];
        if (keep_indices) {
          mpi_ptr[i] = local_mpi[i];
        }
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
    IntegerVector mpi;
    if (idxs) {
      mpi = IntegerVector(profile_len, -1);
    }

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

    MatrixProfileP matrix_profile(data_ref, window_size, compute_order, df, dg, mu, sig, ww, mp, mpi, idxs);

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

  if (mpx_finite_fast_path_safe(data_ref, window_size)) {
    List result = mpx_rcpp_parallel(data_ref, window_size, ez, s_size, idxs, euclidean, progress);
    result["valid_window"] = LogicalVector(data_ref.length() - window_size + 1, true);
    return result;
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
  IntegerVector mpi;
  if (idxs) {
    mpi = IntegerVector(profile_len, NA_INTEGER);
  }

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
                                  mpi, idxs);
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
  int *const mpi_ptr = idxs ? mpi.begin() : nullptr;
  int const *const valid_window_ptr = valid_window.begin();
  for (uint32_t i = 0; i < profile_len; i++) {
    if (!valid_window_ptr[i] || !std::isfinite(mp_ptr[i])) {
      mp_ptr[i] = NA_REAL;
      if (idxs) {
        mpi_ptr[i] = NA_INTEGER;
      }
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
  const bool keep_indices;
  MatrixProfileWorkerDiagnostics *const diagnostics;

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
                   const IntegerVector &mpi_b, bool idxs, MatrixProfileWorkerDiagnostics *diagnostics_)
      : data_ref(data_ref), query_ref(query_ref), window_size(window_size), df_a(df_a), df_b(df_b), dg_a(dg_a),
        dg_b(dg_b), mu_a(mu_a), mu_b(mu_b), sig_a(sig_a), sig_b(sig_b), ww_a(ww_a), ww_b(ww_b),
        compute_ordera(compute_ordera), compute_orderb(compute_orderb), mp_a(mp_a), mp_b(mp_b), mpi_a(mpi_a),
        mpi_b(mpi_b), keep_indices(idxs), diagnostics(diagnostics_) {}

  void set_ab() { this->ab_ba = 0; }
  void set_ba() { this->ab_ba = 1; }

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) override { // exclusion_zone:profile_len
    const bool instrument = diagnostics != nullptr && diagnostics->enabled;
    MatrixProfileWorkerTaskStats task_stats;
    MatrixProfileClock::time_point task_started;
    MatrixProfileClock::time_point scratch_done;
    MatrixProfileClock::time_point compute_done;
    MatrixProfileClock::time_point lock_requested;
    MatrixProfileClock::time_point lock_acquired;
    MatrixProfileClock::time_point merge_done;
    MatrixProfileClock::time_point task_done;
    if (instrument) {
      task_started = MatrixProfileClock::now();
    }

    double c = 0.0, c_cmp = 0.0;
    uint32_t off_max = 0UL, off_diag = 0UL, offset = 0UL;
    uint32_t const a_len = data_ref.length();
    uint32_t const b_len = query_ref.length();
    uint32_t const profile_len_a = mp_a.size();
    uint32_t const profile_len_b = mp_b.size();

    // TBB may invoke this operator many times on the same executor thread.
    // Keep one set of scratch buffers per thread instead of allocating and
    // initializing full profiles for every scheduler task. Only the ranges
    // touched by the previous task are reset before reuse.
    thread_local std::vector<double> inn;
    thread_local std::vector<double> mpp_a;
    thread_local std::vector<int> mpip_a;
    thread_local std::vector<double> mpp_b;
    thread_local std::vector<int> mpip_b;
    thread_local uint32_t previous_a_begin = 0;
    thread_local uint32_t previous_a_end = 0;
    thread_local uint32_t previous_b_begin = 0;
    thread_local uint32_t previous_b_end = 0;
    thread_local bool previous_keep_indices = false;

    if (inn.size() != window_size) {
      inn.resize(window_size);
      if (instrument) {
        task_stats.scratch_initialized += window_size;
      }
    }
    if (mpp_a.size() != profile_len_a || mpp_b.size() != profile_len_b) {
      mpp_a.assign(profile_len_a, -1.0);
      mpp_b.assign(profile_len_b, -1.0);
      if (instrument) {
        task_stats.scratch_initialized += profile_len_a + profile_len_b;
      }
      if (keep_indices) {
        mpip_a.assign(profile_len_a, -1);
        mpip_b.assign(profile_len_b, -1);
        if (instrument) {
          task_stats.scratch_initialized += profile_len_a + profile_len_b;
        }
      }
      previous_a_begin = previous_a_end = previous_b_begin = previous_b_end = 0;
    } else {
      if (keep_indices && (mpip_a.size() != profile_len_a || mpip_b.size() != profile_len_b)) {
        mpip_a.assign(profile_len_a, -1);
        mpip_b.assign(profile_len_b, -1);
        if (instrument) {
          task_stats.scratch_initialized += profile_len_a + profile_len_b;
        }
        previous_a_begin = previous_a_end = previous_b_begin = previous_b_end = 0;
      } else if (keep_indices && !previous_keep_indices) {
        std::fill(mpip_a.begin(), mpip_a.end(), -1);
        std::fill(mpip_b.begin(), mpip_b.end(), -1);
        if (instrument) {
          task_stats.scratch_reset += profile_len_a + profile_len_b;
        }
      }
      if (previous_a_begin < previous_a_end) {
        std::fill(mpp_a.begin() + previous_a_begin, mpp_a.begin() + previous_a_end, -1.0);
        if (instrument) {
          task_stats.scratch_reset += previous_a_end - previous_a_begin;
        }
        if (keep_indices) {
          std::fill(mpip_a.begin() + previous_a_begin, mpip_a.begin() + previous_a_end, -1);
          if (instrument) {
            task_stats.scratch_reset += previous_a_end - previous_a_begin;
          }
        }
      }
      if (previous_b_begin < previous_b_end) {
        std::fill(mpp_b.begin() + previous_b_begin, mpp_b.begin() + previous_b_end, -1.0);
        if (instrument) {
          task_stats.scratch_reset += previous_b_end - previous_b_begin;
        }
        if (keep_indices) {
          std::fill(mpip_b.begin() + previous_b_begin, mpip_b.begin() + previous_b_end, -1);
          if (instrument) {
            task_stats.scratch_reset += previous_b_end - previous_b_begin;
          }
        }
      }
    }

    if (instrument) {
      scratch_done = MatrixProfileClock::now();
    }

    uint32_t current_a_begin = profile_len_a;
    uint32_t current_a_end = 0;
    uint32_t current_b_begin = profile_len_b;
    uint32_t current_b_end = 0;
    if (ab_ba == 0) {
      for (std::size_t order_index = begin; order_index < end; order_index++) {
        // The legacy normal worker processes this range in diagonal order;
        // keep the same order while deriving the touched output ranges.
        uint32_t const diag = order_index;
        uint32_t const diagonal_length = MIN(a_len - window_size - diag + 1, profile_len_b);
        current_a_begin = std::min(current_a_begin, diag);
        current_a_end = std::max(current_a_end, diag + diagonal_length);
        current_b_begin = 0;
        current_b_end = std::max(current_b_end, diagonal_length);
      }
    } else {
      for (std::size_t order_index = begin; order_index < end; order_index++) {
        uint32_t const diag = order_index;
        uint32_t const diagonal_length = MIN(b_len - window_size - diag + 1, profile_len_a);
        current_a_begin = 0;
        current_a_end = std::max(current_a_end, diagonal_length);
        current_b_begin = std::min(current_b_begin, diag);
        current_b_end = std::max(current_b_end, diag + diagonal_length);
      }
    }

    if (ab_ba == 0) {

        for (std::size_t diag = begin; diag < end; diag++) {

          for (uint64_t i = 0; i < window_size; i++) {
            inn[i] = data_ref[diag + i] - mu_a[diag];
          }

          off_max = MIN(a_len - window_size - diag + 1, b_len - window_size + 1);
          if (instrument) {
            task_stats.diagonals++;
            task_stats.pairs_visited += off_max;
          }
          c = std::inner_product(inn.begin(), inn.end(), ww_b.begin(), 0.0);

          // off_max = (a_len - window_size - diag + 1);

          for (offset = 0; offset < off_max; offset++) {
            off_diag = offset + diag;
            c = c + df_a[off_diag] * dg_b[offset] + dg_a[off_diag] * df_b[offset];
            c_cmp = c * sig_b[offset] * sig_a[off_diag];

            if (c_cmp > mpp_b[offset]) {
              mpp_b[offset] = c_cmp;
              if (keep_indices) {
                mpip_b[offset] = off_diag + 1;
              }
            }
            if (c_cmp > mpp_a[off_diag]) {
              mpp_a[off_diag] = c_cmp;
              if (keep_indices) {
                mpip_a[off_diag] = offset + 1;
              }
            }
          }
        }
      } else {
        for (std::size_t diag = begin; diag < end; diag++) {

          for (uint64_t i = 0; i < window_size; i++) {
            inn[i] = query_ref[diag + i] - mu_b[diag];
          }

          off_max = MIN(b_len - window_size - diag + 1, a_len - window_size + 1);
          if (instrument) {
            task_stats.diagonals++;
            task_stats.pairs_visited += off_max;
          }

          c = std::inner_product(inn.begin(), inn.end(), ww_a.begin(), 0.0);

          // off_max = (b_len - window_size - diag + 1);

          for (offset = 0; offset < off_max; offset++) {
            off_diag = offset + diag;
            c = c + df_b[off_diag] * dg_a[offset] + dg_b[off_diag] * df_a[offset];
            c_cmp = c * sig_a[offset] * sig_b[off_diag];

            if (c_cmp > mpp_a[offset]) {
              mpp_a[offset] = c_cmp;
              if (keep_indices) {
                mpip_a[offset] = off_diag + 1;
              }
            }
            if (c_cmp > mpp_b[off_diag]) {
              mpp_b[off_diag] = c_cmp;
              if (keep_indices) {
                mpip_b[off_diag] = offset + 1;
              }
            }
          }
        }
      }

      if (instrument) {
        compute_done = MatrixProfileClock::now();
        lock_requested = compute_done;
      }
      {
        std::lock_guard<decltype(m1)> lock(m1);
        if (instrument) {
          lock_acquired = MatrixProfileClock::now();
        }
        for (uint32_t i = current_a_begin; i < current_a_end; i++) {
          if (mpp_a[i] > mp_a[i]) {
            mp_a[i] = mpp_a[i];
            if (keep_indices) {
              mpi_a[i] = mpip_a[i];
            }
          }
        }
        for (uint32_t i = current_b_begin; i < current_b_end; i++) {
          if (mpp_b[i] > mp_b[i]) {
            mp_b[i] = mpp_b[i];
            if (keep_indices) {
              mpi_b[i] = mpip_b[i];
            }
          }
        }
        if (instrument) {
          if (current_a_end > current_a_begin) {
            task_stats.merge_values += current_a_end - current_a_begin;
          }
          if (current_b_end > current_b_begin) {
            task_stats.merge_values += current_b_end - current_b_begin;
          }
          merge_done = MatrixProfileClock::now();
        }
    }

    previous_a_begin = current_a_begin;
    previous_a_end = current_a_end;
    previous_b_begin = current_b_begin;
    previous_b_end = current_b_end;
    previous_keep_indices = keep_indices;

    if (instrument) {
      task_done = MatrixProfileClock::now();
      task_stats.scratch_ns =
          static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(scratch_done - task_started).count());
      task_stats.compute_ns =
          static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(compute_done - scratch_done).count());
      task_stats.mutex_wait_ns = static_cast<uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(lock_acquired - lock_requested).count());
      task_stats.merge_ns = static_cast<uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(merge_done - lock_acquired).count());
      task_stats.task_ns =
          static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(task_done - task_started).count());
      diagnostics->record(ab_ba != 0, task_stats);
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
  const bool keep_indices;
  const bool skip_invalid_blocks;
  MatrixProfileWorkerDiagnostics *const diagnostics;

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
                     const NumericVector &mp_b, const IntegerVector &mpi_a, const IntegerVector &mpi_b, bool idxs,
                     bool skip_invalid_blocks_, MatrixProfileWorkerDiagnostics *diagnostics_)
      : data(data), query(query), window_size(window_size), df_a(df_a), df_b(df_b), dg_a(dg_a), dg_b(dg_b),
        mu_a(mu_a), mu_b(mu_b), sig_a(sig_a), sig_b(sig_b), first_a(first_a), first_b(first_b), valid_a(valid_a),
        valid_b(valid_b), order_a(order_a), order_b(order_b), mp_a(mp_a), mp_b(mp_b), mpi_a(mpi_a), mpi_b(mpi_b),
        keep_indices(idxs), skip_invalid_blocks(skip_invalid_blocks_), diagnostics(diagnostics_) {}

  void set_ba() { ba = true; }

  void operator()(std::size_t begin, std::size_t end) override {
    const bool instrument = diagnostics != nullptr && diagnostics->enabled;
    MatrixProfileWorkerTaskStats task_stats;
    MatrixProfileClock::time_point task_started;
    MatrixProfileClock::time_point scratch_done;
    MatrixProfileClock::time_point compute_done;
    MatrixProfileClock::time_point lock_requested;
    MatrixProfileClock::time_point lock_acquired;
    MatrixProfileClock::time_point merge_done;
    MatrixProfileClock::time_point task_done;
    if (instrument) {
      task_started = MatrixProfileClock::now();
    }

    uint32_t const profile_len_a = mp_a.size();
    uint32_t const profile_len_b = mp_b.size();
    thread_local std::vector<double> centered_window;
    thread_local std::vector<double> local_mp_a;
    thread_local std::vector<double> local_mp_b;
    thread_local std::vector<int> local_mpi_a;
    thread_local std::vector<int> local_mpi_b;
    thread_local uint32_t previous_a_begin = 0;
    thread_local uint32_t previous_a_end = 0;
    thread_local uint32_t previous_b_begin = 0;
    thread_local uint32_t previous_b_end = 0;
    thread_local bool previous_keep_indices = false;

    if (centered_window.size() != window_size) {
      centered_window.resize(window_size);
      if (instrument) {
        task_stats.scratch_initialized += window_size;
      }
    }
    if (local_mp_a.size() != profile_len_a || local_mp_b.size() != profile_len_b) {
      local_mp_a.assign(profile_len_a, R_NegInf);
      local_mp_b.assign(profile_len_b, R_NegInf);
      if (instrument) {
        task_stats.scratch_initialized += profile_len_a + profile_len_b;
      }
      if (keep_indices) {
        local_mpi_a.assign(profile_len_a, NA_INTEGER);
        local_mpi_b.assign(profile_len_b, NA_INTEGER);
        if (instrument) {
          task_stats.scratch_initialized += profile_len_a + profile_len_b;
        }
      }
      previous_a_begin = previous_a_end = previous_b_begin = previous_b_end = 0;
    } else {
      if (keep_indices && (local_mpi_a.size() != profile_len_a || local_mpi_b.size() != profile_len_b)) {
        local_mpi_a.assign(profile_len_a, NA_INTEGER);
        local_mpi_b.assign(profile_len_b, NA_INTEGER);
        if (instrument) {
          task_stats.scratch_initialized += profile_len_a + profile_len_b;
        }
        previous_a_begin = previous_a_end = previous_b_begin = previous_b_end = 0;
      } else if (keep_indices && !previous_keep_indices) {
        std::fill(local_mpi_a.begin(), local_mpi_a.end(), NA_INTEGER);
        std::fill(local_mpi_b.begin(), local_mpi_b.end(), NA_INTEGER);
        if (instrument) {
          task_stats.scratch_reset += profile_len_a + profile_len_b;
        }
      }
      if (previous_a_begin < previous_a_end) {
        std::fill(local_mp_a.begin() + previous_a_begin, local_mp_a.begin() + previous_a_end, R_NegInf);
        if (instrument) {
          task_stats.scratch_reset += previous_a_end - previous_a_begin;
        }
        if (keep_indices) {
          std::fill(local_mpi_a.begin() + previous_a_begin, local_mpi_a.begin() + previous_a_end, NA_INTEGER);
          if (instrument) {
            task_stats.scratch_reset += previous_a_end - previous_a_begin;
          }
        }
      }
      if (previous_b_begin < previous_b_end) {
        std::fill(local_mp_b.begin() + previous_b_begin, local_mp_b.begin() + previous_b_end, R_NegInf);
        if (instrument) {
          task_stats.scratch_reset += previous_b_end - previous_b_begin;
        }
        if (keep_indices) {
          std::fill(local_mpi_b.begin() + previous_b_begin, local_mpi_b.begin() + previous_b_end, NA_INTEGER);
          if (instrument) {
            task_stats.scratch_reset += previous_b_end - previous_b_begin;
          }
        }
      }
    }

    if (instrument) {
      scratch_done = MatrixProfileClock::now();
    }

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

    uint32_t current_a_begin = profile_len_a;
    uint32_t current_a_end = 0;
    uint32_t current_b_begin = profile_len_b;
    uint32_t current_b_end = 0;
    if (!ba) {
      for (std::size_t order_index = begin; order_index < end; order_index++) {
        uint32_t const diag = order_a_ptr[order_index];
        uint32_t const diagonal_length = MIN(profile_len_a - diag, profile_len_b);
        current_a_begin = std::min(current_a_begin, diag);
        current_a_end = std::max(current_a_end, diag + diagonal_length);
        current_b_begin = 0;
        current_b_end = std::max(current_b_end, diagonal_length);
      }
    } else {
      for (std::size_t order_index = begin; order_index < end; order_index++) {
        uint32_t const diag = order_b_ptr[order_index];
        uint32_t const diagonal_length = MIN(profile_len_b - diag, profile_len_a);
        current_a_begin = 0;
        current_a_end = std::max(current_a_end, diagonal_length);
        current_b_begin = std::min(current_b_begin, diag);
        current_b_end = std::max(current_b_end, diag + diagonal_length);
      }
    }

    if (!ba) {
      for (std::size_t order_index = begin; order_index < end; order_index++) {
        uint32_t const diag = order_a_ptr[order_index];
        for (uint64_t i = 0; i < window_size; i++) {
          centered_window[i] = data_ptr[diag + i] - mu_a_ptr[diag];
        }
        double covariance = std::inner_product(centered_window.begin(), centered_window.end(), first_b_ptr, 0.0);
        uint32_t const diagonal_length = MIN(profile_len_a - diag, profile_len_b);
        if (instrument) {
          task_stats.diagonals++;
          task_stats.pairs_visited += diagonal_length;
        }

        if (!skip_invalid_blocks) {
          uint32_t offset = 0;
          // Process valid runs without a mask branch for every pair. The
          // covariance recurrence still advances through invalid runs, so
          // this preserves the legacy sparse-barrier semantics exactly.
          while (offset < diagonal_length) {
            while (offset < diagonal_length && valid_a_ptr[offset + diag] && valid_b_ptr[offset]) {
              uint32_t const off_diag = offset + diag;
              covariance =
                  covariance + df_a_ptr[off_diag] * dg_b_ptr[offset] + dg_a_ptr[off_diag] * df_b_ptr[offset];
              double const correlation = covariance * sig_a_ptr[off_diag] * sig_b_ptr[offset];
              if (!std::isfinite(correlation)) {
                if (instrument) {
                  task_stats.nonfinite_pairs++;
                }
              } else {
                if (correlation > local_mp_a[off_diag]) {
                  local_mp_a[off_diag] = correlation;
                  if (keep_indices) {
                    local_mpi_a[off_diag] = offset + 1;
                  }
                }
                if (correlation > local_mp_b[offset]) {
                  local_mp_b[offset] = correlation;
                  if (keep_indices) {
                    local_mpi_b[offset] = off_diag + 1;
                  }
                }
              }
              ++offset;
            }
            while (offset < diagonal_length && (!valid_a_ptr[offset + diag] || !valid_b_ptr[offset])) {
              uint32_t const off_diag = offset + diag;
              covariance =
                  covariance + df_a_ptr[off_diag] * dg_b_ptr[offset] + dg_a_ptr[off_diag] * df_b_ptr[offset];
              if (instrument) {
                task_stats.invalid_pairs++;
              }
              ++offset;
            }
          }
          continue;
        }

        uint32_t offset = 0;
        while (offset < diagonal_length) {
          uint32_t off_diag = offset + diag;
          if (!valid_a_ptr[off_diag] || !valid_b_ptr[offset]) {
            do {
              if (instrument) {
                task_stats.invalid_pairs++;
              }
              ++offset;
            } while (offset < diagonal_length &&
                     (!valid_a_ptr[offset + diag] || !valid_b_ptr[offset]));
            if (offset >= diagonal_length) {
              break;
            }

            // The recurrence is not advanced over the invalid gap. Rebuild
            // the covariance at the first valid offset so the next block is
            // independent of samples that touched a non-finite value.
            off_diag = offset + diag;
            covariance = 0.0;
            for (uint64_t i = 0; i < window_size; i++) {
              covariance += (data_ptr[off_diag + i] - mu_a_ptr[off_diag]) *
                            (query_ptr[offset + i] - mu_b_ptr[offset]);
            }
          }

          // Process one contiguous valid block. The first pair uses either
          // the initial dot product or the restarted covariance above;
          // subsequent pairs use the O(1) recurrence.
          while (true) {
            off_diag = offset + diag;
            double const correlation = covariance * sig_a_ptr[off_diag] * sig_b_ptr[offset];
            if (!std::isfinite(correlation)) {
              if (instrument) {
                task_stats.nonfinite_pairs++;
              }
            } else {
              if (correlation > local_mp_a[off_diag]) {
                local_mp_a[off_diag] = correlation;
                if (keep_indices) {
                  local_mpi_a[off_diag] = offset + 1;
                }
              }
              if (correlation > local_mp_b[offset]) {
                local_mp_b[offset] = correlation;
                if (keep_indices) {
                  local_mpi_b[offset] = off_diag + 1;
                }
              }
            }

            ++offset;
            if (offset >= diagonal_length || !valid_a_ptr[offset + diag] || !valid_b_ptr[offset]) {
              break;
            }
            off_diag = offset + diag;
            covariance =
                covariance + df_a_ptr[off_diag] * dg_b_ptr[offset] + dg_a_ptr[off_diag] * df_b_ptr[offset];
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
        if (instrument) {
          task_stats.diagonals++;
          task_stats.pairs_visited += diagonal_length;
        }

        if (!skip_invalid_blocks) {
          uint32_t offset = 0;
          while (offset < diagonal_length) {
            while (offset < diagonal_length && valid_b_ptr[offset + diag] && valid_a_ptr[offset]) {
              uint32_t const off_diag = offset + diag;
              covariance =
                  covariance + df_b_ptr[off_diag] * dg_a_ptr[offset] + dg_b_ptr[off_diag] * df_a_ptr[offset];
              double const correlation = covariance * sig_b_ptr[off_diag] * sig_a_ptr[offset];
              if (!std::isfinite(correlation)) {
                if (instrument) {
                  task_stats.nonfinite_pairs++;
                }
              } else {
                if (correlation > local_mp_b[off_diag]) {
                  local_mp_b[off_diag] = correlation;
                  if (keep_indices) {
                    local_mpi_b[off_diag] = offset + 1;
                  }
                }
                if (correlation > local_mp_a[offset]) {
                  local_mp_a[offset] = correlation;
                  if (keep_indices) {
                    local_mpi_a[offset] = off_diag + 1;
                  }
                }
              }
              ++offset;
            }
            while (offset < diagonal_length && (!valid_b_ptr[offset + diag] || !valid_a_ptr[offset])) {
              uint32_t const off_diag = offset + diag;
              covariance =
                  covariance + df_b_ptr[off_diag] * dg_a_ptr[offset] + dg_b_ptr[off_diag] * df_a_ptr[offset];
              if (instrument) {
                task_stats.invalid_pairs++;
              }
              ++offset;
            }
          }
          continue;
        }

        uint32_t offset = 0;
        while (offset < diagonal_length) {
          uint32_t off_diag = offset + diag;
          if (!valid_b_ptr[off_diag] || !valid_a_ptr[offset]) {
            do {
              if (instrument) {
                task_stats.invalid_pairs++;
              }
              ++offset;
            } while (offset < diagonal_length &&
                     (!valid_b_ptr[offset + diag] || !valid_a_ptr[offset]));
            if (offset >= diagonal_length) {
              break;
            }

            off_diag = offset + diag;
            covariance = 0.0;
            for (uint64_t i = 0; i < window_size; i++) {
              covariance += (query_ptr[off_diag + i] - mu_b_ptr[off_diag]) *
                            (data_ptr[offset + i] - mu_a_ptr[offset]);
            }
          }

          while (true) {
            off_diag = offset + diag;
            double const correlation = covariance * sig_b_ptr[off_diag] * sig_a_ptr[offset];
            if (!std::isfinite(correlation)) {
              if (instrument) {
                task_stats.nonfinite_pairs++;
              }
            } else {
              if (correlation > local_mp_b[off_diag]) {
                local_mp_b[off_diag] = correlation;
                if (keep_indices) {
                  local_mpi_b[off_diag] = offset + 1;
                }
              }
              if (correlation > local_mp_a[offset]) {
                local_mp_a[offset] = correlation;
                if (keep_indices) {
                  local_mpi_a[offset] = off_diag + 1;
                }
              }
            }

            ++offset;
            if (offset >= diagonal_length || !valid_b_ptr[offset + diag] || !valid_a_ptr[offset]) {
              break;
            }
            off_diag = offset + diag;
            covariance =
                covariance + df_b_ptr[off_diag] * dg_a_ptr[offset] + dg_b_ptr[off_diag] * df_a_ptr[offset];
          }
        }
      }
    }

    if (instrument) {
      compute_done = MatrixProfileClock::now();
      lock_requested = compute_done;
    }
    {
      std::lock_guard<decltype(mutex)> lock(mutex);
      if (instrument) {
        lock_acquired = MatrixProfileClock::now();
      }
      for (uint32_t i = current_a_begin; i < current_a_end; i++) {
        if (local_mp_a[i] > mp_a_ptr[i]) {
          mp_a_ptr[i] = local_mp_a[i];
          if (keep_indices) {
            mpi_a_ptr[i] = local_mpi_a[i];
          }
        }
      }
      for (uint32_t i = current_b_begin; i < current_b_end; i++) {
        if (local_mp_b[i] > mp_b_ptr[i]) {
          mp_b_ptr[i] = local_mp_b[i];
          if (keep_indices) {
            mpi_b_ptr[i] = local_mpi_b[i];
          }
        }
      }
      if (instrument) {
        if (current_a_end > current_a_begin) {
          task_stats.merge_values += current_a_end - current_a_begin;
        }
        if (current_b_end > current_b_begin) {
          task_stats.merge_values += current_b_end - current_b_begin;
        }
        merge_done = MatrixProfileClock::now();
      }
    }

    previous_a_begin = current_a_begin;
    previous_a_end = current_a_end;
    previous_b_begin = current_b_begin;
    previous_b_end = current_b_end;
    previous_keep_indices = keep_indices;

    if (instrument) {
      task_done = MatrixProfileClock::now();
      task_stats.scratch_ns =
          static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(scratch_done - task_started).count());
      task_stats.compute_ns =
          static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(compute_done - scratch_done).count());
      task_stats.mutex_wait_ns = static_cast<uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(lock_acquired - lock_requested).count());
      task_stats.merge_ns = static_cast<uint64_t>(
          std::chrono::duration_cast<std::chrono::nanoseconds>(merge_done - lock_acquired).count());
      task_stats.task_ns =
          static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(task_done - task_started).count());
      diagnostics->record(ba, task_stats);
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

    IntegerVector mpi_a;
    IntegerVector mpi_b;
    if (idxs) {
      mpi_a = IntegerVector(profile_len_a, -1);
      mpi_b = IntegerVector(profile_len_b, -1);
    }

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

    MatrixProfileWorkerDiagnostics diagnostics(matrix_profile_worker_diagnostics_enabled());
    MatrixProfilePAB matrix_profile(data_ref, query_ref, window_size, df_a, df_b, dg_a, dg_b, mu_a, mu_b, sig_a, sig_b,
                                    ww_a, ww_b, compute_ordera, compute_orderb, mp_a, mp_b, mpi_a, mpi_b, idxs,
                                    &diagnostics);

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

    List result;
    if (idxs) {
      result = List::create(Rcpp::Named("matrix_profile") = mp_a, Rcpp::Named("profile_index") = mpi_a,
                            Rcpp::Named("mpb") = mp_b, Rcpp::Named("pib") = mpi_b,
                            Rcpp::Named("partial") = partial);
    } else {
      result = List::create(Rcpp::Named("matrix_profile") = mp_a, Rcpp::Named("mpb") = mp_b,
                            Rcpp::Named("partial") = partial);
    }
    if (diagnostics.enabled) {
      result["worker_diagnostics"] = diagnostics.as_list();
    }
    return result;
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

  if (mpx_finite_fast_path_safe(data_ref, window_size) &&
      mpx_finite_fast_path_safe(query_ref, window_size)) {
    List result = mpxab_rcpp_parallel(data_ref, query_ref, window_size, s_size, idxs, euclidean, progress);
    result["valid_window_a"] = LogicalVector(data_ref.length() - window_size + 1, true);
    result["valid_window_b"] = LogicalVector(query_ref.length() - window_size + 1, true);
    return result;
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

  int const *const valid_a_mask_ptr = valid_a.begin();
  int const *const valid_b_mask_ptr = valid_b.begin();
  uint64_t invalid_a_count = 0;
  uint64_t invalid_b_count = 0;
  for (uint32_t i = 0; i < profile_len_a; i++) {
    invalid_a_count += valid_a_mask_ptr[i] ? 0 : 1;
  }
  for (uint32_t i = 0; i < profile_len_b; i++) {
    invalid_b_count += valid_b_mask_ptr[i] ? 0 : 1;
  }
  // Restarting the covariance at valid blocks pays for a fresh O(m) dot
  // product. Use it only when invalid windows are common enough for skipping
  // the invalid runs to offset that cost; sparse barriers retain the cheaper
  // original recurrence path. The default threshold is 1%; an environment
  // override is intentionally kept for benchmark calibration only.
  double const invalid_block_threshold = matrix_profile_invalid_block_threshold();
  bool const skip_invalid_blocks =
      invalid_block_threshold >= 0.0 &&
      ((invalid_a_count / static_cast<double>(profile_len_a)) >= invalid_block_threshold ||
       (invalid_b_count / static_cast<double>(profile_len_b)) >= invalid_block_threshold);

  NumericVector mp_a(profile_len_a, R_NegInf);
  NumericVector mp_b(profile_len_b, R_NegInf);
  IntegerVector mpi_a;
  IntegerVector mpi_b;
  if (idxs) {
    mpi_a = IntegerVector(profile_len_a, NA_INTEGER);
    mpi_b = IntegerVector(profile_len_b, NA_INTEGER);
  }

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
  // A complete join benefits from monotonically ordered diagonals: workers
  // then touch narrow output ranges, reducing scratch resets and merge work.
  // Preserve random order for sampled joins, where it defines partial work.
  if (s_size > 0.0 && s_size < 1.0) {
    order_a = sample(order_a, order_a.size());
    order_b = sample(order_b, order_b.size());
  }

  uint64_t work_a = order_a.size();
  uint64_t work_b = order_b.size();
  if (s_size > 0.0 && s_size < 1.0) {
    work_a = static_cast<uint64_t>(round(work_a * s_size + DBL_EPSILON));
    work_b = static_cast<uint64_t>(round(work_b * s_size + DBL_EPSILON));
    partial = work_a < static_cast<uint64_t>(order_a.size()) || work_b < static_cast<uint64_t>(order_b.size());
  }

  Progress progress_bar(100, progress);
  MatrixProfileWorkerDiagnostics diagnostics(matrix_profile_worker_diagnostics_enabled());
  MatrixProfilePABNA matrix_profile(data, query, window_size, df_a, df_b, dg_a, dg_b, mu_a, mu_b, sig_a, sig_b,
                                    first_a, first_b, valid_a, valid_b, order_a, order_b, mp_a, mp_b, mpi_a, mpi_b,
                                    idxs, skip_invalid_blocks, &diagnostics);

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
  int *const mpi_a_ptr = idxs ? mpi_a.begin() : nullptr;
  int *const mpi_b_ptr = idxs ? mpi_b.begin() : nullptr;
  int const *const valid_a_ptr = valid_a.begin();
  int const *const valid_b_ptr = valid_b.begin();
  for (uint32_t i = 0; i < profile_len_a; i++) {
    if (!valid_a_ptr[i] || !std::isfinite(mp_a_ptr[i])) {
      mp_a_ptr[i] = NA_REAL;
      if (idxs) {
        mpi_a_ptr[i] = NA_INTEGER;
      }
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
      if (idxs) {
        mpi_b_ptr[i] = NA_INTEGER;
      }
      continue;
    }
    mp_b_ptr[i] = std::max(-1.0, std::min(1.0, mp_b_ptr[i]));
    if (euclidean) {
      mp_b_ptr[i] = sqrt(std::max(0.0, 2.0 * window_size * (1.0 - mp_b_ptr[i])));
    }
  }

  List result;
  if (idxs) {
    result = List::create(Rcpp::Named("matrix_profile") = mp_a, Rcpp::Named("profile_index") = mpi_a,
                          Rcpp::Named("mpb") = mp_b, Rcpp::Named("pib") = mpi_b,
                          Rcpp::Named("valid_window_a") = valid_a, Rcpp::Named("valid_window_b") = valid_b,
                          Rcpp::Named("partial") = partial);
  } else {
    result = List::create(Rcpp::Named("matrix_profile") = mp_a, Rcpp::Named("mpb") = mp_b,
                          Rcpp::Named("valid_window_a") = valid_a, Rcpp::Named("valid_window_b") = valid_b,
                          Rcpp::Named("partial") = partial);
  }
  if (diagnostics.enabled) {
    List worker_diagnostics = diagnostics.as_list();
    worker_diagnostics["skip_invalid_blocks"] = skip_invalid_blocks;
    worker_diagnostics["invalid_block_threshold"] = invalid_block_threshold;
    result["worker_diagnostics"] = worker_diagnostics;
  }
  return result;
}
