// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#include "mathtools.h" // math first to fix OSX error
#include "mpx_parallel.h"
#include "windowfunc.h"
#include <cfloat> // DBL_EPSILON when STRICT_R_HEADERS
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <limits>
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

    // Keep scratch storage relative to the output range touched by this
    // scheduler task. In an unequal AB join, this avoids retaining a full
    // profile-B buffer per executor thread when only a much smaller slice is
    // ever read or merged. Reusing vector capacity preserves the allocation
    // benefit of the previous thread-local implementation, while initializing
    // the current range removes the need to reset the preceding range.
    thread_local std::vector<double> inn;
    thread_local std::vector<double> mpp_a;
    thread_local std::vector<int> mpip_a;
    thread_local std::vector<double> mpp_b;
    thread_local std::vector<int> mpip_b;

    if (inn.size() != window_size) {
      inn.resize(window_size);
      if (instrument) {
        task_stats.scratch_initialized += window_size;
      }
    }

    uint32_t const local_a_size = current_a_end - current_a_begin;
    uint32_t const local_b_size = current_b_end - current_b_begin;
    if (mpp_a.capacity() < local_a_size) {
      if (instrument) {
        task_stats.scratch_initialized += local_a_size;
      }
    }
    if (mpp_b.capacity() < local_b_size) {
      if (instrument) {
        task_stats.scratch_initialized += local_b_size;
      }
    }
    mpp_a.assign(local_a_size, -1.0);
    mpp_b.assign(local_b_size, -1.0);
    if (instrument) {
      task_stats.scratch_reset += local_a_size + local_b_size;
    }
    if (keep_indices) {
      if (mpip_a.capacity() < local_a_size) {
        if (instrument) {
          task_stats.scratch_initialized += local_a_size;
        }
      }
      if (mpip_b.capacity() < local_b_size) {
        if (instrument) {
          task_stats.scratch_initialized += local_b_size;
        }
      }
      mpip_a.assign(local_a_size, -1);
      mpip_b.assign(local_b_size, -1);
      if (instrument) {
        task_stats.scratch_reset += local_a_size + local_b_size;
      }
    }

    if (instrument) {
      scratch_done = MatrixProfileClock::now();
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

            if (c_cmp > mpp_b[offset - current_b_begin]) {
              mpp_b[offset - current_b_begin] = c_cmp;
              if (keep_indices) {
                mpip_b[offset - current_b_begin] = off_diag + 1;
              }
            }
            if (c_cmp > mpp_a[off_diag - current_a_begin]) {
              mpp_a[off_diag - current_a_begin] = c_cmp;
              if (keep_indices) {
                mpip_a[off_diag - current_a_begin] = offset + 1;
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

            if (c_cmp > mpp_a[offset - current_a_begin]) {
              mpp_a[offset - current_a_begin] = c_cmp;
              if (keep_indices) {
                mpip_a[offset - current_a_begin] = off_diag + 1;
              }
            }
            if (c_cmp > mpp_b[off_diag - current_b_begin]) {
              mpp_b[off_diag - current_b_begin] = c_cmp;
              if (keep_indices) {
                mpip_b[off_diag - current_b_begin] = offset + 1;
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
          if (mpp_a[i - current_a_begin] > mp_a[i]) {
            mp_a[i] = mpp_a[i - current_a_begin];
            if (keep_indices) {
              mpi_a[i] = mpip_a[i - current_a_begin];
            }
          }
        }
        for (uint32_t i = current_b_begin; i < current_b_end; i++) {
          if (mpp_b[i - current_b_begin] > mp_b[i]) {
            mp_b[i] = mpp_b[i - current_b_begin];
            if (keep_indices) {
              mpi_b[i] = mpip_b[i - current_b_begin];
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

    // Allocate only the output spans touched by this task. The vector capacity
    // remains thread-local and is reused by TBB, but the active working set is
    // bounded by the task's diagonal range rather than by both full profiles.
    thread_local std::vector<double> centered_window;
    thread_local std::vector<double> local_mp_a;
    thread_local std::vector<double> local_mp_b;
    thread_local std::vector<int> local_mpi_a;
    thread_local std::vector<int> local_mpi_b;

    if (centered_window.size() != window_size) {
      centered_window.resize(window_size);
      if (instrument) {
        task_stats.scratch_initialized += window_size;
      }
    }

    uint32_t const local_a_size = current_a_end - current_a_begin;
    uint32_t const local_b_size = current_b_end - current_b_begin;
    if (local_mp_a.capacity() < local_a_size) {
      if (instrument) {
        task_stats.scratch_initialized += local_a_size;
      }
    }
    if (local_mp_b.capacity() < local_b_size) {
      if (instrument) {
        task_stats.scratch_initialized += local_b_size;
      }
    }
    local_mp_a.assign(local_a_size, R_NegInf);
    local_mp_b.assign(local_b_size, R_NegInf);
    if (instrument) {
      task_stats.scratch_reset += local_a_size + local_b_size;
    }
    if (keep_indices) {
      if (local_mpi_a.capacity() < local_a_size) {
        if (instrument) {
          task_stats.scratch_initialized += local_a_size;
        }
      }
      if (local_mpi_b.capacity() < local_b_size) {
        if (instrument) {
          task_stats.scratch_initialized += local_b_size;
        }
      }
      local_mpi_a.assign(local_a_size, NA_INTEGER);
      local_mpi_b.assign(local_b_size, NA_INTEGER);
      if (instrument) {
        task_stats.scratch_reset += local_a_size + local_b_size;
      }
    }

    if (instrument) {
      scratch_done = MatrixProfileClock::now();
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
                if (correlation > local_mp_a[off_diag - current_a_begin]) {
                  local_mp_a[off_diag - current_a_begin] = correlation;
                  if (keep_indices) {
                    local_mpi_a[off_diag - current_a_begin] = offset + 1;
                  }
                }
                if (correlation > local_mp_b[offset - current_b_begin]) {
                  local_mp_b[offset - current_b_begin] = correlation;
                  if (keep_indices) {
                    local_mpi_b[offset - current_b_begin] = off_diag + 1;
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
              if (correlation > local_mp_a[off_diag - current_a_begin]) {
                local_mp_a[off_diag - current_a_begin] = correlation;
                if (keep_indices) {
                  local_mpi_a[off_diag - current_a_begin] = offset + 1;
                }
              }
              if (correlation > local_mp_b[offset - current_b_begin]) {
                local_mp_b[offset - current_b_begin] = correlation;
                if (keep_indices) {
                  local_mpi_b[offset - current_b_begin] = off_diag + 1;
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
                if (correlation > local_mp_b[off_diag - current_b_begin]) {
                  local_mp_b[off_diag - current_b_begin] = correlation;
                  if (keep_indices) {
                    local_mpi_b[off_diag - current_b_begin] = offset + 1;
                  }
                }
                if (correlation > local_mp_a[offset - current_a_begin]) {
                  local_mp_a[offset - current_a_begin] = correlation;
                  if (keep_indices) {
                    local_mpi_a[offset - current_a_begin] = off_diag + 1;
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
              if (correlation > local_mp_b[off_diag - current_b_begin]) {
                local_mp_b[off_diag - current_b_begin] = correlation;
                if (keep_indices) {
                  local_mpi_b[off_diag - current_b_begin] = offset + 1;
                }
              }
              if (correlation > local_mp_a[offset - current_a_begin]) {
                local_mp_a[offset - current_a_begin] = correlation;
                if (keep_indices) {
                  local_mpi_a[offset - current_a_begin] = off_diag + 1;
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
        if (local_mp_a[i - current_a_begin] > mp_a_ptr[i]) {
          mp_a_ptr[i] = local_mp_a[i - current_a_begin];
          if (keep_indices) {
            mpi_a_ptr[i] = local_mpi_a[i - current_a_begin];
          }
        }
      }
      for (uint32_t i = current_b_begin; i < current_b_end; i++) {
        if (local_mp_b[i - current_b_begin] > mp_b_ptr[i]) {
          mp_b_ptr[i] = local_mp_b[i - current_b_begin];
          if (keep_indices) {
            mpi_b_ptr[i] = local_mpi_b[i - current_b_begin];
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

struct MatrixProfileFiniteSegment {
  uint32_t sample_begin;
  uint32_t sample_end;

  uint32_t profile_length(uint32_t window_size) const { return sample_end - sample_begin - window_size + 1; }
};

// The segmented implementation is deliberately restricted to data for which
// invalid subsequences are caused solely by non-finite samples.  This check
// keeps its contract explicit: unlike the monolithic NA-aware implementation,
// it does not also model constant or otherwise non-normalizable windows.
static bool matrix_profile_validity_is_only_nonfinite(const NumericVector &series, uint32_t window_size,
                                                       const LogicalVector &valid_window) {
  uint32_t nonfinite_count = 0;
  for (uint32_t i = 0; i < series.length(); i++) {
    nonfinite_count += std::isfinite(series[i]) ? 0 : 1;
    if (i >= window_size) {
      nonfinite_count -= std::isfinite(series[i - window_size]) ? 0 : 1;
    }
    if (i >= window_size - 1) {
      uint32_t const profile_index = i - window_size + 1;
      if ((nonfinite_count == 0) != static_cast<bool>(valid_window[profile_index])) {
        return false;
      }
    }
  }
  return true;
}

static std::vector<MatrixProfileFiniteSegment> matrix_profile_finite_segments(const NumericVector &series,
                                                                                 uint32_t window_size) {
  std::vector<MatrixProfileFiniteSegment> segments;
  uint32_t position = 0;
  uint32_t const size = series.length();
  while (position < size) {
    while (position < size && !std::isfinite(series[position])) {
      position++;
    }
    uint32_t const begin = position;
    while (position < size && std::isfinite(series[position])) {
      position++;
    }
    if (position - begin >= window_size) {
      segments.push_back(MatrixProfileFiniteSegment{begin, position});
    }
  }
  return segments;
}

static double matrix_profile_segment_correlation(const NumericVector &a, uint32_t a_start, const NumericVector &b,
                                                 uint32_t b_start, uint32_t window_size) {
  double mean_a = 0.0;
  double mean_b = 0.0;
  for (uint32_t k = 0; k < window_size; k++) {
    mean_a += a[a_start + k];
    mean_b += b[b_start + k];
  }
  mean_a /= window_size;
  mean_b /= window_size;

  double dot = 0.0;
  double norm_a = 0.0;
  double norm_b = 0.0;
  for (uint32_t k = 0; k < window_size; k++) {
    double const centered_a = a[a_start + k] - mean_a;
    double const centered_b = b[b_start + k] - mean_b;
    dot += centered_a * centered_b;
    norm_a += centered_a * centered_a;
    norm_b += centered_b * centered_b;
  }
  return dot / sqrt(norm_a * norm_b);
}

struct MatrixProfileSegmentPairTask {
  MatrixProfileFiniteSegment a;
  MatrixProfileFiniteSegment b;
  uint32_t diagonal_begin = 0;
  uint32_t diagonal_end = 0;
};

static uint32_t matrix_profile_native_thread_hint() {
  const char *value = std::getenv("RCPP_PARALLEL_NUM_THREADS");
  if (value != nullptr && value[0] != '\0') {
    char *end = nullptr;
    long parsed = std::strtol(value, &end, 10);
    if (end != value && parsed > 0 && parsed <= std::numeric_limits<uint32_t>::max()) {
      return static_cast<uint32_t>(parsed);
    }
  }
  return 4;
}

static void matrix_profile_add_native_ab_tasks(const std::vector<MatrixProfileFiniteSegment> &segments_a,
                                               const std::vector<MatrixProfileFiniteSegment> &segments_b,
                                               uint32_t window_size, double s_size,
                                               std::vector<MatrixProfileSegmentPairTask> &tasks, bool &partial) {
  uint64_t const pair_count = static_cast<uint64_t>(segments_a.size()) * segments_b.size();
  uint64_t const target_tasks = std::max<uint64_t>(1, 4ULL * matrix_profile_native_thread_hint());
  uint64_t const split_factor = pair_count < target_tasks ? (target_tasks + pair_count - 1) / pair_count : 1;
  for (MatrixProfileFiniteSegment const &a : segments_a) {
    for (MatrixProfileFiniteSegment const &b : segments_b) {
      uint32_t const diagonal_count = a.profile_length(window_size) + b.profile_length(window_size) - 1;
      uint32_t const parts = static_cast<uint32_t>(std::min<uint64_t>(split_factor, diagonal_count));
      uint32_t const tile_size = (diagonal_count + parts - 1) / parts;
      uint32_t const tile_count = (diagonal_count + tile_size - 1) / tile_size;
      uint32_t const keep_tiles = static_cast<uint32_t>(round(tile_count * s_size + DBL_EPSILON));
      if (keep_tiles < tile_count) partial = true;
      IntegerVector order = Range(0, tile_count - 1);
      if (keep_tiles < tile_count && keep_tiles > 0) {
        order = sample(order, order.size());
      }
      std::vector<unsigned char> selected(tile_count, 0);
      for (uint32_t i = 0; i < keep_tiles; i++) selected[order[i]] = 1;
      for (uint32_t tile = 0; tile < tile_count; tile++) {
        if (selected[tile]) {
          uint32_t const begin = tile * tile_size;
          uint32_t const end = std::min(diagonal_count, begin + tile_size);
          tasks.push_back(MatrixProfileSegmentPairTask{a, b, begin, end});
        }
      }
    }
  }
}

static void matrix_profile_add_native_aa_tasks(const std::vector<MatrixProfileFiniteSegment> &segments,
                                               uint32_t window_size, uint32_t exclusion, double s_size,
                                               std::vector<MatrixProfileSegmentPairTask> &tasks, bool &partial) {
  uint64_t const pair_count = static_cast<uint64_t>(segments.size()) * (segments.size() + 1) / 2;
  uint64_t const target_tasks = std::max<uint64_t>(1, 4ULL * matrix_profile_native_thread_hint());
  uint64_t const split_factor = pair_count < target_tasks ? (target_tasks + pair_count - 1) / pair_count : 1;
  for (std::size_t a_index = 0; a_index < segments.size(); a_index++) {
    for (std::size_t b_index = a_index; b_index < segments.size(); b_index++) {
      MatrixProfileFiniteSegment const &a = segments[a_index];
      MatrixProfileFiniteSegment const &b = segments[b_index];
      uint32_t const a_profiles = a.profile_length(window_size);
      uint32_t const b_profiles = b.profile_length(window_size);
      bool const same_segment = a.sample_begin == b.sample_begin;
      uint32_t const diagonal_begin = same_segment ? std::min(exclusion + 1, a_profiles) : 0;
      uint32_t const diagonal_count = same_segment
                                          ? a_profiles - diagonal_begin
                                          : a_profiles + b_profiles - 1;
      if (diagonal_count == 0) continue;
      uint32_t const parts = static_cast<uint32_t>(std::min<uint64_t>(split_factor, diagonal_count));
      uint32_t const tile_size = (diagonal_count + parts - 1) / parts;
      uint32_t const tile_count = (diagonal_count + tile_size - 1) / tile_size;
      uint32_t const keep_tiles = static_cast<uint32_t>(round(tile_count * s_size + DBL_EPSILON));
      if (keep_tiles < tile_count) partial = true;
      IntegerVector order = Range(0, tile_count - 1);
      if (keep_tiles < tile_count && keep_tiles > 0) order = sample(order, order.size());
      std::vector<unsigned char> selected(tile_count, 0);
      for (uint32_t i = 0; i < keep_tiles; i++) selected[order[i]] = 1;
      for (uint32_t tile = 0; tile < tile_count; tile++) {
        if (selected[tile]) {
          uint32_t const begin = diagonal_begin + tile * tile_size;
          uint32_t const end = std::min(diagonal_begin + diagonal_count, begin + tile_size);
          tasks.push_back(MatrixProfileSegmentPairTask{a, b, begin, end});
        }
      }
    }
  }
}

// This worker is the native counterpart to the R-level segment reduction.
// Statistics are calculated once for the complete input, and every task
// evaluates one finite A x B block directly against those global statistics.
// It never calls an exported MPX routine from inside C++, avoiding repeated
// normalization, R allocations, and nested RcppParallel scheduling.
class MatrixProfileSegmentedNativeAB : public Worker {
private:
  const RVector<double> data;
  const RVector<double> query;
  const RVector<double> mu_a;
  const RVector<double> mu_b;
  const RVector<double> sig_a;
  const RVector<double> sig_b;
  const RVector<double> df_a;
  const RVector<double> df_b;
  const RVector<double> dg_a;
  const RVector<double> dg_b;
  const std::vector<MatrixProfileSegmentPairTask> &tasks;
  const uint32_t window_size;
  RVector<double> mp_a;
  RVector<double> mp_b;
  RVector<int> mpi_a;
  RVector<int> mpi_b;
  const bool keep_indices;

#if RCPP_PARALLEL_USE_TBB
  tbb::spin_mutex mutex;
#else
  tthread::mutex mutex;
#endif

  static double covariance_at(const double *data_ptr, const double *query_ptr, const double *mu_a_ptr,
                              const double *mu_b_ptr, uint32_t a_start, uint32_t b_start, uint32_t window_size) {
    double covariance = 0.0;
    for (uint32_t k = 0; k < window_size; k++) {
      covariance += (data_ptr[a_start + k] - mu_a_ptr[a_start]) * (query_ptr[b_start + k] - mu_b_ptr[b_start]);
    }
    return covariance;
  }

public:
  MatrixProfileSegmentedNativeAB(const NumericVector &data, const NumericVector &query, const NumericVector &mu_a,
                                 const NumericVector &mu_b, const NumericVector &sig_a, const NumericVector &sig_b,
                                 const NumericVector &df_a, const NumericVector &df_b, const NumericVector &dg_a,
                                 const NumericVector &dg_b, const std::vector<MatrixProfileSegmentPairTask> &tasks,
                                 uint32_t window_size,
                                 const NumericVector &mp_a, const NumericVector &mp_b, const IntegerVector &mpi_a,
                                 const IntegerVector &mpi_b, bool idxs)
      : data(data), query(query), mu_a(mu_a), mu_b(mu_b), sig_a(sig_a), sig_b(sig_b), df_a(df_a), df_b(df_b),
        dg_a(dg_a), dg_b(dg_b), tasks(tasks), window_size(window_size), mp_a(mp_a), mp_b(mp_b), mpi_a(mpi_a),
        mpi_b(mpi_b), keep_indices(idxs) {}

  void operator()(std::size_t begin, std::size_t end) override {
    double const *const data_ptr = data.begin();
    double const *const query_ptr = query.begin();
    double const *const mu_a_ptr = mu_a.begin();
    double const *const mu_b_ptr = mu_b.begin();
    double const *const sig_a_ptr = sig_a.begin();
    double const *const sig_b_ptr = sig_b.begin();
    double const *const df_a_ptr = df_a.begin();
    double const *const df_b_ptr = df_b.begin();
    double const *const dg_a_ptr = dg_a.begin();
    double const *const dg_b_ptr = dg_b.begin();
    double *const mp_a_ptr = mp_a.begin();
    double *const mp_b_ptr = mp_b.begin();
    int *const mpi_a_ptr = keep_indices ? mpi_a.begin() : nullptr;
    int *const mpi_b_ptr = keep_indices ? mpi_b.begin() : nullptr;

    thread_local std::vector<double> local_mp_a;
    thread_local std::vector<double> local_mp_b;
    thread_local std::vector<int> local_mpi_a;
    thread_local std::vector<int> local_mpi_b;

    for (std::size_t task_index = begin; task_index < end; task_index++) {
      MatrixProfileSegmentPairTask const &task = tasks[task_index];
      uint32_t const a_profiles = task.a.profile_length(window_size);
      uint32_t const b_profiles = task.b.profile_length(window_size);
      uint32_t const total_diagonals = a_profiles + b_profiles - 1;
      uint32_t const diagonal_begin = std::min(task.diagonal_begin, total_diagonals);
      uint32_t const diagonal_end = task.diagonal_end == 0
                                        ? total_diagonals
                                        : std::min(task.diagonal_end, total_diagonals);
      local_mp_a.assign(a_profiles, R_NegInf);
      local_mp_b.assign(b_profiles, R_NegInf);
      if (keep_indices) {
        local_mpi_a.assign(a_profiles, NA_INTEGER);
        local_mpi_b.assign(b_profiles, NA_INTEGER);
      }

      auto update_diagonal = [&](uint32_t a_local_start, uint32_t b_local_start, uint32_t diagonal_length) {
        uint32_t const a_start = task.a.sample_begin + a_local_start;
        uint32_t const b_start = task.b.sample_begin + b_local_start;
        double covariance = covariance_at(data_ptr, query_ptr, mu_a_ptr, mu_b_ptr, a_start, b_start, window_size);
        for (uint32_t offset = 0; offset < diagonal_length; offset++) {
          uint32_t const a_local = a_local_start + offset;
          uint32_t const b_local = b_local_start + offset;
          uint32_t const a_global = task.a.sample_begin + a_local;
          uint32_t const b_global = task.b.sample_begin + b_local;
          double const correlation = covariance * sig_a_ptr[a_global] * sig_b_ptr[b_global];
          if (std::isfinite(correlation)) {
            if (correlation > local_mp_a[a_local]) {
              local_mp_a[a_local] = correlation;
              if (keep_indices) local_mpi_a[a_local] = b_global + 1;
            }
            if (correlation > local_mp_b[b_local]) {
              local_mp_b[b_local] = correlation;
              if (keep_indices) local_mpi_b[b_local] = a_global + 1;
            }
          }
          if (offset + 1 < diagonal_length) {
            uint32_t const next_a = a_global + 1;
            uint32_t const next_b = b_global + 1;
            covariance += df_a_ptr[next_a] * dg_b_ptr[next_b] + dg_a_ptr[next_a] * df_b_ptr[next_b];
          }
        }
      };

      for (uint32_t diagonal = diagonal_begin; diagonal < diagonal_end; diagonal++) {
        if (diagonal < a_profiles) {
          update_diagonal(diagonal, 0, std::min(a_profiles - diagonal, b_profiles));
        } else {
          uint32_t const b_local = diagonal - a_profiles + 1;
          update_diagonal(0, b_local, std::min(a_profiles, b_profiles - b_local));
        }
      }

      std::lock_guard<decltype(mutex)> lock(mutex);
      for (uint32_t a_local = 0; a_local < a_profiles; a_local++) {
        uint32_t const a_global = task.a.sample_begin + a_local;
        if (local_mp_a[a_local] > mp_a_ptr[a_global]) {
          mp_a_ptr[a_global] = local_mp_a[a_local];
          if (keep_indices) mpi_a_ptr[a_global] = local_mpi_a[a_local];
        }
      }
      for (uint32_t b_local = 0; b_local < b_profiles; b_local++) {
        uint32_t const b_global = task.b.sample_begin + b_local;
        if (local_mp_b[b_local] > mp_b_ptr[b_global]) {
          mp_b_ptr[b_global] = local_mp_b[b_local];
          if (keep_indices) mpi_b_ptr[b_global] = local_mpi_b[b_local];
        }
      }
    }
  }
};

// True native finite-block AB segmentation.  Unlike
// mpxab_na_segmented_rcpp_parallel(), this function does not enter the regular
// MPX implementation for each segment pair: statistics and task scheduling
// are shared by the whole join.
// [[Rcpp::export]]
List mpxab_na_segmented_native_rcpp_parallel(NumericVector data_ref, NumericVector query_ref, uint64_t window_size,
                                             double s_size, bool idxs, bool euclidean, bool progress) {
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
  if (!matrix_profile_validity_is_only_nonfinite(data_ref, window_size, valid_a) ||
      !matrix_profile_validity_is_only_nonfinite(query_ref, window_size, valid_b)) {
    Rcpp::stop("the native segmented NA-aware AB join supports only finite windows with non-finite barriers; use mpxab_na_rcpp_parallel for constant or non-normalizable windows");
  }

  std::vector<MatrixProfileFiniteSegment> const segments_a = matrix_profile_finite_segments(data_ref, window_size);
  std::vector<MatrixProfileFiniteSegment> const segments_b = matrix_profile_finite_segments(query_ref, window_size);
  if (segments_a.empty() || segments_b.empty()) {
    Rcpp::stop("the native segmented NA-aware AB join requires at least one finite window in both series");
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
  std::vector<MatrixProfileSegmentPairTask> tasks;
  tasks.reserve(segments_a.size() * segments_b.size());
  bool partial = false;
  matrix_profile_add_native_ab_tasks(segments_a, segments_b, window_size, s_size, tasks, partial);

  NumericVector mp_a(profile_len_a, R_NegInf);
  NumericVector mp_b(profile_len_b, R_NegInf);
  IntegerVector mpi_a;
  IntegerVector mpi_b;
  if (idxs) {
    mpi_a = IntegerVector(profile_len_a, NA_INTEGER);
    mpi_b = IntegerVector(profile_len_b, NA_INTEGER);
  }
  MatrixProfileSegmentedNativeAB worker(data, query, mu_a, mu_b, sig_a, sig_b, df_a, df_b, dg_a, dg_b, tasks,
                                        window_size, mp_a, mp_b, mpi_a, mpi_b, idxs);
  (void)progress;
  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, tasks.size(), worker, 1);
#else
    RcppParallel2::ttParallelFor(0, tasks.size(), worker, 1);
#endif
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Native segmented AB join terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  for (uint32_t i = 0; i < profile_len_a; i++) {
    if (!valid_a[i] || !std::isfinite(mp_a[i])) {
      mp_a[i] = NA_REAL;
      if (idxs) mpi_a[i] = NA_INTEGER;
      continue;
    }
    mp_a[i] = std::max(-1.0, std::min(1.0, mp_a[i]));
    if (euclidean) mp_a[i] = sqrt(std::max(0.0, 2.0 * window_size * (1.0 - mp_a[i])));
  }
  for (uint32_t i = 0; i < profile_len_b; i++) {
    if (!valid_b[i] || !std::isfinite(mp_b[i])) {
      mp_b[i] = NA_REAL;
      if (idxs) mpi_b[i] = NA_INTEGER;
      continue;
    }
    mp_b[i] = std::max(-1.0, std::min(1.0, mp_b[i]));
    if (euclidean) mp_b[i] = sqrt(std::max(0.0, 2.0 * window_size * (1.0 - mp_b[i])));
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

class MatrixProfileSegmentedNativeAA : public Worker {
private:
  const RVector<double> data;
  const RVector<double> mu;
  const RVector<double> sig;
  const RVector<double> df;
  const RVector<double> dg;
  const std::vector<MatrixProfileSegmentPairTask> &tasks;
  const uint32_t window_size;
  const uint32_t exclusion;
  RVector<double> mp;
  RVector<int> mpi;
  const bool keep_indices;

#if RCPP_PARALLEL_USE_TBB
  tbb::spin_mutex mutex;
#else
  tthread::mutex mutex;
#endif

  static double covariance_at(const double *data_ptr, const double *mu_ptr, uint32_t a_start, uint32_t b_start,
                              uint32_t window_size) {
    double covariance = 0.0;
    for (uint32_t k = 0; k < window_size; k++) {
      covariance += (data_ptr[a_start + k] - mu_ptr[a_start]) * (data_ptr[b_start + k] - mu_ptr[b_start]);
    }
    return covariance;
  }

public:
  MatrixProfileSegmentedNativeAA(const NumericVector &data, const NumericVector &mu, const NumericVector &sig,
                                 const NumericVector &df, const NumericVector &dg,
                                 const std::vector<MatrixProfileSegmentPairTask> &tasks, uint32_t window_size,
                                 uint32_t exclusion, const NumericVector &mp, const IntegerVector &mpi, bool idxs)
      : data(data), mu(mu), sig(sig), df(df), dg(dg), tasks(tasks), window_size(window_size), exclusion(exclusion),
        mp(mp), mpi(mpi), keep_indices(idxs) {}

  void operator()(std::size_t begin, std::size_t end) override {
    double const *const data_ptr = data.begin();
    double const *const mu_ptr = mu.begin();
    double const *const sig_ptr = sig.begin();
    double const *const df_ptr = df.begin();
    double const *const dg_ptr = dg.begin();
    double *const mp_ptr = mp.begin();
    int *const mpi_ptr = keep_indices ? mpi.begin() : nullptr;

    thread_local std::vector<double> local_a;
    thread_local std::vector<double> local_b;
    thread_local std::vector<int> local_i_a;
    thread_local std::vector<int> local_i_b;

    for (std::size_t task_index = begin; task_index < end; task_index++) {
      MatrixProfileSegmentPairTask const &task = tasks[task_index];
      uint32_t const a_profiles = task.a.profile_length(window_size);
      uint32_t const b_profiles = task.b.profile_length(window_size);
      bool const same_segment = task.a.sample_begin == task.b.sample_begin;
      uint32_t const total_diagonals = a_profiles + b_profiles - 1;
      local_a.assign(a_profiles, R_NegInf);
      if (!same_segment) local_b.assign(b_profiles, R_NegInf);
      if (keep_indices) {
        local_i_a.assign(a_profiles, NA_INTEGER);
        if (!same_segment) local_i_b.assign(b_profiles, NA_INTEGER);
      }

      auto update_diagonal = [&](uint32_t a_local_start, uint32_t b_local_start, uint32_t diagonal_length,
                                 bool self_diagonal) {
        uint32_t const a_start = task.a.sample_begin + a_local_start;
        uint32_t const b_start = task.b.sample_begin + b_local_start;
        double covariance = covariance_at(data_ptr, mu_ptr, a_start, b_start, window_size);
        for (uint32_t offset = 0; offset < diagonal_length; offset++) {
          uint32_t const a_local = a_local_start + offset;
          uint32_t const b_local = b_local_start + offset;
          uint32_t const a_global = task.a.sample_begin + a_local;
          uint32_t const b_global = task.b.sample_begin + b_local;
          double const correlation = covariance * sig_ptr[a_global] * sig_ptr[b_global];
          if (std::isfinite(correlation)) {
            if (correlation > local_a[a_local]) {
              local_a[a_local] = correlation;
              if (keep_indices) local_i_a[a_local] = b_global + 1;
            }
            if (self_diagonal) {
              if (correlation > local_a[b_local]) {
                local_a[b_local] = correlation;
                if (keep_indices) local_i_a[b_local] = a_global + 1;
              }
            } else if (correlation > local_b[b_local]) {
              local_b[b_local] = correlation;
              if (keep_indices) local_i_b[b_local] = a_global + 1;
            }
          }
          if (offset + 1 < diagonal_length) {
            uint32_t const next_a = a_global + 1;
            uint32_t const next_b = b_global + 1;
            covariance += df_ptr[next_a] * dg_ptr[next_b] + dg_ptr[next_a] * df_ptr[next_b];
          }
        }
      };

      if (same_segment) {
        uint32_t const diagonal_begin = std::max(exclusion + 1, task.diagonal_begin);
        uint32_t const diagonal_end = task.diagonal_end == 0
                                          ? a_profiles
                                          : std::min(task.diagonal_end, a_profiles);
        for (uint32_t diagonal = diagonal_begin; diagonal < diagonal_end; diagonal++) {
          update_diagonal(diagonal, 0, a_profiles - diagonal, true);
        }
      } else {
        uint32_t const diagonal_begin = std::min(task.diagonal_begin, total_diagonals);
        uint32_t const diagonal_end = task.diagonal_end == 0
                                          ? total_diagonals
                                          : std::min(task.diagonal_end, total_diagonals);
        for (uint32_t diagonal = diagonal_begin; diagonal < diagonal_end; diagonal++) {
          if (diagonal < a_profiles) {
            update_diagonal(diagonal, 0, std::min(a_profiles - diagonal, b_profiles), false);
          } else {
            uint32_t const b_local = diagonal - a_profiles + 1;
            update_diagonal(0, b_local, std::min(a_profiles, b_profiles - b_local), false);
          }
        }
      }

      std::lock_guard<decltype(mutex)> lock(mutex);
      for (uint32_t a_local = 0; a_local < a_profiles; a_local++) {
        uint32_t const a_global = task.a.sample_begin + a_local;
        if (local_a[a_local] > mp_ptr[a_global]) {
          mp_ptr[a_global] = local_a[a_local];
          if (keep_indices) mpi_ptr[a_global] = local_i_a[a_local];
        }
      }
      if (!same_segment) {
        for (uint32_t b_local = 0; b_local < b_profiles; b_local++) {
          uint32_t const b_global = task.b.sample_begin + b_local;
          if (local_b[b_local] > mp_ptr[b_global]) {
            mp_ptr[b_global] = local_b[b_local];
            if (keep_indices) mpi_ptr[b_global] = local_i_b[b_local];
          }
        }
      }
    }
  }
};

// Native finite-block AA segmentation.  Only the non-trivial upper half of a
// within-segment block is traversed; cross-segment blocks are evaluated once
// and update both profile positions.
// [[Rcpp::export]]
List mpx_na_segmented_native_rcpp_parallel(NumericVector data_ref, uint64_t window_size, double ez, double s_size,
                                           bool idxs, bool euclidean, bool progress) {
  uint64_t const data_size = data_ref.length();
  if (window_size < 2 || window_size >= data_size) {
    Rcpp::stop("window_size must leave at least two subsequences");
  }
  if (!std::isfinite(ez) || ez < 0.0) {
    Rcpp::stop("exclusion_zone must be non-negative and finite");
  }
  if (!std::isfinite(s_size) || s_size < 0.0 || s_size > 1.0) {
    Rcpp::stop("s_size must be between 0 and 1");
  }

  uint32_t const profile_len = data_size - window_size + 1;
  List const stats = muinvn_na_parallel(data_ref, window_size);
  NumericVector data = stats["data"];
  NumericVector mu = stats["avg"];
  NumericVector sig = stats["sig"];
  LogicalVector valid_window = stats["valid_window"];
  if (!matrix_profile_validity_is_only_nonfinite(data_ref, window_size, valid_window)) {
    Rcpp::stop("the native segmented NA-aware self join supports only finite windows with non-finite barriers; use mpx_na_rcpp_parallel for constant or non-normalizable windows");
  }
  std::vector<MatrixProfileFiniteSegment> const segments = matrix_profile_finite_segments(data_ref, window_size);
  if (segments.empty()) {
    Rcpp::stop("the native segmented NA-aware self join requires at least one finite window");
  }
  NumericVector df = 0.5 * (data[Range(window_size, data_size - 1)] -
                            data[Range(0, data_size - window_size - 1)]);
  df.push_front(0);
  NumericVector dg = (data[Range(window_size, data_size - 1)] - mu[Range(1, profile_len - 1)]) +
                     (data[Range(0, data_size - window_size - 1)] - mu[Range(0, profile_len - 2)]);
  dg.push_front(0);

  std::vector<MatrixProfileSegmentPairTask> tasks;
  tasks.reserve(segments.size() * (segments.size() + 1) / 2);
  bool partial = false;
  matrix_profile_add_native_aa_tasks(segments, window_size, static_cast<uint32_t>(round(window_size * ez)), s_size,
                                     tasks, partial);

  NumericVector mp(profile_len, R_NegInf);
  IntegerVector mpi;
  if (idxs) mpi = IntegerVector(profile_len, NA_INTEGER);
  uint32_t const exclusion = static_cast<uint32_t>(round(window_size * ez));
  MatrixProfileSegmentedNativeAA worker(data, mu, sig, df, dg, tasks, window_size, exclusion, mp, mpi, idxs);
  (void)progress;
  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, tasks.size(), worker, 1);
#else
    RcppParallel2::ttParallelFor(0, tasks.size(), worker, 1);
#endif
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Native segmented self join terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  for (uint32_t i = 0; i < profile_len; i++) {
    if (!valid_window[i] || !std::isfinite(mp[i])) {
      mp[i] = NA_REAL;
      if (idxs) mpi[i] = NA_INTEGER;
      continue;
    }
    mp[i] = std::max(-1.0, std::min(1.0, mp[i]));
    if (euclidean) mp[i] = sqrt(std::max(0.0, 2.0 * window_size * (1.0 - mp[i])));
  }
  if (idxs) {
    return List::create(Rcpp::Named("matrix_profile") = mp, Rcpp::Named("profile_index") = mpi,
                        Rcpp::Named("valid_window") = valid_window, Rcpp::Named("partial") = partial);
  }
  return List::create(Rcpp::Named("matrix_profile") = mp, Rcpp::Named("valid_window") = valid_window,
                      Rcpp::Named("partial") = partial);
}

// Parallel AB-join that decomposes non-finite barriers into finite segments.
// It is intentionally a separate entry point from mpxab_na_rcpp_parallel:
// callers choose the strategy, and unsupported input is reported rather than
// silently switching algorithms.
// [[Rcpp::export]]
List mpxab_na_segmented_rcpp_parallel(NumericVector data_ref, NumericVector query_ref, uint64_t window_size,
                                      double s_size, bool idxs, bool euclidean, bool progress) {
  uint64_t const data_size = data_ref.length();
  uint64_t const query_size = query_ref.length();
  if (window_size < 2 || window_size >= data_size || window_size >= query_size) {
    Rcpp::stop("window_size must leave at least two subsequences in each series");
  }
  if (s_size != 1.0) {
    Rcpp::stop("the segmented NA-aware AB join currently requires s_size = 1; use mpxab_na_rcpp_parallel for sampled joins");
  }

  uint32_t const profile_len_a = data_size - window_size + 1;
  uint32_t const profile_len_b = query_size - window_size + 1;
  List const stats_a = muinvn_na_parallel(data_ref, window_size);
  List const stats_b = muinvn_na_parallel(query_ref, window_size);
  LogicalVector valid_a = stats_a["valid_window"];
  LogicalVector valid_b = stats_b["valid_window"];
  if (!matrix_profile_validity_is_only_nonfinite(data_ref, window_size, valid_a) ||
      !matrix_profile_validity_is_only_nonfinite(query_ref, window_size, valid_b)) {
    Rcpp::stop("the segmented NA-aware AB join supports only finite windows with non-finite barriers; use mpxab_na_rcpp_parallel for constant or non-normalizable windows");
  }

  std::vector<MatrixProfileFiniteSegment> const segments_a = matrix_profile_finite_segments(data_ref, window_size);
  std::vector<MatrixProfileFiniteSegment> const segments_b = matrix_profile_finite_segments(query_ref, window_size);
  if (segments_a.empty() || segments_b.empty()) {
    Rcpp::stop("the segmented NA-aware AB join requires at least one finite window in both series");
  }

  double const initial_profile = euclidean ? R_PosInf : R_NegInf;
  NumericVector mp_a(profile_len_a, initial_profile);
  NumericVector mp_b(profile_len_b, initial_profile);
  IntegerVector mpi_a;
  IntegerVector mpi_b;
  if (idxs) {
    mpi_a = IntegerVector(profile_len_a, NA_INTEGER);
    mpi_b = IntegerVector(profile_len_b, NA_INTEGER);
  }

  bool partial = false;
  // Segment calls delegate to the regular parallel kernel, whose own progress
  // object is intentionally disabled here.  A single outer progress object is
  // not safe while repeatedly entering RcppParallel from this native entry.
  (void)progress;
  for (std::size_t a_segment_index = 0; a_segment_index < segments_a.size() && !partial; a_segment_index++) {
    MatrixProfileFiniteSegment const &a_segment = segments_a[a_segment_index];
    NumericVector const a = data_ref[Range(a_segment.sample_begin, a_segment.sample_end - 1)];
    uint32_t const a_profiles = a_segment.profile_length(window_size);
    for (std::size_t b_segment_index = 0; b_segment_index < segments_b.size() && !partial; b_segment_index++) {
      MatrixProfileFiniteSegment const &b_segment = segments_b[b_segment_index];
      NumericVector const b = query_ref[Range(b_segment.sample_begin, b_segment.sample_end - 1)];
      uint32_t const b_profiles = b_segment.profile_length(window_size);

      if (a_profiles == 1 || b_profiles == 1) {
        for (uint32_t a_local = 0; a_local < a_profiles; a_local++) {
          for (uint32_t b_local = 0; b_local < b_profiles; b_local++) {
            double correlation = matrix_profile_segment_correlation(a, a_local, b, b_local, window_size);
            correlation = std::max(-1.0, std::min(1.0, correlation));
            double const candidate = euclidean ? sqrt(std::max(0.0, 2.0 * window_size * (1.0 - correlation))) : correlation;
            uint32_t const a_global = a_segment.sample_begin + a_local;
            uint32_t const b_global = b_segment.sample_begin + b_local;
            bool const improve_a = euclidean ? candidate < mp_a[a_global] : candidate > mp_a[a_global];
            bool const improve_b = euclidean ? candidate < mp_b[b_global] : candidate > mp_b[b_global];
            if (improve_a) {
              mp_a[a_global] = candidate;
              if (idxs) mpi_a[a_global] = b_global + 1;
            }
            if (improve_b) {
              mp_b[b_global] = candidate;
              if (idxs) mpi_b[b_global] = a_global + 1;
            }
          }
        }
      } else {
        List const local = mpxab_rcpp_parallel(a, b, window_size, 1.0, idxs, euclidean, false);
        NumericVector const local_mp_a = local["matrix_profile"];
        NumericVector const local_mp_b = local["mpb"];
        IntegerVector local_mpi_a;
        IntegerVector local_mpi_b;
        if (idxs) {
          local_mpi_a = local["profile_index"];
          local_mpi_b = local["pib"];
        }
        for (uint32_t a_local = 0; a_local < a_profiles; a_local++) {
          uint32_t const a_global = a_segment.sample_begin + a_local;
          double const candidate = local_mp_a[a_local];
          if ((euclidean ? candidate < mp_a[a_global] : candidate > mp_a[a_global])) {
            mp_a[a_global] = candidate;
            if (idxs) mpi_a[a_global] = b_segment.sample_begin + local_mpi_a[a_local];
          }
        }
        for (uint32_t b_local = 0; b_local < b_profiles; b_local++) {
          uint32_t const b_global = b_segment.sample_begin + b_local;
          double const candidate = local_mp_b[b_local];
          if ((euclidean ? candidate < mp_b[b_global] : candidate > mp_b[b_global])) {
            mp_b[b_global] = candidate;
            if (idxs) mpi_b[b_global] = a_segment.sample_begin + local_mpi_b[b_local];
          }
        }
        partial = as<bool>(local["partial"]);
      }
    }
  }

  for (uint32_t i = 0; i < profile_len_a; i++) {
    if (!valid_a[i] || !std::isfinite(mp_a[i])) {
      mp_a[i] = NA_REAL;
      if (idxs) mpi_a[i] = NA_INTEGER;
    }
  }
  for (uint32_t i = 0; i < profile_len_b; i++) {
    if (!valid_b[i] || !std::isfinite(mp_b[i])) {
      mp_b[i] = NA_REAL;
      if (idxs) mpi_b[i] = NA_INTEGER;
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
