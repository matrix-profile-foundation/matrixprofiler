// Exact, NA-aware Relative Frequency Contrast Profile.

#include "mathtools.h"
#include "mass.h"
#include "rfcp.h"
#include "windowfunc.h"

#include <RcppParallel.h>
#include <RcppThread.h>

#if !RCPP_PARALLEL_USE_TBB
#include "rcpp_parallel_fix.h"
#include "tthread/tinythread.h"
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <vector>

using namespace RcppParallel;

namespace {

constexpr uint64_t rfcp_max_profile_cells = 5000000ULL;

struct RFCPCandidate {
  double distance_squared;
  uint32_t index;
};

// std::make_heap puts the element for which comparator(top, other) is false
// at the front.  This comparator consequently yields the smallest distance,
// then the smallest global index, at the top of the heap.
struct RFCPMinCandidate {
  bool operator()(RFCPCandidate const &left, RFCPCandidate const &right) const {
    if (left.distance_squared != right.distance_squared) {
      return left.distance_squared > right.distance_squared;
    }
    return left.index > right.index;
  }
};

struct RFCPQueryTask {
  uint32_t anchor_begin;
  uint32_t output_begin;
  uint32_t output_end;
};

struct RFCPTaskBest {
  bool available = false;
  uint32_t query_index = 0;
  double rms = R_NegInf;
  std::vector<double> aa_distances;
  std::vector<double> ab_distances;
  std::vector<double> rfcp;
  std::vector<int> aa_indices;
  std::vector<int> ab_indices;
};

static uint32_t rfcp_query_rows_per_task(uint32_t window_size) {
  // A direct covariance profile initializes each fixed task.  Keeping at
  // least sixteen times m rows amortizes that restart to roughly six percent
  // of the pairwise work, while the bounds retain enough tasks for small and
  // very large inputs.
  uint64_t const scaled = 16ULL * window_size;
  return static_cast<uint32_t>(std::max<uint64_t>(1024ULL, std::min<uint64_t>(32768ULL, scaled)));
}

static void rfcp_fill_differences(const NumericVector &data, const NumericVector &mean,
                                  uint32_t window_size, NumericVector &df, NumericVector &dg) {
  uint32_t const profile_length = mean.length();
  double const *const data_ptr = data.begin();
  double const *const mean_ptr = mean.begin();
  double *const df_ptr = df.begin();
  double *const dg_ptr = dg.begin();
  df_ptr[0] = 0.0;
  dg_ptr[0] = 0.0;
  for (uint32_t i = 1; i < profile_length; i++) {
    uint32_t const incoming = i + window_size - 1;
    uint32_t const outgoing = i - 1;
    df_ptr[i] = 0.5 * (data_ptr[incoming] - data_ptr[outgoing]);
    dg_ptr[i] = (data_ptr[incoming] - mean_ptr[i]) + (data_ptr[outgoing] - mean_ptr[i - 1]);
  }
}

class RFCPWorker : public Worker {
private:
  const RVector<double> positive;
  const RVector<double> negative;
  const RVector<double> mean_positive;
  const RVector<double> mean_negative;
  const RVector<double> sig_positive;
  const RVector<double> sig_negative;
  const RVector<double> df_positive;
  const RVector<double> df_negative;
  const RVector<double> dg_positive;
  const RVector<double> dg_negative;
  const RVector<int> valid_positive;
  const RVector<int> valid_negative;
  const std::vector<RFCPQueryTask> &tasks;
  const uint32_t window_size;
  const uint32_t max_freq;
  const uint32_t exclusion_radius;
  const uint32_t request_begin;
  RVector<double> rms_profile;
  RVector<int> completed;
  RVector<int> paired_frequency;
  RVector<double> profile_aa;
  RVector<double> profile_ab;
  RVector<double> profile_rfcp;
  RVector<int> profile_index_aa;
  RVector<int> profile_index_ab;
  const bool keep_profiles;
  std::vector<RFCPTaskBest> &task_best;

  static double centered_covariance(const double *query, const double *reference,
                                    const double *query_mean, const double *reference_mean,
                                    uint32_t query_start, uint32_t reference_start,
                                    uint32_t window_size) {
    double covariance = 0.0;
    double const query_mu = query_mean[query_start];
    double const reference_mu = reference_mean[reference_start];
    for (uint32_t k = 0; k < window_size; k++) {
      covariance += (query[query_start + k] - query_mu) *
                    (reference[reference_start + k] - reference_mu);
    }
    return covariance;
  }

  static void initialize_covariances(const double *query, const double *reference,
                                     const double *query_mean, const double *reference_mean,
                                     uint32_t query_start, uint32_t reference_size,
                                     uint32_t window_size, std::vector<double> &covariance) {
    uint32_t const reference_profiles = reference_size - window_size + 1;
    covariance.resize(reference_profiles);
    uint64_t const direct_work =
        static_cast<uint64_t>(reference_profiles) * window_size;
    if (direct_work >= 1000000ULL) {
      // MASS-style overlap-save convolution. Centering only the query keeps
      // the covariance numerically stable without needing to center every
      // overlapping reference window: sum((q-mu_q) * r) is corrected by the
      // (small) residual centered-query sum below.
      uint32_t const grain = set_k_cpp(8192, reference_size, window_size);
      uint32_t const jump = grain - window_size + 1;
      std::vector<double> reversed_query(grain, 0.0);
      long double centered_sum_long = 0.0L;
      for (uint32_t k = 0; k < window_size; k++) {
        double const centered = query[query_start + k] - query_mean[query_start];
        reversed_query[window_size - 1 - k] = centered;
        centered_sum_long += static_cast<long double>(centered);
      }
      double const centered_sum = static_cast<double>(centered_sum_long);
      std::vector<std::complex<double>> query_fft = fft_rcpp(reversed_query, false);
      std::vector<std::complex<double>> product(grain);

      uint64_t reference_offset = 0;
      uint64_t const full_chunk_end = reference_size - grain;
      for (; reference_offset <= full_chunk_end; reference_offset += jump) {
        std::vector<double> reference_chunk(
            reference + reference_offset, reference + reference_offset + grain);
        std::vector<std::complex<double>> reference_fft =
            fft_rcpp(reference_chunk, false);
        std::transform(reference_fft.begin(), reference_fft.end(),
                       query_fft.begin(), product.begin(),
                       std::multiplies<std::complex<double>>());
        std::vector<double> convolution = fft_rcpp_real(product, true);
        for (uint32_t i = 0; i < jump; i++) {
          uint64_t const profile_index = reference_offset + i;
          covariance[profile_index] =
              convolution[window_size - 1 + i] -
              reference_mean[profile_index] * centered_sum;
        }
      }

      uint32_t const remaining = reference_size - reference_offset;
      if (remaining >= window_size) {
        std::vector<double> reference_chunk(
            reference + reference_offset, reference + reference_size);
        std::vector<double> final_query(remaining, 0.0);
        std::copy(reversed_query.begin(), reversed_query.begin() + window_size,
                  final_query.begin());
        std::vector<std::complex<double>> reference_fft =
            fft_rcpp(reference_chunk, false);
        query_fft = fft_rcpp(final_query, false);
        product.resize(remaining);
        std::transform(reference_fft.begin(), reference_fft.end(),
                       query_fft.begin(), product.begin(),
                       std::multiplies<std::complex<double>>());
        std::vector<double> convolution = fft_rcpp_real(product, true);
        uint32_t const remaining_profiles = remaining - window_size + 1;
        for (uint32_t i = 0; i < remaining_profiles; i++) {
          uint64_t const profile_index = reference_offset + i;
          covariance[profile_index] =
              convolution[window_size - 1 + i] -
              reference_mean[profile_index] * centered_sum;
        }
      }
      return;
    }

    for (uint32_t reference_start = 0; reference_start < reference_profiles; reference_start++) {
      covariance[reference_start] = centered_covariance(
          query, reference, query_mean, reference_mean, query_start, reference_start, window_size);
    }
  }

  static void advance_covariances(const double *query, const double *reference,
                                  const double *query_mean, const double *reference_mean,
                                  const double *query_df, const double *reference_df,
                                  const double *query_dg, const double *reference_dg,
                                  uint32_t query_start, uint32_t reference_profiles,
                                  uint32_t window_size, std::vector<double> &covariance) {
    for (uint32_t reference_start = reference_profiles - 1; reference_start > 0; reference_start--) {
      covariance[reference_start] = covariance[reference_start - 1] +
                                    query_df[query_start] * reference_dg[reference_start] +
                                    query_dg[query_start] * reference_df[reference_start];
    }
    covariance[0] = centered_covariance(
        query, reference, query_mean, reference_mean, query_start, 0, window_size);
  }

  static uint32_t next_stamp(std::vector<uint32_t> &marks, uint32_t &stamp) {
    stamp++;
    if (stamp == 0) {
      std::fill(marks.begin(), marks.end(), 0);
      stamp = 1;
    }
    return stamp;
  }

  static uint32_t select_neighbors(const std::vector<double> &covariance,
                                   const double *reference_sig, const int *reference_valid,
                                   double query_sig, uint32_t query_index, bool self_join,
                                   uint32_t window_size, uint32_t exclusion_radius,
                                   uint32_t max_freq, std::vector<RFCPCandidate> &candidates,
                                   std::vector<uint32_t> &exclusion_marks, uint32_t &stamp,
                                   std::vector<double> &distances, std::vector<int> &indices) {
    uint32_t const reference_profiles = covariance.size();
    candidates.clear();
    uint32_t const self_left = self_join && query_index > exclusion_radius
                                   ? query_index - exclusion_radius
                                   : 0;
    uint32_t const self_right = self_join
                                    ? static_cast<uint32_t>(std::min<uint64_t>(
                                          reference_profiles - 1,
                                          static_cast<uint64_t>(query_index) + exclusion_radius))
                                    : 0;
    double const distance_scale = 2.0 * window_size;
    for (uint32_t reference_index = 0; reference_index < reference_profiles; reference_index++) {
      if (!reference_valid[reference_index] ||
          (self_join && reference_index >= self_left && reference_index <= self_right)) {
        continue;
      }
      double correlation = covariance[reference_index] * query_sig * reference_sig[reference_index];
      if (!std::isfinite(correlation)) {
        continue;
      }
      correlation = std::max(-1.0, std::min(1.0, correlation));
      candidates.push_back(RFCPCandidate{
          std::max(0.0, distance_scale * (1.0 - correlation)), reference_index});
    }

    std::fill(distances.begin(), distances.end(), NA_REAL);
    std::fill(indices.begin(), indices.end(), NA_INTEGER);
    std::make_heap(candidates.begin(), candidates.end(), RFCPMinCandidate());
    uint32_t const current_stamp = next_stamp(exclusion_marks, stamp);
    uint32_t accepted = 0;
    while (!candidates.empty() && accepted < max_freq) {
      std::pop_heap(candidates.begin(), candidates.end(), RFCPMinCandidate());
      RFCPCandidate const candidate = candidates.back();
      candidates.pop_back();
      if (exclusion_marks[candidate.index] == current_stamp) {
        continue;
      }
      distances[accepted] = std::sqrt(candidate.distance_squared);
      indices[accepted] = static_cast<int>(candidate.index + 1);
      accepted++;

      uint32_t const left = candidate.index > exclusion_radius
                                ? candidate.index - exclusion_radius
                                : 0;
      uint32_t const right = static_cast<uint32_t>(std::min<uint64_t>(
          reference_profiles - 1,
          static_cast<uint64_t>(candidate.index) + exclusion_radius));
      std::fill(exclusion_marks.begin() + left, exclusion_marks.begin() + right + 1,
                current_stamp);
    }
    return accepted;
  }

public:
  RFCPWorker(const NumericVector &positive, const NumericVector &negative,
             const NumericVector &mean_positive, const NumericVector &mean_negative,
             const NumericVector &sig_positive, const NumericVector &sig_negative,
             const NumericVector &df_positive, const NumericVector &df_negative,
             const NumericVector &dg_positive, const NumericVector &dg_negative,
             const LogicalVector &valid_positive, const LogicalVector &valid_negative,
             const std::vector<RFCPQueryTask> &tasks, uint32_t window_size,
             uint32_t max_freq, uint32_t exclusion_radius, uint32_t request_begin,
             const NumericVector &rms_profile, const LogicalVector &completed,
             const IntegerVector &paired_frequency, const NumericVector &profile_aa,
             const NumericVector &profile_ab, const NumericVector &profile_rfcp,
             const IntegerVector &profile_index_aa, const IntegerVector &profile_index_ab,
             bool keep_profiles, std::vector<RFCPTaskBest> &task_best)
      : positive(positive), negative(negative), mean_positive(mean_positive),
        mean_negative(mean_negative), sig_positive(sig_positive), sig_negative(sig_negative),
        df_positive(df_positive), df_negative(df_negative), dg_positive(dg_positive),
        dg_negative(dg_negative), valid_positive(valid_positive), valid_negative(valid_negative),
        tasks(tasks), window_size(window_size), max_freq(max_freq),
        exclusion_radius(exclusion_radius), request_begin(request_begin),
        rms_profile(rms_profile), completed(completed), paired_frequency(paired_frequency),
        profile_aa(profile_aa), profile_ab(profile_ab), profile_rfcp(profile_rfcp),
        profile_index_aa(profile_index_aa), profile_index_ab(profile_index_ab),
        keep_profiles(keep_profiles), task_best(task_best) {}

  void operator()(std::size_t begin, std::size_t end) override {
    double const *const positive_ptr = positive.begin();
    double const *const negative_ptr = negative.begin();
    double const *const mean_positive_ptr = mean_positive.begin();
    double const *const mean_negative_ptr = mean_negative.begin();
    double const *const sig_positive_ptr = sig_positive.begin();
    double const *const sig_negative_ptr = sig_negative.begin();
    double const *const df_positive_ptr = df_positive.begin();
    double const *const df_negative_ptr = df_negative.begin();
    double const *const dg_positive_ptr = dg_positive.begin();
    double const *const dg_negative_ptr = dg_negative.begin();
    int const *const valid_positive_ptr = valid_positive.begin();
    int const *const valid_negative_ptr = valid_negative.begin();
    double *const rms_ptr = rms_profile.begin();
    int *const completed_ptr = completed.begin();
    int *const paired_frequency_ptr = paired_frequency.begin();
    double *const profile_aa_ptr = keep_profiles ? profile_aa.begin() : nullptr;
    double *const profile_ab_ptr = keep_profiles ? profile_ab.begin() : nullptr;
    double *const profile_rfcp_ptr = keep_profiles ? profile_rfcp.begin() : nullptr;
    int *const profile_index_aa_ptr = keep_profiles ? profile_index_aa.begin() : nullptr;
    int *const profile_index_ab_ptr = keep_profiles ? profile_index_ab.begin() : nullptr;
    uint32_t const positive_profiles = mean_positive.size();
    uint32_t const negative_profiles = mean_negative.size();
    double const clip = std::sqrt(2.0 * window_size);
    uint32_t const max_reference_profiles = std::max(positive_profiles, negative_profiles);

    // These buffers belong to one executor invocation and are deliberately
    // reused across its query blocks.  Reallocating them per block is costly
    // for large segmented streams, while the exclusion stamps let us avoid
    // clearing the full reference-profile-sized mask between blocks.
    std::vector<double> covariance_aa;
    std::vector<double> covariance_ab;
    std::vector<RFCPCandidate> candidates;
    candidates.reserve(max_reference_profiles);
    std::vector<uint32_t> exclusion_marks(max_reference_profiles, 0);
    uint32_t exclusion_stamp = 0;
    std::vector<double> aa_distances(max_freq, NA_REAL);
    std::vector<double> ab_distances(max_freq, NA_REAL);
    std::vector<double> rfcp_values(max_freq, NA_REAL);
    std::vector<int> aa_indices(max_freq, NA_INTEGER);
    std::vector<int> ab_indices(max_freq, NA_INTEGER);
    RFCPTaskBest best;

    for (std::size_t task_index = begin; task_index < end; task_index++) {
      RFCPQueryTask const &task = tasks[task_index];
      initialize_covariances(positive_ptr, positive_ptr, mean_positive_ptr,
                             mean_positive_ptr, task.anchor_begin, positive.size(),
                             window_size, covariance_aa);
      initialize_covariances(positive_ptr, negative_ptr, mean_positive_ptr,
                             mean_negative_ptr, task.anchor_begin, negative.size(),
                             window_size, covariance_ab);
      best.available = false;
      best.query_index = 0;
      best.rms = R_NegInf;

      for (uint32_t query_index = task.anchor_begin; query_index < task.output_end;
           query_index++) {
        if (query_index > task.anchor_begin) {
          advance_covariances(positive_ptr, positive_ptr, mean_positive_ptr,
                              mean_positive_ptr, df_positive_ptr, df_positive_ptr,
                              dg_positive_ptr, dg_positive_ptr, query_index,
                              positive_profiles, window_size, covariance_aa);
          advance_covariances(positive_ptr, negative_ptr, mean_positive_ptr,
                              mean_negative_ptr, df_positive_ptr, df_negative_ptr,
                              dg_positive_ptr, dg_negative_ptr, query_index,
                              negative_profiles, window_size, covariance_ab);
        }
        if (query_index < task.output_begin) {
          continue;
        }

        uint32_t const output_index = query_index - request_begin;
        if (!valid_positive_ptr[query_index]) {
          completed_ptr[output_index] = true;
          continue;
        }

        uint32_t const aa_count = select_neighbors(
            covariance_aa, sig_positive_ptr, valid_positive_ptr,
            sig_positive_ptr[query_index], query_index, true, window_size,
            exclusion_radius, max_freq, candidates, exclusion_marks,
            exclusion_stamp, aa_distances, aa_indices);
        uint32_t const ab_count = select_neighbors(
            covariance_ab, sig_negative_ptr, valid_negative_ptr,
            sig_positive_ptr[query_index], query_index, false, window_size,
            exclusion_radius, max_freq, candidates, exclusion_marks,
            exclusion_stamp, ab_distances, ab_indices);
        uint32_t const paired = std::min(aa_count, ab_count);
        paired_frequency_ptr[output_index] = static_cast<int>(paired);
        double squared_sum = 0.0;
        std::fill(rfcp_values.begin(), rfcp_values.end(), NA_REAL);
        for (uint32_t rank = 0; rank < paired; rank++) {
          double const aa = std::min(aa_distances[rank], clip);
          double const ab = std::min(ab_distances[rank], clip);
          double const contrast = std::max(0.0, (ab - aa) / clip);
          rfcp_values[rank] = contrast;
          squared_sum += contrast * contrast;
        }
        double const rms = std::sqrt(squared_sum / max_freq);
        rms_ptr[output_index] = rms;

        if (keep_profiles) {
          uint64_t const column_offset = static_cast<uint64_t>(output_index) * max_freq;
          for (uint32_t rank = 0; rank < max_freq; rank++) {
            uint64_t const cell = column_offset + rank;
            profile_aa_ptr[cell] = aa_distances[rank];
            profile_ab_ptr[cell] = ab_distances[rank];
            profile_rfcp_ptr[cell] = rfcp_values[rank];
            profile_index_aa_ptr[cell] = aa_indices[rank];
            profile_index_ab_ptr[cell] = ab_indices[rank];
          }
        }

        if (!best.available || rms > best.rms) {
          best.available = true;
          best.query_index = query_index;
          best.rms = rms;
          best.aa_distances = aa_distances;
          best.ab_distances = ab_distances;
          best.rfcp = rfcp_values;
          best.aa_indices = aa_indices;
          best.ab_indices = ab_indices;
        }
        completed_ptr[output_index] = true;
      }
      task_best[task_index] = best;
    }
  }
};

} // namespace

// query_begin is zero-based and query_end is exclusive.  The public R wrapper
// exposes one-based inclusive coordinates and performs that normalization.
// [[Rcpp::export]]
List rfcp_na_segmented_native_rcpp_parallel(NumericVector positive_ref,
                                             NumericVector negative_ref,
                                             uint64_t window_size,
                                             uint32_t max_freq,
                                             uint32_t exclusion_radius,
                                             uint64_t query_begin,
                                             uint64_t query_end,
                                             bool return_profiles,
                                             bool progress) {
  uint64_t const positive_size = positive_ref.length();
  uint64_t const negative_size = negative_ref.length();
  if (window_size < 2 || window_size > positive_size || window_size > negative_size) {
    Rcpp::stop("window_size must fit both input series and be at least 2");
  }
  if (positive_size > std::numeric_limits<uint32_t>::max() ||
      negative_size > std::numeric_limits<uint32_t>::max()) {
    Rcpp::stop("RFCP currently supports input series shorter than 2^32 samples");
  }
  if (max_freq < 1) {
    Rcpp::stop("max_freq must be at least 1");
  }

  uint32_t const positive_profiles = positive_size - window_size + 1;
  uint32_t const negative_profiles = negative_size - window_size + 1;
  if (positive_profiles > static_cast<uint64_t>(std::numeric_limits<int>::max()) ||
      negative_profiles > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
    Rcpp::stop("RFCP window indices must fit in an R integer");
  }
  if (query_begin >= query_end || query_end > positive_profiles) {
    Rcpp::stop("query range must be a non-empty subset of positive window starts");
  }
  uint32_t const request_begin = query_begin;
  uint32_t const request_end = query_end;
  uint32_t const query_count = request_end - request_begin;
  uint64_t const profile_cells = static_cast<uint64_t>(query_count) * max_freq;
  if (return_profiles && profile_cells > rfcp_max_profile_cells) {
    Rcpp::stop("return_profiles exceeds the limit of 5,000,000 rank-query cells");
  }
  (void)progress;

  List const positive_stats = muinvn_na_parallel(positive_ref, window_size);
  List const negative_stats = muinvn_na_parallel(negative_ref, window_size);
  NumericVector positive = positive_stats["data"];
  NumericVector negative = negative_stats["data"];
  NumericVector mean_positive = positive_stats["avg"];
  NumericVector mean_negative = negative_stats["avg"];
  NumericVector sig_positive = positive_stats["sig"];
  NumericVector sig_negative = negative_stats["sig"];
  LogicalVector valid_positive = positive_stats["valid_window"];
  LogicalVector valid_negative = negative_stats["valid_window"];

  NumericVector df_positive(positive_profiles);
  NumericVector dg_positive(positive_profiles);
  NumericVector df_negative(negative_profiles);
  NumericVector dg_negative(negative_profiles);
  rfcp_fill_differences(positive, mean_positive, window_size, df_positive, dg_positive);
  rfcp_fill_differences(negative, mean_negative, window_size, df_negative, dg_negative);

  uint32_t const rows_per_task = rfcp_query_rows_per_task(window_size);
  uint32_t const first_anchor = (request_begin / rows_per_task) * rows_per_task;
  std::vector<RFCPQueryTask> tasks;
  for (uint64_t anchor64 = first_anchor; anchor64 < request_end;
       anchor64 += rows_per_task) {
    uint32_t const anchor = static_cast<uint32_t>(anchor64);
    uint32_t const block_end = static_cast<uint32_t>(std::min<uint64_t>(
        positive_profiles, anchor64 + rows_per_task));
    uint32_t const output_begin = std::max(request_begin, anchor);
    uint32_t const output_end = std::min(request_end, block_end);
    if (output_begin < output_end) {
      tasks.push_back(RFCPQueryTask{anchor, output_begin, output_end});
    }
  }

  NumericVector rms_profile(query_count, NA_REAL);
  LogicalVector completed(query_count, false);
  IntegerVector paired_frequency(query_count, 0);
  NumericVector profile_aa;
  NumericVector profile_ab;
  NumericVector profile_rfcp;
  IntegerVector profile_index_aa;
  IntegerVector profile_index_ab;
  if (return_profiles) {
    profile_aa = NumericVector(profile_cells, NA_REAL);
    profile_ab = NumericVector(profile_cells, NA_REAL);
    profile_rfcp = NumericVector(profile_cells, NA_REAL);
    profile_index_aa = IntegerVector(profile_cells, NA_INTEGER);
    profile_index_ab = IntegerVector(profile_cells, NA_INTEGER);
  }

  std::vector<RFCPTaskBest> task_best(tasks.size());
  RFCPWorker worker(positive, negative, mean_positive, mean_negative,
                    sig_positive, sig_negative, df_positive, df_negative,
                    dg_positive, dg_negative, valid_positive, valid_negative,
                    tasks, window_size, max_freq, exclusion_radius, request_begin,
                    rms_profile, completed, paired_frequency, profile_aa, profile_ab,
                    profile_rfcp, profile_index_aa, profile_index_ab,
                    return_profiles, task_best);

  bool partial = false;
  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(0, tasks.size(), worker, 1);
#else
    RcppParallel2::ttParallelFor(0, tasks.size(), worker, 1);
#endif
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
  } catch (...) {
    Rcpp::stop("c++ exception while computing RFCP");
  }

  uint32_t completed_count = 0;
  while (completed_count < query_count && completed[completed_count]) {
    completed_count++;
  }
  if (completed_count < query_count) {
    partial = true;
    for (uint32_t i = completed_count; i < query_count; i++) {
      rms_profile[i] = NA_REAL;
      paired_frequency[i] = 0;
      completed[i] = false;
    }
    if (return_profiles) {
      for (uint64_t cell = static_cast<uint64_t>(completed_count) * max_freq;
           cell < profile_cells; cell++) {
        profile_aa[cell] = NA_REAL;
        profile_ab[cell] = NA_REAL;
        profile_rfcp[cell] = NA_REAL;
        profile_index_aa[cell] = NA_INTEGER;
        profile_index_ab[cell] = NA_INTEGER;
      }
    }
  }

  uint32_t effective_max_freq = 0;
  bool plato_available = false;
  uint32_t plato_query = 0;
  double plato_rms = R_NegInf;
  for (uint32_t i = 0; i < completed_count; i++) {
    effective_max_freq = std::max<uint32_t>(effective_max_freq, paired_frequency[i]);
    if (std::isfinite(rms_profile[i]) &&
        (!plato_available || rms_profile[i] > plato_rms)) {
      plato_available = true;
      plato_query = request_begin + i;
      plato_rms = rms_profile[i];
    }
  }

  NumericVector rfmp_aa_at_plato(max_freq, NA_REAL);
  NumericVector rfmp_ab_at_plato(max_freq, NA_REAL);
  NumericVector rfcp_at_plato(max_freq, NA_REAL);
  IntegerVector rfmpi_aa_at_plato(max_freq, NA_INTEGER);
  IntegerVector rfmpi_ab_at_plato(max_freq, NA_INTEGER);
  if (plato_available) {
    for (RFCPTaskBest const &best : task_best) {
      if (!best.available || best.query_index != plato_query) {
        continue;
      }
      std::copy(best.aa_distances.begin(), best.aa_distances.end(),
                rfmp_aa_at_plato.begin());
      std::copy(best.ab_distances.begin(), best.ab_distances.end(),
                rfmp_ab_at_plato.begin());
      std::copy(best.rfcp.begin(), best.rfcp.end(), rfcp_at_plato.begin());
      std::copy(best.aa_indices.begin(), best.aa_indices.end(),
                rfmpi_aa_at_plato.begin());
      std::copy(best.ab_indices.begin(), best.ab_indices.end(),
                rfmpi_ab_at_plato.begin());
      break;
    }
  }

  List result = List::create(
      Rcpp::Named("rms_profile") = rms_profile,
      Rcpp::Named("valid_window_positive") = valid_positive,
      Rcpp::Named("valid_window_negative") = valid_negative,
      Rcpp::Named("plato_index") = plato_available
                                             ? Rcpp::wrap(static_cast<double>(plato_query + 1))
                                             : Rcpp::wrap(NA_REAL),
      Rcpp::Named("plato_rms") = plato_available ? plato_rms : NA_REAL,
      Rcpp::Named("plato_twin_index") = plato_available && rfmpi_aa_at_plato[0] != NA_INTEGER
                                                  ? Rcpp::wrap(rfmpi_aa_at_plato[0])
                                                  : Rcpp::wrap(NA_INTEGER),
      Rcpp::Named("rfcp_at_plato") = rfcp_at_plato,
      Rcpp::Named("rfmp_aa_at_plato") = rfmp_aa_at_plato,
      Rcpp::Named("rfmp_ab_at_plato") = rfmp_ab_at_plato,
      Rcpp::Named("rfmpi_aa_at_plato") = rfmpi_aa_at_plato,
      Rcpp::Named("rfmpi_ab_at_plato") = rfmpi_ab_at_plato,
      Rcpp::Named("requested_max_freq") = static_cast<int>(max_freq),
      Rcpp::Named("effective_max_freq") = static_cast<int>(effective_max_freq),
      Rcpp::Named("query_begin") = static_cast<double>(request_begin + 1),
      Rcpp::Named("query_end") = static_cast<double>(request_end),
      Rcpp::Named("completed_query_end") = completed_count > 0
                                                       ? Rcpp::wrap(static_cast<double>(request_begin + completed_count))
                                                       : Rcpp::wrap(NA_REAL),
      Rcpp::Named("partial") = partial,
      Rcpp::Named("query_rows_per_task") = static_cast<int>(rows_per_task));

  if (return_profiles) {
    IntegerVector const dimensions = IntegerVector::create(max_freq, query_count);
    profile_aa.attr("dim") = dimensions;
    profile_ab.attr("dim") = dimensions;
    profile_rfcp.attr("dim") = dimensions;
    profile_index_aa.attr("dim") = dimensions;
    profile_index_ab.attr("dim") = dimensions;
    result["rfmp_aa"] = profile_aa;
    result["rfmp_ab"] = profile_ab;
    result["rfcp"] = profile_rfcp;
    result["rfmpi_aa"] = profile_index_aa;
    result["rfmpi_ab"] = profile_index_ab;
  }

  return result;
}
