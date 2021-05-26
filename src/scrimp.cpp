#include "math.h" // math first to fix OSX error
#include "scrimp.h"
#include "mass.h"
#include "windowfunc.h"
#include <numeric>
#include <cfloat> // DBL_EPSILON when STRICT_R_HEADERS
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

// [[Rcpp::export]]
List scrimp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                 double s_size, double pre_scrimp, bool progress) {

  bool partial = false;
  uint32_t exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint32_t data_size = data_ref.length();

  uint32_t matrix_profile_size = data_size - window_size + 1;

  // TODO: check skip position (DBL_EPSILON, etc)
  LogicalVector skip_location(matrix_profile_size, 0);

  for (uint64_t i = 0; i < matrix_profile_size; i++) {
    NumericVector range = data_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_location[i] = TRUE;
    }
  }

  NumericVector data = data_ref;
  NumericVector query = query_ref;

  data[is_na(data)] = 0;
  data[is_infinite(data)] = 0;
  query[is_na(query)] = 0;
  query[is_infinite(query)] = 0;

  NumericVector matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector profile_index(matrix_profile_size, -1);

  IntegerVector orig_index = Range(0, matrix_profile_size - 1);

  uint32_t k = set_k_rcpp(window_size, data_size, window_size);

  List pre = mass_pre_rcpp(data, query, window_size);

  NumericVector data_mean = pre["data_mean"];
  NumericVector data_sd = pre["data_sd"];
  NumericVector query_mean = pre["query_mean"];
  NumericVector query_sd = pre["query_sd"];

  try {
    //// PRE-SCRIMP ----

    if (pre_scrimp > 0) {
      // initialization
      int64_t current_step = floor(window_size * pre_scrimp + DBL_EPSILON);
      IntegerVector pre_scrimp_idxs = seq_by(0, matrix_profile_size - 1, current_step);
      Progress ps(pre_scrimp_idxs.size(), progress);
      // compute the matrix profile
      NumericVector dotproduct(matrix_profile_size);
      NumericVector refine_distance(matrix_profile_size, R_PosInf);

      int64_t j = 1;
      for (int64_t &&i : pre_scrimp_idxs) {

        RcppThread::checkUserInterrupt();
        ps.increment();

        // compute the distance profile
        List nn = mass3_rcpp(query[Range(i, i + window_size - 1)], data, pre["data_size"], pre["window_size"],
                             data_mean, data_sd, query_mean[i], query_sd[i], k);

        NumericVector distance_profile = as<NumericVector>(nn["distance_profile"]);

        int64_t exc_st = 0;
        int64_t exc_ed = 0;

        // apply exclusion zone
        if (exclusion_zone > 0) {
          exc_st = MAX(0, i > exclusion_zone ? (i - exclusion_zone) : 0);
          exc_ed = MIN(matrix_profile_size - 1, i + exclusion_zone);
        }

        IntegerVector exc_idxs = Range(exc_st, exc_ed);
        distance_profile[exc_idxs] = R_PosInf;

        distance_profile[data_sd < DBL_EPSILON] = R_PosInf;
        if (skip_location[i] || query_sd[i] < DBL_EPSILON) {
          distance_profile.fill(R_PosInf);
        }
        distance_profile[skip_location] = R_PosInf;

        // figure out and store the neareest neighbor
        if (j == 1) {
          matrix_profile = distance_profile;
          profile_index.fill(i);
          int64_t min_idx = which_min(distance_profile);
          profile_index[i] = min_idx;
          matrix_profile[i] = distance_profile[min_idx];
        } else {
          LogicalVector update_pos = distance_profile < matrix_profile;
          profile_index[update_pos] = (int32_t)i;
          matrix_profile[update_pos] = distance_profile[update_pos];
          int64_t min_idx = which_min(distance_profile);
          profile_index[i] = min_idx;
          matrix_profile[i] = distance_profile[min_idx];
        }

        int64_t idx_nn = profile_index[i];
        int64_t idx_diff = (idx_nn - i);
        dotproduct[i] = (window_size - (matrix_profile[i] / 2)) * data_sd[i] * data_sd[idx_nn] +
                        window_size * data_mean[i] * data_mean[idx_nn];

        int64_t endidx = MIN(matrix_profile_size - 1, (int64_t)i + current_step - 1);
        endidx = MIN((int64_t)endidx, matrix_profile_size - idx_diff - 1);

        if (i < endidx) {
          // confirmed, sequences asc
          Range dot_idxs1 = Range((i + 1), endidx);
          Range dot_idxs2 = Range((i + window_size), (endidx + window_size - 1));
          Range dot_idxs3 = Range((idx_nn + window_size), (endidx + window_size - 1 + idx_diff));
          Range dot_idxs4 = Range(i, (endidx - 1));
          Range dot_idxs5 = Range(idx_nn, (endidx - 1 + idx_diff));
          dotproduct[dot_idxs1] = (NumericVector)(
              (NumericVector)(cumsum(data[dot_idxs2] * data[dot_idxs3] - data[dot_idxs4] * data[dot_idxs5])) +
              dotproduct[i]);

          // confirmed, sequences asc
          Range ref_idxs1 = Range((idx_nn + 1), (endidx + idx_diff));
          refine_distance[dot_idxs1] =
              2 * (window_size - (dotproduct[dot_idxs1] - data_mean[dot_idxs1] * data_mean[ref_idxs1] * window_size) /
                                     (data_sd[dot_idxs1] * data_sd[ref_idxs1]));
        }

        int64_t beginidx = ((i + 1) <= current_step) ? 0 : (i + 1 - current_step);
        if (idx_diff < 0) {
          beginidx = MAX(beginidx, abs(idx_diff));
        }

        if (i > 0 && i > beginidx) {
          // sequences reversed
          IntegerVector dot_rev_idxs1 = ::seq((i - 1), beginidx);
          IntegerVector dot_rev_idxs2 = ::seq((idx_nn - 1), (beginidx + idx_diff));
          IntegerVector dot_rev_idxs3 = ::seq((i - 1 + window_size), (beginidx + window_size));
          IntegerVector dot_rev_idxs4 = ::seq((idx_nn - 1 + window_size), (beginidx + idx_diff + window_size));
          dotproduct[dot_rev_idxs1] = (NumericVector)((NumericVector)cumsum(data[dot_rev_idxs1] * data[dot_rev_idxs2] -
                                                                            data[dot_rev_idxs3] * data[dot_rev_idxs4]) +
                                                      dotproduct[i]);

          Range ref_idxs2 = Range(beginidx, (i - 1));
          Range ref_idxs3 = Range((beginidx + idx_diff), (idx_nn - 1));

          refine_distance[ref_idxs2] =
              2 * (window_size - (dotproduct[ref_idxs2] - data_mean[ref_idxs2] * data_mean[ref_idxs3] * window_size) /
                                     (data_sd[ref_idxs2] * data_sd[ref_idxs3]));
        }

        LogicalVector rd = (refine_distance < 0);

        if (sum(as<IntegerVector>(rd)) > 0) {
          Rcout << "Debug: refine_distance < 0" << std::endl;
        }

        refine_distance[rd] = 0;

        Range upd_idxs1 = Range(beginidx, endidx);
        Range upd_idxs2 = Range((beginidx + idx_diff), (endidx + idx_diff));
        IntegerVector update_pos1 = which_cpp(refine_distance[upd_idxs1] < matrix_profile[upd_idxs1]);
        matrix_profile[(update_pos1 + beginidx)] = refine_distance[(update_pos1 + beginidx)];
        IntegerVector new_idxs = as<IntegerVector>(orig_index[(update_pos1 + beginidx)]) + idx_diff;
        profile_index[(update_pos1 + beginidx)] = new_idxs;
        IntegerVector update_pos2 = which_cpp(refine_distance[upd_idxs1] < matrix_profile[upd_idxs2]);
        matrix_profile[(update_pos2 + beginidx + idx_diff)] = refine_distance[(update_pos2 + beginidx)];
        new_idxs = as<IntegerVector>(orig_index[(update_pos2 + beginidx + idx_diff)]) - idx_diff;
        profile_index[(update_pos2 + beginidx + idx_diff)] = new_idxs;

        j++;
      }
    }

    //// SCRIMP ----

    IntegerVector compute_order = orig_index[orig_index > exclusion_zone];

    NumericVector curlastz(matrix_profile_size);
    NumericVector curdistance(matrix_profile_size);
    NumericVector dist1(matrix_profile_size, R_PosInf);
    NumericVector dist2(matrix_profile_size, R_PosInf);

    Progress p(compute_order.size(), progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    if (s_size < 1.0) {
      stop = round(compute_order.size() * s_size + DBL_EPSILON);
    }

    uint64_t j = 1;
    for (int64_t &&i : compute_order) {

      RcppThread::checkUserInterrupt();
      p.increment();

      curlastz[i] = sum(data[Range(0, window_size - 1)] * data[Range(i, i + window_size - 1)]);

      if (i < (matrix_profile_size - 1)) {
        curlastz[Range(i + 1, matrix_profile_size - 1)] =
            (NumericVector)cumsum(
                data[Range(window_size, data_size - i - 1)] * data[Range(i + window_size, data_size - 1)] -
                data[Range(0, matrix_profile_size - i - 2)] * data[Range(i, matrix_profile_size - 2)]) +
            curlastz[i];
      }
      curdistance[Range(i, matrix_profile_size - 1)] =
          2 *
          (window_size -
           (curlastz[Range(i, matrix_profile_size - 1)] - window_size * data_mean[Range(i, matrix_profile_size - 1)] *
                                                              data_mean[Range(0, matrix_profile_size - i - 1)]) /
               (data_sd[Range(i, matrix_profile_size - 1)] * data_sd[Range(0, matrix_profile_size - i - 1)]));

      LogicalVector cd = (curdistance < 0);

      if (sum(as<IntegerVector>(cd)) > 0) {
        Rcout << "Debug: curdistance < 0" << std::endl;
      }

      curdistance[cd] = 0;

      dist1[::seq(0, i - 1)] = R_PosInf;
      dist1[::seq(i, matrix_profile_size - 1)] = curdistance[::seq(i, matrix_profile_size - 1)];

      dist2[::seq(0, matrix_profile_size - i - 1)] = curdistance[::seq(i, matrix_profile_size - 1)];
      dist2[::seq(matrix_profile_size - i, matrix_profile_size - 1)] = R_PosInf;

      LogicalVector loc1 = dist1 < matrix_profile;
      matrix_profile[loc1] = dist1[loc1];
      profile_index[loc1] = (IntegerVector)(as<IntegerVector>(orig_index[loc1]) - i);

      LogicalVector loc2 = dist2 < matrix_profile;
      matrix_profile[loc2] = dist2[loc2];
      profile_index[loc2] = (IntegerVector)(as<IntegerVector>(orig_index[loc2]) + i);

      if (stop > 0 && j >= stop) {
        partial = true;
        break;
      }
      j++;
    }

  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  matrix_profile = sqrt(matrix_profile);
  profile_index = profile_index + 1;

  return (List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("profile_index") = profile_index,
                       Rcpp::Named("partial") = partial, Rcpp::Named("ez") = ez));
}

struct ScrimpWorker : public Worker {
  // input
  const RVector<double> data_ref;
  const RVector<double> window_ref;
  const uint64_t w_size;
  const uint64_t d_size;
  const RVector<double> d_mean;
  const RVector<double> d_std;
  const RVector<int> skip_location;

  Progress *p;

  RVector<double> mp;
  RVector<int> pi;

#if RCPP_PARALLEL_USE_TBB
  tbb::spin_mutex m;
#else
  tthread::mutex m;
#endif

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  ScrimpWorker(const NumericVector data_ref, const NumericVector window_ref, const uint64_t w_size,
               const uint64_t d_size, const NumericVector d_mean, const NumericVector d_std,
               const IntegerVector skip_location, Progress *p, NumericVector mp, IntegerVector pi)
      : data_ref(data_ref), window_ref(window_ref), w_size(w_size), d_size(d_size), d_mean(d_mean), d_std(d_std),
        skip_location(skip_location), p(p), mp(mp), pi(pi) {}

  ~ScrimpWorker() {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    // begin and end are the compute_order
    //// SCRIMP ----
    uint64_t dp_size = d_size - w_size + 1;

    std::vector<double> curlastz(dp_size);
    std::vector<double> curdistance(dp_size);
    std::vector<double> dist1(dp_size, R_PosInf);
    std::vector<double> dist2(dp_size, R_PosInf);

    // index of sliding window
    try {
      for (std::size_t i = begin; i < end; i++) {

        if (i % 10 == 0) {
          RcppThread::checkUserInterrupt();
          m.lock();
          p->increment();
          m.unlock();
        }

        std::vector<double> lastz(w_size);
        std::transform(data_ref.begin(), data_ref.begin() + w_size, data_ref.begin() + i, lastz.begin(),
                       std::multiplies<double>());

        curlastz[i] = std::accumulate(lastz.begin(), lastz.end(), decltype(lastz)::value_type(0)); // just a sum(vector)

        if (i < (dp_size - 1)) {

          uint64_t j = i + 1;
          double cum = 0.0;
          for (uint64_t k = w_size; k <= d_size - i - 1; k++) {
            cum = cum + data_ref[k] * data_ref[i + k] - data_ref[k - w_size] * data_ref[i + k - w_size];
            curlastz[j++] = cum + curlastz[i];
          }
        }

        for (uint64_t j = i; j <= dp_size - 1; j++) {
          curdistance[j] =
              2 * (w_size - (curlastz[j] - w_size * d_mean[j] * d_mean[j - i]) / (d_std[j] * d_std[j - i]));

          if (curdistance[j] < 0) {
            curdistance[j] = 0;
            Rcout << "Debug: curdistance < 0" << std::endl;
          }
        }

        for (uint64_t j = 0; j <= i - 1; j++) {
          dist1[j] = R_PosInf;
        }

        for (uint64_t j = i; j <= dp_size - 1; j++) {
          dist1[j] = curdistance[j];
        }

        for (uint64_t j = 0; j <= dp_size - i - 1; j++) {
          dist2[j] = curdistance[j + i];
        }

        for (uint64_t j = dp_size - i; j <= dp_size - 1; j++) {
          dist2[j] = R_PosInf;
        }

        for (uint64_t j = 0; j < dp_size; j++) {
          m.lock();
          if (dist1[j] < mp[j]) {
            mp[j] = dist1[j];
            pi[j] = j - i;
          }
          if (dist2[j] < mp[j]) {
            mp[j] = dist2[j];
            pi[j] = j + i;
          }
          m.unlock();
        }
      }

    } catch (RcppThread::UserInterruptException &e) {
      Rcout << "Computation interrupted by the user." << std::endl;
      Rcout << "Please wait for other threads to stop." << std::endl;
      throw;
    }
  }
};

// [[Rcpp::export]]
List scrimp_rcpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                          double s_size, bool progress) {
  uint64_t exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint64_t data_size = data_ref.length();
  // uint64_t query_size = query_ref.length();
  uint64_t matrix_profile_size = data_size - window_size + 1;
  bool partial = false;

  // TODO: check skip position (DBL_EPSILON, etc)
  IntegerVector skip_location(matrix_profile_size, 0);

  for (uint64_t i = 0; i < matrix_profile_size; i++) {
    NumericVector range = data_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_location[i] = 1;
    }
  }

  NumericVector data = data_ref;
  NumericVector query = query_ref;

  data[is_na(data)] = 0;
  data[is_infinite(data)] = 0;
  query[is_na(query)] = 0;
  query[is_infinite(query)] = 0;

  NumericVector matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector profile_index(matrix_profile_size, -1);

  // NumericVector first_product = rnn["last_product"];

  List pre = mass_pre_rcpp(data, query, window_size);

  Progress p((matrix_profile_size - exclusion_zone - 1) / 10, progress);

  ScrimpWorker scrimp_worker(data, query, window_size, data_size, pre["data_mean"], pre["data_sd"], skip_location, &p,
                             matrix_profile, profile_index);

  // call parallelFor to do the work
  try {
#if RCPP_PARALLEL_USE_TBB
    RcppParallel::parallelFor(exclusion_zone + 1, matrix_profile_size, scrimp_worker);
#else
    RcppParallel2::ttParallelFor(exclusion_zone + 1, matrix_profile_size, scrimp_worker);
#endif

  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results "
             "were returned.";
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("matrix_profile") = sqrt(matrix_profile),
                       Rcpp::Named("profile_index") = profile_index + 1, Rcpp::Named("partial") = partial,
                       Rcpp::Named("ez") = ez));
}

// [[Rcpp::export]]
List scrimpab_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double s_size,
                   bool progress) {
  // double s_size = R_PosInf;
  bool partial = false;
  uint32_t data_size = data_ref.length();
  uint32_t query_size = query_ref.length();

  uint32_t mmpa_size = data_size - window_size + 1;
  uint32_t mmpb_size = query_size - window_size + 1;

  // TODO: check skip position (DBL_EPSILON, etc)
  LogicalVector skip_location(mmpa_size, 0);

  for (uint64_t i = 0; i < mmpa_size; i++) {
    NumericVector range = data_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_location[i] = TRUE;
    }
  }

  NumericVector data = data_ref;
  NumericVector query = query_ref;

  data[is_na(data)] = 0;
  data[is_infinite(data)] = 0;
  query[is_na(query)] = 0;
  query[is_infinite(query)] = 0;

  NumericVector mmpa(mmpa_size, R_PosInf);
  IntegerVector mpia(mmpa_size, -1);
  NumericVector mmpb(mmpb_size, R_PosInf);
  IntegerVector mpib(mmpb_size, -1);

  IntegerVector orig_index = Range(0, mmpa_size - 1);

  List dd = movmean_std_rcpp(data, window_size);
  List qq = movmean_std_rcpp(query, window_size);

  NumericVector data_mean = dd["avg"];
  NumericVector data_sd = dd["sd"];
  NumericVector query_mean = qq["avg"];
  NumericVector query_sd = qq["sd"];

  Progress p(orig_index.size() * 2 - 2, progress);

  try {
    //// SCRIMP ----

    IntegerVector compute_order = orig_index[orig_index > 1];

    NumericVector curlastz(mmpb_size);
    NumericVector curdistance(mmpb_size);
    NumericVector dist1(mmpb_size, R_PosInf);
    NumericVector dist2(mmpb_size, R_PosInf);

    int64_t j = 1;
    for (int64_t &&i : compute_order) {

      RcppThread::checkUserInterrupt();
      p.increment();

      curlastz[i] = sum(data[Range(0, window_size - 1)] * query[Range(i, i + window_size - 1)]);

      if (i < (mmpb_size - 1)) {
        curlastz[Range(i + 1, mmpb_size - 1)] =
            (NumericVector)cumsum(data[Range(window_size, data_size - i - 1)] *
                                      query[Range(i + window_size, query_size - 1)] -
                                  data[Range(0, mmpb_size - i - 2)] * query[Range(i, mmpb_size - 2)]) +
            curlastz[i];
      }
      curdistance[Range(i, mmpb_size - 1)] =
          2 * (window_size - (curlastz[Range(i, mmpb_size - 1)] - window_size * query_mean[Range(i, mmpb_size - 1)] *
                                                                      data_mean[Range(0, mmpb_size - i - 1)]) /
                                 (query_sd[Range(i, mmpb_size - 1)] * data_sd[Range(0, mmpb_size - i - 1)]));

      LogicalVector cd = (curdistance < 0);

      if (sum(as<IntegerVector>(cd)) > 0) {
        Rcout << "Debug: curdistance < 0" << std::endl;
      }

      curdistance[cd] = 0;

      dist1[::seq(0, i - 1)] = R_PosInf;
      dist1[::seq(i, mmpb_size - 1)] = curdistance[::seq(i, mmpb_size - 1)];
      LogicalVector loc1 = dist1 < mmpb;
      mmpb[loc1] = dist1[loc1];
      mpib[loc1] = (IntegerVector)(as<IntegerVector>(orig_index[loc1]) - i);

      dist2[::seq(0, mmpb_size - i - 1)] = curdistance[::seq(i, mmpb_size - 1)];
      dist2[::seq(mmpb_size - i, mmpb_size - 1)] = R_PosInf;
      LogicalVector loc2 = dist2 < mmpa;
      mmpa[loc2] = dist2[loc2];
      mpia[loc2] = (IntegerVector)(as<IntegerVector>(orig_index[loc2]) + i);

      j++;
    }

  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  try {
    //// SCRIMP ----

    IntegerVector compute_order = orig_index[orig_index > 1];

    NumericVector curlastz(mmpa_size);
    NumericVector curdistance(mmpa_size);
    NumericVector dist1(mmpa_size, R_PosInf);
    NumericVector dist2(mmpa_size, R_PosInf);

    int64_t j = 1;
    for (int64_t &&i : compute_order) {

      RcppThread::checkUserInterrupt();
      p.increment();

      curlastz[i] = sum(query[Range(0, window_size - 1)] * data[Range(i, i + window_size - 1)]);

      if (i < (mmpa_size - 1)) {
        curlastz[Range(i + 1, mmpa_size - 1)] =
            (NumericVector)cumsum(query[Range(window_size, data_size - i - 1)] *
                                      data[Range(i + window_size, query_size - 1)] -
                                  query[Range(0, mmpa_size - i - 2)] * data[Range(i, mmpa_size - 2)]) +
            curlastz[i];
      }
      curdistance[Range(i, mmpa_size - 1)] =
          2 * (window_size - (curlastz[Range(i, mmpa_size - 1)] - window_size * data_mean[Range(i, mmpa_size - 1)] *
                                                                      query_mean[Range(0, mmpa_size - i - 1)]) /
                                 (data_sd[Range(i, mmpa_size - 1)] * query_sd[Range(0, mmpa_size - i - 1)]));

      LogicalVector cd = (curdistance < 0);

      if (sum(as<IntegerVector>(cd)) > 0) {
        Rcout << "Debug: curdistance < 0" << std::endl;
      }

      curdistance[cd] = 0;

      dist1[::seq(0, i - 1)] = R_PosInf;
      dist1[::seq(i, mmpa_size - 1)] = curdistance[::seq(i, mmpa_size - 1)];
      LogicalVector loc1 = dist1 < mmpa;
      mmpa[loc1] = dist1[loc1];
      mpia[loc1] = (IntegerVector)(as<IntegerVector>(orig_index[loc1]) - i);

      dist2[::seq(0, mmpa_size - i - 1)] = curdistance[::seq(i, mmpa_size - 1)];
      dist2[::seq(mmpa_size - i, mmpa_size - 1)] = R_PosInf;
      LogicalVector loc2 = dist2 < mmpb;
      mmpb[loc2] = dist2[loc2];
      mpib[loc2] = (IntegerVector)(as<IntegerVector>(orig_index[loc2]) + i);

      j++;
    }

  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  mmpa = sqrt(mmpa);
  mmpb = sqrt(mmpb);
  mpia = mpia + 1;
  mpib = mpib + 1;

  return (List::create(Rcpp::Named("matrix_profile") = mmpa, Rcpp::Named("profile_index") = mpia,
                       Rcpp::Named("mpb") = mmpb, Rcpp::Named("pib") = mpib, Rcpp::Named("partial") = partial,
                       Rcpp::Named("ez") = 0));
}
