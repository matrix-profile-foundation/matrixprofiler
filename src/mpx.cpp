#include "math.h" // math first to fix OSX error
#include "mpx.h"
#include "windowfunc.h"
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

// MPX stream version for the right side, getting stats from input, in review

List mpxis_rcpp(NumericVector data_ref, uint64_t batch_size, List object, List stats, uint64_t history,
                uint64_t mp_time_constraint, bool progress, float threshold) {

  try {
    double c, c_cmp, ez;
    uint32_t off_min, off_diag, offset, off_start;
    uint32_t n, window_size, exclusion_zone, profile_len;
    uint64_t mp_offset;
    NumericVector data, ddf_t, ddg_t;
    bool initial = false;
    bool partial = false;

    mp_offset = (uint64_t)object["offset"] + batch_size;
    window_size = (uint32_t)object["w"];
    ez = (double)object["ez"];
    exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;
    n = data_ref.length();
    profile_len = n - window_size + 1;

    if (n == batch_size) {
      initial = true;
    }

    // matrix profile using cross correlation,

    NumericVector mmu = as<NumericVector>(stats["avg"]);  // shallow copy
    NumericVector ssig = as<NumericVector>(stats["sig"]); // shallow copy
    NumericVector ddf = as<NumericVector>(stats["ddf"]);  // shallow copy
    NumericVector ddg = as<NumericVector>(stats["ddg"]);  // shallow copy

    double *mu = &mmu[0];
    double *sig = &ssig[0];
    double *df = &ddf[0];
    double *dg = &ddg[0];

    uint32_t diag_start = 0;
    uint32_t diag_end = n - window_size - exclusion_zone;

    if (mp_time_constraint > 0) {
      diag_start = n - mp_time_constraint - window_size;
    }

    IntegerVector compute_order = Range(diag_start, diag_end);

    NumericVector rmmp(profile_len, -1.0);
    IntegerVector rmmpi(profile_len, -1);

    if (!initial) {
      // copy the current profiles
      rmmp[Range(0, n - window_size - batch_size)] = as<NumericVector>(object["right_matrix_profile"]);
      rmmpi[Range(0, n - window_size - batch_size)] = as<IntegerVector>(object["right_profile_index"]);
    }

    double *rmp = &rmmp[0];
    int *rpi = &rmmpi[0];

    uint64_t num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress p(100, progress);

    // compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    try {
      uint64_t i = 1;
      // first window demeaned
      NumericVector ww = (data_ref[Range(n - window_size, n - 1)] - mmu[n - window_size]);

      for (int32_t diag : compute_order) {

        if ((i % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          p.increment();
        }
        // this always use the first window (ww); c is the sum of element-wise product
        c = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        // off_max goes from profile_len - ez to 0
        if (initial) {
          off_min = n - window_size - diag - 1;
        } else {
          off_min = MAX(n - window_size - batch_size, n - window_size - diag - 1);
        }
        off_start = n - window_size;

        // loop if batch_size > 1
        for (offset = off_start; offset > off_min; offset--) {
          // min is offset + diag; max is (profile_len - 1); each iteration has the size of off_max
          off_diag = offset - (n - window_size - diag);

          c = c + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = c * sig[offset] * sig[off_diag];

          // RMP
          if (threshold > 0) {
            // v1 ~ v2 doesn't matter for FLOSS, since it only cares about the rpi
            // if (c_cmp >= threshold) {
            //   if (c_cmp > rmp[off_diag]) {
            //     rmp[off_diag] = c_cmp;
            //   }
            //   rpi[off_diag] = offset + 1;
            // }

            if (threshold >= 10) {
              threshold /= 100;
              // hard
              // v2
              if (c_cmp >= threshold) {
                rmp[off_diag] = c_cmp;
                rpi[off_diag] = offset + 1;
              }
            } else {
              // soft
              // v3
              if (c_cmp >= threshold) {
                if (c_cmp > rmp[off_diag]) {
                  rmp[off_diag] = c_cmp;
                  rpi[off_diag] = offset + 1;
                }
              }
            }
          } else {
            if (c_cmp > rmp[off_diag]) {
              rmp[off_diag] = c_cmp;
              rpi[off_diag] = offset + 1;
            }
          }
        }

        if (stop > 0 && i >= stop) {
          partial = true;
          break;
        }
        i++;
      }
    } catch (RcppThread::UserInterruptException &ex) {
      partial = true;
      Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
    }

    // to do ed
    rmmp[rmmp > 1.0] = 1.0;

    // if (euclidean) { // correlation to ed
    //   rmmp = sqrt(2 * window_size * (1 - rmmp));
    //   rmmp[rmmpi < 0] = R_PosInf;
    // }

    if (history > 0 && (n > history)) {
      // data_ref = tail(data_ref, history);
      uint64_t mp_new_size = history - window_size + 1;
      uint32_t diff = n - history;
      // ddf = tail(ddf, mp_new_size);
      // ddg = tail(ddg, mp_new_size);
      // mmu = tail(mmu, mp_new_size);
      // ssig = tail(ssig, mp_new_size);
      rmmp = tail(rmmp, mp_new_size);
      rmmpi = tail((rmmpi - diff), mp_new_size);
      rmmpi[rmmpi < -1] = -1;
    }

    return (List::create(Rcpp::Named("right_matrix_profile") = rmmp, Rcpp::Named("right_profile_index") = rmmpi,
                         Rcpp::Named("w") = window_size, Rcpp::Named("ez") = ez,
                         Rcpp::Named("mp_time_constraint") = mp_time_constraint, Rcpp::Named("offset") = mp_offset,
                         //  Rcpp::Named("ddf") = ddf, Rcpp::Named("ddg") = ddg, Rcpp::Named("avg") = mmu,
                         //  Rcpp::Named("sig") = ssig,
                         //  Rcpp::Named("data") = data_ref,
                         Rcpp::Named("partial") = partial));
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

// MPX stream version for the right side, aiming to compute only the necessary. MP is in pearson values, in review

List mpxi_rcpp(NumericVector new_data, List object, uint64_t history, uint64_t mp_time_constraint, bool progress) {

  try {
    double c, c_cmp, ez;
    uint32_t off_min, off_diag, offset, off_start;
    uint32_t upd_size, n, window_size, exclusion_zone;
    uint64_t data_size, mp_offset;
    NumericVector data, ddf_t, ddg_t;
    bool partial = false;

    upd_size = new_data.length();

    if (object.containsElementNamed("offset"))
      mp_offset = (uint64_t)object["offset"] + upd_size;
    else
      mp_offset = upd_size;

    window_size = (uint32_t)object["w"];

    ez = (double)object["ez"];

    data = as<NumericVector>(object["data"]);

    data_size = data.length();

    exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

    NumericVector data_ref(data_size + upd_size);

    n = data_ref.length();

    // data concatenation
    data_ref[Range(0, data_size - 1)] = data;
    data_ref[Range(data_size, n - 1)] = new_data;

    // matrix profile using cross correlation,

    // copy current avg and sig
    NumericVector mmu(n - window_size + 1);
    mmu[Range(0, data_size - window_size)] = as<NumericVector>(object["avg"]);
    NumericVector ssig(n - window_size + 1);
    ssig[Range(0, data_size - window_size)] = as<NumericVector>(object["sig"]);

    // compute the new avg and sig
    List msd = muinvn_rcpp(data_ref[Range(data_size - window_size + 1, n - 1)], window_size);
    mmu[Range(data_size - window_size + 1, n - window_size)] = as<NumericVector>(msd["avg"]);
    ssig[Range(data_size - window_size + 1, n - window_size)] = as<NumericVector>(msd["sig"]);

    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t profile_len = n - window_size + 1;

    uint32_t diag_start = 0;
    uint32_t diag_end = n - window_size - exclusion_zone;

    if (mp_time_constraint > 0) {
      diag_start = n - mp_time_constraint - window_size;
    }

    IntegerVector compute_order = Range(diag_start, diag_end);

    NumericVector rmmp(profile_len, -1.0);
    IntegerVector rmmpi(profile_len, -1);

    double *rmp = &rmmp[0];
    int *rpi = &rmmpi[0];

    // copy the current profiles
    rmmp[Range(0, data_size - window_size)] = as<NumericVector>(object["right_matrix_profile"]);
    rmmpi[Range(0, data_size - window_size)] = as<IntegerVector>(object["right_profile_index"]);

    ddf_t = as<NumericVector>(object["ddf"]); // shallow copy
    ddg_t = as<NumericVector>(object["ddg"]); // shallow copy

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.
    // ddf is the diff(data_ref, lag = window_size) / 2
    NumericVector ddf(n - window_size + 1, 0);
    ddf[Range(0, data_size - window_size)] = ddf_t;
    ddf[Range(data_size - window_size, n - window_size - 1)] =
        0.5 * (data_ref[Range(data_size - window_size, n - window_size - 1)] - data_ref[Range(data_size, n - 1)]);

    // ddg: (data[(w+1):data_len] - mov_avg[2:(data_len - w + 1)]) + (data[1:(data_len - w)] - mov_avg[1:(data_len -
    // w)]) (subtract the mov_mean of all data, but the first window) + (subtract the mov_mean of all data, but the last
    // window)

    NumericVector ddg(n - window_size + 1, 0);
    ddg[Range(0, data_size - window_size)] = ddg_t;
    ddg[Range(data_size - window_size, n - window_size - 1)] =
        (data_ref[Range(data_size, n - 1)] - mmu[Range(data_size - window_size + 1, profile_len - 1)]) +
        (data_ref[Range(data_size - window_size, n - window_size - 1)] -
         mmu[Range(data_size - window_size, n - window_size - 1)]);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    uint64_t num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress p(100, progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    try {
      uint64_t i = 1;
      // first window demeaned
      NumericVector ww = (data_ref[Range(n - window_size, n - 1)] - mmu[n - window_size]);

      for (int32_t diag : compute_order) {

        if ((i % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          p.increment();
        }
        // this always use the first window (ww); c is the sum of element-wise product
        c = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        // off_max goes from profile_len - ez to 0
        off_min = MAX(n - window_size - upd_size, n - window_size - diag - 1);
        off_start = n - window_size;

        for (offset = off_start; offset > off_min; offset--) {
          // min is offset + diag; max is (profile_len - 1); each iteration has the size of off_max
          off_diag = offset - (n - window_size - diag);

          c = c + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = c * sig[offset] * sig[off_diag];

          // RMP
          if (c_cmp > rmp[off_diag]) {
            rmp[off_diag] = c_cmp;
            rpi[off_diag] = offset + 1;
          }
        }

        if (stop > 0 && i >= stop) {
          partial = true;
          break;
        }
        i++;
      }
    } catch (RcppThread::UserInterruptException &ex) {
      partial = true;
      Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
    }

    // to do ed
    rmmp[rmmp > 1.0] = 1.0;

    // if (euclidean) { // correlation to ed
    //   rmmp = sqrt(2 * window_size * (1 - rmmp));
    //   rmmp[rmmpi < 0] = R_PosInf;
    // }

    if (history > 0 && (n > history)) {
      data_ref = tail(data_ref, history);
      uint64_t mp_new_size = history - window_size + 1;
      uint32_t diff = n - history;
      ddf = tail(ddf, mp_new_size);
      ddg = tail(ddg, mp_new_size);
      mmu = tail(mmu, mp_new_size);
      ssig = tail(ssig, mp_new_size);
      rmmp = tail(rmmp, mp_new_size);
      rmmpi = tail((rmmpi - diff), mp_new_size);
      rmmpi[rmmpi < -1] = -1;
    }

    return (List::create(Rcpp::Named("right_matrix_profile") = rmmp, Rcpp::Named("right_profile_index") = rmmpi,
                         Rcpp::Named("data") = data_ref, Rcpp::Named("w") = window_size, Rcpp::Named("ez") = ez,
                         Rcpp::Named("mp_time_constraint") = mp_time_constraint, Rcpp::Named("ddf") = ddf,
                         Rcpp::Named("ddg") = ddg, Rcpp::Named("avg") = mmu, Rcpp::Named("offset") = mp_offset,
                         Rcpp::Named("sig") = ssig, Rcpp::Named("partial") = partial));
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

// MPX stream version for the right side, but recomputes all again, in review

List mpxiright_rcpp(NumericVector data_ref, uint64_t window_size, double ez, uint64_t mp_time_constraint, double s_size,
                    bool idxs, bool euclidean, bool progress, uint64_t start, List old) {

  uint64_t exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

  try {
    double c, c_cmp;
    uint32_t off_min, off_diag, offset, off_start;
    bool partial = false;
    // matrix profile using cross correlation,
    uint32_t n = data_ref.length();

    List msd = muinvn_rcpp(data_ref, window_size);

    NumericVector mmu = msd["avg"];
    NumericVector ssig = msd["sig"];
    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t profile_len = n - window_size + 1;
    uint32_t diag_start = 0;
    uint32_t diag_end = n - window_size - exclusion_zone;

    if (mp_time_constraint > 0) {
      diag_start = n - mp_time_constraint - window_size;
    }

    IntegerVector compute_order = Range(diag_start, diag_end);

    NumericVector rmmp(profile_len, -1.0);
    IntegerVector rmmpi(profile_len, -1);

    double *rmp = &rmmp[0];
    int *rpi = &rmmpi[0];

    if (start > 0) {

      rmmp[Range(0, profile_len - start - 1)] = as<NumericVector>(old["right_matrix_profile"]);
      rmmpi[Range(0, profile_len - start - 1)] = as<IntegerVector>(old["right_profile_index"]);
    }

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.
    // ddf is the diff(data_ref, lag = window_size) / 2
    NumericVector ddf = 0.5 * (data_ref[Range(0, n - window_size - 1)] - data_ref[Range(window_size, n - 1)]);
    ddf.push_back(0);

    // ddg: (data[(w+1):data_len] - mov_avg[2:(data_len - w + 1)]) + (data[1:(data_len - w)] - mov_avg[1:(data_len -
    // w)]) (subtract the mov_mean of all data, but the first window) + (subtract the mov_mean of all data, but the last
    // window)
    NumericVector ddg = (data_ref[Range(window_size, n - 1)] - mmu[Range(1, profile_len - 1)]) +
                        (data_ref[Range(0, n - window_size - 1)] - mmu[Range(0, n - window_size - 1)]);
    ddg.push_back(0);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    uint64_t num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress p(100, progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    if (s_size < 1.0) {
      stop = round(compute_order.size() * s_size + DBL_EPSILON);
    }

    try {
      uint64_t i = 1;
      // first window demeaned
      NumericVector ww = (data_ref[Range(n - window_size, n - 1)] - mmu[n - window_size]);

      for (int32_t diag : compute_order) {

        if ((i % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          p.increment();
        }
        // this always use the first window (ww); c is the sum of element-wise product
        c = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        // off_max goes from profile_len - ez to 0
        if (start > 0)
          off_min = MAX(n - window_size - start, n - window_size - diag - 1);
        else
          off_min = (profile_len - diag);
        off_start = n - window_size;

        for (offset = off_start; offset > off_min; offset--) {
          // min is offset + diag; max is (profile_len - 1); each iteration has the size of off_max
          off_diag = offset - (n - window_size - diag);

          c = c + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = c * sig[offset] * sig[off_diag];

          // RMP
          if (c_cmp > rmp[off_diag]) {
            rmp[off_diag] = c_cmp;
            if (idxs) {
              rpi[off_diag] = offset + 1;
            }
          }
        }

        if (stop > 0 && i >= stop) {
          partial = true;
          break;
        }
        i++;
      }
    } catch (RcppThread::UserInterruptException &ex) {
      partial = true;
      Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
    }

    // to do ed
    rmmp[rmmp > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      rmmp = sqrt(2 * window_size * (1 - rmmp));
      rmmp[rmmpi < 0] = R_PosInf;
    }

    if (idxs) {
      return (List::create(Rcpp::Named("right_matrix_profile") = rmmp, Rcpp::Named("right_profile_index") = rmmpi,
                           //  Rcpp::Named("data") = data_ref,
                           Rcpp::Named("w") = window_size, Rcpp::Named("ez") = ez,
                           Rcpp::Named("mp_time_constraint") = mp_time_constraint,
                           //  Rcpp::Named("ddf") = ddf,
                           //  Rcpp::Named("ddg") = ddg, Rcpp::Named("avg") = mmu, Rcpp::Named("sig") = ssig,
                           Rcpp::Named("partial") = partial));
    } else {
      return (List::create(Rcpp::Named("right_matrix_profile") = rmmp, Rcpp::Named("partial") = partial));
    }
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

// MPX stream version for the left side, but recomputes all again, in review

List mpxileft_rcpp(NumericVector data_ref, uint64_t window_size, double ez, double s_size, bool idxs, bool euclidean,
                   bool progress, uint64_t start, List old) {

  uint64_t exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

  try {
    double c, c_cmp;
    uint32_t off_max, off_diag, offset, off_start;
    bool partial = false;
    // matrix profile using cross correlation,
    uint32_t n = data_ref.length();

    List msd = muinvn_rcpp(data_ref, window_size);

    NumericVector mmu = msd["avg"];
    NumericVector ssig = msd["sig"];
    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t profile_len = n - window_size + 1;

    IntegerVector compute_order = Range(exclusion_zone, profile_len - 1);

    NumericVector lmmp(profile_len, -1.0);
    IntegerVector lmmpi(profile_len, -1);

    double *lmp = &lmmp[0];
    int *lpi = &lmmpi[0];

    if (start > 0) {

      lmmp[Range(start, profile_len - 1)] = as<NumericVector>(old["left_matrix_profile"]);
      lmmpi[Range(start, profile_len - 1)] = as<IntegerVector>(old["left_profile_index"]);

      // add indexes
      for (uint64_t i = (start + exclusion_zone); i < profile_len; i++) {
        lpi[i] = lpi[i] + start;
      }
    }

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    // ddf is the diff(data_ref, lag = window_size) / 2
    NumericVector ddf = 0.5 * (data_ref[Range(window_size, n - 1)] - data_ref[Range(0, n - window_size - 1)]);
    ddf.push_front(0);

    // ddg: (data[(w+1):data_len] - mov_avg[2:(data_len - w + 1)]) + (data[1:(data_len - w)] - mov_avg[1:(data_len -
    // w)]) (subtract the mov_mean of all data, but the first window) + (subtract the mov_mean of all data, but the last
    // window)
    NumericVector ddg = (data_ref[Range(window_size, n - 1)] - mmu[Range(1, profile_len - 1)]) +
                        (data_ref[Range(0, n - window_size - 1)] - mmu[Range(0, n - window_size - 1)]);
    ddg.push_front(0);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    uint64_t num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress p(100, progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    if (s_size < 1.0) {
      stop = round(compute_order.size() * s_size + DBL_EPSILON);
    }

    try {
      uint64_t i = 1;
      // first window demeaned
      NumericVector ww = (data_ref[Range(0, window_size - 1)] - mmu[0]);

      for (int32_t diag : compute_order) {

        if ((i % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          p.increment();
        }

        // this always use the first window (ww); c is the sum of element-wise product
        c = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        // off_max goes from profile_len - ez to 0
        if (start > 0)
          off_max = MIN(start, (profile_len - diag));
        else
          off_max = (profile_len - diag);
        off_start = 0;

        for (offset = off_start; offset < off_max; offset++) {
          // min is offset + diag; max is (profile_len - 1); each iteration has the size of off_max
          off_diag = offset + diag;

          c = c + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = c * sig[offset] * sig[off_diag];

          // LMP
          if (c_cmp > lmp[off_diag]) {
            lmp[off_diag] = c_cmp;
            if (idxs) {
              lpi[off_diag] = offset + 1;
            }
          }
        }

        if (stop > 0 && i >= stop) {
          partial = true;
          break;
        }
        i++;
      }
    } catch (RcppThread::UserInterruptException &ex) {
      partial = true;
      Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
    }

    // to do ed
    lmmp[lmmp > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      lmmp = sqrt(2 * window_size * (1 - lmmp));
      lmmp[lmmpi < 0] = R_PosInf;
    }

    if (idxs) {
      return (List::create(Rcpp::Named("left_matrix_profile") = lmmp, Rcpp::Named("left_profile_index") = lmmpi,
                           Rcpp::Named("partial") = partial));
    } else {
      return (List::create(Rcpp::Named("left_matrix_profile") = lmmp, Rcpp::Named("partial") = partial));
    }
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

// MPX classic version, in review

// IMPROVE: check data for infinity
// IMPROVE: Zero invalid data, since we are aiming on streaming
// CONCEPT: correction for noise: d_corrected = sqrt(d^2 - (2 + 2m) * std_n^2 / max(std_X, std_Y)^2)
// (10.5220/0007314100830093)

List mpx_rcpp_new(NumericVector data_ref, uint64_t window_size, double ez, uint64_t mp_time_constraint, double s_size,
                  bool idxs, bool euclidean, bool progress) {

  uint64_t exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

  try {
    double c, c_cmp;
    uint32_t off_max, off_diag, offset;
    bool partial = false;
    // matrix profile using cross correlation,
    uint32_t n = data_ref.length();

    List msd = muinvn_rcpp(data_ref, window_size);

    NumericVector mmu = msd["avg"];
    NumericVector ssig = msd["sig"];
    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t profile_len = n - window_size + 1;
    uint32_t diag_start = exclusion_zone;
    uint32_t diag_end = profile_len - 1;

    if (mp_time_constraint > 0) {
      diag_end = mp_time_constraint;
    }

    IntegerVector compute_order = Range(diag_start, diag_end);

    NumericVector mmp(profile_len, -1.0);
    IntegerVector mmpi(profile_len, -1);
    NumericVector rmmp(profile_len, -1.0);
    IntegerVector rmmpi(profile_len, -1);
    NumericVector lmmp(profile_len, -1.0);
    IntegerVector lmmpi(profile_len, -1);

    double *mp = &mmp[0];
    int *mpi = &mmpi[0];
    double *rmp = &rmmp[0];
    int *rpi = &rmmpi[0];
    double *lmp = &lmmp[0];
    int *lpi = &lmmpi[0];

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector ddf = 0.5 * (data_ref[Range(window_size, n - 1)] - data_ref[Range(0, n - window_size - 1)]);
    ddf.push_front(0);
    NumericVector ddg = (data_ref[Range(window_size, n - 1)] - mmu[Range(1, profile_len - 1)]) +
                        (data_ref[Range(0, n - window_size - 1)] - mmu[Range(0, n - window_size - 1)]);
    ddg.push_front(0);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    NumericVector ww = (data_ref[Range(0, window_size - 1)] - mmu[0]);

    uint64_t num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress p(100, progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    if (s_size < 1.0) {
      stop = round(compute_order.size() * s_size + DBL_EPSILON);
    }

    try {
      uint64_t i = 1;
      for (int32_t diag : compute_order) {

        if ((i % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          p.increment();
        }

        c = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        off_max = (n - window_size - diag + 1);

        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          c = c + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];

          if ((sig[offset] > 60) || (sig[off_diag] > 60)) { // wild sig, misleading
            continue;
          }

          c_cmp = c * sig[offset] * sig[off_diag];

          // RMP
          if (c_cmp > rmp[offset]) {
            rmp[offset] = c_cmp;
            if (idxs) {
              rpi[offset] = off_diag + 1;
            }
          }

          // FRANZ TEST START

          // MP
          // if (c_cmp > mp[offset]) {
          //   if (mp[offset] > -1)
          //     mp[offset] = (mp[offset] + c_cmp) / 2;
          //   else
          //     mp[offset] = c_cmp;

          //   if (idxs) {
          //     mpi[offset] = off_diag + 1;
          //   }
          // }
          // if (c_cmp > mp[off_diag]) {
          //   if (mp[off_diag] > -1)
          //     mp[off_diag] = (mp[off_diag] + c_cmp) / 2;
          //   else
          //     mp[off_diag] = c_cmp;

          //   if (idxs) {
          //     mpi[off_diag] = offset + 1;
          //   }
          // }

          // MP
          if (c_cmp > mp[offset]) {
            mp[offset] = c_cmp;
            if (idxs) {
              mpi[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mp[off_diag]) {
            mp[off_diag] = c_cmp;
            if (idxs) {
              mpi[off_diag] = offset + 1;
            }
          }

          // FRANZ TEST END

          // LMP
          if (c_cmp > lmp[off_diag]) {
            lmp[off_diag] = c_cmp;
            if (idxs) {
              lpi[off_diag] = offset + 1;
            }
          }
        }

        if (stop > 0 && i >= stop) {
          partial = true;
          break;
        }

        i++;
      }
    } catch (RcppThread::UserInterruptException &ex) {
      partial = true;
      Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
    }

    // to do ed
    mmp[mmp > 1.0] = 1.0;
    rmmp[rmmp > 1.0] = 1.0;
    lmmp[lmmp > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      mmp = sqrt(2 * window_size * (1 - mmp));
      rmmp = sqrt(2 * window_size * (1 - rmmp));
      lmmp = sqrt(2 * window_size * (1 - lmmp));

      if (idxs) {
        mmp[mmpi < 0] = R_PosInf;
        rmmp[rmmpi < 0] = R_PosInf;
        lmmp[lmmpi < 0] = R_PosInf;
      }
    }

    if (idxs) {
      return (List::create(Rcpp::Named("matrix_profile") = mmp, Rcpp::Named("profile_index") = mmpi,
                           Rcpp::Named("right_matrix_profile") = rmmp, Rcpp::Named("right_profile_index") = rmmpi,
                           Rcpp::Named("left_matrix_profile") = lmmp, Rcpp::Named("left_profile_index") = lmmpi,
                           Rcpp::Named("data") = data_ref, Rcpp::Named("w") = window_size, Rcpp::Named("ez") = ez,
                           Rcpp::Named("mp_time_constraint") = mp_time_constraint, Rcpp::Named("ddf") = ddf,
                           Rcpp::Named("ddg") = ddg, Rcpp::Named("avg") = mmu, Rcpp::Named("sig") = ssig,
                           Rcpp::Named("partial") = partial));
    } else {
      return (List::create(Rcpp::Named("matrix_profile") = mmp, Rcpp::Named("right_matrix_profile") = rmmp,
                           Rcpp::Named("left_matrix_profile") = lmmp, Rcpp::Named("partial") = partial));
    }
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

// FIXME: check skip_locations in mpx

// MPX
//
// @param data_ref Time Series
// @return data_ref List
// [[Rcpp::export]]
List mpx_rcpp(NumericVector data_ref, uint64_t window_size, double ez, double s_size, bool idxs, bool euclidean,
              bool progress) {

  uint64_t exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

  try {
    double c, c_cmp;
    uint32_t off_max, off_diag, offset;
    bool partial = false;
    // matrix profile using cross correlation,
    uint32_t n = data_ref.length();

    List msd = muinvn_rcpp(data_ref, window_size);

    NumericVector mmu = msd["avg"];
    NumericVector ssig = msd["sig"];
    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t profile_len = n - window_size + 1;
    IntegerVector compute_order = Range(exclusion_zone, profile_len - 1);

    NumericVector mmp(profile_len, -1.0);
    IntegerVector mmpi(profile_len, -1);

    double *mp = &mmp[0];
    int *mpi = &mmpi[0];

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector ddf = 0.5 * (data_ref[Range(window_size, n - 1)] - data_ref[Range(0, n - window_size - 1)]);
    ddf.push_front(0);
    NumericVector ddg = (data_ref[Range(window_size, n - 1)] - mmu[Range(1, profile_len - 1)]) +
                        (data_ref[Range(0, n - window_size - 1)] - mmu[Range(0, n - window_size - 1)]);
    ddg.push_front(0);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    NumericVector ww = (data_ref[Range(0, window_size - 1)] - mmu[0]);

    uint64_t num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress p(100, progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    if (s_size < 1.0) {
      stop = round(compute_order.size() * s_size + DBL_EPSILON);
    }

    try {
      uint64_t i = 1;
      for (int32_t diag : compute_order) {

        if ((i % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          p.increment();
        }

        c = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        off_max = (n - window_size - diag + 1);

        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          c = c + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = c * sig[offset] * sig[off_diag];
          if (c_cmp > mp[offset]) {
            mp[offset] = c_cmp;
            if (idxs) {
              mpi[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mp[off_diag]) {
            mp[off_diag] = c_cmp;
            if (idxs) {
              mpi[off_diag] = offset + 1;
            }
          }
        }

        if (stop > 0 && i >= stop) {
          partial = true;
          break;
        }

        i++;
      }
    } catch (RcppThread::UserInterruptException &ex) {
      partial = true;
      Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
    }

    // to do ed
    mmp[mmp > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      mmp = sqrt(2 * window_size * (1 - mmp));
    }

    if (idxs) {
      return (List::create(Rcpp::Named("matrix_profile") = mmp, Rcpp::Named("profile_index") = mmpi,
                           Rcpp::Named("partial") = partial));
    } else {
      return (List::create(Rcpp::Named("matrix_profile") = mmp, Rcpp::Named("partial") = partial));
    }
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

// [[Rcpp::export]]
List mpxab_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double s_size, bool idxs,
                bool euclidean, bool progress) {

  try {
    uint64_t exclusion_zone = 0;
    bool partial = false;
    double c, c_cmp;
    uint32_t off_max, off_diag, offset;
    // matrix profile using cross correlation,
    uint32_t a_len = data_ref.length();
    uint32_t b_len = query_ref.length();

    List msd_a = muinvn_rcpp(data_ref, window_size);
    List msd_b = muinvn_rcpp(query_ref, window_size);

    NumericVector mmu_a = msd_a["avg"];
    NumericVector ssig_a = msd_a["sig"];
    NumericVector mmu_b = msd_b["avg"];
    NumericVector ssig_b = msd_b["sig"];
    double *mu_a = &mmu_a[0];
    double *sig_a = &ssig_a[0];
    double *mu_b = &mmu_b[0];
    double *sig_b = &ssig_b[0];

    uint32_t profile_len_a = a_len - window_size + 1;
    NumericVector mmp_a(profile_len_a, -1.0);
    IntegerVector mmpi_a(profile_len_a, -1);

    double *mp_a = &mmp_a[0];
    int *mpi_a = &mmpi_a[0];

    uint32_t profile_len_b = b_len - window_size + 1;
    NumericVector mmp_b(profile_len_b, -1.0);
    IntegerVector mmpi_b(profile_len_b, -1);

    double *mp_b = &mmp_b[0];
    int *mpi_b = &mmpi_b[0];

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector ddf_a = 0.5 * (data_ref[Range(window_size, a_len - 1)] - data_ref[Range(0, a_len - window_size - 1)]);
    ddf_a.push_front(0);
    NumericVector ddg_a = (data_ref[Range(window_size, a_len - 1)] - mmu_a[Range(1, profile_len_a - 1)]) +
                          (data_ref[Range(0, a_len - window_size - 1)] - mmu_a[Range(0, a_len - window_size - 1)]);
    ddg_a.push_front(0);
    NumericVector ddf_b =
        0.5 * (query_ref[Range(window_size, b_len - 1)] - query_ref[Range(0, b_len - window_size - 1)]);
    ddf_b.push_front(0);
    NumericVector ddg_b = (query_ref[Range(window_size, b_len - 1)] - mmu_b[Range(1, profile_len_b - 1)]) +
                          (query_ref[Range(0, b_len - window_size - 1)] - mmu_b[Range(0, b_len - window_size - 1)]);
    ddg_b.push_front(0);

    double *df_a = &ddf_a[0];
    double *dg_a = &ddg_a[0];
    double *df_b = &ddf_b[0];
    double *dg_b = &ddg_b[0];

    NumericVector ww = (query_ref[Range(0, window_size - 1)] - mmu_b[0]);
    IntegerVector compute_order = Range(exclusion_zone, profile_len_a - 1);
    compute_order = sample(compute_order, compute_order.size());

    uint64_t num_progress = ceil((double)(profile_len_a + profile_len_b - 2 * exclusion_zone) /
                                 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress p(100, progress);

    try {
      //// AB ----
      for (int32_t diag : compute_order) {

        if ((diag % num_progress) == 10) {
          p.increment();
          RcppThread::checkUserInterrupt();
        }

        off_max = MIN(a_len - window_size - diag + 1, b_len - window_size + 1);

        // 1. Towsley, A. et al.: Correlation Angles and Inner Products: Application to a Problem from Physics. ISRN
        // Applied Mathematics. 2011, 1-12 (2011). https://doi.org/10.5402/2011/323864.
        // c is the Covar(X,Y); Additionally it follows from the Cauchy-Schwartz inequality[1] that
        // |Covar(X,Y)| <= SD(X) SD(Y). Here the Covar(X,Y) is obtained using the inner_product(X,Y)
        // [1] K. Hoffman and R. Kunze, Linear Algebra, Prentice Hall, Englewood Cliffs, NJ, USA, 2nd edition, 1971.

        c = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu_a[diag]), ww);
        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          c = c + df_a[off_diag] * dg_b[offset] + dg_a[off_diag] * df_b[offset];

          if ((sig_b[offset] > 60) || (sig_a[off_diag] > 60)) { // wild sig, misleading
            continue;
          }

          // or c / (sd_b * sd_a), but fractions are expensive
          // also, c_cmp is cos(Γ(X,Y))
          c_cmp = c * sig_b[offset] * sig_a[off_diag];

          if (c_cmp > mp_b[offset]) { // mpb
            mp_b[offset] = c_cmp;
            if (idxs) {
              mpi_b[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mp_a[off_diag]) { // mpa
            mp_a[off_diag] = c_cmp;
            if (idxs) {
              mpi_a[off_diag] = offset + 1;
            }
          }
        }
      }

      //// BA ----
      ww = (data_ref[Range(0, window_size - 1)] - mmu_a[0]);
      compute_order = Range(exclusion_zone, profile_len_b - 1);
      compute_order = sample(compute_order, compute_order.size());

      for (int32_t diag : compute_order) {

        if ((diag % num_progress) == 10) {
          p.increment();
          RcppThread::checkUserInterrupt();
        }

        off_max = MIN(b_len - window_size - diag + 1, a_len - window_size + 1);
        c = inner_product((query_ref[Range(diag, (diag + window_size - 1))] - mu_b[diag]), ww);
        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          c = c + df_b[off_diag] * dg_a[offset] + dg_b[off_diag] * df_a[offset];

          if ((sig_a[offset] > 60) || (sig_b[off_diag] > 60)) { // wild sig, misleading
            continue;
          }

          c_cmp = c * sig_a[offset] * sig_b[off_diag];

          if (c_cmp > mp_a[offset]) {
            mp_a[offset] = c_cmp;
            if (idxs) {
              mpi_a[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mp_b[off_diag]) {
            mp_b[off_diag] = c_cmp;
            if (idxs) {
              mpi_b[off_diag] = offset + 1;
            }
          }
        }
      }
    } catch (RcppThread::UserInterruptException &ex) {
      partial = true;
      Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
    }

    // to do ed
    mmp_a[mmp_a > 1.0] = 1.0;
    mmp_b[mmp_b > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      mmp_a = sqrt(2 * window_size * (1 - mmp_a));
      mmp_b = sqrt(2 * window_size * (1 - mmp_b));
    }

    if (idxs) {
      return (List::create(Rcpp::Named("matrix_profile") = mmp_a, Rcpp::Named("profile_index") = mmpi_a,
                           Rcpp::Named("mpb") = mmp_b, Rcpp::Named("pib") = mmpi_b, Rcpp::Named("partial") = partial));
    } else {
      return (List::create(Rcpp::Named("matrix_profile") = mmp_a, Rcpp::Named("mpb") = mmp_b,
                           Rcpp::Named("partial") = partial));
    }
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

//##### Parallel version #####

struct MatrixProfileP : public Worker {
  // input
  const RVector<double> data_ref;
  const uint64_t window_size;
  const RVector<int> compute_order;
  const RVector<double> df;
  const RVector<double> dg;
  const RVector<double> mu;
  const RVector<double> sig;
  const RVector<double> ww;

  Progress *p;
  uint64_t num_progress;
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
  MatrixProfileP(const NumericVector data_ref, const uint64_t window_size, const IntegerVector compute_order,
                 const NumericVector df, const NumericVector dg, const NumericVector mmu, const NumericVector sig,
                 const NumericVector ww, Progress *p, uint64_t num_progress, NumericVector mp, IntegerVector mpi)
      : data_ref(data_ref), window_size(window_size), compute_order(compute_order), df(df), dg(dg), mu(mmu), sig(sig),
        ww(ww), p(p), num_progress(num_progress), mp(mp), mpi(mpi) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) { // exclusion_zone:profile_len
    double c, c_cmp;
    uint32_t off_max, off_diag, offset;
    uint32_t n = data_ref.length();
    std::vector<double> aa(window_size);

    std::vector<double> mpp(mp.size(), -1.0);
    std::vector<int> mpip(mp.size(), -1);

    try {
      for (uint64_t dd = begin; dd < end; dd++) {
        int64_t diag = compute_order[dd];

        if ((diag % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          m.lock();
          p->increment();
          m.unlock();
        }

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

      m.lock();
      for (uint64_t i = 0; i < mp.size(); i++) {
        if (mpp[i] > mp[i]) {
          mp[i] = mpp[i];
          mpi[i] = mpip[i];
        }
      }
      m.unlock();
    } catch (RcppThread::UserInterruptException &e) {
      Rcout << "Computation interrupted by the user." << std::endl;
      Rcout << "Please wait for other threads to stop." << std::endl;
      throw;
    }
  }
};

// [[Rcpp::export]]
List mpx_rcpp_parallel(NumericVector data_ref, uint64_t window_size, double ez, double s_size, bool idxs,
                       bool euclidean, bool progress) {

  uint64_t exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

  try {
    // matrix profile using cross correlation,
    bool partial = false;
    uint32_t n = data_ref.length();

    List msd = muinvn_rcpp_parallel(data_ref, window_size);

    NumericVector mu = msd["avg"];
    NumericVector sig = msd["sig"];

    uint32_t profile_len = n - window_size + 1;
    NumericVector mp(profile_len, -1.0);
    IntegerVector mpi(profile_len, -1);

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector df = 0.5 * (data_ref[Range(window_size, n - 1)] - data_ref[Range(0, n - window_size - 1)]);
    df.push_front(0);
    NumericVector dg = (data_ref[Range(window_size, n - 1)] - mu[Range(1, profile_len - 1)]) +
                       (data_ref[Range(0, n - window_size - 1)] - mu[Range(0, n - window_size - 1)]);
    dg.push_front(0);

    NumericVector ww = (data_ref[Range(0, window_size - 1)] - mu[0]);

    IntegerVector compute_order = Range(exclusion_zone, profile_len);
    compute_order = sample(compute_order, compute_order.size());

    uint64_t num_progress = (profile_len - exclusion_zone) / 100;

    Progress p(100, progress);

    MatrixProfileP matrix_profile(data_ref, window_size, compute_order, df, dg, mu, sig, ww, &p, num_progress, mp, mpi);

    try {
#if RCPP_PARALLEL_USE_TBB
      RcppParallel::parallelFor(0, profile_len - exclusion_zone, matrix_profile, 4 * window_size);
#else
      RcppParallel2::ttParallelFor(0, profile_len - exclusion_zone, matrix_profile, 4 * window_size);
#endif
    } catch (RcppThread::UserInterruptException &e) {
      partial = true;
      Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
    } catch (...) {
      ::Rf_error("c++ exception (unknown reason)");
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
    ::Rf_error("c++ exception (unknown reason)");
  }
}

struct MatrixProfilePAB : public Worker {
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
  Progress *p;
  uint64_t num_progress;

  // output
  RVector<double> mp_a;
  RVector<double> mp_b;
  RVector<int> mpi_a;
  RVector<int> mpi_b;

  // AB == 0, BA == 1
  uint8_t ab_ba;

#if RCPP_PARALLEL_USE_TBB
  tbb::spin_mutex m;
  tbb::mutex m1;
#else
  tthread::mutex m;
  tthread::mutex m1;
#endif

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  MatrixProfilePAB(const NumericVector data_ref, const NumericVector query_ref, const uint64_t window_size,
                   const NumericVector df_a, const NumericVector df_b, const NumericVector dg_a,
                   const NumericVector dg_b, const NumericVector mu_a, const NumericVector mu_b,
                   const NumericVector sig_a, const NumericVector sig_b, const NumericVector ww_a,
                   const NumericVector ww_b, IntegerVector compute_ordera, IntegerVector compute_orderb, Progress *p,
                   uint64_t num_progress, NumericVector mp_a, NumericVector mp_b, IntegerVector mpi_a,
                   IntegerVector mpi_b)
      : data_ref(data_ref), query_ref(query_ref), window_size(window_size), df_a(df_a), df_b(df_b), dg_a(dg_a),
        dg_b(dg_b), mu_a(mu_a), mu_b(mu_b), sig_a(sig_a), sig_b(sig_b), ww_a(ww_a), ww_b(ww_b),
        compute_ordera(compute_ordera), compute_orderb(compute_orderb), p(p), num_progress(num_progress), mp_a(mp_a),
        mp_b(mp_b), mpi_a(mpi_a), mpi_b(mpi_b), ab_ba(0) {}

  void set_ab() { this->ab_ba = 0; }
  void set_ba() { this->ab_ba = 1; }

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) { // exclusion_zone:profile_len
    double c, c_cmp;
    uint32_t off_max, off_diag, offset;
    uint32_t a_len = data_ref.length();
    uint32_t b_len = query_ref.length();
    std::vector<double> inn(window_size);
    std::vector<double> mpp_a(mp_a.size(), -1.0);
    std::vector<int> mpip_a(mp_a.size(), -1);
    std::vector<double> mpp_b(mp_b.size(), -1.0);
    std::vector<int> mpip_b(mp_b.size(), -1);

    try {
      if (ab_ba == 0) {

        for (uint32_t diag = begin; diag < end; diag++) {

          if ((diag % num_progress) == 0) {
            RcppThread::checkUserInterrupt();
            m.lock();
            p->increment();
            m.unlock();
          }

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
        for (uint32_t diag = begin; diag < end; diag++) {

          if ((diag % num_progress) == 0) {
            RcppThread::checkUserInterrupt();
            m.lock();
            p->increment();
            m.unlock();
          }

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

      m1.lock();
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
      m1.unlock();
    } catch (RcppThread::UserInterruptException &e) {
      Rcout << "Computation interrupted by the user." << std::endl;
      Rcout << "Please wait for other threads to stop." << std::endl;
      throw;
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

    uint32_t profile_len_a = a_len - window_size + 1;
    uint32_t profile_len_b = b_len - window_size + 1;

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

    uint64_t num_progress =
        ceil((double)(profile_len_a + profile_len_b) / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress p(100, progress);

    IntegerVector compute_ordera = Range(0, profile_len_a);
    compute_ordera = sample(compute_ordera, compute_ordera.size());
    IntegerVector compute_orderb = Range(0, profile_len_b);
    compute_orderb = sample(compute_orderb, compute_orderb.size());

    MatrixProfilePAB matrix_profile(data_ref, query_ref, window_size, df_a, df_b, dg_a, dg_b, mu_a, mu_b, sig_a, sig_b,
                                    ww_a, ww_b, compute_ordera, compute_orderb, &p, num_progress, mp_a, mp_b, mpi_a,
                                    mpi_b);

    try {
#if RCPP_PARALLEL_USE_TBB
      RcppParallel::parallelFor(0, profile_len_a, matrix_profile, 4 * window_size);
#else
      RcppParallel2::ttParallelFor(0, profile_len_a, matrix_profile, 4 * window_size);
#endif
    } catch (RcppThread::UserInterruptException &e) {
      partial = true;
      Rcout << "Process AB terminated by the user successfully, partial results were returned." << std::endl;
    } catch (...) {
      ::Rf_error("c++ exception (unknown reason)");
    }

    // switch worker to BA
    matrix_profile.set_ba();

    try {
#if RCPP_PARALLEL_USE_TBB
      RcppParallel::parallelFor(1, profile_len_b, matrix_profile, 4 * window_size);
#else
      RcppParallel2::ttParallelFor(1, profile_len_b, matrix_profile, 4 * window_size);
#endif
    } catch (RcppThread::UserInterruptException &e) {
      partial = true;
      Rcout << "Process BA terminated by the user successfully, partial results were returned." << std::endl;
    } catch (...) {
      ::Rf_error("c++ exception (unknown reason)");
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
    ::Rf_error("c++ exception (unknown reason)");
  }
}
