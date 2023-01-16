// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#include "mathtools.h" // math first to fix OSX error
#include "mpx.h"
#include "windowfunc.h"
#include <cfloat> // DBL_EPSILON when STRICT_R_HEADERS

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

// MPX stream version for the right side, getting stats from input, in review

List mpxis_rcpp(NumericVector data_ref, uint64_t batch_size, List object, List stats, uint64_t history,
                uint64_t mp_time_constraint, bool progress, float threshold) {

  try {
    double cc = 0.0, c_cmp = 0.0, ez = 0.0;
    uint32_t off_min = 0UL, off_diag = 0UL, offset = 0UL, off_start = 0UL;
    uint32_t data_size = 0UL, window_size = 0UL, exclusion_zone = 0UL, profile_len = 0UL;
    uint64_t mp_offset = 0UL;
    NumericVector const data, ddf_t, ddg_t;
    bool initial = false;
    bool partial = false;

    mp_offset = (uint64_t)object["offset"] + batch_size;
    window_size = (uint32_t)object["w"];
    ez = (double)object["ez"];
    exclusion_zone = static_cast<uint32_t>(round(window_size * ez + DBL_EPSILON) + 1);
    data_size = data_ref.length();
    profile_len = data_size - window_size + 1;

    if (data_size == batch_size) {
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
    uint32_t const diag_end = data_size - window_size - exclusion_zone;

    if (mp_time_constraint > 0) {
      diag_start = data_size - mp_time_constraint - window_size;
    }

    IntegerVector const compute_order = Range(diag_start, diag_end);

    NumericVector rmmp(profile_len, -1.0);
    IntegerVector rmmpi(profile_len, -1);

    if (!initial) {
      // copy the current profiles
      rmmp[Range(0, data_size - window_size - batch_size)] = as<NumericVector>(object["right_matrix_profile"]);
      rmmpi[Range(0, data_size - window_size - batch_size)] = as<IntegerVector>(object["right_profile_index"]);
    }

    double *rmp = &rmmp[0];
    int *rpi = &rmmpi[0];

    // added double inside sqrt to avoid ambiguity on Solaris
    uint64_t const num_progress = ceil(static_cast<double>(compute_order.size()) / 100.0);

    Progress prog(100, progress);

    // compute_order = sample(compute_order, compute_order.size());

    uint64_t const stop = 0;

    try {
      uint64_t it = 1;
      // first window demeaned
      NumericVector const ww = (data_ref[Range(data_size - window_size, data_size - 1)] - mmu[data_size - window_size]);

      for (int32_t const &diag : compute_order) {

        if ((it % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          prog.increment();
        }
        // this always use the first window (ww); cc is the sum of element-wise product
        cc = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        // off_max goes from profile_len - ez to 0
        if (initial) {
          off_min = data_size - window_size - diag - 1;
        } else {
          off_min = MAX(data_size - window_size - batch_size, data_size - window_size - diag - 1);
        }
        off_start = data_size - window_size;

        // loop if batch_size > 1
        for (offset = off_start; offset > off_min; offset--) {
          // min is offset + diag; max is (profile_len - 1); each iteration has the size of off_max
          off_diag = offset - (data_size - window_size - diag);

          cc = cc + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = cc * sig[offset] * sig[off_diag];

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

        if (stop > 0 && it >= stop) {
          partial = true;
          break;
        }
        it++;
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

    if (history > 0 && (data_size > history)) {
      // data_ref = tail(data_ref, history);
      uint64_t const mp_new_size = history - window_size + 1;
      uint32_t const diff = data_size - history;
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
    Rcpp::stop("c++ exception (unknown reason)");
  }
}

// MPX stream version for the right side, aiming to compute only the necessary. MP is in pearson values, in review

List mpxi_rcpp(NumericVector new_data, List object, uint64_t history, uint64_t mp_time_constraint, bool progress) {

  try {
    double cc = 0.0, c_cmp = 0.0, ez = 0.0;
    uint32_t off_min = 0UL, off_diag = 0UL, offset = 0UL, off_start = 0UL;
    uint32_t upd_size = 0UL, data_ref_size = 0UL, window_size = 0UL, exclusion_zone = 0UL;
    uint64_t data_size = 0, mp_offset = 0;
    NumericVector data, ddf_t, ddg_t;
    bool partial = false;

    upd_size = new_data.length();

    if (object.containsElementNamed("offset")) {
      mp_offset = (uint64_t)object["offset"] + upd_size;
    } else {
      mp_offset = upd_size;
    }

    window_size = (uint32_t)object["w"];

    ez = (double)object["ez"];

    data = as<NumericVector>(object["data"]);

    data_size = data.length();

    exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

    NumericVector data_ref(data_size + upd_size);

    data_ref_size = data_ref.length();

    // data concatenation
    data_ref[Range(0, data_size - 1)] = data;
    data_ref[Range(data_size, data_ref_size - 1)] = new_data;

    // matrix profile using cross correlation,

    // copy current avg and sig
    NumericVector mmu(data_ref_size - window_size + 1);
    mmu[Range(0, data_size - window_size)] = as<NumericVector>(object["avg"]);
    NumericVector ssig(data_ref_size - window_size + 1);
    ssig[Range(0, data_size - window_size)] = as<NumericVector>(object["sig"]);

    // compute the new avg and sig
    List msd = muinvn_rcpp(data_ref[Range(data_size - window_size + 1, data_ref_size - 1)], window_size);
    mmu[Range(data_size - window_size + 1, data_ref_size - window_size)] = as<NumericVector>(msd["avg"]);
    ssig[Range(data_size - window_size + 1, data_ref_size - window_size)] = as<NumericVector>(msd["sig"]);

    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t const profile_len = data_ref_size - window_size + 1;

    uint32_t diag_start = 0;
    uint32_t const diag_end = data_size - window_size - exclusion_zone;

    if (mp_time_constraint > 0) {
      diag_start = data_ref_size - mp_time_constraint - window_size;
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
    NumericVector ddf(data_ref_size - window_size + 1, 0);
    ddf[Range(0, data_size - window_size)] = ddf_t;
    ddf[Range(data_size - window_size, data_ref_size - window_size - 1)] =
        0.5 * (data_ref[Range(data_size - window_size, data_ref_size - window_size - 1)] -
               data_ref[Range(data_size, data_ref_size - 1)]);

    // ddg: (data[(w+1):data_len] - mov_avg[2:(data_len - w + 1)]) + (data[1:(data_len - w)] - mov_avg[1:(data_len -
    // w)]) (subtract the mov_mean of all data, but the first window) + (subtract the mov_mean of all data, but the last
    // window)

    NumericVector ddg(data_ref_size - window_size + 1, 0);
    ddg[Range(0, data_size - window_size)] = ddg_t;
    ddg[Range(data_size - window_size, data_ref_size - window_size - 1)] =
        (data_ref[Range(data_size, data_ref_size - 1)] - mmu[Range(data_size - window_size + 1, profile_len - 1)]) +
        (data_ref[Range(data_size - window_size, data_ref_size - window_size - 1)] -
         mmu[Range(data_size - window_size, data_ref_size - window_size - 1)]);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    uint64_t const num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress prog(100, progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t const stop = 0;

    try {
      uint64_t iteration = 1;
      // first window demeaned
      NumericVector const ww =
          (data_ref[Range(data_ref_size - window_size, data_ref_size - 1)] - mmu[data_ref_size - window_size]);

      for (int32_t const diag : compute_order) {

        if ((iteration % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          prog.increment();
        }
        // this always use the first window (ww); cc is the sum of element-wise product
        cc = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        // off_max goes from profile_len - ez to 0
        off_min = MAX(data_ref_size - window_size - upd_size, data_ref_size - window_size - diag - 1);
        off_start = data_ref_size - window_size;

        for (offset = off_start; offset > off_min; offset--) {
          // min is offset + diag; max is (profile_len - 1); each iteration has the size of off_max
          off_diag = offset - (data_ref_size - window_size - diag);

          cc = cc + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = cc * sig[offset] * sig[off_diag];

          // RMP
          if (c_cmp > rmp[off_diag]) {
            rmp[off_diag] = c_cmp;
            rpi[off_diag] = offset + 1;
          }
        }

        if (stop > 0 && iteration >= stop) {
          partial = true;
          break;
        }
        iteration++;
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

    if (history > 0 && (data_ref_size > history)) {
      data_ref = tail(data_ref, history);
      uint64_t const mp_new_size = history - window_size + 1;
      uint32_t const diff = data_ref_size - history;
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
    Rcpp::stop("c++ exception (unknown reason)");
  }
}

// MPX stream version for the right side, but recomputes all again, in review

List mpxiright_rcpp(NumericVector data_ref, uint64_t window_size, double ez, uint64_t mp_time_constraint, double s_size,
                    bool idxs, bool euclidean, bool progress, uint64_t start, List old) {

  uint64_t const exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

  try {
    double cc = 0.0, c_cmp = 0.0;
    uint32_t off_min = 0UL, off_diag = 0UL, offset = 0UL, off_start = 0UL;
    bool partial = false;
    // matrix profile using cross correlation,
    uint32_t const data_size = data_ref.length();

    List msd = muinvn_rcpp(data_ref, window_size);

    NumericVector mmu = msd["avg"];
    NumericVector ssig = msd["sig"];
    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t const profile_len = data_size - window_size + 1;
    uint32_t diag_start = 0;
    uint32_t const diag_end = data_size - window_size - exclusion_zone;

    if (mp_time_constraint > 0) {
      diag_start = data_size - mp_time_constraint - window_size;
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
    NumericVector ddf =
        0.5 * (data_ref[Range(0, data_size - window_size - 1)] - data_ref[Range(window_size, data_size - 1)]);
    ddf.push_back(0);

    // ddg: (data[(w+1):data_len] - mov_avg[2:(data_len - w + 1)]) + (data[1:(data_len - w)] - mov_avg[1:(data_len -
    // w)]) (subtract the mov_mean of all data, but the first window) + (subtract the mov_mean of all data, but the last
    // window)
    NumericVector ddg = (data_ref[Range(window_size, data_size - 1)] - mmu[Range(1, profile_len - 1)]) +
                        (data_ref[Range(0, data_size - window_size - 1)] - mmu[Range(0, data_size - window_size - 1)]);
    ddg.push_back(0);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    uint64_t const num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress prog(100, progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    if (s_size < 1.0) {
      stop = round(compute_order.size() * s_size + DBL_EPSILON);
    }

    try {
      uint64_t it = 1;
      // first window demeaned
      NumericVector const ww = (data_ref[Range(data_size - window_size, data_size - 1)] - mmu[data_size - window_size]);

      for (int32_t const diag : compute_order) {

        if ((it % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          prog.increment();
        }
        // this always use the first window (ww); cc is the sum of element-wise product
        cc = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        // off_max goes from profile_len - ez to 0
        if (start > 0) {
          off_min = MAX(data_size - window_size - start, data_size - window_size - diag - 1);
        } else {
          off_min = (profile_len - diag);
        }
        off_start = data_size - window_size;

        for (offset = off_start; offset > off_min; offset--) {
          // min is offset + diag; max is (profile_len - 1); each iteration has the size of off_max
          off_diag = offset - (data_size - window_size - diag);

          cc = cc + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = cc * sig[offset] * sig[off_diag];

          // RMP
          if (c_cmp > rmp[off_diag]) {
            rmp[off_diag] = c_cmp;
            if (idxs) {
              rpi[off_diag] = offset + 1;
            }
          }
        }

        if (stop > 0 && it >= stop) {
          partial = true;
          break;
        }
        it++;
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
    Rcpp::stop("c++ exception (unknown reason)");
  }
}

// MPX stream version for the left side, but recomputes all again, in review

List mpxileft_rcpp(NumericVector data_ref, uint64_t window_size, double ez, double s_size, bool idxs, bool euclidean,
                   bool progress, uint64_t start, List old) {

  uint64_t const exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

  try {
    double cc = 0.0, c_cmp = 0.0;
    uint32_t off_max = 0UL, off_diag = 0UL, offset = 0UL, off_start = 0UL;
    bool partial = false;
    // matrix profile using cross correlation,
    uint32_t const data_size = data_ref.length();

    List msd = muinvn_rcpp(data_ref, window_size);

    NumericVector mmu = msd["avg"];
    NumericVector ssig = msd["sig"];
    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t const profile_len = data_size - window_size + 1;

    IntegerVector compute_order = Range(exclusion_zone, profile_len - 1);

    NumericVector lmmp(profile_len, -1.0);
    IntegerVector lmmpi(profile_len, -1);

    double *lmp = &lmmp[0];
    int *lpi = &lmmpi[0];

    if (start > 0) {

      lmmp[Range(start, profile_len - 1)] = as<NumericVector>(old["left_matrix_profile"]);
      lmmpi[Range(start, profile_len - 1)] = as<IntegerVector>(old["left_profile_index"]);

      // add indexes
      for (uint64_t it = (start + exclusion_zone); it < profile_len; it++) {
        lpi[it] = lpi[it] + start;
      }
    }

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    // ddf is the diff(data_ref, lag = window_size) / 2
    NumericVector ddf =
        0.5 * (data_ref[Range(window_size, data_size - 1)] - data_ref[Range(0, data_size - window_size - 1)]);
    ddf.push_front(0);

    // ddg: (data[(w+1):data_len] - mov_avg[2:(data_len - w + 1)]) + (data[1:(data_len - w)] - mov_avg[1:(data_len -
    // w)]) (subtract the mov_mean of all data, but the first window) + (subtract the mov_mean of all data, but the last
    // window)
    NumericVector ddg = (data_ref[Range(window_size, data_size - 1)] - mmu[Range(1, profile_len - 1)]) +
                        (data_ref[Range(0, data_size - window_size - 1)] - mmu[Range(0, data_size - window_size - 1)]);
    ddg.push_front(0);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    uint64_t const num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress prog(100, progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    if (s_size < 1.0) {
      stop = round(compute_order.size() * s_size + DBL_EPSILON);
    }

    try {
      uint64_t it = 1;
      // first window demeaned
      NumericVector const ww = (data_ref[Range(0, window_size - 1)] - mmu[0]);

      for (int32_t const diag : compute_order) {

        if ((it % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          prog.increment();
        }

        // this always use the first window (ww); cc is the sum of element-wise product
        cc = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        // off_max goes from profile_len - ez to 0
        if (start > 0) {
          off_max = MIN(start, (profile_len - diag));
        } else {
          off_max = (profile_len - diag);
        }
        off_start = 0;

        for (offset = off_start; offset < off_max; offset++) {
          // min is offset + diag; max is (profile_len - 1); each iteration has the size of off_max
          off_diag = offset + diag;

          cc = cc + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = cc * sig[offset] * sig[off_diag];

          // LMP
          if (c_cmp > lmp[off_diag]) {
            lmp[off_diag] = c_cmp;
            if (idxs) {
              lpi[off_diag] = offset + 1;
            }
          }
        }

        if (stop > 0 && it >= stop) {
          partial = true;
          break;
        }
        it++;
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
    Rcpp::stop("c++ exception (unknown reason)");
  }
}

// MPX classic version, in review

// IMPROVE: check data for infinity
// IMPROVE: Zero invalid data, since we are aiming on streaming
// CONCEPT: correction for noise: d_corrected = sqrt(d^2 - (2 + 2m) * std_n^2 / max(std_X, std_Y)^2)
// (10.5220/0007314100830093)

List mpx_rcpp_new(NumericVector data_ref, uint64_t window_size, double ez, uint64_t mp_time_constraint, double s_size,
                  bool idxs, bool euclidean, bool progress) {

  uint64_t const exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

  try {
    double cc = 0.0, c_cmp = 0.0;
    uint32_t off_max = 0UL, off_diag = 0UL, offset = 0UL;
    bool partial = false;
    // matrix profile using cross correlation,
    uint32_t const data_size = data_ref.length();

    List msd = muinvn_rcpp(data_ref, window_size);

    NumericVector mmu = msd["avg"];
    NumericVector ssig = msd["sig"];
    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t const profile_len = data_size - window_size + 1;
    uint32_t const diag_start = exclusion_zone;
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

    NumericVector ddf =
        0.5 * (data_ref[Range(window_size, data_size - 1)] - data_ref[Range(0, data_size - window_size - 1)]);
    ddf.push_front(0);
    NumericVector ddg = (data_ref[Range(window_size, data_size - 1)] - mmu[Range(1, profile_len - 1)]) +
                        (data_ref[Range(0, data_size - window_size - 1)] - mmu[Range(0, data_size - window_size - 1)]);
    ddg.push_front(0);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    NumericVector const ww = (data_ref[Range(0, window_size - 1)] - mmu[0]);

    uint64_t const num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress prog(100, progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    if (s_size < 1.0) {
      stop = round(compute_order.size() * s_size + DBL_EPSILON);
    }

    try {
      uint64_t it = 1;
      for (int32_t const diag : compute_order) {

        if ((it % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          prog.increment();
        }

        cc = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        off_max = (data_size - window_size - diag + 1);

        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          cc = cc + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];

          if ((sig[offset] > 60) || (sig[off_diag] > 60)) { // wild sig, misleading
            continue;
          }

          c_cmp = cc * sig[offset] * sig[off_diag];

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

        if (stop > 0 && it >= stop) {
          partial = true;
          break;
        }

        it++;
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
    Rcpp::stop("c++ exception (unknown reason)");
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

  uint64_t const exclusion_zone = round(window_size * ez + DBL_EPSILON) + 1;

  try {
    double cc = 0.0, c_cmp = 0.0;
    uint32_t off_max = 0UL, off_diag = 0UL, offset = 0UL;
    bool partial = false;
    // matrix profile using cross correlation,
    uint32_t const data_size = data_ref.length();

    List msd = muinvn_rcpp(data_ref, window_size);

    NumericVector mmu = msd["avg"];
    NumericVector ssig = msd["sig"];
    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t const profile_len = data_size - window_size + 1;
    IntegerVector compute_order = Range(exclusion_zone, profile_len - 1);

    NumericVector mmp(profile_len, -1.0);
    IntegerVector mmpi(profile_len, -1);

    double *mp = &mmp[0];
    int *mpi = &mmpi[0];

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector ddf =
        0.5 * (data_ref[Range(window_size, data_size - 1)] - data_ref[Range(0, data_size - window_size - 1)]);
    ddf.push_front(0);
    NumericVector ddg = (data_ref[Range(window_size, data_size - 1)] - mmu[Range(1, profile_len - 1)]) +
                        (data_ref[Range(0, data_size - window_size - 1)] - mmu[Range(0, data_size - window_size - 1)]);
    ddg.push_front(0);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    NumericVector const ww = (data_ref[Range(0, window_size - 1)] - mmu[0]);

    uint64_t const num_progress =
        ceil((double)compute_order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress prog(100, progress);

    compute_order = sample(compute_order, compute_order.size());

    uint64_t stop = 0;

    if (s_size < 1.0) {
      stop = round(compute_order.size() * s_size + DBL_EPSILON);
    }

    try {
      uint64_t it = 1;
      for (int32_t const diag : compute_order) {

        if ((it % num_progress) == 0) {
          RcppThread::checkUserInterrupt();
          prog.increment();
        }

        cc = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu[diag]), ww);

        off_max = (data_size - window_size - diag + 1);

        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          cc = cc + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = cc * sig[offset] * sig[off_diag];
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

        if (stop > 0 && it >= stop) {
          partial = true;
          break;
        }

        it++;
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
    Rcpp::stop("c++ exception (unknown reason)");
  }
}

// [[Rcpp::export]]
List mpxab_rcpp(NumericVector data_ref, NumericVector query_ref, uint64_t window_size, double s_size, bool idxs,
                bool euclidean, bool progress) {

  try {
    uint64_t const exclusion_zone = 0;
    bool partial = false;
    double cc = 0.0, c_cmp = 0.0;
    uint32_t off_max = 0UL, off_diag = 0UL, offset = 0UL;
    // matrix profile using cross correlation,
    uint32_t const a_len = data_ref.length();
    uint32_t const b_len = query_ref.length();

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

    uint32_t const profile_len_a = a_len - window_size + 1;
    NumericVector mmp_a(profile_len_a, -1.0);
    IntegerVector mmpi_a(profile_len_a, -1);

    double *mp_a = &mmp_a[0];
    int *mpi_a = &mmpi_a[0];

    uint32_t const profile_len_b = b_len - window_size + 1;
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

    uint64_t const num_progress = ceil((double)(profile_len_a + profile_len_b - 2 * exclusion_zone) /
                                       100); // added double inside sqrt to avoid ambiguity on Solaris

    Progress prog(100, progress);

    try {
      //// AB ----
      for (int32_t const diag : compute_order) {

        if ((diag % num_progress) == 10) {
          prog.increment();
          RcppThread::checkUserInterrupt();
        }

        off_max = MIN(a_len - window_size - diag + 1, b_len - window_size + 1);

        // 1. Towsley, A. et al.: Correlation Angles and Inner Products: Application to a Problem from Physics. ISRN
        // Applied Mathematics. 2011, 1-12 (2011). https://doi.org/10.5402/2011/323864.
        // cc is the Covar(X,Y); Additionally it follows from the Cauchy-Schwartz inequality[1] that
        // |Covar(X,Y)| <= SD(X) SD(Y). Here the Covar(X,Y) is obtained using the inner_product(X,Y)
        // [1] K. Hoffman and R. Kunze, Linear Algebra, Prentice Hall, Englewood Cliffs, NJ, USA, 2nd edition, 1971.

        cc = inner_product((data_ref[Range(diag, (diag + window_size - 1))] - mu_a[diag]), ww);
        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          cc = cc + df_a[off_diag] * dg_b[offset] + dg_a[off_diag] * df_b[offset];

          if ((sig_b[offset] > 60) || (sig_a[off_diag] > 60)) { // wild sig, misleading
            continue;
          }

          // or cc / (sd_b * sd_a), but fractions are expensive
          // also, c_cmp is cos(Γ(X,Y))
          c_cmp = cc * sig_b[offset] * sig_a[off_diag];

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

      for (int32_t const diag : compute_order) {

        if ((diag % num_progress) == 10) {
          prog.increment();
          RcppThread::checkUserInterrupt();
        }

        off_max = MIN(b_len - window_size - diag + 1, a_len - window_size + 1);
        cc = inner_product((query_ref[Range(diag, (diag + window_size - 1))] - mu_b[diag]), ww);
        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          cc = cc + df_b[off_diag] * dg_a[offset] + dg_b[off_diag] * df_a[offset];

          if ((sig_a[offset] > 60) || (sig_b[off_diag] > 60)) { // wild sig, misleading
            continue;
          }

          c_cmp = cc * sig_a[offset] * sig_b[off_diag];

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
    Rcpp::stop("c++ exception (unknown reason)");
  }
}
