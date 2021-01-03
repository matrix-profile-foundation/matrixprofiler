#include "scrimp.h"
#include "mass.h"
#include "math.h"
#include <numeric>
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
List scrimp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez = 0.5,
                 double pre_scrimp = 0.25, bool progress = false) {
  // double s_size = R_PosInf;
  bool partial = false;
  double exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint32_t data_size = data_ref.length();
  // uint32_t query_size = query_ref.length();

  uint32_t matrix_profile_size = data_size - window_size + 1;
  // uint32_t num_queries = query_size - window_size + 1;

  // check skip position
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

  IntegerVector orig_index = Range(0, matrix_profile_size - 1);

  uint32_t k = set_k_rcpp(window_size, data_size, window_size);

  List pre = mass_pre_rcpp(data, query, window_size);

  NumericVector data_mean = pre["data_mean"];
  NumericVector data_sd = pre["data_sd"];
  NumericVector query_mean = pre["query_mean"];
  NumericVector query_sd = pre["query_sd"];

  List nn = mass3_rcpp(query[Range(0, window_size - 1)], data, pre["data_size"], pre["window_size"], data_mean, data_sd,
                       query_mean[0], query_sd[0], k);

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

      uint64_t j = 1;
      for (uint64_t &&i : pre_scrimp_idxs) {

        RcppThread::checkUserInterrupt();
        ps.increment();

        // compute the distance profile
        nn = mass3_rcpp(query[Range(i, i + window_size - 1)], data, pre["data_size"], pre["window_size"], data_mean,
                        data_sd, query_mean[i], query_sd[i], k);

        NumericVector distance_profile = sqrt(as<NumericVector>(nn["distance_profile"]));

        uint64_t exc_st = 0;
        uint64_t exc_ed = 0;
        // apply exclusion zone
        if (exclusion_zone > 0) {
          exc_st = MAX(0, i > exclusion_zone ? (i - exclusion_zone) : 0);
          exc_ed = MIN(matrix_profile_size - 1, i + exclusion_zone);
        }

        IntegerVector exc_idxs = Range(exc_st, exc_ed);
        distance_profile[exc_idxs] = R_PosInf;

        // figure out and store the neareest neighbor
        if (j == 1) {
          matrix_profile = distance_profile;
          profile_index.fill(i);
          uint64_t min_idx = which_min(distance_profile);
          profile_index[i] = min_idx;
          matrix_profile[i] = distance_profile[min_idx];
        } else {
          LogicalVector update_pos = distance_profile < matrix_profile;
          profile_index[update_pos] = (int32_t)i;
          matrix_profile[update_pos] = distance_profile[update_pos];
          uint64_t min_idx = which_min(distance_profile);
          profile_index[i] = min_idx;
          matrix_profile[i] = distance_profile[min_idx];
        }

        uint64_t idx_nn = profile_index[i];
        int64_t idx_diff = idx_nn - i;
        dotproduct[i] = (window_size - (pow(matrix_profile[i], 2) / 2)) * data_sd[i] * data_sd[idx_nn] +
                        window_size * data_mean[i] * data_mean[idx_nn];

        uint64_t endidx = MIN(matrix_profile_size - 1, (int64_t)i + current_step - 1);
        endidx = MIN((int64_t)endidx, matrix_profile_size - idx_diff - 1);

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
        refine_distance[dot_idxs1] = sqrt(
            abs(2 * (window_size - (dotproduct[dot_idxs1] - data_mean[dot_idxs1] * data_mean[ref_idxs1] * window_size) /
                                       (data_sd[dot_idxs1] * data_sd[ref_idxs1]))));

        int64_t beginidx = MAX(0, (int64_t)i - current_step + 1);
        beginidx = MAX(beginidx, 1 - idx_diff);

        if (i > 0) {
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

          refine_distance[ref_idxs2] = sqrt(abs(
              2 * (window_size - (dotproduct[ref_idxs2] - data_mean[ref_idxs2] * data_mean[ref_idxs3] * window_size) /
                                     (data_sd[ref_idxs2] * data_sd[ref_idxs3]))));
        }

        Range upd_idxs1 = Range(beginidx, endidx);
        Range upd_idxs2 = Range((beginidx + idx_diff), (endidx + idx_diff));
        IntegerVector update_pos1 = which(refine_distance[upd_idxs1] < matrix_profile[upd_idxs1]);
        matrix_profile[(update_pos1 + beginidx)] = refine_distance[(update_pos1 + beginidx)];
        IntegerVector new_idxs = as<IntegerVector>(orig_index[(update_pos1 + beginidx)]) + idx_diff;
        profile_index[(update_pos1 + beginidx)] = new_idxs;
        IntegerVector update_pos2 = which(refine_distance[upd_idxs1] < matrix_profile[upd_idxs2]);
        matrix_profile[(update_pos2 + beginidx + idx_diff)] = refine_distance[(update_pos2 + beginidx)];
        new_idxs = as<IntegerVector>(orig_index[(update_pos2 + beginidx + idx_diff)]) - idx_diff;
        profile_index[(update_pos2 + beginidx + idx_diff)] = new_idxs;

        j = j + 1;
      }
    }

    //// SCRIMP ----

    IntegerVector compute_order = orig_index[orig_index > exclusion_zone];
    // uint64_t ssize = MIN(s_size, order.size());
    // order = sample(order, ssize);
    // compute_order = sample(compute_order, compute_order.size());

    NumericVector curlastz(matrix_profile_size);
    NumericVector curdistance(matrix_profile_size);
    NumericVector dist1(matrix_profile_size, R_PosInf);
    NumericVector dist2(matrix_profile_size, R_PosInf);

    Progress p(compute_order.size(), progress);

    uint64_t j = 1;
    for (uint64_t &&i : compute_order) {

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
      curdistance[Range(i, matrix_profile_size - 1)] = sqrt(abs(
          2 *
          (window_size -
           (curlastz[Range(i, matrix_profile_size - 1)] - window_size * data_mean[Range(i, matrix_profile_size - 1)] *
                                                              data_mean[Range(0, matrix_profile_size - i - 1)]) /
               (data_sd[Range(i, matrix_profile_size - 1)] * data_sd[Range(0, matrix_profile_size - i - 1)]))));

      dist1[::seq(0, i - 1)] = R_PosInf;
      dist1[::seq(i, matrix_profile_size - 1)] = curdistance[::seq(i, matrix_profile_size - 1)];

      if (i < (matrix_profile_size - 1)) {
        dist2[::seq(0, matrix_profile_size - i - 1)] = curdistance[::seq(i, matrix_profile_size - 1)];
      }

      dist2[::seq(matrix_profile_size - i, matrix_profile_size - 1)] = R_PosInf;

      LogicalVector loc1 = dist1 < matrix_profile;
      matrix_profile[loc1] = dist1[loc1];
      profile_index[loc1] = (IntegerVector)(as<IntegerVector>(orig_index[loc1]) - i);

      LogicalVector loc2 = dist2 < matrix_profile;
      matrix_profile[loc2] = dist2[loc2];
      profile_index[loc2] = (IntegerVector)(as<IntegerVector>(orig_index[loc2]) + i - 1);

      j++;
    }

  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    std::cout << "Process terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  // profile_index = profile_index + 1;

  return (List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("profile_index") = profile_index,
                       Rcpp::Named("partial") = partial, Rcpp::Named("ez") = ez));
}

// struct ScrimpWorker : public Worker {
//   // input
//   const RVector<double> data_ref;
//   const RVector<double> window_ref;
//   const uint64_t w_size;
//   const uint64_t d_size;
//   const RVector<double> d_mean;
//   const RVector<double> d_std;
//   const RVector<double> q_mean;
//   const RVector<double> q_std;
//   const RVector<int> skip_location;
//   const RVector<double> first_product;
//   const uint64_t ez;

//   Progress *p;

//   RVector<double> mp;
//   RVector<int> pi;

// #if RCPP_PARALLEL_USE_TBB
//   tbb::mutex m;
// #else
//   tthread::mutex m;
// #endif

//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   ScrimpWorker(const NumericVector data_ref, const NumericVector window_ref,
//                const uint64_t w_size, const uint64_t d_size,
//                const NumericVector d_mean, const NumericVector d_std,
//                const NumericVector q_mean, const NumericVector q_std,
//                const IntegerVector skip_location,
//                const NumericVector first_product, const uint64_t ez,
//                Progress *p, NumericVector mp, IntegerVector pi)
//       : data_ref(data_ref), window_ref(window_ref), w_size(w_size),
//         d_size(d_size), d_mean(d_mean), d_std(d_std), q_mean(q_mean),
//         q_std(q_std), skip_location(skip_location),
//         first_product(first_product), ez(ez), p(p), mp(mp), pi(pi) {}

//   ~ScrimpWorker() {}

//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     // begin and end are the query window

//     // index of sliding window
//     try {
//       RcppThread::checkUserInterrupt();

//       uint64_t chunk = (end - begin);

//       if (chunk <= w_size) {
//         std::cout << "Chunk size is too small (" << chunk
//                   << ") for a window size of " << w_size << std::endl;
//         return;
//       }

//       uint64_t k = set_k_rcpp(w_size * 2, chunk, w_size);

//       m.lock();
//       List nn = mass3_cpp(window_ref.begin() + begin, data_ref.begin(),
//       d_size,
//                           w_size, d_mean.begin(), d_std.begin(),
//                           q_mean[begin], q_std[begin], k);
//       m.unlock();
//       std::vector<double> distance_profile(
//           as<NumericVector>(nn["distance_profile"]).begin(),
//           as<NumericVector>(nn["distance_profile"]).end());
//       std::vector<double> last_product(
//           as<NumericVector>(nn["last_product"]).begin(),
//           as<NumericVector>(nn["last_product"]).end());

//       std::vector<double> matrix_profile(d_mean.size(), R_PosInf);
//       std::vector<int> profile_index(d_mean.size(), -1);
//       double drop_value = 0;
//       uint32_t exc_st = 0;
//       uint32_t exc_ed = 0;

//       for (uint64_t i = begin; i < end; i++) {

//         if (i % 100 == 0) {
//           RcppThread::checkUserInterrupt();
//           m.lock();
//           p->increment();
//           m.unlock();
//         }

//         // compute the distance profile
//         if (i > begin) {

//           for (uint64_t j = d_mean.size() - 1; j > 0; j--) {
//             last_product[j] =
//                 last_product[j - 1] - data_ref[j - 1] * drop_value +
//                 data_ref[w_size + j - 1] * window_ref[i + w_size - 1];
//           }
//           last_product[0] = first_product[i];
//         }

//         if (ez > 0) {
//           exc_st = MAX(0, i > ez ? (i - ez) : 0);
//           exc_ed = MIN(d_std.size() - 1, i + ez);
//         }

//         for (uint64_t j = 0; j < d_mean.size(); j++) {
//           double dp = R_PosInf;

//           if (skip_location[j] == 0) {
//             if (ez == 0 || j < exc_st || j > exc_ed) {

//               dp = 2 * (w_size -
//                         (last_product[j] - w_size * d_mean[j] * q_mean[i]) /
//                             (d_std[j] * q_std[i]));
//             } else if (i == begin) {
//               distance_profile[j] = R_PosInf;
//               continue;
//             }
//           }

//           distance_profile[j] = (dp > 0) ? dp : 0;
//         }

//         drop_value = window_ref[i];

//         for (uint64_t j = 0; j < d_mean.size(); j++) {
//           if (distance_profile[j] < matrix_profile[j]) {
//             matrix_profile[j] = distance_profile[j];
//             profile_index[j] = i + 1;
//           }
//         }
//       }

//       m.lock();
//       for (uint64_t j = 0; j < d_mean.size(); j++) {
//         if (matrix_profile[j] < mp[j]) {
//           mp[j] = matrix_profile[j];
//           pi[j] = profile_index[j];
//         }
//       }
//       m.unlock();

//     } catch (RcppThread::UserInterruptException &e) {
//       Rcout << "Computation interrupted by the user." << std::endl;
//       Rcout << "Please wait for other threads to stop." << std::endl;
//       throw;
//     }
//   }
// };

// // [[Rcpp::export]]
// List scrimp_rcpp_parallel(const NumericVector data_ref,
//                           const NumericVector query_ref, uint32_t
//                           window_size, double ez = 0.5, bool progress =
//                           false) {
//   uint64_t exclusion_zone = round(window_size * ez + DBL_EPSILON);
//   uint64_t data_size = data_ref.length();
//   uint64_t query_size = query_ref.length();
//   uint64_t matrix_profile_size = data_size - window_size + 1;
//   uint64_t num_queries = query_size - window_size + 1;
//   bool partial = false;

//   // check skip position
//   IntegerVector skip_location(matrix_profile_size, 0);

//   for (uint64_t i = 0; i < matrix_profile_size; i++) {
//     NumericVector range = data_ref[Range(i, (i + window_size - 1))];
//     if (any(is_na(range) | is_infinite(range))) {
//       skip_location[i] = 1;
//     }
//   }

//   NumericVector data = data_ref;
//   NumericVector query = query_ref;

//   data[is_na(data)] = 0;
//   data[is_infinite(data)] = 0;
//   query[is_na(query)] = 0;
//   query[is_infinite(query)] = 0;

//   NumericVector matrix_profile(matrix_profile_size, R_PosInf);
//   IntegerVector profile_index(matrix_profile_size, -1);

//   uint64_t k = set_k_rcpp(256, data_size, window_size);

//   ///// This is needed for JOIN similarity
//   List rpre = mass_pre_rcpp(query, data, window_size);
//   List rnn = mass3_rcpp(data[Range(0, window_size - 1)], query, query_size,
//                         rpre["window_size"], rpre["data_mean"],
//                         rpre["data_sd"],
//                         as<NumericVector>(rpre["query_mean"])[0],
//                         as<NumericVector>(rpre["query_sd"])[0], k);

//   NumericVector first_product = rnn["last_product"];

//   List pre = mass_pre_rcpp(data, query, window_size);

//   Progress p(num_queries / 100, progress);

//   ScrimpWorker scrimp_worker(
//       data, query, pre["window_size"], data_size, pre["data_mean"],
//       pre["data_sd"], pre["query_mean"], pre["query_sd"], skip_location,
//       first_product, exclusion_zone, &p, matrix_profile, profile_index);

//   k = set_k_rcpp(1024, num_queries, window_size);

//   // call parallelFor to do the work
//   try {
// #if RCPP_PARALLEL_USE_TBB
//     RcppParallel::parallelFor(0, num_queries, scrimp_worker, 2 * k);
// #else
//     RcppParallel2::ttParallelFor(0, num_queries, scrimp_worker, 2 * k);
// #endif

//   } catch (RcppThread::UserInterruptException &e) {
//     partial = true;
//     Rcout << "Process terminated by the user successfully, partial results "
//              "were returned.";
//   } catch (...) {
//     ::Rf_error("c++ exception (unknown reason)");
//   }

//   return (List::create(Rcpp::Named("matrix_profile") = sqrt(matrix_profile),
//                        Rcpp::Named("profile_index") = profile_index,
//                        Rcpp::Named("partial") = partial,
//                        Rcpp::Named("ez") = ez));
// }
