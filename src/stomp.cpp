// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#include "mathtools.h" // math first to fix OSX error
#include "stomp.h"
#include "mass.h"
#include <numeric>
#include <mutex>
#include <vector>
#include <utility>
#include <algorithm>
#include <cfloat> // DBL_EPSILON when STRICT_R_HEADERS
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

// [[Rcpp::export]]
List stomp_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                bool progress, bool left_right_profile) {
  bool partial = false;
  double const exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint32_t const data_size = data_ref.length();
  uint32_t const query_size = query_ref.length();
  uint32_t const matrix_profile_size = data_size - window_size + 1;
  uint32_t const num_queries = query_size - window_size + 1;

  // check skip position
  LogicalVector skip_data(matrix_profile_size, FALSE);
  for (uint64_t i = 0; i < matrix_profile_size; i++) {
    NumericVector const range = data_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_data[i] = TRUE;
    }
  }

  LogicalVector skip_query(num_queries, FALSE);
  for (uint64_t i = 0; i < num_queries; i++) {
    NumericVector const range = query_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_query[i] = TRUE;
    }
  }

  NumericVector data = clone(data_ref);
  NumericVector query = clone(query_ref);

  data[is_na(data)] = 0;
  data[is_infinite(data)] = 0;
  query[is_na(query)] = 0;
  query[is_infinite(query)] = 0;

  NumericVector matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector profile_index(matrix_profile_size, NA_INTEGER);
  NumericVector right_matrix_profile;
  IntegerVector right_profile_index;
  NumericVector left_matrix_profile;
  IntegerVector left_profile_index;

  if (left_right_profile) {
    right_matrix_profile = NumericVector(matrix_profile_size, R_PosInf);
    right_profile_index = IntegerVector(matrix_profile_size, NA_INTEGER);
    left_matrix_profile = NumericVector(matrix_profile_size, R_PosInf);
    left_profile_index = IntegerVector(matrix_profile_size, NA_INTEGER);
  }

  uint32_t const grain = find_best_k_rcpp(data, query, window_size);

  List pre = mass_pre_rcpp(data, query, window_size);
  // NOLINTNEXTLINE(readability-suspicious-call-argument): it is really swapped
  List rpre = mass_pre_rcpp(query, data, window_size); //-V764
  List nn =
      mass3_rcpp(query[Range(0, window_size - 1)], data, pre["data_size"], rpre["window_size"], pre["data_mean"],
                 pre["data_sd"], as<NumericVector>(pre["query_mean"])[0], as<NumericVector>(pre["query_sd"])[0], grain);

  ///// This is needed for JOIN similarity
  List rnn = mass3_rcpp(data[Range(0, window_size - 1)], query, query_size, rpre["window_size"], rpre["data_mean"],
                        rpre["data_sd"], as<NumericVector>(rpre["query_mean"])[0],
                        as<NumericVector>(rpre["query_sd"])[0], grain);

  NumericVector first_product = rnn["last_product"];

  ///// This is needed for JOIN similarity
  IntegerVector const lp_range = Range(1, (data_size - window_size));
  IntegerVector const lp2_range = Range(0, (data_size - window_size - 1));
  IntegerVector const dt_range = Range(window_size, data_size - 1);
  NumericVector distance_profile;
  NumericVector last_product;
  double drop_value = query[0];

  IntegerVector const order = Range(0, num_queries - 1);

  RcppThread::ProgressBar bar(100, 1);

  uint32_t const num_progress =
      ceil((double)order.size() / 100); // added double inside sqrt to avoid ambiguity on Solaris

  try {
    for (int32_t const &i : order) {

      if ((i % num_progress) == 0) {
        RcppThread::checkUserInterrupt();
        if (progress) {
          bar++;
        }
      }
      // compute the distance profile
      NumericVector query_window = query[Range(i, (i + window_size - 1))];
      if (i == 0) {
        distance_profile = nn["distance_profile"];
        last_product = nn["last_product"];
      } else {
        last_product[lp_range] = (NumericVector)(as<NumericVector>(last_product[lp2_range]) -
                                                 as<NumericVector>(data[lp2_range]) * drop_value +
                                                 as<NumericVector>(data[dt_range]) * query_window[window_size - 1]);
        last_product[0] = first_product[i];
        distance_profile =
            2 * (window_size - (last_product - window_size * as<NumericVector>(pre["data_mean"]) *
                                                   as<NumericVector>(pre["query_mean"])[i]) /
                                   (as<NumericVector>(pre["data_sd"]) * as<NumericVector>(pre["query_sd"])[i]));
        distance_profile[distance_profile < 0] = 0;
      }

      // distance_profile = sqrt(distance_profile);
      drop_value = query_window[0];

      // apply exclusion zone
      if (exclusion_zone > 0) {
        uint32_t const exc_st = MAX(0, i - exclusion_zone);
        uint32_t const exc_ed = MIN(matrix_profile_size - 1, i + exclusion_zone);
        IntegerVector const dp_range = Range(exc_st, exc_ed);
        distance_profile[dp_range] = R_PosInf;
      }

      distance_profile[as<NumericVector>(pre["data_sd"]) < DBL_EPSILON] = R_PosInf;
      if (skip_query[i] == TRUE || as<NumericVector>(pre["query_sd"])[i] < DBL_EPSILON) {
        distance_profile.fill(R_PosInf);
      }

      distance_profile[skip_data] = R_PosInf;

      if (left_right_profile) {
        // left matrix_profile
        LogicalVector ldx(matrix_profile_size, FALSE);
        ldx[Range(i, matrix_profile_size - 1)] = (distance_profile[Range(i, matrix_profile_size - 1)] <
                                                  left_matrix_profile[Range(i, matrix_profile_size - 1)]);
        left_matrix_profile[ldx] = distance_profile[ldx];
        left_profile_index[which_cpp(ldx)] = i + 1;

        // right matrix_profile
        LogicalVector rdx(matrix_profile_size, FALSE);
        rdx[Range(0, i)] = (distance_profile[Range(0, i)] < right_matrix_profile[Range(0, i)]);
        right_matrix_profile[rdx] = distance_profile[rdx];
        right_profile_index[which_cpp(rdx)] = i + 1;
      }

      LogicalVector const idx = (distance_profile < matrix_profile);
      matrix_profile[idx] = distance_profile[idx];
      profile_index[which_cpp(idx)] = i + 1;
    }

    if (left_right_profile) {
      left_matrix_profile = sqrt(left_matrix_profile);
      right_matrix_profile = sqrt(right_matrix_profile);
    }

    matrix_profile = sqrt(matrix_profile);
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  if (left_right_profile) {
    return (List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("profile_index") = profile_index,
                         Rcpp::Named("right_matrix_profile") = right_matrix_profile,
                         Rcpp::Named("right_profile_index") = right_profile_index,
                         Rcpp::Named("left_matrix_profile") = left_matrix_profile,
                         Rcpp::Named("left_profile_index") = left_profile_index, Rcpp::Named("partial") = partial,
                         Rcpp::Named("ez") = ez));
  } else {
    return (List::create(Rcpp::Named("matrix_profile") = matrix_profile, Rcpp::Named("profile_index") = profile_index,
                         Rcpp::Named("partial") = partial, Rcpp::Named("ez") = ez));
  }
}

// [[Rcpp::export]]
List stomp_rcpp_parallel(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez,
                         uint32_t n_workers = 1, bool progress = false, bool left_right_profile = false) {
  uint64_t const exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint64_t const data_size = data_ref.length();
  uint64_t const query_size = query_ref.length();
  uint64_t const matrix_profile_size = data_size - window_size + 1;
  uint64_t const num_queries = query_size - window_size + 1;
  bool partial = false;

  // check skip position
  std::vector<int> skip_data(matrix_profile_size, 0);
  for (uint64_t i = 0; i < matrix_profile_size; i++) {
    NumericVector const range = data_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_data[i] = 1;
    }
  }

  std::vector<int> skip_query(num_queries, 0);
  for (uint64_t i = 0; i < num_queries; i++) {
    NumericVector const range = query_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_query[i] = 1;
    }
  }

  NumericVector data = clone(data_ref);
  NumericVector query = clone(query_ref);

  data[is_na(data)] = 0;
  data[is_infinite(data)] = 0;
  query[is_na(query)] = 0;
  query[is_infinite(query)] = 0;

  NumericVector matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector profile_index(matrix_profile_size, NA_INTEGER);
  NumericVector right_matrix_profile;
  IntegerVector right_profile_index;
  NumericVector left_matrix_profile;
  IntegerVector left_profile_index;

  if (left_right_profile) {
    right_matrix_profile = NumericVector(matrix_profile_size, R_PosInf);
    right_profile_index = IntegerVector(matrix_profile_size, NA_INTEGER);
    left_matrix_profile = NumericVector(matrix_profile_size, R_PosInf);
    left_profile_index = IntegerVector(matrix_profile_size, NA_INTEGER);
  }

  uint64_t grain = set_k_rcpp(256, data_size, window_size);

  ///// This is needed for JOIN similarity
  // NOLINTNEXTLINE(readability-suspicious-call-argument): it is really swapped
  List rpre = mass_pre_rcpp(query, data, window_size); //-V764
  List rnn = mass3_rcpp(data[Range(0, window_size - 1)], query, query_size, rpre["window_size"], rpre["data_mean"],
                        rpre["data_sd"], as<NumericVector>(rpre["query_mean"])[0],
                        as<NumericVector>(rpre["query_sd"])[0], grain);

  NumericVector const first_product = rnn["last_product"];
  List pre = mass_pre_rcpp(data, query, window_size);

  const double* data_ptr = data.begin();
  const double* query_ptr = query.begin();
  const double* first_prod_ptr = first_product.begin();
  NumericVector d_mean_vec = pre["data_mean"];
  NumericVector d_std_vec = pre["data_sd"];
  NumericVector q_mean_vec = pre["query_mean"];
  NumericVector q_std_vec = pre["query_sd"];
  const double* d_mean_ptr = d_mean_vec.begin();
  const double* d_std_ptr = d_std_vec.begin();
  const double* q_mean_ptr = q_mean_vec.begin();
  const double* q_std_ptr = q_std_vec.begin();

  double* mp_ptr = matrix_profile.begin();
  int* pi_ptr = profile_index.begin();
  double* rmp_ptr = left_right_profile ? right_matrix_profile.begin() : nullptr;
  int* rpi_ptr = left_right_profile ? right_profile_index.begin() : nullptr;
  double* lmp_ptr = left_right_profile ? left_matrix_profile.begin() : nullptr;
  int* lpi_ptr = left_right_profile ? left_profile_index.begin() : nullptr;

  if (n_workers < 1) {
    n_workers = 1;
  }

  grain = set_k_rcpp(1024, num_queries, window_size);
  size_t chunk_size = 2 * grain;
  if (chunk_size < 1) {
    chunk_size = 1;
  }
  if (num_queries / chunk_size < n_workers) {
    chunk_size = std::max<size_t>(1, (num_queries + n_workers - 1) / n_workers);
  }

  std::vector<std::pair<size_t, size_t>> chunks;
  for (size_t cur = 0; cur < num_queries; cur += chunk_size) {
    size_t cur_end = std::min<size_t>(cur + chunk_size, static_cast<size_t>(num_queries));
    chunks.emplace_back(cur, cur_end);
  }

  std::mutex red_mutex;
  RcppThread::ProgressBar bar(chunks.size(), 1);

  try {
    RcppThread::parallelFor(0, chunks.size(), [&](size_t task_idx) {
      size_t begin = chunks[task_idx].first;
      size_t end = chunks[task_idx].second;

      uint64_t const k_grain = set_k_cpp(window_size * 2, data_size, window_size);
      Mass3Result nn = mass3_cpp_result(query_ptr + begin, data_ptr, data_size, window_size,
                                        d_mean_ptr, d_std_ptr, q_mean_ptr[begin], q_std_ptr[begin], k_grain);
      std::vector<double> distance_profile = std::move(nn.distance_profile);
      std::vector<double> last_product = std::move(nn.last_product);

      std::vector<double> local_mp(matrix_profile_size, R_PosInf);
      std::vector<int> local_pi(matrix_profile_size, NA_INTEGER);
      std::vector<double> local_rmp(left_right_profile ? matrix_profile_size : 0, R_PosInf);
      std::vector<int> local_rpi(left_right_profile ? matrix_profile_size : 0, NA_INTEGER);
      std::vector<double> local_lmp(left_right_profile ? matrix_profile_size : 0, R_PosInf);
      std::vector<int> local_lpi(left_right_profile ? matrix_profile_size : 0, NA_INTEGER);

      double drop_value = 0;
      uint32_t exc_st = 0;
      uint32_t exc_ed = 0;

      for (uint64_t i = begin; i < end; i++) {
        if (i > begin) {
          for (uint64_t j = matrix_profile_size - 1; j > 0; j--) {
            last_product[j] = last_product[j - 1] - data_ptr[j - 1] * drop_value +
                              data_ptr[window_size + j - 1] * query_ptr[i + window_size - 1];
          }
          last_product[0] = first_prod_ptr[i];
        }

        if (exclusion_zone > 0) {
          exc_st = MAX(0, i > exclusion_zone ? (i - exclusion_zone) : 0);
          exc_ed = MIN(matrix_profile_size - 1, i + exclusion_zone);
        }

        for (uint64_t j = 0; j < matrix_profile_size; j++) {
          if (skip_data[j] == 1 || skip_query[i] == 1 || d_std_ptr[j] < DBL_EPSILON || q_std_ptr[i] < DBL_EPSILON) {
            distance_profile[j] = R_PosInf;
          } else if (exclusion_zone == 0 || j < exc_st || j > exc_ed) {
            double const dp = 2 * (window_size - (last_product[j] - window_size * d_mean_ptr[j] * q_mean_ptr[i]) /
                                                     (d_std_ptr[j] * q_std_ptr[i]));
            distance_profile[j] = (dp > 0) ? dp : 0;
          } else {
            distance_profile[j] = R_PosInf;
          }
        }

        drop_value = query_ptr[i];

        for (uint64_t j = 0; j < matrix_profile_size; j++) {
          if (distance_profile[j] < local_mp[j]) {
            local_mp[j] = distance_profile[j];
            local_pi[j] = static_cast<int>(i + 1);
          }
          if (left_right_profile && j >= i && distance_profile[j] < local_lmp[j]) {
            local_lmp[j] = distance_profile[j];
            local_lpi[j] = static_cast<int>(i + 1);
          }
          if (left_right_profile && j <= i && distance_profile[j] < local_rmp[j]) {
            local_rmp[j] = distance_profile[j];
            local_rpi[j] = static_cast<int>(i + 1);
          }
        }
      }

      {
        std::lock_guard<std::mutex> lock(red_mutex);
        for (uint64_t j = 0; j < matrix_profile_size; j++) {
          if (local_mp[j] < mp_ptr[j]) {
            mp_ptr[j] = local_mp[j];
            pi_ptr[j] = local_pi[j];
          }
          if (left_right_profile) {
            if (local_rmp[j] < rmp_ptr[j]) {
              rmp_ptr[j] = local_rmp[j];
              rpi_ptr[j] = local_rpi[j];
            }
            if (local_lmp[j] < lmp_ptr[j]) {
              lmp_ptr[j] = local_lmp[j];
              lpi_ptr[j] = local_lpi[j];
            }
          }
        }
      }

      if (progress) {
        bar++;
      }
    }, n_workers);
  } catch (RcppThread::UserInterruptException &e) {
    partial = true;
    Rcout << "Process terminated by the user successfully, partial results were returned." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  if (left_right_profile) {
    return List::create(Rcpp::Named("matrix_profile") = sqrt(matrix_profile),
                        Rcpp::Named("profile_index") = profile_index,
                        Rcpp::Named("right_matrix_profile") = sqrt(right_matrix_profile),
                        Rcpp::Named("right_profile_index") = right_profile_index,
                        Rcpp::Named("left_matrix_profile") = sqrt(left_matrix_profile),
                        Rcpp::Named("left_profile_index") = left_profile_index,
                        Rcpp::Named("partial") = partial,
                        Rcpp::Named("ez") = ez);
  }

  return List::create(Rcpp::Named("matrix_profile") = sqrt(matrix_profile),
                      Rcpp::Named("profile_index") = profile_index,
                      Rcpp::Named("partial") = partial,
                      Rcpp::Named("ez") = ez);
}
