#include "mathtools.h"
#include "mstomp.h"
#include "mass.h"
#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cfloat>
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

// [[Rcpp::export]]
List mstomp_rcpp(const NumericMatrix data_ref, uint32_t window_size, double ez,
                 Nullable<IntegerVector> must_dim,
                 Nullable<IntegerVector> exc_dim,
                 uint32_t n_workers, bool progress) {
  bool partial = false;
  uint64_t const data_size = data_ref.nrow();
  uint64_t const n_dim = data_ref.ncol();
  uint64_t const matrix_profile_size = data_size - window_size + 1;
  uint64_t const exclusion_zone = static_cast<uint64_t>(round(window_size * ez + DBL_EPSILON));

  // Parse must_dim and exc_dim
  std::vector<int> must_dims;
  std::vector<int> exc_dims;
  std::vector<bool> mask_must(n_dim, false);
  std::vector<bool> mask_exc(n_dim, false);

  if (must_dim.isNotNull()) {
    IntegerVector mv(must_dim);
    for (int idx_1based : mv) {
      if (idx_1based < 1 || static_cast<uint64_t>(idx_1based) > n_dim) {
        Rcpp::stop("`must_dim` must be less then the total dimension.");
      }
      int idx_0based = idx_1based - 1;
      if (!mask_must[idx_0based]) {
        mask_must[idx_0based] = true;
        must_dims.push_back(idx_0based);
      }
    }
  }

  if (exc_dim.isNotNull()) {
    IntegerVector ev(exc_dim);
    for (int idx_1based : ev) {
      if (idx_1based < 1 || static_cast<uint64_t>(idx_1based) > n_dim) {
        Rcpp::stop("`exc_dim` must be less then the total dimension.");
      }
      int idx_0based = idx_1based - 1;
      if (mask_must[idx_0based]) {
        Rcpp::stop("The same dimension is presented in both the exclusion dimension and must have dimension.");
      }
      if (!mask_exc[idx_0based]) {
        mask_exc[idx_0based] = true;
        exc_dims.push_back(idx_0based);
      }
    }
  }

  std::vector<int> other_dims;
  for (uint64_t d = 0; d < n_dim; ++d) {
    if (!mask_must[d]) {
      other_dims.push_back(static_cast<int>(d));
    }
  }

  // Check skip positions
  std::vector<uint8_t> skip_location(matrix_profile_size, 0);
  for (uint64_t i = 0; i < matrix_profile_size; ++i) {
    for (uint64_t d = 0; d < n_dim; ++d) {
      if (!mask_exc[d]) {
        for (uint64_t t = 0; t < window_size; ++t) {
          double val = data_ref(i + t, d);
          if (ISNA(val) || ISNAN(val) || !R_FINITE(val)) {
            skip_location[i] = 1;
            break;
          }
        }
        if (skip_location[i] == 1) {
          break;
        }
      }
    }
  }

  // Clean data copy (replacing NA/Inf with 0.0)
  std::vector<double> clean_data(data_size * n_dim);
  for (uint64_t d = 0; d < n_dim; ++d) {
    for (uint64_t i = 0; i < data_size; ++i) {
      double val = data_ref(i, d);
      if (ISNA(val) || ISNAN(val) || !R_FINITE(val)) {
        val = 0.0;
      }
      clean_data[d * data_size + i] = val;
    }
  }

  // Precomputation of means, standard deviations, and first products
  std::vector<double> all_means(n_dim * matrix_profile_size);
  std::vector<double> all_sds(n_dim * matrix_profile_size);
  std::vector<double> all_first_products(n_dim * matrix_profile_size);

  for (uint64_t d = 0; d < n_dim; ++d) {
    const double *col_ptr = clean_data.data() + d * data_size;
    NumericVector col_v(col_ptr, col_ptr + data_size);
    List pre = mass_pre_rcpp(col_v, col_v, window_size);
    NumericVector m_v = pre["data_mean"];
    NumericVector s_v = pre["data_sd"];
    for (uint64_t i = 0; i < matrix_profile_size; ++i) {
      all_means[d * matrix_profile_size + i] = m_v[i];
      all_sds[d * matrix_profile_size + i] = s_v[i];
    }
    uint64_t grain = set_k_rcpp(256, data_size, window_size);
    List rnn = mass3_rcpp(col_v[Range(0, window_size - 1)], col_v, data_size, window_size,
                          m_v, s_v, m_v[0], s_v[0], grain);
    NumericVector fp_v = rnn["last_product"];
    for (uint64_t i = 0; i < matrix_profile_size; ++i) {
      all_first_products[d * matrix_profile_size + i] = fp_v[i];
    }
  }

  // Output matrices
  NumericMatrix matrix_profile(matrix_profile_size, n_dim);
  NumericMatrix profile_index(matrix_profile_size, n_dim);
  NumericMatrix left_matrix_profile(matrix_profile_size, n_dim);
  NumericMatrix left_profile_index(matrix_profile_size, n_dim);
  NumericMatrix right_matrix_profile(matrix_profile_size, n_dim);
  NumericMatrix right_profile_index(matrix_profile_size, n_dim);

  std::fill(matrix_profile.begin(), matrix_profile.end(), R_PosInf);
  std::fill(profile_index.begin(), profile_index.end(), R_NegInf);
  std::fill(left_matrix_profile.begin(), left_matrix_profile.end(), R_PosInf);
  std::fill(left_profile_index.begin(), left_profile_index.end(), R_NegInf);
  std::fill(right_matrix_profile.begin(), right_matrix_profile.end(), R_PosInf);
  std::fill(right_profile_index.begin(), right_profile_index.end(), R_NegInf);

  if (n_workers < 1) {
    n_workers = 1;
  }

  size_t chunk_size = std::max<size_t>(1, (matrix_profile_size + n_workers - 1) / n_workers);
  if (matrix_profile_size / chunk_size < n_workers) {
    chunk_size = std::max<size_t>(1, matrix_profile_size / (n_workers * 2));
  }

  std::vector<std::pair<size_t, size_t>> chunks;
  for (size_t cur = 0; cur < matrix_profile_size; cur += chunk_size) {
    size_t cur_end = std::min<size_t>(cur + chunk_size, static_cast<size_t>(matrix_profile_size));
    chunks.emplace_back(cur, cur_end);
  }

  RcppThread::ProgressBar bar(chunks.size(), 1);

  try {
    RcppThread::parallelFor(0, chunks.size(), [&](size_t task_idx) {
      size_t begin = chunks[task_idx].first;
      size_t end = chunks[task_idx].second;

      std::vector<double> last_product(n_dim * matrix_profile_size);
      std::vector<double> sort_buf(other_dims.size());
      std::vector<double> dist_merged(matrix_profile_size * n_dim, R_PosInf);
      std::vector<double> row_dist(n_dim, R_PosInf);

      // Initialize last_product for query `begin`
      for (uint64_t d = 0; d < n_dim; ++d) {
        if (begin == 0) {
          for (uint64_t j = 0; j < matrix_profile_size; ++j) {
            last_product[d * matrix_profile_size + j] = all_first_products[d * matrix_profile_size + j];
          }
        } else {
          const double *col_ptr = clean_data.data() + d * data_size;
          const double *m_ptr = all_means.data() + d * matrix_profile_size;
          const double *s_ptr = all_sds.data() + d * matrix_profile_size;
          uint64_t const k_grain = set_k_cpp(window_size * 2, data_size, window_size);
          Mass3Result nn = mass3_cpp_result(col_ptr + begin, col_ptr, data_size, window_size,
                                            m_ptr, s_ptr, m_ptr[begin], s_ptr[begin], k_grain);
          for (uint64_t j = 0; j < matrix_profile_size; ++j) {
            last_product[d * matrix_profile_size + j] = nn.last_product[j];
          }
        }
      }

      for (size_t i = begin; i < end; ++i) {
        if (i > begin) {
          for (uint64_t d = 0; d < n_dim; ++d) {
            const double *col_ptr = clean_data.data() + d * data_size;
            double drop_val = col_ptr[i - 1];
            double query_val = col_ptr[i + window_size - 1];
            double *lp = last_product.data() + d * matrix_profile_size;
            for (uint64_t j = matrix_profile_size - 1; j > 0; --j) {
              lp[j] = lp[j - 1] - col_ptr[j - 1] * drop_val + col_ptr[window_size + j - 1] * query_val;
            }
            lp[0] = all_first_products[d * matrix_profile_size + i];
          }
        }

        bool query_invalid = (skip_location[i] == 1);
        if (!query_invalid) {
          for (uint64_t d = 0; d < n_dim; ++d) {
            if (!mask_exc[d] && all_sds[d * matrix_profile_size + i] < DBL_EPSILON) {
              query_invalid = true;
              break;
            }
          }
        }

        if (query_invalid) {
          for (uint64_t p = 0; p < n_dim; ++p) {
            matrix_profile(i, p) = R_PosInf;
            profile_index(i, p) = 1.0;
            left_matrix_profile(i, p) = R_PosInf;
            left_profile_index(i, p) = R_NegInf;
            right_matrix_profile(i, p) = R_PosInf;
            right_profile_index(i, p) = R_NegInf;
          }
          continue;
        }

        int64_t exc_st = std::max<int64_t>(0, static_cast<int64_t>(i) - static_cast<int64_t>(exclusion_zone));
        int64_t exc_ed = std::min<int64_t>(static_cast<int64_t>(matrix_profile_size) - 1, static_cast<int64_t>(i) + static_cast<int64_t>(exclusion_zone));

        for (uint64_t k = 0; k < matrix_profile_size; ++k) {
          bool cand_invalid = (skip_location[k] == 1) || (static_cast<int64_t>(k) >= exc_st && static_cast<int64_t>(k) <= exc_ed);
          if (cand_invalid) {
            for (uint64_t p = 0; p < n_dim; ++p) {
              dist_merged[k * n_dim + p] = R_PosInf;
            }
          } else {
            // Must dimensions
            for (size_t m_idx = 0; m_idx < must_dims.size(); ++m_idx) {
              int d = must_dims[m_idx];
              double sd_k = all_sds[d * matrix_profile_size + k];
              if (sd_k < DBL_EPSILON) {
                row_dist[m_idx] = R_PosInf;
              } else {
                double lp = last_product[d * matrix_profile_size + k];
                double mean_k = all_means[d * matrix_profile_size + k];
                double mean_i = all_means[d * matrix_profile_size + i];
                double sd_i = all_sds[d * matrix_profile_size + i];
                double dp = 2.0 * (window_size - (lp - window_size * mean_k * mean_i) / (sd_k * sd_i));
                row_dist[m_idx] = (dp > 0.0) ? dp : 0.0;
              }
            }

            // Other dimensions
            for (size_t o_idx = 0; o_idx < other_dims.size(); ++o_idx) {
              int d = other_dims[o_idx];
              if (mask_exc[d]) {
                sort_buf[o_idx] = R_PosInf;
              } else {
                double sd_k = all_sds[d * matrix_profile_size + k];
                if (sd_k < DBL_EPSILON) {
                  sort_buf[o_idx] = R_PosInf;
                } else {
                  double lp = last_product[d * matrix_profile_size + k];
                  double mean_k = all_means[d * matrix_profile_size + k];
                  double mean_i = all_means[d * matrix_profile_size + i];
                  double sd_i = all_sds[d * matrix_profile_size + i];
                  double dp = 2.0 * (window_size - (lp - window_size * mean_k * mean_i) / (sd_k * sd_i));
                  sort_buf[o_idx] = (dp > 0.0) ? dp : 0.0;
                }
              }
            }

            std::sort(sort_buf.begin(), sort_buf.end());
            for (size_t o_idx = 0; o_idx < other_dims.size(); ++o_idx) {
              row_dist[must_dims.size() + o_idx] = sort_buf[o_idx];
            }

            double cum_sum = 0.0;
            for (size_t p = 0; p < n_dim; ++p) {
              cum_sum += row_dist[p];
              dist_merged[k * n_dim + p] = cum_sum / (p + 1);
            }
          }
        }

        // Nearest neighbors extraction for query i
        size_t start_p = must_dims.empty() ? 0 : (must_dims.size() - 1);
        size_t end_p = n_dim - exc_dims.size();

        for (size_t p = start_p; p < end_p; ++p) {
          // Overall nearest neighbor
          double min_val = R_PosInf;
          int min_idx = 0;
          for (uint64_t k = 0; k < matrix_profile_size; ++k) {
            double v = dist_merged[k * n_dim + p];
            if (v < min_val) {
              min_val = v;
              min_idx = static_cast<int>(k);
            }
          }
          if (std::isinf(min_val)) {
            min_idx = 0;
          }
          matrix_profile(i, p) = min_val;
          profile_index(i, p) = min_idx + 1; // 1-based

          // Left nearest neighbor
          if (static_cast<int64_t>(i) > static_cast<int64_t>(exclusion_zone)) {
            double l_min_val = R_PosInf;
            int l_min_idx = 0;
            size_t limit = i - exclusion_zone;
            for (size_t k = 0; k < limit; ++k) {
              double v = dist_merged[k * n_dim + p];
              if (v < l_min_val) {
                l_min_val = v;
                l_min_idx = static_cast<int>(k);
              }
            }
            left_matrix_profile(i, p) = l_min_val;
            left_profile_index(i, p) = l_min_idx + 1;
          } else {
            left_matrix_profile(i, p) = R_PosInf;
            left_profile_index(i, p) = R_NegInf;
          }

          // Right nearest neighbor
          if (i + exclusion_zone + 1 < matrix_profile_size) {
            double r_min_val = R_PosInf;
            size_t start_k = i + exclusion_zone + 1;
            int r_min_idx = static_cast<int>(start_k);
            for (size_t k = start_k; k < matrix_profile_size; ++k) {
              double v = dist_merged[k * n_dim + p];
              if (v < r_min_val) {
                r_min_val = v;
                r_min_idx = static_cast<int>(k);
              }
            }
            right_matrix_profile(i, p) = r_min_val;
            right_profile_index(i, p) = r_min_idx + 1;
          } else {
            right_matrix_profile(i, p) = R_PosInf;
            right_profile_index(i, p) = R_NegInf;
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
  } catch (const std::exception &e) {
    Rcpp::stop(e.what());
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  // Sqrt on distance matrices
  for (uint64_t cell = 0; cell < matrix_profile_size * n_dim; ++cell) {
    if (matrix_profile[cell] < R_PosInf) {
      matrix_profile[cell] = std::sqrt(matrix_profile[cell]);
    }
    if (left_matrix_profile[cell] < R_PosInf) {
      left_matrix_profile[cell] = std::sqrt(left_matrix_profile[cell]);
    }
    if (right_matrix_profile[cell] < R_PosInf) {
      right_matrix_profile[cell] = std::sqrt(right_matrix_profile[cell]);
    }
  }

  // Mask invalid k settings (must_dims and exc_dims)
  if (must_dims.size() > 1) {
    for (size_t p = 0; p < must_dims.size() - 1; ++p) {
      for (uint64_t i = 0; i < matrix_profile_size; ++i) {
        matrix_profile(i, p) = NA_REAL;
        left_matrix_profile(i, p) = NA_REAL;
        right_matrix_profile(i, p) = NA_REAL;
        profile_index(i, p) = NA_REAL;
        left_profile_index(i, p) = NA_REAL;
        right_profile_index(i, p) = NA_REAL;
      }
    }
  }

  if (!exc_dims.empty()) {
    for (size_t p = n_dim - exc_dims.size(); p < n_dim; ++p) {
      for (uint64_t i = 0; i < matrix_profile_size; ++i) {
        matrix_profile(i, p) = NA_REAL;
        left_matrix_profile(i, p) = NA_REAL;
        right_matrix_profile(i, p) = NA_REAL;
        profile_index(i, p) = NA_REAL;
        left_profile_index(i, p) = NA_REAL;
        right_profile_index(i, p) = NA_REAL;
      }
    }
  }

  return List::create(
    Rcpp::Named("matrix_profile") = matrix_profile,
    Rcpp::Named("profile_index") = profile_index,
    Rcpp::Named("right_matrix_profile") = right_matrix_profile,
    Rcpp::Named("right_profile_index") = right_profile_index,
    Rcpp::Named("left_matrix_profile") = left_matrix_profile,
    Rcpp::Named("left_profile_index") = left_profile_index,
    Rcpp::Named("ez") = ez,
    Rcpp::Named("w") = window_size,
    Rcpp::Named("n_dim") = static_cast<int>(n_dim),
    Rcpp::Named("must") = must_dim,
    Rcpp::Named("exc") = exc_dim,
    Rcpp::Named("partial") = partial
  );
}
