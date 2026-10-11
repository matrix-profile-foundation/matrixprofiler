// This is a personal academic project. Dear PVS-Studio, please check it.

// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#include "mathtools.h" // math first to fix OSX error
#include "mass.h"
#include "fft.h"
#include "windowfunc.h"
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
#include <memory>
#include <Rcpp/Benchmark/Timer.h>
using std::multiplies;

inline void clamp_negative(NumericVector &v) {
  for (R_xlen_t i = 0; i < v.length(); ++i) {
    if (!NumericVector::is_na(v[i]) && !std::isnan(v[i]) && v[i] < 0.0) {
      v[i] = 0.0;
    }
  }
}

// #### MASS functions ####
//[[Rcpp::export]]
List mass_weighted_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint32_t data_size,
                        uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd,
                        double query_mean, double query_sd, const NumericVector data_pre, const NumericVector weight,
                        const bool normalized, uint32_t n_workers = 1) {
  const uint32_t profile_len = data_size - window_size + 1;
  NumericVector distance_profile(profile_len);
  NumericVector last_product(profile_len);

  try {
    const R_xlen_t fft_len = data_fft.length();
    NumericVector query;
    NumericVector rev_query(fft_len);
    NumericVector rev_weight(fft_len);

    // normalized query
    if (normalized) {
      query = (query_window - query_mean) / query_sd;
    } else {
      query = query_window;
    }

    double const sumwy = sum(query * weight);
    double const sumwy2 = sum(weight * query * query);

    std::reverse_copy(query.begin(), query.end(), rev_query.begin());
    std::reverse_copy(weight.begin(), weight.end(), rev_weight.begin());

    ComplexVector const q_fft = fft_rcpp(rev_weight * rev_query);
    ComplexVector prod(fft_len);

    uint32_t actual_workers = std::max<uint32_t>(1, n_workers);

    if (actual_workers > 1 && fft_len >= 1024) {
      const Rcomplex *d_ptr = reinterpret_cast<const Rcomplex *>(data_fft.begin());
      const Rcomplex *q_ptr = reinterpret_cast<const Rcomplex *>(q_fft.begin());
      Rcomplex *p_ptr = reinterpret_cast<Rcomplex *>(prod.begin());

      RcppThread::parallelFor(0, static_cast<std::size_t>(fft_len), [&](std::size_t i) {
        p_ptr[i].r = d_ptr[i].r * q_ptr[i].r - d_ptr[i].i * q_ptr[i].i;
        p_ptr[i].i = d_ptr[i].r * q_ptr[i].i + d_ptr[i].i * q_ptr[i].r;
      }, actual_workers);
    } else {
      prod = data_fft * q_fft;
    }

    NumericVector data_queryw = Re(fft_rcpp(prod, true));

    double *dp = distance_profile.begin();
    double *lp = last_product.begin();
    const double *dq_w = data_queryw.begin();
    const double *d_pre = data_pre.begin();
    const double *d_mean = data_mean.begin();
    const double *d_sd = data_sd.begin();

    if (actual_workers > 1 && profile_len >= 1024) {
      RcppThread::parallelFor(0, static_cast<std::size_t>(profile_len), [&](std::size_t i) {
        double prod_val = dq_w[window_size - 1 + i];
        lp[i] = prod_val;
        double dist = d_pre[i] - 2.0 * (prod_val - sumwy * d_mean[i]) / d_sd[i] + sumwy2;
        if (!NumericVector::is_na(dist) && !std::isnan(dist) && dist < 0.0) {
          dist = 0.0;
        }
        dp[i] = dist;
      }, actual_workers);
    } else {
      for (std::size_t i = 0; i < profile_len; ++i) {
        double prod_val = dq_w[window_size - 1 + i];
        lp[i] = prod_val;
        double dist = d_pre[i] - 2.0 * (prod_val - sumwy * d_mean[i]) / d_sd[i] + sumwy2;
        if (!NumericVector::is_na(dist) && !std::isnan(dist) && dist < 0.0) {
          dist = 0.0;
        }
        dp[i] = dist;
      }
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (const std::exception &ex) {
    Rcpp::stop(ex.what());
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product);
}

//[[Rcpp::export]]
List mass_absolute_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint32_t data_size,
                        uint32_t window_size, const NumericVector sumx2, double sumy2, uint32_t n_workers = 1) {
  const uint32_t profile_len = data_size - window_size + 1;
  NumericVector distance_profile(profile_len);
  NumericVector last_product(profile_len);

  try {
    const R_xlen_t fft_len = data_fft.length();
    NumericVector rev_query(fft_len); // use data_fft because is padded to a power of 2
    std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());

    ComplexVector const q_fft = fft_rcpp(rev_query);
    ComplexVector prod(fft_len);

    uint32_t actual_workers = std::max<uint32_t>(1, n_workers);

    if (actual_workers > 1 && fft_len >= 1024) {
      const Rcomplex *d_ptr = reinterpret_cast<const Rcomplex *>(data_fft.begin());
      const Rcomplex *q_ptr = reinterpret_cast<const Rcomplex *>(q_fft.begin());
      Rcomplex *p_ptr = reinterpret_cast<Rcomplex *>(prod.begin());

      RcppThread::parallelFor(0, static_cast<std::size_t>(fft_len), [&](std::size_t i) {
        p_ptr[i].r = d_ptr[i].r * q_ptr[i].r - d_ptr[i].i * q_ptr[i].i;
        p_ptr[i].i = d_ptr[i].r * q_ptr[i].i + d_ptr[i].i * q_ptr[i].r;
      }, actual_workers);
    } else {
      prod = data_fft * q_fft;
    }

    NumericVector z = Re(fft_rcpp(prod, true));

    double *dp = distance_profile.begin();
    double *lp = last_product.begin();
    const double *z_ptr = z.begin();
    const double *sx2_ptr = sumx2.begin();

    if (actual_workers > 1 && profile_len >= 1024) {
      RcppThread::parallelFor(0, static_cast<std::size_t>(profile_len), [&](std::size_t i) {
        double prod_val = z_ptr[window_size - 1 + i];
        lp[i] = prod_val;
        double dist = sx2_ptr[i] - 2.0 * prod_val + sumy2;
        if (!NumericVector::is_na(dist) && !std::isnan(dist) && dist < 0.0) {
          dist = 0.0;
        }
        dp[i] = dist;
      }, actual_workers);
    } else {
      for (std::size_t i = 0; i < profile_len; ++i) {
        double prod_val = z_ptr[window_size - 1 + i];
        lp[i] = prod_val;
        double dist = sx2_ptr[i] - 2.0 * prod_val + sumy2;
        if (!NumericVector::is_na(dist) && !std::isnan(dist) && dist < 0.0) {
          dist = 0.0;
        }
        dp[i] = dist;
      }
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (const std::exception &ex) {
    Rcpp::stop(ex.what());
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product);
}

//[[Rcpp::export]]
List mass2_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint64_t data_size,
                uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd, double query_mean,
                double query_sd, uint32_t n_workers = 1) {
  const uint64_t profile_len = data_size - window_size + 1;
  NumericVector distance_profile(profile_len);
  NumericVector last_product(profile_len);
  double q_mean = query_mean;
  double q_std = query_sd;

  try {
    const R_xlen_t fft_len = data_fft.length();
    NumericVector rev_query(fft_len); // use data_fft because is padded to a power of 2
    std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());

    ComplexVector const q_fft = fft_rcpp(rev_query);
    ComplexVector prod(fft_len);

    uint32_t actual_workers = std::max<uint32_t>(1, n_workers);

    if (actual_workers > 1 && fft_len >= 1024) {
      const Rcomplex *d_ptr = reinterpret_cast<const Rcomplex *>(data_fft.begin());
      const Rcomplex *q_ptr = reinterpret_cast<const Rcomplex *>(q_fft.begin());
      Rcomplex *p_ptr = reinterpret_cast<Rcomplex *>(prod.begin());

      RcppThread::parallelFor(0, static_cast<std::size_t>(fft_len), [&](std::size_t i) {
        p_ptr[i].r = d_ptr[i].r * q_ptr[i].r - d_ptr[i].i * q_ptr[i].i;
        p_ptr[i].i = d_ptr[i].r * q_ptr[i].i + d_ptr[i].i * q_ptr[i].r;
      }, actual_workers);
    } else {
      prod = data_fft * q_fft;
    }

    NumericVector z = Re(fft_rcpp(prod, true));

    double *dp = distance_profile.begin();
    double *lp = last_product.begin();
    const double *z_ptr = z.begin();
    const double *d_mean = data_mean.begin();
    const double *d_std = data_sd.begin();

    if (actual_workers > 1 && profile_len >= 1024) {
      RcppThread::parallelFor(0, static_cast<std::size_t>(profile_len), [&](std::size_t i) {
        double z_val = z_ptr[window_size - 1 + i];
        lp[i] = z_val;
        double dist = 2.0 * (static_cast<double>(window_size) -
                             (z_val - static_cast<double>(window_size) * d_mean[i] * q_mean) /
                             (d_std[i] * q_std));
        if (!NumericVector::is_na(dist) && !std::isnan(dist) && dist < 0.0) {
          dist = 0.0;
        }
        dp[i] = dist;
      }, actual_workers);
    } else {
      for (std::size_t i = 0; i < profile_len; ++i) {
        double z_val = z_ptr[window_size - 1 + i];
        lp[i] = z_val;
        double dist = 2.0 * (static_cast<double>(window_size) -
                             (z_val - static_cast<double>(window_size) * d_mean[i] * q_mean) /
                             (d_std[i] * q_std));
        if (!NumericVector::is_na(dist) && !std::isnan(dist) && dist < 0.0) {
          dist = 0.0;
        }
        dp[i] = dist;
      }
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (const std::exception &ex) {
    Rcpp::stop(ex.what());
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product);
}

//[[Rcpp::export]]
List mass3_rcpp(const NumericVector &query_window, const NumericVector &data_ref, uint64_t data_size,
                uint32_t window_size, const NumericVector &data_mean, const NumericVector &data_sd, double query_mean,
                double query_sd, uint32_t grain) {
  // data_ref is the long time series
  // query_window is the query
  // grain is the size of pieces, preferably a power of two
  // data_ref is the data, query_window is the query
  //

  // FIXME: mass3 is throwing error when data std is zero; currently a hack is being used on find_k_neighbors.R

  uint32_t const w_size = window_size;
  uint64_t const d_size = data_size;
  NumericVector dist(data_mean.length());
  NumericVector::iterator dist_it = dist.begin();
  NumericVector last(data_mean.length());
  NumericVector::iterator last_it = last.begin();

  grain = set_k_rcpp(grain, d_size, w_size);
  double q_mean = 0.0, q_std = 0.0;

  // compute query_window stats -- O(d_size)
  q_mean = query_mean;
  q_std = query_sd;

  // compute data stats -- O(d_size)
  NumericVector d_mean = data_mean;
  NumericVector d_std = data_sd;

  NumericVector rev_query(grain);
  std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());
  ComplexVector Y = fft_rcpp(rev_query);

  uint64_t j = 0;
  uint64_t jump = grain - w_size + 1;
  uint64_t const seq_end = d_size - grain;
  ComplexVector Z;
  NumericVector z;
  NumericVector d;

  try {
    for (j = 0; j <= seq_end; j = j + jump) {
      // The main trick of getting dot products in O(d_size log d_size) time
      uint64_t const idx_begin = j;
      uint64_t const idx_end = j + grain - w_size;                         // idx_begin + jump - 1
      ComplexVector const X = fft_rcpp(data_ref[Range(j, j + grain - 1)]); // idx_begin:(idx_begin + jump + w_size - 2)

      Z = X * Y;

      z = Re(fft_rcpp(Z, true));

      d = 2 * (w_size - (z[Range(w_size - 1, grain - 1)] - w_size * d_mean[Range(idx_begin, idx_end)] * q_mean) /
                            (d_std[Range(idx_begin, idx_end)] * q_std));

      // noise correction
      // d = d - (2 + 2 * w_size) * pow(std_n, 2) / pow(max(std_x, std_y), 2);

      std::copy(d.begin(), d.end(), dist_it + j);
      std::copy(z.begin() + w_size - 1, z.begin() + grain, last_it + j);
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (const std::exception &ex) {
    Rcout << "C++ std::exception in mass3 loop 1: " << ex.what() << std::endl;
    Rcpp::stop(ex.what());
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  jump = d_size - j;

  try {
    if (jump >= w_size) {
      uint64_t const idx_begin = j;
      uint64_t const idx_end = d_size - w_size;

      if ((jump - (w_size - 1) + j) > (uint64_t)data_mean.length()) {
        Rcout << "DEBUG: error." << std::endl;
      } else {
        ComplexVector const X = fft_rcpp(data_ref[Range(j, d_size - 1)]);
        Y = fft_rcpp(rev_query[Range(0, jump - 1)]);

        Z = X * Y;
        z = Re(fft_rcpp(Z, true));

        d = 2 * (w_size - (z[Range(w_size - 1, jump - 1)] - w_size * d_mean[Range(idx_begin, idx_end)] * q_mean) /
                              (d_std[Range(idx_begin, idx_end)] * q_std));

        std::copy(d.begin(), d.end(), dist_it + j);
        std::copy(z.begin() + w_size - 1, z.begin() + jump, last_it + j);
      }
    }

    clamp_negative(dist);
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (const std::exception &ex) {
    Rcout << "C++ std::exception in mass3 loop 2: " << ex.what() << std::endl;
    Rcpp::stop(ex.what());
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("distance_profile") = dist, Rcpp::Named("last_product") = last));
}

//[[Rcpp::export]]
List mass3_rcpp_parallel(const NumericVector query_window, const NumericVector data_ref, uint64_t data_size,
                         uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd,
                         double query_mean, double query_sd, uint16_t grain, uint32_t n_workers) {
  try {
    grain = set_k_rcpp(grain, data_size, window_size);

    // allocate the output matrix
    uint64_t const profile_len = data_mean.length();
    NumericVector distance_profile(profile_len);
    NumericVector last_product(profile_len);

    double *dp = distance_profile.begin();
    double *lp = last_product.begin();

    const double *d_mean = data_mean.begin();
    const double *d_std = data_sd.begin();
    const double *data_ptr = data_ref.begin();
    const double *win_ptr = query_window.begin();

    uint64_t const w_size = window_size;
    uint64_t const d_size = data_size;
    uint64_t const p_size = profile_len;

    std::vector<double> rev_query(grain, 0.0);
    std::reverse_copy(win_ptr, win_ptr + w_size, rev_query.begin());
    std::vector<std::complex<double>> const Y_full = fft_rcpp(rev_query, false);

    uint64_t const jump = grain - w_size + 1;
    uint64_t const seq_end = (d_size >= grain) ? (d_size - grain) : 0;

    struct Mass3Task {
      uint64_t j;
      uint64_t data_len;
      uint64_t n_pts;
      bool is_tail;
    };

    std::vector<Mass3Task> chunks;
    uint64_t j = 0;
    if (d_size >= grain) {
      for (j = 0; j <= seq_end; j += jump) {
        chunks.push_back({j, grain, jump, false});
      }
    }

    uint64_t const tail_jump = d_size - j;
    if (tail_jump >= w_size && (tail_jump - (w_size - 1) + j) <= p_size) {
      chunks.push_back({j, tail_jump, tail_jump - w_size + 1, true});
    }

    uint32_t actual_workers = std::max<uint32_t>(1, n_workers);
    actual_workers = std::min<uint32_t>(actual_workers, static_cast<uint32_t>(chunks.size()));

    try {
      RcppThread::parallelFor(0, chunks.size(), [&](std::size_t chunk_idx) {
        const Mass3Task &task = chunks[chunk_idx];
        uint64_t const c_j = task.j;
        uint64_t const c_data_len = task.data_len;
        uint64_t const c_n_pts = task.n_pts;

        std::vector<double> data_chunk(data_ptr + c_j, data_ptr + c_j + c_data_len);
        std::vector<std::complex<double>> X = fft_rcpp(data_chunk, false);

        std::vector<std::complex<double>> Z(c_data_len);
        if (!task.is_tail) {
          std::transform(X.begin(), X.end(), Y_full.begin(), Z.begin(), std::multiplies<std::complex<double>>());
        } else {
          std::vector<double> const rev_query_tail(rev_query.begin(), rev_query.begin() + static_cast<int64_t>(c_data_len));
          std::vector<std::complex<double>> const Y_tail = fft_rcpp(rev_query_tail, false);
          std::transform(X.begin(), X.end(), Y_tail.begin(), Z.begin(), std::multiplies<std::complex<double>>());
        }

        std::vector<double> z = fft_rcpp_real(Z, true);

        for (uint64_t i = 0; i < c_n_pts; ++i) {
          dp[c_j + i] =
              2.0 * (w_size - (z[w_size - 1 + i] - w_size * d_mean[c_j + i] * query_mean) /
                                       (d_std[c_j + i] * query_sd));
          lp[c_j + i] = z[w_size - 1 + i];
        }
      }, actual_workers);

      clamp_negative(distance_profile);

    } catch (RcppThread::UserInterruptException &ex) {
      Rcout << "Process terminated.\n";
    } catch (...) {
      Rcpp::stop("c++ exception (unknown reason)");
    }

    return (
        List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product));
  } catch (std::out_of_range &ex) {
    stop("out_of_range\n");
  } catch (...) {
    stop("Interrupted\n");
  }
}

//[[Rcpp::export]]
uint32_t set_k_rcpp(uint32_t grain, uint64_t data_size, uint64_t window_size) {
  try {
    grain = set_k_cpp(grain, data_size, window_size);
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return (grain);
}

//[[Rcpp::export]]
uint32_t find_best_k_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size) {
  uint64_t const data_size = data_ref.length();
  uint32_t grain = set_k_rcpp(window_size, data_size, window_size); // Set baseline
  uint64_t best_time = pow((double)2, (double)50); // added double inside sqrt to avoid ambiguity on Solaris
  uint32_t best_k = grain;
  List pre = mass_pre_rcpp(data_ref, query_ref, window_size);
  Timer const timer;

  try {
    for (uint16_t j = 0; j < 10; j++) {
      uint64_t const tictoc = timer.now();
      for (uint16_t i = 0; i < 10; i++) {
        List const nn = mass3_rcpp(
            query_ref[Range(i, i + window_size - 1)], data_ref, pre["data_size"], pre["window_size"], pre["data_mean"],
            pre["data_sd"], as<NumericVector>(pre["query_mean"])[i], as<NumericVector>(pre["query_sd"])[i], grain);
      }
      uint64_t const time_res = timer.now() - tictoc;
      if (time_res < best_time) {
        best_time = time_res;
        best_k = grain;
        grain = grain * 2;
        if (grain > data_size) {
          break;
        }
      } else {
        break;
      }
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return best_k;
}

// #### MASS PRE functions ####

//[[Rcpp::export]]
List mass_pre_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size) {
  uint64_t const data_size = data_ref.length();
  uint64_t const query_size = query_ref.length();

  List data_avgsd = movmean_std_rcpp(data_ref, window_size); // precompute moving average and SD

  uint32_t const pad_size = pow(2, ceil(log2(data_size))); // padded to a power of 2
  NumericVector data_padded(pad_size);

  std::copy(data_ref.begin(), data_ref.end(), data_padded.begin());

  ComplexVector const data_fft = fft_rcpp(data_padded); // precompute fft of data

  NumericVector query_mean;
  NumericVector query_sd;

  try {
    if (query_size > 0) {
      List query_avgsd = movmean_std_rcpp(query_ref, window_size); // precompute moving average and SD
      query_mean = query_avgsd["avg"];
      query_sd = query_avgsd["sd"];
    } else {
      query_mean = data_avgsd["avg"];
      query_sd = data_avgsd["sd"];
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("data_fft") = data_fft, Rcpp::Named("data_size") = data_size,
                       Rcpp::Named("window_size") = window_size, Rcpp::Named("data_mean") = data_avgsd["avg"],
                       Rcpp::Named("data_sd") = data_avgsd["sd"], Rcpp::Named("query_mean") = query_mean,
                       Rcpp::Named("query_sd") = query_sd));
}

//[[Rcpp::export]]
List mass_pre_abs_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size) {
  uint64_t const data_size = data_ref.length();
  uint64_t const query_size = query_ref.length();

  uint32_t const pad_size = pow(2, ceil(log2(data_size))); // padded to a power of 2
  NumericVector data_padded(pad_size);

  std::copy(data_ref.begin(), data_ref.end(), data_padded.begin());

  ComplexVector const data_fft = fft_rcpp(data_padded); // precompute fft of data
  NumericVector const sumx2 = movsum_ogita_rcpp(data_ref * data_ref, window_size);
  NumericVector sumy2;
  try {
    if (query_size > 0) {
      sumy2 = movsum_ogita_rcpp(query_ref * query_ref, window_size);
    } else {
      sumy2 = sumx2;
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("data_fft") = data_fft, Rcpp::Named("window_size") = window_size,
                       Rcpp::Named("data_size") = data_size, Rcpp::Named("sumx2") = sumx2,
                       Rcpp::Named("sumy2") = sumy2));
}

//[[Rcpp::export]]
List mass_pre_weighted_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
                            const NumericVector weight) {
  uint64_t const data_size = data_ref.length();
  uint64_t const query_size = query_ref.length();

  List data_avgsd = movmean_std_rcpp(data_ref, window_size); // precompute moving average and SD

  uint32_t const pad_size = pow(2, ceil(log2(data_size))); // padded to a power of 2
  NumericVector data_padded(pad_size);
  NumericVector rev_weight(pad_size);
  double const sumw = sum(weight);

  std::reverse_copy(weight.begin(), weight.end(), rev_weight.begin());
  std::copy(data_ref.begin(), data_ref.end(), data_padded.begin());

  ComplexVector const data_fft = fft_rcpp(data_padded); // precompute fft of data
  ComplexVector const w_fft = fft_rcpp(rev_weight);
  ComplexVector const data_fft_w_fft = data_fft * w_fft;

  NumericVector query_mean;
  NumericVector query_sd;

  try {
    if (query_size > 0) {
      List query_avgsd = movmean_std_rcpp(query_ref, window_size); // precompute moving average and SD
      query_mean = query_avgsd["avg"];
      query_sd = query_avgsd["sd"];
    } else {
      query_mean = data_avgsd["avg"];
      query_sd = data_avgsd["sd"];
    }
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  IntegerVector const range_s = Range(window_size - 1, data_size - 1);

  NumericVector data_w = Re(fft_rcpp(data_fft_w_fft, true));
  ComplexVector const data2_fft = fft_rcpp(data_padded * data_padded);
  ComplexVector const data2_fft_w_fft = data2_fft * w_fft;
  NumericVector data2_w = Re(fft_rcpp(data2_fft_w_fft, true));
  NumericVector const sumxw2 = data2_w[range_s];
  NumericVector const sumxw = data_w[range_s];

  NumericVector const data_mean = data_avgsd["avg"];
  NumericVector const data_sd = data_avgsd["sd"];

  NumericVector const data_pre =
      (sumxw2 - 2 * sumxw * data_mean + sumw * (data_mean * data_mean)) / (data_sd * data_sd);

  return (List::create(
      Rcpp::Named("data_fft") = data_fft, Rcpp::Named("data_pre") = data_pre, Rcpp::Named("data_size") = data_size,
      Rcpp::Named("window_size") = window_size, Rcpp::Named("data_mean") = data_mean, Rcpp::Named("data_sd") = data_sd,
      Rcpp::Named("query_mean") = query_mean, Rcpp::Named("query_sd") = query_sd, Rcpp::Named("weight") = weight));
}
