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
                        const bool normalized) {
  NumericVector distance_profile;
  NumericVector last_product;

  NumericVector query;
  NumericVector rev_query(data_fft.length());
  NumericVector rev_weight(data_fft.length());

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

  ComplexVector const prod = data_fft * fft_rcpp(rev_weight * rev_query);
  NumericVector data_queryw = Re(fft_rcpp(prod, true));

  IntegerVector const range_d = Range(window_size - 1, data_size - 1);

  last_product = data_queryw[range_d];

  distance_profile = data_pre - 2 * (last_product - sumwy * data_mean) / data_sd + sumwy2;
  clamp_negative(distance_profile);

  return (List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product));
}

//[[Rcpp::export]]
List mass_absolute_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint32_t data_size,
                        uint32_t window_size, const NumericVector sumx2, double sumy2) {
  NumericVector distance_profile;
  NumericVector last_product;

  // pre-process query for fft
  try {
    NumericVector rev_query(data_fft.length()); // use data_fft because is padded to a power of 2
    std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());

    // compute the product
    ComplexVector const prod = data_fft * fft_rcpp(rev_query);
    NumericVector z = Re(fft_rcpp(prod, true));
    // compute the distance profile
    IntegerVector const range_z = Range(window_size - 1, data_size - 1);
    IntegerVector const range_d = Range(0, data_size - window_size);
    last_product = z[range_z];
    distance_profile = as<NumericVector>(sumx2[range_d]) - 2 * last_product + sumy2;
    clamp_negative(distance_profile);
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product));
}

//[[Rcpp::export]]
List mass2_rcpp(const ComplexVector data_fft, const NumericVector query_window, uint64_t data_size,
                uint32_t window_size, const NumericVector data_mean, const NumericVector data_sd, double query_mean,
                double query_sd) {
  NumericVector distance_profile;
  NumericVector last_product;
  double q_mean = 0.0, q_std = 0.0;

  // compute query_window stats -- O(d_size)
  q_mean = query_mean;
  q_std = query_sd;

  // pre-process query for fft
  try {
    NumericVector rev_query(data_fft.length()); // use data_fft because is padded to a power of 2
    std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());

    // compute the product
    ComplexVector const prod = data_fft * fft_rcpp(rev_query);
    NumericVector z = Re(fft_rcpp(prod, true));
    // compute the distance profile
    last_product = z[Range(window_size - 1, data_size - 1)];

    distance_profile = 2 * (window_size - (last_product - window_size * data_mean * q_mean) / (data_sd * q_std));
    clamp_negative(distance_profile);
  } catch (RcppThread::UserInterruptException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    Rcpp::stop("c++ exception (unknown reason)");
  }

  return (List::create(Rcpp::Named("distance_profile") = distance_profile, Rcpp::Named("last_product") = last_product));
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

    std::vector<std::pair<std::size_t, std::size_t>> chunks;
    for (std::size_t cur = 0; cur < data_size; cur += grain) {
      std::size_t end = std::min(cur + grain, static_cast<std::size_t>(data_size));
      chunks.emplace_back(cur, end);
    }

    uint32_t actual_workers = std::max<uint32_t>(1, n_workers);
    actual_workers = std::min<uint32_t>(actual_workers, static_cast<uint32_t>(chunks.size()));

    try {
      RcppThread::parallelFor(0, chunks.size(), [&](std::size_t chunk_idx) {
        std::size_t begin = chunks[chunk_idx].first;
        std::size_t end = chunks[chunk_idx].second;
        if (begin >= data_size - window_size + 1) return;

        uint64_t jump = end - begin;
        uint64_t k = jump + window_size - 1;
        if (end > data_size - window_size) {
          jump = (data_size - window_size + 1) - begin;
          k = end - begin;
        }
        if (jump == 0) return;

        uint64_t const s = static_cast<uint64_t>(std::pow(2.0, std::ceil(std::log2(static_cast<double>(k)))));
        std::unique_ptr<FFT::fftw> fft(new FFT::fftw());

        std::vector<std::complex<double>> rev_query(s);
        for (uint64_t i = window_size, j = 0; i > 0; --i, ++j) {
          rev_query[i - 1] = std::complex<double>(win_ptr[j], 0.0);
        }
        std::vector<std::complex<double>> Y = fft->fft(rev_query, false);

        std::vector<std::complex<double>> data(s);
        for (uint64_t i = 0; i < k; ++i) {
          data[i] = std::complex<double>(data_ptr[begin + i], 0.0);
        }
        std::vector<std::complex<double>> X = fft->fft(data, false);

        std::vector<std::complex<double>> Z(s);
        std::transform(X.begin(), X.end(), Y.begin(), Z.begin(), std::multiplies<std::complex<double>>());
        std::vector<std::complex<double>> z = fft->fft(Z, true);

        for (uint64_t i = 0; i < jump; ++i) {
          dp[begin + i] =
              2.0 * (window_size - (z[k - jump + i].real() - window_size * d_mean[begin + i] * query_mean) /
                                       (d_std[begin + i] * query_sd));
          lp[begin + i] = z[k - jump + i].real();
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
