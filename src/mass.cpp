#include "mass.h"
#include "math.h"
#include "windowfunc.h"
#include <Rcpp/Benchmark/Timer.h>

//[[Rcpp::export]]
List mass_weighted_rcpp(ComplexVector data_fft, NumericVector query_window, uint32_t data_size,
                        uint32_t window_size, NumericVector data_mean, NumericVector data_sd,
                        double query_mean, double query_sd,
                        NumericVector data_pre, NumericVector weight, Function fft) {

  NumericVector distance_profile;
  NumericVector last_product;

  NumericVector query;
  NumericVector rev_query(data_fft.length());
  NumericVector rev_weight(data_fft.length());

  // normalized query
  query = (query_window - query_mean) / query_sd;

  double sumwy = sum(query * weight);
  double sumwy2 = sum(weight * pow(query, 2));

  std::reverse_copy(query.begin(), query.end(), rev_query.begin());
  std::reverse_copy(weight.begin(), weight.end(), rev_weight.begin());

  ComplexVector prod = data_fft * as<ComplexVector>(fft(rev_weight * rev_query));
  NumericVector data_queryw = Re(as<ComplexVector>(fft(prod, true))) / data_fft.length();

  IntegerVector range_d = Range(window_size - 1, data_size - 1);

  last_product = data_queryw[range_d];

  distance_profile = data_pre - 2 * (last_product - sumwy * data_mean) / data_sd + sumwy2;
  distance_profile[distance_profile < 0] = 0;

  return (List::create(
      Rcpp::Named("distance_profile") = distance_profile,
      Rcpp::Named("last_product") = last_product
  ));
}

//[[Rcpp::export]]
List mass_absolute_rcpp(ComplexVector data_fft, NumericVector query_window, uint32_t data_size,
                        uint32_t window_size, NumericVector sumx2, double sumy2, Function fft) {

  NumericVector distance_profile;
  NumericVector last_product;
  // pre-process query for fft
  try {
    NumericVector rev_query(data_fft.length()); // use data_fft because is padded to a power of 2
    std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());

    // compute the product
    ComplexVector prod = data_fft * as<ComplexVector>(fft(rev_query));
    NumericVector z = Re(as<ComplexVector>(fft(prod, true))) / prod.length();
    // compute the distance profile
    IntegerVector range_z = Range(window_size - 1, data_size - 1);
    IntegerVector range_d = Range(0, data_size - window_size);
    last_product = z[range_z];
    distance_profile = as<NumericVector>(sumx2[range_d]) - 2 * last_product + sumy2;
    distance_profile[distance_profile < 0] = 0;
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(
            Rcpp::Named("distance_profile") = distance_profile,
            Rcpp::Named("last_product") = last_product
          ));
}

//[[Rcpp::export]]
List mass2_rcpp(ComplexVector data_fft, NumericVector query_window, uint64_t data_size,
                uint32_t window_size, NumericVector data_mean, NumericVector data_sd,
                double query_mean, double query_sd, Function fft) {

  NumericVector distance_profile;
  NumericVector last_product;
  // pre-process query for fft
  try {
    NumericVector rev_query(data_fft.length()); // use data_fft because is padded to a power of 2
    std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());

    // compute the product
    ComplexVector prod = data_fft * as<ComplexVector>(fft(rev_query));
    NumericVector z = Re(as<ComplexVector>(fft(prod, true))) / prod.length();
    // compute the distance profile
    last_product = z[Range(window_size - 1, data_size - 1)];
    distance_profile = 2 * (window_size - (last_product - window_size * data_mean * query_mean) / (data_sd * query_sd));
    distance_profile[distance_profile < 0] = 0;
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(
            Rcpp::Named("distance_profile") = distance_profile,
            Rcpp::Named("last_product") = last_product
          ));
}

//[[Rcpp::export]]
List mass3_rcpp(NumericVector query_window, NumericVector data_ref,
                uint64_t data_size, uint32_t window_size, NumericVector data_mean,
                NumericVector data_sd, double query_mean, double query_sd, Function fft,
                uint32_t k) {

  // data_ref is the long time series
  // query_window is the query
  // k is the size of pieces, preferably a power of two
  // data_ref is the data, query_window is the query
  //

  uint32_t w_size = window_size;
  uint64_t d_size = data_size;
  NumericVector dist(data_mean.length());
  NumericVector::iterator dist_it = dist.begin();
  NumericVector last(data_mean.length());
  NumericVector::iterator last_it = last.begin();

  k = set_k(k, d_size, w_size);

  // compute query_window stats -- O(d_size)
  double q_mean = query_mean;
  double q_std = query_sd;

  // compute data stats -- O(d_size)
  NumericVector d_mean = data_mean;
  NumericVector d_std = data_sd;

  // k = 4096; // assume k > w_size
  // k = pow2(nextpow2(sqrt(d_size)));

  NumericVector rev_query(k);
  std::reverse_copy(query_window.begin(), query_window.end(), rev_query.begin());
  ComplexVector Y = fft(rev_query);

  uint64_t j = 0;
  uint64_t jump = k - w_size + 1;
  uint64_t seq_end = d_size - k;
  ComplexVector Z;
  NumericVector z;
  NumericVector d;

  try {
    for (j = 0; j <= seq_end; j = j + jump) {
      // The main trick of getting dot products in O(d_size log d_size) time
      uint64_t idx_begin = j;
      uint64_t idx_end = j + k - w_size;
      ComplexVector X = fft(data_ref[Range(j, j + k - 1)]);
      Z = X * Y;
      z = Re(as<ComplexVector>(fft(Z, true))) / Z.length();
      d = 2 * (w_size - (z[Range(w_size - 1, k - 1)] - w_size * d_mean[Range(idx_begin, idx_end)]
                         * q_mean) / (d_std[Range(idx_begin, idx_end)] * q_std));
      std::copy(d.begin(), d.end(), dist_it + j);
      std::copy(z.begin() + w_size - 1, z.begin() + k, last_it + j);
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  // j = j + jump;
  jump = d_size - j;

  try {
    if (jump >= w_size) {
      uint64_t idx_begin = j;
      uint64_t idx_end = d_size - w_size;

      if ((jump - (w_size - 1) + j) > (uint64_t)data_mean.length()) {
        Rcout << "DEBUG: error." << std::endl;
      } else {
        ComplexVector X = fft(data_ref[Range(j, d_size - 1)]);
        Y = fft(rev_query[Range(0, jump - 1)]);
        Z = X * Y;
        z = Re(as<ComplexVector>(fft(Z, true))) / Z.length();
        d = 2 * (w_size - (z[Range(w_size - 1, jump - 1)] - w_size * d_mean[Range(idx_begin, idx_end)] * q_mean) / (d_std[Range(idx_begin, idx_end)] * q_std));
        std::copy(d.begin(), d.end(), dist_it + j);
        std::copy(z.begin() + w_size - 1, z.begin() + jump, last_it + j);
      }
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(
            Rcpp::Named("distance_profile") = dist,
            Rcpp::Named("last_product") = last
          ));
}

//[[Rcpp::export]]
uint32_t set_k(uint32_t k, uint64_t data_size, uint64_t window_size) {

  try {
    if (k > data_size) {
      k = pow(2, ceil(log2(sqrt(data_size))));
    }

    if (k <= window_size) {
      k = pow(2, (ceil(log2(window_size)) + 1));
      if (k > data_size) {
        k = data_size;
      }
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (k);
}

//[[Rcpp::export]]
uint32_t find_best_k(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size) {
  uint64_t data_size = data_ref.length();
  uint32_t k = set_k(window_size, data_size, window_size); // Set baseline
  uint64_t best_time = pow(2, 50);
  uint32_t best_k = k;
  Environment stats = Environment::namespace_env("stats");
  Function fft = stats["fft"];
  List pre = mass_pre_rcpp(data_ref, query_ref, window_size);
  Timer timer;

  try {
    for (uint16_t j = 0; j < 10; j++) {
      uint64_t tictoc = timer.now();
      for (uint16_t i = 0; i < 10; i++) {
        List nn = mass3_rcpp(query_ref[Range(i, i + window_size - 1)], data_ref, pre["data_size"], pre["window_size"],
                             pre["data_mean"], pre["data_sd"], as<NumericVector>(pre["query_mean"])[i],
                             as<NumericVector>(pre["query_sd"])[i], fft, k);
      }
      uint64_t time_res = timer.now() - tictoc;
      if (time_res < best_time) {
        best_time = time_res;
        best_k = k;
        k = k * 2;
        if (k > data_size) {
          break;
        }
      } else {
        break;
      }
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return best_k;
}

//[[Rcpp::export]]
List mass_pre_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size) {

  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();

  List data_avgsd = fast_avg_sd_rcpp(data_ref, window_size); // precompute moving average and SD

  uint32_t pad_size = pow(2, ceil(log2(data_size))); // padded to a power of 2
  NumericVector data_padded(pad_size);

  std::copy(data_ref.begin(), data_ref.end(), data_padded.begin());

  Environment stats = Environment::namespace_env("stats");
  Function fft = stats["fft"];

  ComplexVector data_fft = as<ComplexVector>(fft(data_padded)); // precompute fft of data

  NumericVector query_mean;
  NumericVector query_sd;

  try {
    if (query_size > 0) {
      List query_avgsd = fast_avg_sd_rcpp(query_ref, window_size); // precompute moving average and SD
      query_mean = query_avgsd["avg"];
      query_sd = query_avgsd["sd"];
    } else {
      query_mean = data_avgsd["avg"];
      query_sd = data_avgsd["sd"];
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(
            Rcpp::Named("data_fft") = data_fft,
            Rcpp::Named("data_size") = data_size,
            Rcpp::Named("window_size") = window_size,
            Rcpp::Named("data_mean") = data_avgsd["avg"],
            Rcpp::Named("data_sd") = data_avgsd["sd"],
            Rcpp::Named("query_mean") = query_mean,
            Rcpp::Named("query_sd") = query_sd
          ));
}

//[[Rcpp::export]]
List mass_pre_abs_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size) {

  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();

  uint32_t pad_size = pow(2, ceil(log2(data_size))); // padded to a power of 2
  NumericVector data_padded(pad_size);

  std::copy(data_ref.begin(), data_ref.end(), data_padded.begin());

  Environment stats = Environment::namespace_env("stats");
  Function fft = stats["fft"];

  ComplexVector data_fft = as<ComplexVector>(fft(data_padded)); // precompute fft of data
  NumericVector sumx2 = movsum(pow(data_ref, 2), window_size);
  NumericVector sumy2;
  try {
    if (query_size > 0) {
      sumy2 = movsum(pow(query_ref, 2), window_size);
    } else {
      sumy2 = sumx2;
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(
            Rcpp::Named("data_fft") = data_fft,
            Rcpp::Named("window_size") = window_size,
            Rcpp::Named("data_size") = data_size,
            Rcpp::Named("sumx2") = sumx2,
            Rcpp::Named("sumy2") = sumy2
          ));
}

//[[Rcpp::export]]
List mass_pre_weighted_rcpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
                            const NumericVector weight) {

  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();

  List data_avgsd = fast_avg_sd_rcpp(data_ref, window_size); // precompute moving average and SD

  uint32_t pad_size = pow(2, ceil(log2(data_size))); // padded to a power of 2
  NumericVector data_padded(pad_size);
  NumericVector rev_weight(pad_size);
  double sumw = sum(weight);

  std::reverse_copy(weight.begin(), weight.end(), rev_weight.begin());
  std::copy(data_ref.begin(), data_ref.end(), data_padded.begin());

  Environment stats = Environment::namespace_env("stats");
  Function fft = stats["fft"];

  ComplexVector data_fft = as<ComplexVector>(fft(data_padded)); // precompute fft of data
  ComplexVector w_fft = as<ComplexVector>(fft(rev_weight));
  ComplexVector data_fft_w_fft = data_fft * w_fft;

  NumericVector query_mean;
  NumericVector query_sd;

  try {
    if (query_size > 0) {
      List query_avgsd = fast_avg_sd_rcpp(query_ref, window_size); // precompute moving average and SD
      query_mean = query_avgsd["avg"];
      query_sd = query_avgsd["sd"];
    } else {
      query_mean = data_avgsd["avg"];
      query_sd = data_avgsd["sd"];
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  IntegerVector range_s = Range(window_size - 1, data_size - 1);

  NumericVector data_w = Re(as<ComplexVector>(fft(data_fft_w_fft, true))) / data_fft_w_fft.length();
  ComplexVector data2_fft = fft(pow(data_padded, 2));
  ComplexVector data2_fft_w_fft = data2_fft * w_fft;
  NumericVector data2_w = Re(as<ComplexVector>(fft(data2_fft_w_fft, true))) / data2_fft_w_fft.length();
  NumericVector sumxw2 = data2_w[range_s];
  NumericVector sumxw = data_w[range_s];

  NumericVector data_mean = data_avgsd["avg"];
  NumericVector data_sd = data_avgsd["sd"];

  NumericVector data_pre = (sumxw2 - 2 * sumxw * data_mean + sumw * pow(data_mean, 2)) / pow(data_sd, 2);

  return (List::create(
      Rcpp::Named("data_fft") = data_fft,
      Rcpp::Named("data_pre") = data_pre,
      Rcpp::Named("data_size") = data_size,
      Rcpp::Named("window_size") = window_size,
      Rcpp::Named("data_mean") = data_mean,
      Rcpp::Named("data_sd") = data_sd,
      Rcpp::Named("query_mean") = query_mean,
      Rcpp::Named("query_sd") = query_sd,
      Rcpp::Named("weight") = weight
  ));
}

