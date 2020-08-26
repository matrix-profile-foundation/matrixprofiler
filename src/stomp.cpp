#include "stomp.h"
#include "mass.h"
#include "math.h"

// [[Rcpp::export]]
List stomp_cpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size, double ez) {

  double exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint32_t data_size = data_ref.length();
  uint32_t query_size = query_ref.length();
  uint32_t matrix_profile_size = data_size - window_size + 1;
  uint32_t num_queries = query_size - window_size + 1;

// check skip position
  LogicalVector skip_location(matrix_profile_size);

  for (uint32_t i = 0; i < matrix_profile_size; i++) {
    NumericVector range = data_ref[Range(i, (i + window_size - 1))];
    if (any(is_na(range) | is_infinite(range))) {
      skip_location[i] = true;
    }
  }

  NumericVector data = data_ref;
  NumericVector query = query_ref;

  data[is_na(data)] = 0;
  data[is_infinite(data)] = 0;
  query[is_na(query)] = 0;
  query[is_infinite(query)] = 0;

  NumericVector matrix_profile(matrix_profile_size, R_PosInf);
  IntegerVector profile_index(matrix_profile_size, R_NegInf);
  Environment stats = Environment::namespace_env("stats");
  Function fft = stats["fft"];

  List pre = mass_pre_rcpp(data, query, window_size, fft);
  List rpre = mass_pre_rcpp(query, data, window_size, fft);
  List nn = mass3_rcpp(query[Range(0, window_size - 1)], data, rpre["window_size"],
                        pre["data_size"], pre["data_mean"], pre["data_sd"], as<NumericVector>(pre["query_mean"])[0],
                        as<NumericVector>(pre["query_sd"])[0], fft);


  ///// This is needed for JOIN similarity
  List rnn = mass3_rcpp(data[Range(0, window_size - 1)], query, rpre["window_size"],
                        query_size, rpre["data_mean"], rpre["data_sd"], as<NumericVector>(rpre["query_mean"])[0],
                        as<NumericVector>(rpre["query_sd"])[0], fft);

  NumericVector first_product = rnn["last_product"];


  ///// This is needed for JOIN similarity
  IntegerVector lp_range = Range(1, (data_size - window_size));
  IntegerVector lp2_range = Range(0, (data_size - window_size - 1));
  IntegerVector dt_range = Range(window_size, data_size - 1);
  NumericVector distance_profile;
  NumericVector last_product;
  double drop_value = query[0];

  IntegerVector order = Range(0, num_queries - 1);

  for (int32_t i : order) {
// compute the distance profile
    NumericVector query_window = query[Range(i, (i + window_size - 1))];
    if (i == 0) {
      distance_profile = nn["distance_profile"];
      last_product = nn["last_product"];
    } else {
      last_product[lp_range] = (NumericVector)(as<NumericVector>(last_product[lp2_range]) - as<NumericVector>(data[lp2_range]) * drop_value +
        as<NumericVector>(data[dt_range]) * query_window[window_size - 1]);
      last_product[0] = first_product[i];
      distance_profile = 2 * (window_size - (last_product - window_size * as<NumericVector>(pre["data_mean"]) * as<NumericVector>(pre["query_mean"])[i]) /
                              (as<NumericVector>(pre["data_sd"]) * as<NumericVector>(pre["query_sd"])[i]));
    }
    distance_profile[distance_profile < 0] = 0;
    // distance_profile = sqrt(distance_profile);
    drop_value = query_window[0];

// apply exclusion zone
    if (exclusion_zone > 0) {
      uint32_t exc_st = MAX(0, i - exclusion_zone);
      uint32_t exc_ed = MIN(matrix_profile_size - 1, i + exclusion_zone);
      IntegerVector dp_range = Range(exc_st, exc_ed);
      distance_profile[dp_range] = R_PosInf;
    }

    distance_profile[as<NumericVector>(pre["data_sd"]) < DBL_EPSILON] = R_PosInf;
    if (skip_location[i] || as<NumericVector>(pre["query_sd"])[i] < DBL_EPSILON) {
      distance_profile.fill(R_PosInf);
    }
    distance_profile[skip_location] = R_PosInf;

    LogicalVector idx = (distance_profile < matrix_profile);
    matrix_profile[idx] = distance_profile[idx];
    profile_index[which(idx)] = i + 1;
  }

  return (List::create(
            Rcpp::Named("matrix_profile") = sqrt(matrix_profile),
            Rcpp::Named("profile_index") = profile_index
          ));
}
