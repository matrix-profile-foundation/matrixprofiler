#include "stamp.h"
#include "mass.h"
#include "math.h"

// [[Rcpp::export]]
List stamp_cpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
               double ez) {

  double exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();
  uint64_t matrix_profile_size = data_size - window_size + 1;
  uint64_t num_queries = query_size - window_size + 1;


// check skip position
  LogicalVector skip_location(matrix_profile_size);

  for (uint64_t i = 0; i < matrix_profile_size; i++) {
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
  List pre = mass_pre_rcpp(data, query, window_size);

  IntegerVector order = Range(0, num_queries - 1);
  order = sample(order, num_queries);

  uint32_t k = find_best_k_rcpp(data, query, window_size);

  try {
    for (int32_t i : order) {
      List nn = mass3_rcpp(query[Range(i, i + window_size - 1)], data, pre["data_size"], pre["window_size"],
                           pre["data_mean"], pre["data_sd"], as<NumericVector>(pre["query_mean"])[i],
                           as<NumericVector>(pre["query_sd"])[i], k);

      NumericVector distance_profile = sqrt(as<NumericVector>(nn["distance_profile"]));

      // apply exclusion zone
      if (exclusion_zone > 0) {
        uint64_t exc_st = MAX(0, i - exclusion_zone);
        uint64_t exc_ed = MIN(matrix_profile_size - 1, i + exclusion_zone);
        IntegerVector dp_range = Range(exc_st, exc_ed);
        distance_profile[dp_range] = R_PosInf;
      }

      distance_profile[as<NumericVector>(pre["data_sd"]) < DBL_EPSILON] = R_PosInf;
      if (skip_location[i] || as<NumericVector>(pre["query_sd"])[i] < DBL_EPSILON) {
        distance_profile.fill(R_PosInf);
      }
      distance_profile[skip_location] = R_PosInf;

      // normal matrix_profile
      LogicalVector idx = (distance_profile < matrix_profile);
      matrix_profile[idx] = distance_profile[idx];
      profile_index[which(idx)] = i + 1;
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(
            Rcpp::Named("matrix_profile") = matrix_profile,
            Rcpp::Named("profile_index") = profile_index
          ));
}

// [[Rcpp::export]]
List stamp_par_cpp(const NumericVector data_ref, const NumericVector query_ref, uint32_t window_size,
               double ez) {



  double exclusion_zone = round(window_size * ez + DBL_EPSILON);
  uint64_t data_size = data_ref.length();
  uint64_t query_size = query_ref.length();
  uint64_t matrix_profile_size = data_size - window_size + 1;
  uint64_t num_queries = query_size - window_size + 1;


  // check skip position
  LogicalVector skip_location(matrix_profile_size);

  for (uint64_t i = 0; i < matrix_profile_size; i++) {
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
  List pre = mass_pre_rcpp(data, query, window_size);

  IntegerVector order = Range(0, num_queries - 1);
  order = sample(order, num_queries);

  uint32_t k = find_best_k_rcpp(data, query, window_size);

  try {
    for (int32_t i : order) {
      List nn = mass3_rcpp(query[Range(i, i + window_size - 1)], data, pre["data_size"], pre["window_size"],
                           pre["data_mean"], pre["data_sd"], as<NumericVector>(pre["query_mean"])[i],
                                                                                                 as<NumericVector>(pre["query_sd"])[i], k);

      NumericVector distance_profile = sqrt(as<NumericVector>(nn["distance_profile"]));

      // apply exclusion zone
      if (exclusion_zone > 0) {
        uint64_t exc_st = MAX(0, i - exclusion_zone);
        uint64_t exc_ed = MIN(matrix_profile_size - 1, i + exclusion_zone);
        IntegerVector dp_range = Range(exc_st, exc_ed);
        distance_profile[dp_range] = R_PosInf;
      }

      distance_profile[as<NumericVector>(pre["data_sd"]) < DBL_EPSILON] = R_PosInf;
      if (skip_location[i] || as<NumericVector>(pre["query_sd"])[i] < DBL_EPSILON) {
        distance_profile.fill(R_PosInf);
      }
      distance_profile[skip_location] = R_PosInf;

      // normal matrix_profile
      LogicalVector idx = (distance_profile < matrix_profile);
      matrix_profile[idx] = distance_profile[idx];
      profile_index[which(idx)] = i + 1;
    }
  } catch (Rcpp::internal::InterruptedException &ex) {
    Rcout << "Process terminated." << std::endl;
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }

  return (List::create(
      Rcpp::Named("matrix_profile") = matrix_profile,
      Rcpp::Named("profile_index") = profile_index
  ));
}
