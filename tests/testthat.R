library(testthat)
library(matrixprofiler)

if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
}

test_check("matrixprofiler")
