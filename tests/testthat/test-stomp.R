if (!testthat:::on_cran()) {
  context("Testing STOMP")
  library(matrixprofiler)

  test_that("Stomp Rcpp", {
    expect_silent(matrixprofiler:::stomp_rcpp(
      data_ref = motifs_discords_small,
      query_ref = motifs_discords_small,
      window_size = 200,
      ez = 0.5, progress = FALSE
    ))
  })

  test_that("Stomp Parallel Rcpp", {
    expect_silent(matrixprofiler:::stomp_rcpp_parallel(
      data_ref = motifs_discords_small,
      query_ref = motifs_discords_small,
      window_size = 200,
      ez = 0.5, progress = FALSE
    ))
  })
}
