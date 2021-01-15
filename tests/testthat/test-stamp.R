if (!testthat:::on_cran()) {
  context("Testing STAMP")
  library(matrixprofiler)

  test_that("Stamp Rcpp", {
    expect_silent(matrixprofiler:::stamp_rcpp(
      data_ref = motifs_discords_small,
      query_ref = motifs_discords_small,
      window_size = 200,
      ez = 0.5, progress = FALSE
    ))
  })

  test_that("Stamp Parallel Rcpp", {
    expect_silent(matrixprofiler:::stamp_rcpp_parallel(
      data_ref = motifs_discords_small,
      query_ref = motifs_discords_small,
      window_size = 200,
      ez = 0.5, progress = FALSE
    ))
  })
}
