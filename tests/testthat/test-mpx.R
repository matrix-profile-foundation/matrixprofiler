if (!testthat:::on_cran()) {
  context("Testing MPX")
  library(matrixprofiler)

  test_that("MPX Rcpp", {
    expect_silent(matrixprofiler:::mpx_rcpp(motifs_discords_small,
                                             window_size = 200, ez = 0.5,
                                            idxs = TRUE,
                                            euclidean = TRUE,
                                            progress = FALSE))
  })

  test_that("MPX Rcpp Parallel", {
    expect_silent(matrixprofiler:::mpx_rcpp_parallel(motifs_discords_small,
                                                     window_size = 200, ez = 0.5,
                                                     idxs = TRUE,
                                                     euclidean = TRUE,
                                                     progress = FALSE))
  })

  test_that("MPXAB Rcpp", {
    expect_silent(matrixprofiler:::mpxab_rcpp(motifs_discords_small, rev(motifs_discords_small),
                                            window_size = 200,
                                            idxs = TRUE,
                                            euclidean = TRUE,
                                            progress = FALSE))
  })

  test_that("MPXAB Rcpp Parallel", {
    expect_silent(matrixprofiler:::mpxab_rcpp_parallel(motifs_discords_small, rev(motifs_discords_small),
                                              window_size = 200,
                                              idxs = TRUE,
                                              euclidean = TRUE,
                                              progress = FALSE))
  })

  test_that("MPX R", {
    expect_silent(mpx(data = motifs_discords_small, window_size = 200, progress = FALSE))
  })

  test_that("MPXAB R", {
    expect_silent(mpx(data = motifs_discords_small, window_size = 200,
                      query = rev(motifs_discords_small), progress = FALSE))
  })

  test_that("MPX Parallel R", {
    expect_silent(mpx(data = motifs_discords_small, window_size = 200, progress = FALSE, n_workers = 2))
  })

  test_that("MPXAB Parallel R", {
    expect_silent(mpx(data = motifs_discords_small, window_size = 200,
                      query = rev(motifs_discords_small), progress = FALSE, n_workers = 2))
  })
}
