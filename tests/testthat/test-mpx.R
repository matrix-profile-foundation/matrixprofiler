if (!testthat:::on_cran()) {
  context("Testing MPX")
  library(matrixprofiler)

  test_that("it works", {
    expect_silent(mpx(motifs_discords_small, 200))
  })
}
