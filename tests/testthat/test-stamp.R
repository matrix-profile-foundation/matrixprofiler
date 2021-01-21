if (!testthat:::on_cran()) {
  test_that("Stamp", {
    expect_silent(stamp(
      data = motifs_discords_small, window_size = 200, exclusion_zone = 0.5,
      progress = FALSE
    ))
  })

  test_that("Stamp Parallel", {
    expect_silent(stamp(
      data = motifs_discords_small, window_size = 200, exclusion_zone = 0.5,
      n_workers = 2, progress = FALSE
    ))
  })
}
