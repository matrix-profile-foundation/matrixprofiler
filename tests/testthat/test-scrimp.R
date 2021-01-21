if (!testthat:::on_cran()) {
  test_that("Scrimp", {
    expect_silent(scrimp(
      data = motifs_discords_small, window_size = 200, exclusion_zone = 0.5,
      progress = FALSE
    ))
  })

  test_that("Scrimp Parallel", {
    expect_silent(scrimp(
      data = motifs_discords_small, window_size = 200, exclusion_zone = 0.5,
      n_workers = 2, progress = FALSE
    ))
  })
}
