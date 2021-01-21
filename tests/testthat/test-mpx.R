if (!testthat:::on_cran()) {
  test_that("MPX", {
    expect_silent(mpx(
      data = motifs_discords_small,
      window_size = 200, exclusion_zone = 0.5,
      idxs = TRUE, distance = "euclidean", progress = FALSE
    ))
  })

  test_that("MPX Parallel", {
    expect_silent(mpx(
      data = motifs_discords_small,
      window_size = 200, exclusion_zone = 0.5,
      idxs = TRUE, distance = "euclidean", n_workers = 2, progress = FALSE
    ))
  })

  test_that("MPXAB", {
    expect_silent(mpx(
      data = motifs_discords_small, query = rev(motifs_discords_small),
      window_size = 200, exclusion_zone = 0.5,
      idxs = TRUE, distance = "euclidean", progress = FALSE
    ))
  })

  test_that("MPXAB Parallel", {
    expect_silent(mpx(
      data = motifs_discords_small, query = rev(motifs_discords_small),
      window_size = 200, exclusion_zone = 0.5,
      idxs = TRUE, distance = "euclidean", n_workers = 2, progress = FALSE
    ))
  })
}
