if (!testthat:::on_cran()) {
  mpx_res <- NULL
  mpx_res_par <- NA
  mpxab_res <- NULL
  mpxab_res_par <- NA

  test_that("MPX", {
    expect_silent(mpx_res <<- mpx(
      data = motifs_discords_small,
      window_size = 150, exclusion_zone = 0.5,
      idxs = TRUE, distance = "euclidean", progress = FALSE
    ))
    expect_type(mpx_res, "list")
    expect_snapshot_value(mpx_res, style = "serialize")
  })

  test_that("MPX Parallel", {
    expect_silent(mpx_res_par <<- mpx(
      data = motifs_discords_small,
      window_size = 150, exclusion_zone = 0.5,
      idxs = TRUE, distance = "euclidean", n_workers = 2, progress = FALSE
    ))
    expect_type(mpx_res_par, "list")
    expect_snapshot_value(mpx_res_par, style = "serialize")
  })

  test_that("MPXAB", {
    expect_silent(mpxab_res <<- mpx(
      data = motifs_discords_small, query = rev(motifs_discords_small),
      window_size = 150, exclusion_zone = 0.5,
      idxs = TRUE, distance = "euclidean", progress = FALSE
    ))
    expect_type(mpxab_res, "list")
    expect_snapshot_value(mpxab_res, style = "serialize")
  })

  test_that("MPXAB Parallel", {
    expect_silent(mpxab_res_par <<- mpx(
      data = motifs_discords_small, query = rev(motifs_discords_small),
      window_size = 150, exclusion_zone = 0.5,
      idxs = TRUE, distance = "euclidean", n_workers = 2, progress = FALSE
    ))
    expect_type(mpxab_res_par, "list")
    expect_snapshot_value(mpxab_res_par, style = "serialize")
  })

  test_that("MPXs are equal", {
    expect_identical(mpx_res, mpx_res_par)
  })

  test_that("MPXABs are equal", {
    expect_identical(mpxab_res, mpxab_res_par)
  })
}
