stamp_res <- NULL
stamp_res_par <- NA

test_that("Stamp", {
  set.seed(2021)
  expect_silent(stamp_res <<- stamp(
    data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
    progress = FALSE
  ))
  expect_type(stamp_res, "list")
  expect_snapshot_value(stamp_res, style = "serialize")
})

test_that("Stamp Parallel", {
  set.seed(2021)
  expect_silent(stamp_res_par <<- stamp(
    data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
    n_workers = 2, progress = FALSE
  ))
  expect_type(stamp_res_par, "list")
  # expect_snapshot_value(stamp_res_par, style = "serialize")
})

test_that("Stamps are equal", {
  expect_equal(stamp_res, stamp_res_par)
})
