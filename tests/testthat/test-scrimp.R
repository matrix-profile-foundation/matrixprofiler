scrimp_res <- NULL
scrimp_res_par <- NA

test_that("Scrimp", {
  set.seed(2021)
  expect_silent(scrimp_res <<- scrimp(
    data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
    progress = FALSE
  ))
  expect_type(scrimp_res, "list")
  expect_snapshot_value(scrimp_res, style = "serialize")
})

test_that("Scrimp Parallel", {
  set.seed(2021)
  expect_silent(scrimp_res_par <<- scrimp(
    data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
    n_workers = 2, progress = FALSE
  ))
  expect_type(scrimp_res_par, "list")
  expect_snapshot_value(scrimp_res_par, style = "serialize")
})

test_that("Scrimps are equal", {
  expect_equal(scrimp_res, scrimp_res_par)
})
