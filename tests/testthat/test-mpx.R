mpx_res <- NULL
mpx_res_par <- NA
mpxab_res <- NULL
mpxab_res_par <- NA

test_that("MPX", {
  set.seed(2021)
  expect_silent(mpx_res <<- mpx(
    data = motifs_discords_small,
    window_size = 150, exclusion_zone = 0.5,
    idxs = TRUE, distance = "euclidean", progress = FALSE
  ))
  expect_type(mpx_res, "list")
  expect_snapshot_value(mpx_res, style = "serialize")
})

test_that("MPX Parallel", {
  set.seed(2021)
  expect_silent(mpx_res_par <<- mpx(
    data = motifs_discords_small,
    window_size = 150, exclusion_zone = 0.5,
    idxs = TRUE, distance = "euclidean", n_workers = 2, progress = FALSE
  ))
  expect_type(mpx_res_par, "list")
  expect_snapshot_value(mpx_res_par, style = "serialize")
})

test_that("MPXAB", {
  set.seed(2021)
  expect_silent(mpxab_res <<- mpx(
    data = motifs_discords_small, query = rev(motifs_discords_small),
    window_size = 150, exclusion_zone = 0.5,
    idxs = TRUE, distance = "euclidean", progress = FALSE
  ))
  expect_type(mpxab_res, "list")
  expect_snapshot_value(mpxab_res, style = "serialize")
})

test_that("MPXAB Parallel", {
  set.seed(2021)
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

test_that("parallel MPX handles fewer than 100 diagonals", {
  set.seed(2001)
  data <- rnorm(100)

  set.seed(2002)
  serial <- matrixprofiler:::mpx_rcpp(data, 40L, 0.5, 1, TRUE, TRUE, FALSE, 60)
  set.seed(2002)
  parallel <- matrixprofiler:::mpx_rcpp_parallel(data, 40L, 0.5, 1, TRUE, TRUE, FALSE)

  expect_equal(parallel, serial, tolerance = 1e-10)
})

test_that("parallel MPX honors partial sample size", {
  set.seed(2003)
  data <- rnorm(1000)

  set.seed(2004)
  serial <- matrixprofiler:::mpx_rcpp(data, 100L, 0.5, 0.25, TRUE, TRUE, FALSE, 60)
  set.seed(2004)
  parallel <- matrixprofiler:::mpx_rcpp_parallel(data, 100L, 0.5, 0.25, TRUE, TRUE, FALSE)

  expect_true(parallel$partial)
  expect_equal(parallel, serial, tolerance = 1e-10)
})

test_that("mpx validates sample size", {
  expect_error(mpx(rnorm(100), 40L, s_size = -0.1), "between 0 and 1", fixed = TRUE)
  expect_error(mpx(rnorm(100), 40L, s_size = 1.1), "between 0 and 1", fixed = TRUE)
})
