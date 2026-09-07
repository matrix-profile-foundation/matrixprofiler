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

test_that("parallel MPXAB omits indices when idxs is FALSE", {
  with_mpx_test_threads <- function(code, n = 2L) {
    previous <- RcppParallel::defaultNumThreads()
    on.exit(RcppParallel::setThreadOptions(numThreads = previous), add = TRUE)
    RcppParallel::setThreadOptions(numThreads = n)
    force(code)
  }

  with_mpx_test_threads({
    set.seed(2022)
    data <- rnorm(180)
    query <- rnorm(140)
    set.seed(2023)
    serial <- matrixprofiler:::mpxab_rcpp(data, query, 15L, 1, FALSE, TRUE, FALSE)
    set.seed(2023)
    parallel <- matrixprofiler:::mpxab_rcpp_parallel(data, query, 15L, 1, FALSE, TRUE, FALSE)

    expect_identical(names(parallel), names(serial))
    expect_equal(parallel, serial, tolerance = 1e-10)
  })
})

test_that("parallel MPXAB worker diagnostics are opt-in", {
  previous <- Sys.getenv("MATRIXPROFILER_PROFILE_WORKERS", unset = NA_character_)
  on.exit({
    if (is.na(previous)) Sys.unsetenv("MATRIXPROFILER_PROFILE_WORKERS") else Sys.setenv(MATRIXPROFILER_PROFILE_WORKERS = previous)
  }, add = TRUE)

  set.seed(2024)
  data <- rnorm(180)
  query <- rnorm(140)
  Sys.unsetenv("MATRIXPROFILER_PROFILE_WORKERS")
  set.seed(2025)
  normal <- matrixprofiler:::mpxab_rcpp_parallel(data, query, 15L, 1, FALSE, TRUE, FALSE)
  expect_false("worker_diagnostics" %in% names(normal))

  Sys.setenv(MATRIXPROFILER_PROFILE_WORKERS = "1")
  set.seed(2025)
  profiled <- matrixprofiler:::mpxab_rcpp_parallel(data, query, 15L, 1, FALSE, TRUE, FALSE)
  expect_true("worker_diagnostics" %in% names(profiled))
  expect_true(profiled$worker_diagnostics$enabled)
  expect_gt(profiled$worker_diagnostics$ab$tasks, 0)
  expect_gt(profiled$worker_diagnostics$ba$tasks, 0)
  expect_gt(profiled$worker_diagnostics$ab$pairs_visited, 0)
  expect_gt(profiled$worker_diagnostics$ba$pairs_visited, 0)
  expect_equal(profiled$matrix_profile, normal$matrix_profile, tolerance = 1e-10)
  expect_equal(profiled$mpb, normal$mpb, tolerance = 1e-10)
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
