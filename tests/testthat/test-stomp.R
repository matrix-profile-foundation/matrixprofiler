test_that("parallel STOMP processes joins with fewer query windows than the window size", {
  data <- sin(seq(0, 16 * pi, length.out = 500))
  query <- cos(seq(0, 12 * pi, length.out = 500))

  serial <- matrixprofiler:::stomp_rcpp(data, query, 308L, 0, FALSE, FALSE)
  parallel <- expect_silent(
    matrixprofiler:::stomp_rcpp_parallel(data, query, 308L, 0, FALSE, FALSE)
  )

  expect_true(all(is.finite(parallel$matrix_profile)))
  expect_equal(parallel, serial, tolerance = 1e-10)
})

test_that("parallel STOMP handles fewer than 100 diagonals", {
  set.seed(1001)
  data <- rnorm(100)

  serial <- matrixprofiler:::stomp_rcpp(data, data, 40L, 0.5, FALSE, FALSE)
  parallel <- matrixprofiler:::stomp_rcpp_parallel(data, data, 40L, 0.5, FALSE, FALSE)

  expect_equal(parallel, serial, tolerance = 1e-10)
})

test_that("public parallel STOMP AB join supplies all native arguments", {
  set.seed(1002)
  data <- rnorm(800)
  query <- rnorm(700)

  serial <- stomp(data, 100L, query = query, n_workers = 1L, progress = FALSE)
  parallel <- expect_silent(
    stomp(data, 100L, query = query, n_workers = 2L, progress = FALSE)
  )

  expect_equal(parallel, serial, tolerance = 1e-10)
})

test_that("parallel STOMP returns left and right profiles", {
  set.seed(1003)
  data <- rnorm(800)

  serial <- stomp(data, 100L, n_workers = 1L, progress = FALSE, left_right_profile = TRUE)
  parallel <- stomp(data, 100L, n_workers = 2L, progress = FALSE, left_right_profile = TRUE)

  expect_named(parallel, names(serial))
  expect_equal(parallel, serial, tolerance = 1e-10)
  expect_equal(
    parallel$matrix_profile,
    pmin(parallel$left_matrix_profile, parallel$right_matrix_profile),
    tolerance = 1e-10
  )
})
