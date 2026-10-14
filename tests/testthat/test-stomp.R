test_that("parallel STOMP processes joins with fewer query windows than the window size", {
  data <- sin(seq(0, 16 * pi, length.out = 500))
  query <- cos(seq(0, 12 * pi, length.out = 500))

  serial <- matrixprofiler:::stomp_rcpp(data, query, 308L, 0, FALSE, FALSE)
  parallel <- expect_silent(
    matrixprofiler:::stomp_rcpp_parallel(data, query, 308L, 0, 2L, FALSE, FALSE)
  )

  expect_true(all(is.finite(parallel$matrix_profile)))
  expect_equal(parallel, serial, tolerance = 1e-10)
})

test_that("parallel STOMP handles fewer than 100 diagonals", {
  set.seed(1001)
  data <- rnorm(100)

  serial <- matrixprofiler:::stomp_rcpp(data, data, 40L, 0.5, FALSE, FALSE)
  parallel <- matrixprofiler:::stomp_rcpp_parallel(data, data, 40L, 0.5, 2L, FALSE, FALSE)

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

test_that("STOMP warns and uses the serial path in a forked child", {
  skip_on_os("windows")
  set.seed(1004)
  data <- rnorm(400)

  # Do not fork a test process that has already initialized native thread pools:
  # the children inherit their allocations, but not their worker threads.
  # Simulate child detection in a local copy without modifying any namespace.
  child_stomp <- stomp
  child_parallel <- list2env(list(isChild = function() TRUE), parent = baseenv())
  parallel_called <- FALSE
  environment(child_stomp) <- list2env(list(
    asNamespace = function(ns) {
      if (identical(ns, "parallel")) child_parallel else base::asNamespace(ns)
    },
    stomp_rcpp_parallel = function(...) {
      parallel_called <<- TRUE
      stop("A forked child must not start native workers")
    }
  ), parent = environment(stomp))

  expect_warning(
    child_result <- child_stomp(data, 50L, n_workers = 2L, progress = FALSE),
    "will use 1 worker", fixed = TRUE
  )
  expect_false(parallel_called)
  serial <- stomp(data, 50L, n_workers = 1L, progress = FALSE)
  expect_equal(child_result, serial, tolerance = 1e-10)
})

test_that("STOMP directional profile indices are 1-based and use NA_integer_ for invalid matches", {
  set.seed(1005)
  data <- rnorm(300)

  serial <- stomp(data, 30L, n_workers = 1L, progress = FALSE, left_right_profile = TRUE)
  parallel <- stomp(data, 30L, n_workers = 2L, progress = FALSE, left_right_profile = TRUE)

  # First window has no left match
  expect_identical(serial$left_profile_index[1], NA_integer_)
  expect_identical(parallel$left_profile_index[1], NA_integer_)

  # Last window has no right match
  n_windows <- length(serial$profile_index)
  expect_identical(serial$right_profile_index[n_windows], NA_integer_)
  expect_identical(parallel$right_profile_index[n_windows], NA_integer_)

  # All valid indices must be in 1..n_windows
  valid_lpi <- serial$left_profile_index[!is.na(serial$left_profile_index)]
  valid_rpi <- serial$right_profile_index[!is.na(serial$right_profile_index)]
  expect_true(all(valid_lpi >= 1L & valid_lpi <= n_windows))
  expect_true(all(valid_rpi >= 1L & valid_rpi <= n_windows))

  expect_equal(serial$left_profile_index, parallel$left_profile_index)
  expect_equal(serial$right_profile_index, parallel$right_profile_index)
})

test_that("STOMP supports NA and infinite values with skip_location", {
  set.seed(1006)
  data <- rnorm(200)
  data[25:30] <- NA_real_
  data[120] <- Inf

  w <- 20L
  original_data <- data

  serial <- stomp(data, w, n_workers = 1L, progress = FALSE)
  parallel <- stomp(data, w, n_workers = 2L, progress = FALSE)

  # Inputs must remain untouched
  expect_identical(data, original_data)

  # Serial and parallel must match
  expect_equal(serial$matrix_profile, parallel$matrix_profile, tolerance = 1e-10)
  expect_identical(serial$profile_index, parallel$profile_index)

  # Windows overlapping NAs/Infs must be infinite distance and NA index
  invalid_windows <- vapply(seq_len(length(data) - w + 1L), function(i) {
    any(!is.finite(data[i:(i + w - 1L)]))
  }, logical(1))

  expect_true(all(is.infinite(serial$matrix_profile[invalid_windows])))
  expect_true(all(is.na(serial$profile_index[invalid_windows])))

  # Valid windows with matches should have finite distance
  valid_windows <- which(!invalid_windows)
  expect_true(any(is.finite(serial$matrix_profile[valid_windows])))
})
