# ── mpx_stream_right() and mpx_stream_s_right() ───────────────────────────────

init_stream_obj <- function(data, w, ez = 0.5) {
  n <- length(data)
  stats <- muinvn(data, w)
  ddf <- c(0.5 * (data[1:(n - w)] - data[(w + 1):n]), 0)
  ddg <- c((data[(w + 1):n] - stats$avg[2:(n - w + 1)]) + (data[1:(n - w)] - stats$avg[1:(n - w)]), 0)
  list(
    data = as.numeric(data),
    w = as.integer(w),
    ez = as.numeric(ez),
    avg = as.numeric(stats$avg),
    sig = as.numeric(stats$sig),
    right_matrix_profile = rep(-1.0, n - w + 1L),
    right_profile_index = rep(-1L, n - w + 1L),
    ddf = as.numeric(ddf),
    ddg = as.numeric(ddg),
    offset = 0L
  )
}

test_that("mpx_stream_right updates a profile with new data", {
  set.seed(200L)
  data0 <- rnorm(100L)
  w <- 10L

  # Build initial stream object
  start_obj <- init_stream_obj(data0, w)

  # Feed 5 new observations
  new_obs <- rnorm(5L)
  res <- expect_silent(
    mpx_stream_right(new_obs, start_obj, history = 0L, mp_time_constraint = 0L, progress = FALSE)
  )
  expect_type(res, "list")
  expect_true("right_matrix_profile" %in% names(res))
  expect_equal(length(res$data), length(data0) + length(new_obs))

  # Feed another observation iteratively
  res2 <- expect_silent(
    mpx_stream_right(rnorm(3L), res, progress = FALSE)
  )
  expect_equal(length(res2$data), length(res$data) + 3L)
})

test_that("mpx_stream_right with time constraint and history", {
  set.seed(201L)
  data0 <- rnorm(100L)
  w <- 10L
  start_obj <- init_stream_obj(data0, w)

  new_obs <- rnorm(5L)
  res <- expect_silent(
    mpx_stream_right(new_obs, start_obj, history = 50L, mp_time_constraint = 20L, progress = FALSE)
  )
  expect_type(res, "list")
  expect_true("right_matrix_profile" %in% names(res))
})

test_that("mpx_stream_s_right processes a batch of data", {
  set.seed(202L)
  full_data <- rnorm(150L)
  w <- 10L

  stats_obj <- list(
    avg = muinvn(full_data, w)$avg,
    sig = muinvn(full_data, w)$sig,
    ddf = c(0.5 * (full_data[1:(150L - w)] - full_data[(w + 1L):150L]), 0),
    ddg = c((full_data[(w + 1L):150L] - muinvn(full_data, w)$avg[2:(150L - w + 1L)]) +
      (full_data[1:(150L - w)] - muinvn(full_data, w)$avg[1:(150L - w)]), 0)
  )

  start_obj <- list(offset = 0L, w = w, ez = 0.5)

  res <- expect_silent(
    mpx_stream_s_right(
      data               = full_data,
      batch_size         = length(full_data),
      start_obj          = start_obj,
      stats              = stats_obj,
      history            = 0L,
      mp_time_constraint = 0L,
      progress           = FALSE,
      threshold          = -1.0
    )
  )
  expect_type(res, "list")
  expect_true("right_matrix_profile" %in% names(res))
})

test_that("mpx_stream_s_right respects history and threshold parameters", {
  set.seed(203L)
  full_data <- rnorm(120L)
  w <- 10L

  stats_obj <- list(
    avg = muinvn(full_data, w)$avg,
    sig = muinvn(full_data, w)$sig,
    ddf = c(0.5 * (full_data[1:(120L - w)] - full_data[(w + 1L):120L]), 0),
    ddg = c((full_data[(w + 1L):120L] - muinvn(full_data, w)$avg[2:(120L - w + 1L)]) +
      (full_data[1:(120L - w)] - muinvn(full_data, w)$avg[1:(120L - w)]), 0)
  )

  start_obj <- list(offset = 0L, w = w, ez = 0.5)

  res <- expect_silent(
    mpx_stream_s_right(
      data               = full_data,
      batch_size         = length(full_data),
      start_obj          = start_obj,
      stats              = stats_obj,
      history            = 50L,
      mp_time_constraint = 30L,
      progress           = FALSE,
      threshold          = 0.8
    )
  )
  expect_type(res, "list")
  expect_true("right_matrix_profile" %in% names(res))
})
