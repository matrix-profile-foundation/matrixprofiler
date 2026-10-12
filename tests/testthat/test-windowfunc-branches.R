# ── parse_na_strategy() and windowfunc branches not yet hit ───────────────────

# Internal helper -------------------------------------------------------

test_that("parse_na_strategy errors on non-character na_strategy", {
  # Passing a number instead of a string must error
  expect_error(
    matrixprofiler:::parse_na_strategy(42L),
    "'na_strategy' must be one of"
  )
})

test_that("parse_na_strategy na.rm = TRUE overrides na_strategy", {
  expect_equal(matrixprofiler:::parse_na_strategy(na.rm = TRUE), 1L)
  expect_equal(matrixprofiler:::parse_na_strategy(na.rm = FALSE), 0L)
})

# mov_mean – na_strategy = "error" branch ---------------------------------

test_that("mov_mean errors when na_strategy = 'error' and data contains NA", {
  expect_error(mov_mean(c(1, 2, NA, 4, 5), window_size = 3L, na_strategy = "error"))
})

test_that("mov_mean respects na.rm = TRUE shortcut", {
  x <- c(1, 2, NA, 4, 5)
  res <- mov_mean(x, 3L, na.rm = TRUE)
  expect_false(is.na(res[2L])) # window [2,NA] → mean(2) = 2
})

# mov_var – na_strategy = "error" -----------------------------------------

test_that("mov_var errors when na_strategy = 'error' and data contains NA", {
  expect_error(mov_var(c(1, 2, NA, 4, 5), window_size = 3L, na_strategy = "error"))
})

# mov_sum – na_strategy = "error" -----------------------------------------

test_that("mov_sum errors when na_strategy = 'error' and data contains NA", {
  expect_error(mov_sum(c(1, 2, NA, 4, 5), window_size = 3L, na_strategy = "error"))
})

# mov_max / mov_min – error branch ----------------------------------------

test_that("mov_max errors when na_strategy = 'error' and data contains NA", {
  expect_error(mov_max(c(1, 2, NA, 4, 5), window_size = 3L, na_strategy = "error"))
})

test_that("mov_min errors when na_strategy = 'error' and data contains NA", {
  expect_error(mov_min(c(1, 2, NA, 4, 5), window_size = 3L, na_strategy = "error"))
})

# mov_std – error branch and pure-R path ----------------------------------

test_that("mov_std errors when na_strategy = 'error' and data contains NA", {
  expect_error(mov_std(c(1, 2, NA, 4, 5), window_size = 3L, na_strategy = "error"))
})

test_that("mov_std pure-R path (rcpp = FALSE) with finite data equals Rcpp path", {
  data <- rnorm(200L)
  r_res <- mov_std(data, 30L, rcpp = FALSE)
  cpp_res <- mov_std(data, 30L, rcpp = TRUE)
  expect_equal(r_res, cpp_res, tolerance = 1e-10)
})

test_that("mov_std na.rm = TRUE shortcut works", {
  x <- c(1, 2, NA, 4, 5, 6)
  res <- mov_std(x, 3L, na.rm = TRUE)
  expect_false(is.na(res[2L]))
})

# movmean_std – error and pure-R branches ---------------------------------

test_that("movmean_std errors when na_strategy = 'error' and data contains NA", {
  expect_error(movmean_std(c(1, 2, NA, 4, 5), window_size = 3L, na_strategy = "error"))
})

test_that("movmean_std pure-R path (rcpp = FALSE) matches Rcpp path", {
  data <- rnorm(200L)
  r_res <- movmean_std(data, 30L, rcpp = FALSE)
  cpp_res <- movmean_std(data, 30L, rcpp = TRUE)
  expect_equal(r_res$avg, cpp_res$avg, tolerance = 1e-10)
  expect_equal(r_res$sd, cpp_res$sd, tolerance = 1e-10)
})

# muinvn – non-numeric n_workers fallback ---------------------------------

test_that("muinvn falls back to 1 worker when n_workers is non-numeric", {
  data <- rnorm(200L)
  # Passing a character should not error — it falls back to 1L internally
  expect_silent(muinvn(data, 30L, n_workers = "bad"))
})

test_that("muinvn na_strategy = 'error' errors on NA data", {
  expect_error(muinvn(c(1, 2, NA, 4, 5), window_size = 3L, na_strategy = "error"))
})

# window_size < 2 guards (all functions) ----------------------------------

test_that("all window functions error when window_size < 2", {
  x <- rnorm(20L)
  expect_error(mov_mean(x, 1L), "'window_size' must be at least 2")
  expect_error(mov_var(x, 1L), "'window_size' must be at least 2")
  expect_error(mov_sum(x, 1L), "'window_size' must be at least 2")
  expect_error(mov_max(x, 1L), "'window_size' must be at least 2")
  expect_error(mov_min(x, 1L), "'window_size' must be at least 2")
  expect_error(mov_std(x, 1L), "'window_size' must be at least 2")
  expect_error(movmean_std(x, 1L), "'window_size' must be at least 2")
  expect_error(muinvn(x, 1L), "'window_size' must be at least 2")
})
