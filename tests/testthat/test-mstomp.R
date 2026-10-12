# ── mstomp() ──────────────────────────────────────────────────────────────────

## shared multivariate fixture (3 dimensions, 200 observations)
set.seed(42L)
mv_mat <- matrix(rnorm(600L), nrow = 200L, ncol = 3L)

test_that("mstomp works with a matrix input", {
  res <- expect_silent(mstomp(mv_mat, window_size = 30L, progress = FALSE))
  expect_type(res, "list")
  expect_true("matrix_profile" %in% names(res))
  expect_true("profile_index" %in% names(res))
  expect_equal(res$n_dim, 3L)
  expect_equal(res$w, 30L)
})

test_that("mstomp works with a data.frame input", {
  df <- as.data.frame(mv_mat)
  res <- expect_silent(mstomp(df, window_size = 30L, progress = FALSE))
  expect_type(res, "list")
  expect_equal(res$n_dim, 3L)
})

test_that("mstomp transposes matrix when ncol > nrow", {
  # 3 rows × 200 cols  →  should be transposed to 200 × 3
  wide <- t(mv_mat)
  expect_equal(nrow(wide), 3L)
  expect_equal(ncol(wide), 200L)
  res <- expect_silent(mstomp(wide, window_size = 30L, progress = FALSE))
  expect_equal(res$n_dim, 3L)
})

test_that("mstomp works with a list input (equal lengths)", {
  lst <- lapply(seq_len(3L), function(i) mv_mat[, i])
  res <- expect_silent(mstomp(lst, window_size = 30L, progress = FALSE))
  expect_type(res, "list")
  expect_equal(res$n_dim, 3L)
})

test_that("mstomp handles list with unequal lengths (pads with NA)", {
  lst <- list(
    rnorm(200L),
    rnorm(190L), # shorter — should be padded
    rnorm(200L)
  )
  res <- expect_silent(mstomp(lst, window_size = 30L, progress = FALSE))
  expect_type(res, "list")
  expect_equal(res$n_dim, 3L)
})

test_that("mstomp works with a univariate vector input", {
  vec <- rnorm(150L)
  res <- expect_silent(mstomp(vec, window_size = 30L, progress = FALSE))
  expect_equal(res$n_dim, 1L)
})

test_that("mstomp errors on unknown data type", {
  # An environment is not a matrix / data.frame / list / vector, so it hits
  # the stop("Unknown type of data...") branch.
  expect_error(mstomp(new.env(), 20L), "Unknown type")
})

test_that("mstomp errors when window_size > data_size / 2", {
  # 200 rows: window_size must be <= 100
  expect_error(mstomp(mv_mat, window_size = 101L, progress = FALSE), "too short")
})

test_that("mstomp errors when window_size < 4", {
  expect_error(mstomp(mv_mat, window_size = 3L, progress = FALSE), "at least 4")
})

test_that("mstomp respects must_dim", {
  res <- expect_silent(
    mstomp(mv_mat, window_size = 30L, must_dim = 1L, progress = FALSE)
  )
  expect_equal(res$must, 1L)
})

test_that("mstomp respects exc_dim", {
  res <- expect_silent(
    mstomp(mv_mat, window_size = 30L, exc_dim = 3L, progress = FALSE)
  )
  expect_equal(res$exc, 3L)
})

test_that("mstomp errors when must_dim is out of range", {
  expect_error(mstomp(mv_mat, window_size = 30L, must_dim = 5L, progress = FALSE), "must_dim")
})

test_that("mstomp errors when exc_dim is out of range", {
  expect_error(mstomp(mv_mat, window_size = 30L, exc_dim = 0L, progress = FALSE), "exc_dim")
})

test_that("mstomp errors when must_dim and exc_dim intersect", {
  expect_error(
    mstomp(mv_mat, window_size = 30L, must_dim = 2L, exc_dim = 2L, progress = FALSE),
    "same dimension"
  )
})

test_that("mstomp runs in parallel (n_workers = 2)", {
  res_serial <- mstomp(mv_mat, window_size = 30L, n_workers = 1L, progress = FALSE)
  res_par <- mstomp(mv_mat, window_size = 30L, n_workers = 2L, progress = FALSE)
  expect_equal(res_serial$matrix_profile, res_par$matrix_profile, tolerance = 1e-10)
})
