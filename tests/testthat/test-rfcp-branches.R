# ── rfcp() – branches not yet exercised ───────────────────────────────────────

test_that("rfcp errors when positive_data is empty", {
  expect_error(
    rfcp(numeric(0), rnorm(100L), window_size = 10L, max_freq = 2L),
    "non-empty"
  )
})

test_that("rfcp errors when negative_data is empty", {
  expect_error(
    rfcp(rnorm(100L), numeric(0), window_size = 10L, max_freq = 2L),
    "non-empty"
  )
})

test_that("rfcp errors when window_size exceeds either input length", {
  pos <- rnorm(50L)
  neg <- rnorm(80L)
  expect_error(
    rfcp(pos, neg, window_size = 51L, max_freq = 2L),
    "must not exceed"
  )
})

test_that("rfcp errors when query range is invalid (begin > end)", {
  pos <- rnorm(100L)
  neg <- rnorm(150L)
  expect_error(
    rfcp(pos, neg,
      window_size = 10L, max_freq = 2L,
      query_begin = 80L, query_end = 10L
    ),
    "query range"
  )
})

test_that("rfcp errors when query_end exceeds window count", {
  pos <- rnorm(100L)
  neg <- rnorm(150L)
  # n_positive_windows = 100 - 10 + 1 = 91
  expect_error(
    rfcp(pos, neg,
      window_size = 10L, max_freq = 2L,
      query_begin = 1L, query_end = 200L
    ),
    "query range"
  )
})

test_that("rfcp runs with explicit exclusion_radius", {
  set.seed(301L)
  pos <- sin(seq(0, 8 * pi, length.out = 200L))
  neg <- cos(seq(0, 12 * pi, length.out = 250L))
  res <- expect_silent(
    rfcp(pos, neg,
      window_size = 20L, max_freq = 3L,
      exclusion_radius = 5L, progress = FALSE
    )
  )
  expect_type(res, "list")
  expect_true("rms_profile" %in% names(res))
})

test_that("rfcp runs with return_profiles = TRUE", {
  set.seed(302L)
  pos <- rnorm(100L)
  neg <- rnorm(120L)
  res <- expect_silent(
    rfcp(pos, neg,
      window_size = 10L, max_freq = 2L,
      return_profiles = TRUE, progress = FALSE
    )
  )
  expect_true("rfcp" %in% names(res))
  expect_true("rfmp_aa" %in% names(res))
  expect_true("rfmp_ab" %in% names(res))
})

test_that("rfcp runs with parallel workers (n_workers = 2)", {
  set.seed(303L)
  pos <- rnorm(200L)
  neg <- rnorm(250L)
  res_s <- rfcp(pos, neg, window_size = 20L, max_freq = 2L, n_workers = 1L, progress = FALSE)
  res_p <- rfcp(pos, neg, window_size = 20L, max_freq = 2L, n_workers = 2L, progress = FALSE)
  expect_equal(res_s$rms_profile, res_p$rms_profile, tolerance = 1e-10)
})

test_that("rfcp default exclusion_radius equals window_size", {
  set.seed(304L)
  pos <- rnorm(100L)
  neg <- rnorm(120L)
  w <- 10L
  # Run once with NULL (default) and once with explicit window_size
  res_default <- rfcp(pos, neg, w, max_freq = 2L, exclusion_radius = NULL, progress = FALSE)
  res_explicit <- rfcp(pos, neg, w, max_freq = 2L, exclusion_radius = w, progress = FALSE)
  expect_equal(res_default$rms_profile, res_explicit$rms_profile, tolerance = 1e-10)
})

test_that("rfcp with a partial query range respects query_begin and query_end", {
  set.seed(305L)
  pos <- rnorm(200L)
  neg <- rnorm(250L)
  w <- 20L
  n_wins <- length(pos) - w + 1L # 181

  res_full <- rfcp(pos, neg, w,
    max_freq = 2L,
    query_begin = 1L, query_end = n_wins, progress = FALSE
  )
  res_partial <- rfcp(pos, neg, w,
    max_freq = 2L,
    query_begin = 10L, query_end = 50L, progress = FALSE
  )

  # Partial run covers a sub-range
  expect_equal(res_partial$query_begin, 10L)
  expect_equal(res_partial$query_end, 50L)
  # Full run's profile at positions 10:50 should match the 41-element partial profile
  expect_equal(
    res_full$rms_profile[10:50],
    as.numeric(res_partial$rms_profile),
    tolerance = 1e-10
  )
})
