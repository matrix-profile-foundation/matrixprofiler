# ── zero_crossing() ───────────────────────────────────────────────────────────
# This function had zero test coverage.

test_that("zero_crossing returns integer counts per window", {
  # Known signal: alternating +1 / -1  → every window crosses zero
  data <- rep(c(1, -1), 100L)
  res <- expect_silent(zero_crossing(data, window_size = 4L))
  expect_type(res, "integer")
  expect_length(res, length(data) - 4L + 1L)
  # For window_size = 4, the algorithm checks w - 2 transitions, yielding 2L
  expect_true(all(res == 2L))
})

test_that("zero_crossing returns zero counts for constant positive signal", {
  data <- rep(2.5, 100L)
  res <- expect_silent(zero_crossing(data, window_size = 10L))
  expect_true(all(res == 0L))
})

test_that("zero_crossing works with a sinusoidal signal", {
  data <- sin(seq(0, 4 * pi, length.out = 200L))
  res <- expect_silent(zero_crossing(data, window_size = 20L))
  expect_type(res, "integer")
  expect_length(res, 181L)
  # Values must be non-negative counts
  expect_true(all(res >= 0L))
})

test_that("zero_crossing uses the real-data fixture", {
  res <- expect_silent(zero_crossing(motifs_discords_small, window_size = 50L))
  expect_type(res, "integer")
  expect_length(res, length(motifs_discords_small) - 50L + 1L)
})
