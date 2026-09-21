test_that("rfcp_top_k validates inputs appropriately", {
  expect_error(
    rfcp_top_k(c(1, 2, 3), c(4, 5, 6), window_size = 2L, k_platos = 0L),
    "k_platos"
  )
  expect_error(
    rfcp_top_k(c(1, 2, 3), c(4, 5, 6), window_size = 2L, k_platos = 2L, max_freq = c(2L, 3L, 4L)),
    "length `k_platos`"
  )
  expect_error(
    rfcp_top_k(numeric(0), c(1, 2, 3), window_size = 2L),
    "non-empty"
  )
})

test_that("rfcp_top_k discovers multiple distinct morphologies with frequency adaptation", {
  set.seed(12345)
  m <- 16L
  t_val <- seq(0, 2 * pi, length.out = m)
  sine_shape <- sin(t_val)
  triangle_shape <- c(seq(0, 1, length.out = 8), seq(1, 0, length.out = 8))

  # Create a positive series with 5 sine pulses and 2 triangle pulses
  n_pos <- 350L
  positive <- rnorm(n_pos, sd = 0.02)

  sine_starts <- c(25L, 75L, 125L, 175L, 225L)
  for (st in sine_starts) {
    positive[st:(st + m - 1L)] <- positive[st:(st + m - 1L)] + sine_shape
  }

  triangle_starts <- c(280L, 320L)
  for (st in triangle_starts) {
    positive[st:(st + m - 1L)] <- positive[st:(st + m - 1L)] + triangle_shape * 1.5
  }

  # Negative series is random Gaussian noise
  negative <- rnorm(400L, sd = 0.05)

  # Run rfcp_top_k for 2 Platos with max_freq = c(4, 2)
  res <- rfcp_top_k(
    positive_data = positive,
    negative_data = negative,
    window_size = m,
    k_platos = 2L,
    max_freq = c(4L, 2L),
    exclusion_radius = m,
    progress = FALSE
  )

  expect_s3_class(res, "rfcp_top_k")
  expect_identical(res$k_found, 2L)
  expect_identical(res$k_requested, 2L)
  expect_equal(res$max_freq, c(4L, 2L))

  # Plato 1 should correspond to one of the sine pulse starts (peak offset within window)
  p1_idx <- res$platos[[1L]]$plato_index
  p1_match <- any(abs(sine_starts - p1_idx) <= 4L)
  expect_true(p1_match)
  expect_gt(res$platos[[1L]]$plato_rms, 0.1)

  # All 4 positive self-join neighbors of Plato 1 should match the remaining sine starts
  p1_neighbors <- res$platos[[1L]]$rfmpi_aa_at_plato
  for (nb in p1_neighbors) {
    expect_true(any(abs(sine_starts - nb) <= 4L))
  }

  # Plato 2 should correspond to one of the triangle pulse starts
  p2_idx <- res$platos[[2L]]$plato_index
  p2_match <- any(abs(triangle_starts - p2_idx) <= 15L)
  expect_true(p2_match)
  expect_gt(res$platos[[2L]]$plato_rms, 0.1)

  # The twin of Plato 2 should match the other triangle pulse
  p2_twin <- res$platos[[2L]]$plato_twin_index
  expect_true(any(abs(triangle_starts - p2_twin) <= 15L))

  # Verify that Plato 2 did NOT re-discover any sine pulse
  p2_sine_match <- any(abs(sine_starts - p2_idx) <= 4L)
  expect_false(p2_sine_match)

  # Test print method runs without error
  expect_output(print(res), "Top-2 Relative Frequency Contrast Profile")
})

test_that("rfcp_top_k terminates early when min_rms threshold is not met", {
  set.seed(42)
  m <- 10L
  pos <- rnorm(100L)
  neg <- rnorm(100L)

  # With very high min_rms, no random pattern should qualify
  res <- rfcp_top_k(
    positive_data = pos,
    negative_data = neg,
    window_size = m,
    k_platos = 5L,
    max_freq = 2L,
    min_rms = 0.99,
    progress = FALSE
  )

  expect_identical(res$k_found, 0L)
  expect_identical(length(res$platos), 0L)
  expect_output(print(res), "No Platos found")
})

test_that("rfcp_top_k adapts max_freq when adaptive_max_freq = TRUE", {
  # Synthetic stream with 1 dominant pattern (5 repeats) and 1 secondary pattern (2 repeats)
  set.seed(42)
  pos <- rnorm(500, sd = 0.05)
  sine_pulse <- sin(seq(0, pi, length.out = 20))
  cos_pulse <- cos(seq(0, pi, length.out = 20)) * 1.5

  # Dominant morphology: 5 occurrences
  for (p in c(30, 80, 130, 180, 230)) pos[p:(p + 19)] <- pos[p:(p + 19)] + sine_pulse
  # Secondary morphology: 2 occurrences
  for (p in c(330, 400)) pos[p:(p + 19)] <- pos[p:(p + 19)] + cos_pulse

  neg <- rnorm(500, sd = 0.05)

  res <- rfcp_top_k(
    positive_data = pos,
    negative_data = neg,
    window_size = 20L,
    k_platos = 2L,
    max_freq = 5L,
    adaptive_max_freq = TRUE,
    min_freq = 2L,
    backend = "cpu"
  )

  expect_s3_class(res, "rfcp_top_k")
  expect_equal(res$k_found, 2L)
  expect_true(res$max_freq[[1L]] == 5L)
  # After suppressing dominant pulse (5 occurrences), the next frequency adapts down
  expect_true(res$max_freq[[2L]] <= res$max_freq[[1L]])
  expect_true(res$adaptive_max_freq)
  expect_equal(res$min_freq, 2L)
  expect_output(print(res), "Adaptive: TRUE")
})
