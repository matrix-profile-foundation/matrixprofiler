test_that("mpdist simple and vector work correctly", {
  ref_data <- motifs_discords_small[1:200]
  qe_data <- motifs_discords_small[201:400]
  qd_data <- motifs_discords_small[150:180]
  w <- 20

  # Simple mode
  d_simple <- mpdist(ref_data, qe_data, w, type = "simple")
  expect_type(d_simple, "double")
  expect_length(d_simple, 1L)
  expect_gt(d_simple, 0.0)

  # Vector mode - serial and parallel
  d_vect1 <- mpdist(ref_data, qd_data, w, type = "vector", n_workers = 1L)
  d_vect2 <- mpdist(ref_data, qd_data, w, type = "vector", n_workers = 2L)

  expect_type(d_vect1, "double")
  expect_equal(length(d_vect1), length(ref_data) - length(qd_data) + 1L)
  expect_equal(d_vect1, d_vect2)

  # Argument validation
  expect_error(mpdist(ref_data, qd_data, window_size = 2L), "must be at least 4")
  expect_error(mpdist(ref_data, qd_data[1:10], window_size = 15L), "cannot be larger than query length")

  # Query longer than reference (swapping branch)
  d_swap <- mpdist(qd_data, ref_data, w, type = "simple")
  expect_type(d_swap, "double")
  expect_gt(d_swap, 0.0)
})
