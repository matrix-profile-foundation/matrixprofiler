cp <- NULL
cp_par <- NA

set.seed(2021)
test_that("Contrast Profile", {
  data1 <- motifs_discords_small[1:400]
  data2 <- motifs_discords_small[450:750]
  cp <<- contrast(data1, data2, 50, progress = FALSE)
  expect_snapshot_value(cp, style = "serialize")
})

set.seed(2021)
test_that("Contrast Profile Parallel", {
  data1 <- motifs_discords_small[1:400]
  data2 <- motifs_discords_small[450:750]
  cp_par <<- contrast(data1, data2, 50, n_workers = 2L, progress = FALSE)
  expect_snapshot_value(cp_par, style = "serialize")
})

set.seed(2021)
test_that("Results are equal", {
  expect_equal(cp, cp_par)
})

test_that("Contrast Profile branches: pearson, precomputed matrix, error checks", {
  data1 <- motifs_discords_small[1:200]
  data2 <- motifs_discords_small[250:450]
  w <- 30L

  # Pearson distance
  cp_pearson <- expect_silent(contrast(data1, data2, w, distance = "pearson", progress = FALSE))
  expect_type(cp_pearson, "list")
  expect_false(cp_pearson$euclidean)

  # Precomputed positive_matrix
  pos_mp <- mpx(data2, w, progress = FALSE)
  cp_precomp <- expect_silent(contrast(data1, data2, w, positive_matrix = pos_mp, progress = FALSE))
  expect_type(cp_precomp, "list")

  # Validation errors: window_size >= smaller_size
  expect_error(contrast(data1[1:20], data2, window_size = 25L), "too short relative to desired window size")

  # Validation errors: exclusion zone too large for positive_data
  expect_error(contrast(data1, data2[1:40], window_size = 30L, exclusion_zone = 0.5), "too short relative to desired window size and exclusion zone")
})
