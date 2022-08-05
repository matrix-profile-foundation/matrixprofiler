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
