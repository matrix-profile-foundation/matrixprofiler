test_that("Conversion between Euclidean Distance and Correlation are equal.", {
  set.seed(2021)
  corr <- runif(1000)

  ed <- corr_ed(corr, 100)
  newcorr <- ed_corr(ed, 100)

  expect_equal(corr, newcorr)
})

test_that("Mode", {
  set.seed(2021)
  data_mode <- c(round(rnorm(500, mean = 10, sd = 2)), round(rnorm(500, mean = 5, sd = 2.5)))
  test_mode <- mode(data_mode)
  expect_snapshot_value(test_mode, style = "json")
})
test_that("Std", {
  data_sd <- std(motifs_discords_small)
  expect_snapshot_value(data_sd, style = "serialize")
})
test_that("znorm", {
  data_znorm <- znorm(motifs_discords_small)
  expect_snapshot_value(data_znorm, style = "serialize")
})
test_that("normalize", {
  data_normalized <- normalize(motifs_discords_small, 1, 5)
  expect_snapshot_value(data_normalized, style = "serialize")
})
test_that("complexity", {
  data_complex <- complexity(motifs_discords_small) # Computes the complexity index of the data
  expect_snapshot_value(data_complex, style = "serialize")
})
test_that("binary_split", {
  split <- binary_split(50)
  expect_snapshot_value(split, style = "json2")
})
