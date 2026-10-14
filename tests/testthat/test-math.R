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

test_that("Math functions R branches and edge cases", {
  x <- motifs_discords_small[1:50]

  # znorm pure-R path and low-variance branch
  expect_equal(znorm(x, rcpp = FALSE), znorm(x, rcpp = TRUE), tolerance = 1e-6)
  const_data <- rep(5, 20)
  expect_equal(znorm(const_data, rcpp = FALSE), rep(0, 20))

  # ed_corr and corr_ed pure-R paths
  w <- 10L
  corr <- c(0.1, 0.5, 0.9, 1.2)
  ed <- corr_ed(corr, w, rcpp = FALSE)
  expect_type(ed, "double")
  corr_back <- ed_corr(ed, w, rcpp = FALSE)
  expect_type(corr_back, "double")

  # mode Rcpp path
  m_rcpp <- mode(c(1L, 2L, 2L, 3L), rcpp = TRUE)
  m_r <- mode(c(1L, 2L, 2L, 3L), rcpp = FALSE)
  expect_equal(m_rcpp, m_r)

  # std pure-R path and NA handling
  expect_equal(std(x, rcpp = FALSE), std(x, rcpp = TRUE), tolerance = 1e-6)
  expect_true(is.na(std(c(1, 2, NA), na.rm = FALSE, rcpp = FALSE)))
  expect_false(is.na(std(c(1, 2, NA), na.rm = TRUE, rcpp = FALSE)))

  # normalize Rcpp path with and without NAs
  norm_rcpp <- normalize(x, 0, 1, rcpp = TRUE)
  norm_r <- normalize(x, 0, 1, rcpp = FALSE)
  expect_equal(norm_rcpp, norm_r, tolerance = 1e-6)
  x_na <- c(x[1:10], NA, x[11:20])
  norm_na <- normalize(x_na, 0, 1, rcpp = TRUE)
  expect_true(is.na(norm_na[11]))

  # binary_split pure-R path and edge cases
  expect_equal(binary_split(1, rcpp = FALSE), 1)
  expect_equal(binary_split(20, rcpp = FALSE), binary_split(20, rcpp = TRUE))
})
