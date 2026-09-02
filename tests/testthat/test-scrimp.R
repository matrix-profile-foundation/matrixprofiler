scrimp_res <- NULL
scrimp_res_par <- NA

test_that("Scrimp", {
  set.seed(2021)
  expect_silent(scrimp_res <<- scrimp(
    data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
    progress = FALSE
  ))
  expect_type(scrimp_res, "list")
  expect_snapshot_value(scrimp_res, style = "serialize")
})

test_that("Scrimp Parallel", {
  set.seed(2021)
  expect_silent(scrimp_res_par <<- scrimp(
    data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
    n_workers = 2, progress = FALSE
  ))
  expect_type(scrimp_res_par, "list")
  expect_snapshot_value(scrimp_res_par, style = "serialize")
})

test_that("Scrimps are equal", {
  expect_equal(scrimp_res, scrimp_res_par)
})

test_that("SCRIMP agrees with MPX on finite self-joins", {
  set.seed(52)
  data <- rnorm(170)
  truth <- mpx(data, 24, exclusion_zone = 0.5, progress = FALSE)

  for (n_workers in c(1, 2)) {
    result <- scrimp(data, 24, exclusion_zone = 0.5, n_workers = n_workers, progress = FALSE)
    expect_equal(result$matrix_profile, truth$matrix_profile, tolerance = 1e-10)
    expect_equal(result$profile_index, truth$profile_index)
  }
})

test_that("SCRIMP handles non-finite and constant windows safely", {
  set.seed(53)
  data <- c(rnorm(30), rep(1, 21), rnorm(69))
  data[c(15, 100)] <- c(NA_real_, Inf)
  original <- data
  truth <- mpx(data, 10, progress = FALSE)

  for (n_workers in c(1, 2)) {
    result <- scrimp(data, 10, n_workers = n_workers, progress = FALSE)
    expect_identical(data, original)
    expect_equal(result$matrix_profile, truth$matrix_profile, tolerance = 1e-10)
    expect_equal(result$profile_index, truth$profile_index)
  }
})

test_that("SCRIMP AB-join agrees with MPXAB and accepts unequal lengths", {
  set.seed(54)
  data <- rnorm(100)
  query <- rnorm(130)
  data[40] <- NA_real_
  query[75] <- Inf
  original_data <- data
  original_query <- query
  truth <- mpx(data, 12, query = query, progress = FALSE)
  result <- scrimp(data, 12, query = query, progress = FALSE)

  expect_identical(data, original_data)
  expect_identical(query, original_query)
  expect_equal(result$matrix_profile, truth$matrix_profile, tolerance = 1e-10)
  expect_equal(result$profile_index, truth$profile_index)
  expect_equal(result$mpb, truth$mpb, tolerance = 1e-10)
  expect_equal(result$pib, truth$pib)
})

test_that("parallel SCRIMP respects anytime sampling", {
  set.seed(55)
  result <- scrimp(rnorm(140), 20, s_size = 0.2, n_workers = 2, progress = FALSE)
  expect_true(result$partial)
})

test_that("SCRIMP validates anytime parameters", {
  expect_error(scrimp(1:20, 5, s_size = 1.1, progress = FALSE), "between 0 and 1")
  expect_error(scrimp(1:20, 5, pre_scrimp = 1.1, progress = FALSE), "between 0 and 1")
})
