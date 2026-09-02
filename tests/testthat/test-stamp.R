stamp_res <- NULL
stamp_res_par <- NA

test_that("Stamp", {
  set.seed(2021)
  expect_silent(stamp_res <<- stamp(
    data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
    progress = FALSE
  ))
  expect_type(stamp_res, "list")
  expect_snapshot_value(stamp_res, style = "serialize")
})

test_that("Stamp Parallel", {
  set.seed(2021)
  expect_silent(stamp_res_par <<- stamp(
    data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
    n_workers = 2, progress = FALSE
  ))
  expect_type(stamp_res_par, "list")
  # expect_snapshot_value(stamp_res_par, style = "serialize")
})

test_that("Stamps are equal", {
  expect_equal(stamp_res, stamp_res_par)
})

test_that("STAMP agrees with MPX on finite self-joins", {
  set.seed(42)
  data <- rnorm(160)
  truth <- mpx(data, 24, exclusion_zone = 0.5, progress = FALSE)

  for (n_workers in c(1, 2)) {
    set.seed(7)
    result <- stamp(data, 24, exclusion_zone = 0.5, n_workers = n_workers, progress = FALSE)
    expect_equal(result$matrix_profile, truth$matrix_profile, tolerance = 1e-10)
    expect_equal(result$profile_index, truth$profile_index)
  }
})

test_that("STAMP handles non-finite data without mutating its inputs", {
  set.seed(43)
  data <- rnorm(120)
  query <- rnorm(135)
  data[c(35, 80)] <- c(NA_real_, Inf)
  query[c(20, 100)] <- c(NaN, -Inf)
  original_data <- data
  original_query <- query

  serial <- stamp(data, 16, query = query, n_workers = 1, progress = FALSE)
  set.seed(44)
  parallel <- stamp(data, 16, query = query, n_workers = 2, progress = FALSE)

  expect_identical(data, original_data)
  expect_identical(query, original_query)
  expect_equal(parallel, serial, tolerance = 1e-10)

  invalid <- vapply(seq_len(length(data) - 16 + 1), function(i) {
    any(!is.finite(data[i:(i + 15)]))
  }, logical(1))
  expect_true(all(is.infinite(serial$matrix_profile[invalid])))
  expect_true(all(serial$profile_index[invalid] == -1L))
})

test_that("parallel STAMP respects anytime sampling", {
  set.seed(45)
  result <- stamp(rnorm(140), 20, s_size = 0.2, n_workers = 2, progress = FALSE)
  expect_true(result$partial)
})

test_that("STAMP accepts the smallest supported self-join", {
  set.seed(46)
  expect_type(stamp(rnorm(10), 5, progress = FALSE), "list")
})
