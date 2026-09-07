brute_force_mpx_na <- function(data, window_size, exclusion_zone = 0.5) {
  profile_length <- length(data) - window_size + 1L
  # C++ std::round rounds positive half values away from zero, unlike R's
  # round-to-even behavior.
  exclusion <- floor(window_size * exclusion_zone + 0.5)
  matrix_profile <- rep(NA_real_, profile_length)
  valid_window <- rep(FALSE, profile_length)
  normalized <- vector("list", profile_length)

  for (i in seq_len(profile_length)) {
    window <- data[i:(i + window_size - 1L)]
    if (!all(is.finite(window))) {
      next
    }

    centered <- window - mean(window)
    window_norm <- sqrt(sum(centered^2))
    if (!is.finite(window_norm) || window_norm == 0) {
      next
    }

    valid_window[i] <- TRUE
    normalized[[i]] <- centered / window_norm
  }

  for (i in which(valid_window)) {
    candidates <- which(valid_window & abs(seq_len(profile_length) - i) > exclusion)
    if (length(candidates) == 0L) {
      next
    }

    distances <- vapply(candidates, function(j) {
      sqrt(window_size) * sqrt(sum((normalized[[i]] - normalized[[j]])^2))
    }, numeric(1))
    matrix_profile[i] <- min(distances)
  }

  list(matrix_profile = matrix_profile, valid_window = valid_window)
}

brute_force_mpxab_na <- function(data, query, window_size) {
  normalize_windows <- function(series) {
    profile_length <- length(series) - window_size + 1L
    valid <- rep(FALSE, profile_length)
    normalized <- vector("list", profile_length)
    for (i in seq_len(profile_length)) {
      window <- series[i:(i + window_size - 1L)]
      if (!all(is.finite(window))) {
        next
      }
      centered <- window - mean(window)
      window_norm <- sqrt(sum(centered^2))
      if (!is.finite(window_norm) || window_norm == 0) {
        next
      }
      valid[i] <- TRUE
      normalized[[i]] <- centered / window_norm
    }
    list(valid = valid, normalized = normalized)
  }

  windows_a <- normalize_windows(data)
  windows_b <- normalize_windows(query)
  profile_a <- rep(NA_real_, length(windows_a$valid))
  profile_b <- rep(NA_real_, length(windows_b$valid))

  for (i in which(windows_a$valid)) {
    distances <- vapply(which(windows_b$valid), function(j) {
      sqrt(window_size) * sqrt(sum((windows_a$normalized[[i]] - windows_b$normalized[[j]])^2))
    }, numeric(1))
    if (length(distances) > 0L) {
      profile_a[i] <- min(distances)
    }
  }
  for (j in which(windows_b$valid)) {
    distances <- vapply(which(windows_a$valid), function(i) {
      sqrt(window_size) * sqrt(sum((windows_b$normalized[[j]] - windows_a$normalized[[i]])^2))
    }, numeric(1))
    if (length(distances) > 0L) {
      profile_b[j] <- min(distances)
    }
  }

  list(
    matrix_profile = profile_a, mpb = profile_b,
    valid_window_a = windows_a$valid, valid_window_b = windows_b$valid
  )
}

with_mpx_test_threads <- function(code, n = 4L) {
  previous <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = previous), add = TRUE)
  RcppParallel::setThreadOptions(numThreads = n)
  force(code)
}

test_that("public mpx routes non-finite self-joins to NA implementations", {
  set.seed(2100)
  data <- rnorm(120)
  data[c(4, 35, 80, 119)] <- c(NA_real_, NaN, Inf, -Inf)

  set.seed(2101)
  serial <- mpx(data, 12L, n_workers = 1, progress = FALSE)
  set.seed(2101)
  parallel <- mpx(data, 12L, n_workers = 2, progress = FALSE)

  expect_identical(serial$valid_window, parallel$valid_window)
  expect_equal(serial$matrix_profile, parallel$matrix_profile, tolerance = 1e-10)
  expect_identical(serial$profile_index, parallel$profile_index)
  expect_true(all(is.na(serial$matrix_profile[!serial$valid_window])))
})

test_that("public mpx routes non-finite AB-joins to NA implementations", {
  set.seed(2102)
  data <- rnorm(120)
  query <- rnorm(90)
  query[c(3, 40, 89)] <- c(NA_real_, NaN, Inf)

  set.seed(2103)
  serial <- mpx(data, 12L, query = query, n_workers = 1, progress = FALSE)
  set.seed(2103)
  parallel <- mpx(data, 12L, query = query, n_workers = 2, progress = FALSE)

  expect_identical(serial$valid_window_a, parallel$valid_window_a)
  expect_identical(serial$valid_window_b, parallel$valid_window_b)
  expect_equal(serial$matrix_profile, parallel$matrix_profile, tolerance = 1e-10)
  expect_equal(serial$mpb, parallel$mpb, tolerance = 1e-10)
  expect_identical(serial$profile_index, parallel$profile_index)
  expect_identical(serial$pib, parallel$pib)
})

test_that("public mpx keeps finite inputs on regular implementations", {
  set.seed(2104)
  data <- rnorm(100)
  query <- rnorm(80)

  self_join <- mpx(data, 10L, n_workers = 1, progress = FALSE)
  ab_join <- mpx(data, 10L, query = query, n_workers = 2, progress = FALSE)

  expect_null(self_join$valid_window)
  expect_null(ab_join$valid_window_a)
  expect_null(ab_join$valid_window_b)
})

test_that("mpx_na_rcpp masks non-finite windows and recovers afterwards", {
  data <- sin(seq(0, 12, length.out = 60))
  data[c(18, 36, 50)] <- c(NA_real_, NaN, Inf)
  window_size <- 8L

  result <- matrixprofiler:::mpx_na_rcpp(data, window_size, 0.5, 1, TRUE, TRUE, FALSE)
  expected_valid <- vapply(seq_len(length(data) - window_size + 1L), function(i) {
    all(is.finite(data[i:(i + window_size - 1L)]))
  }, logical(1))

  expect_identical(result$valid_window, expected_valid)
  expect_true(all(is.na(result$matrix_profile[!expected_valid])))
  expect_true(all(is.na(result$profile_index[!expected_valid])))
  expect_true(all(is.finite(result$matrix_profile[expected_valid])))
})

test_that("mpx_na_rcpp agrees with serial MPX when all samples are finite", {
  set.seed(2101)
  data <- rnorm(120)

  set.seed(2102)
  expected <- matrixprofiler:::mpx_rcpp(data, 16L, 0.5, 1, TRUE, TRUE, FALSE, 60)
  set.seed(2102)
  result <- matrixprofiler:::mpx_na_rcpp(data, 16L, 0.5, 1, TRUE, TRUE, FALSE)

  expect_true(all(result$valid_window))
  expect_equal(result$matrix_profile, expected$matrix_profile, tolerance = 1e-10)
  expect_identical(result$profile_index, expected$profile_index)
})

test_that("parallel mpx_na agrees with serial mpx_na", {
  with_mpx_test_threads({
    set.seed(2111)
    data <- rnorm(1000)
    data[c(4, 191, 507, 998)] <- c(NA_real_, NaN, Inf, -Inf)

    set.seed(2112)
    serial <- matrixprofiler:::mpx_na_rcpp(data, 20L, 0.5, 1, TRUE, TRUE, FALSE)
    set.seed(2112)
    parallel <- matrixprofiler:::mpx_na_rcpp_parallel(data, 20L, 0.5, 1, TRUE, TRUE, FALSE)

    expect_identical(parallel$valid_window, serial$valid_window)
    expect_equal(parallel$matrix_profile, serial$matrix_profile, tolerance = 1e-10)
    expect_identical(parallel$profile_index, serial$profile_index)
    expect_identical(parallel$partial, serial$partial)
  })
})

test_that("parallel mpx_na agrees with regular parallel MPX for finite data", {
  with_mpx_test_threads({
    set.seed(2113)
    data <- rnorm(300)

    set.seed(2114)
    expected <- matrixprofiler:::mpx_rcpp_parallel(data, 25L, 0.5, 1, TRUE, TRUE, FALSE)
    set.seed(2114)
    result <- matrixprofiler:::mpx_na_rcpp_parallel(data, 25L, 0.5, 1, TRUE, TRUE, FALSE)

    expect_true(all(result$valid_window))
    expect_equal(result$matrix_profile, expected$matrix_profile, tolerance = 1e-10)
    expect_identical(result$profile_index, expected$profile_index)
  })
})

test_that("parallel mpx_na honors partial sample size", {
  with_mpx_test_threads({
    set.seed(2115)
    data <- rnorm(500)
    data[c(150, 400)] <- c(NA_real_, Inf)

    set.seed(2116)
    serial <- matrixprofiler:::mpx_na_rcpp(data, 30L, 0.5, 0.25, TRUE, TRUE, FALSE)
    set.seed(2116)
    parallel <- matrixprofiler:::mpx_na_rcpp_parallel(data, 30L, 0.5, 0.25, TRUE, TRUE, FALSE)

    expect_true(parallel$partial)
    expect_identical(parallel$valid_window, serial$valid_window)
    expect_equal(parallel$matrix_profile, serial$matrix_profile, tolerance = 1e-10)
  })
})

test_that("mpxab_na agrees with regular MPXAB for finite data", {
  set.seed(2121)
  data <- rnorm(120)
  query <- rnorm(95)

  set.seed(2122)
  expected <- matrixprofiler:::mpxab_rcpp(data, query, 15L, 1, TRUE, TRUE, FALSE)
  set.seed(2122)
  result <- matrixprofiler:::mpxab_na_rcpp(data, query, 15L, 1, TRUE, TRUE, FALSE)

  expect_true(all(result$valid_window_a))
  expect_true(all(result$valid_window_b))
  expect_equal(result$matrix_profile, expected$matrix_profile, tolerance = 1e-10)
  expect_equal(result$mpb, expected$mpb, tolerance = 1e-10)
  expect_identical(result$profile_index, expected$profile_index)
  expect_identical(result$pib, expected$pib)
})

test_that("mpxab_na agrees with STAMP in both directions with non-finite data", {
  set.seed(2123)
  data <- rnorm(90)
  query <- rnorm(75)
  data[c(2, 40)] <- c(NA_real_, Inf)
  query[c(20, 74)] <- c(NaN, -Inf)
  window_size <- 10L

  result <- matrixprofiler:::mpxab_na_rcpp(data, query, window_size, 1, TRUE, TRUE, FALSE)
  stamp_a <- matrixprofiler:::stamp_rcpp(data, query, window_size, 0, 1, FALSE)
  stamp_b <- matrixprofiler:::stamp_rcpp(query, data, window_size, 0, 1, FALSE)

  expect_equal(result$matrix_profile[result$valid_window_a], stamp_a$matrix_profile[result$valid_window_a],
    tolerance = 1e-10
  )
  expect_equal(result$mpb[result$valid_window_b], stamp_b$matrix_profile[result$valid_window_b],
    tolerance = 1e-10
  )
  expect_true(all(is.na(result$matrix_profile[!result$valid_window_a])))
  expect_true(all(is.na(result$mpb[!result$valid_window_b])))
})

test_that("mpxab_na honors partial sample size", {
  set.seed(2124)
  data <- rnorm(100)
  query <- rnorm(80)
  data[30] <- NA_real_
  query[50] <- Inf

  result <- matrixprofiler:::mpxab_na_rcpp(data, query, 12L, 0.25, TRUE, TRUE, FALSE)

  expect_true(result$partial)
  expect_length(result$matrix_profile, length(data) - 12L + 1L)
  expect_length(result$mpb, length(query) - 12L + 1L)
  expect_true(all(is.na(result$matrix_profile[!result$valid_window_a])))
  expect_true(all(is.na(result$mpb[!result$valid_window_b])))
})

test_that("parallel mpxab_na agrees with serial mpxab_na", {
  with_mpx_test_threads({
    set.seed(2125)
    data <- rnorm(1000)
    query <- rnorm(800)
    data[c(3, 301, 999)] <- c(NA_real_, Inf, -Inf)
    query[c(150, 799)] <- c(NaN, NA_real_)

    set.seed(2126)
    serial <- matrixprofiler:::mpxab_na_rcpp(data, query, 20L, 1, TRUE, TRUE, FALSE)
    set.seed(2126)
    parallel <- matrixprofiler:::mpxab_na_rcpp_parallel(data, query, 20L, 1, TRUE, TRUE, FALSE)

    expect_identical(parallel$valid_window_a, serial$valid_window_a)
    expect_identical(parallel$valid_window_b, serial$valid_window_b)
    expect_equal(parallel$matrix_profile, serial$matrix_profile, tolerance = 1e-10)
    expect_equal(parallel$mpb, serial$mpb, tolerance = 1e-10)
    expect_identical(parallel$profile_index, serial$profile_index)
    expect_identical(parallel$pib, serial$pib)
  })
})

test_that("parallel mpxab_na agrees with regular parallel MPXAB for finite data", {
  with_mpx_test_threads({
    set.seed(2127)
    data <- rnorm(300)
    query <- rnorm(250)

    set.seed(2128)
    expected <- matrixprofiler:::mpxab_rcpp_parallel(data, query, 20L, 1, TRUE, TRUE, FALSE)
    set.seed(2128)
    result <- matrixprofiler:::mpxab_na_rcpp_parallel(data, query, 20L, 1, TRUE, TRUE, FALSE)

    expect_true(all(result$valid_window_a))
    expect_true(all(result$valid_window_b))
    expect_equal(result$matrix_profile, expected$matrix_profile, tolerance = 1e-10)
    expect_equal(result$mpb, expected$mpb, tolerance = 1e-10)
  })
})

test_that("parallel mpxab_na restarts covariance for dense invalid blocks", {
  with_mpx_test_threads({
    set.seed(2131)
    data <- rnorm(900)
    query <- rnorm(1000)
    data[seq(80, length(data), by = 80)] <- NA_real_
    query[seq(100, length(query), by = 100)] <- NaN

    set.seed(2132)
    serial <- matrixprofiler:::mpxab_na_rcpp(data, query, 30L, 1, TRUE, TRUE, FALSE)
    set.seed(2132)
    parallel <- matrixprofiler:::mpxab_na_rcpp_parallel(data, query, 30L, 1, TRUE, TRUE, FALSE)

    expect_equal(parallel$matrix_profile, serial$matrix_profile, tolerance = 1e-10)
    expect_equal(parallel$mpb, serial$mpb, tolerance = 1e-10)
    expect_identical(parallel$profile_index, serial$profile_index)
    expect_identical(parallel$pib, serial$pib)
  })
})

test_that("parallel mpxab_na threshold override preserves profiles", {
  with_mpx_test_threads({
    old_threshold <- Sys.getenv("MATRIXPROFILER_INVALID_BLOCK_THRESHOLD", unset = NA_character_)
    on.exit({
      if (is.na(old_threshold)) Sys.unsetenv("MATRIXPROFILER_INVALID_BLOCK_THRESHOLD")
      else Sys.setenv(MATRIXPROFILER_INVALID_BLOCK_THRESHOLD = old_threshold)
    }, add = TRUE)

    set.seed(2133)
    data <- rnorm(700)
    query <- rnorm(800)
    data[seq(90, length(data), by = 90)] <- NA_real_
    query[seq(110, length(query), by = 110)] <- Inf

    Sys.setenv(MATRIXPROFILER_INVALID_BLOCK_THRESHOLD = "legacy")
    set.seed(2134)
    legacy <- matrixprofiler:::mpxab_na_rcpp_parallel(data, query, 30L, 1, TRUE, TRUE, FALSE)
    Sys.setenv(MATRIXPROFILER_INVALID_BLOCK_THRESHOLD = "force")
    set.seed(2134)
    forced <- matrixprofiler:::mpxab_na_rcpp_parallel(data, query, 30L, 1, TRUE, TRUE, FALSE)

    expect_identical(forced$valid_window_a, legacy$valid_window_a)
    expect_identical(forced$valid_window_b, legacy$valid_window_b)
    expect_equal(forced$matrix_profile, legacy$matrix_profile, tolerance = 1e-10)
    expect_equal(forced$mpb, legacy$mpb, tolerance = 1e-10)
    expect_identical(forced$profile_index, legacy$profile_index)
    expect_identical(forced$pib, legacy$pib)
  })
})

test_that("parallel NA-aware joins omit indices when idxs is FALSE", {
  with_mpx_test_threads({
    set.seed(2128)
    data <- rnorm(300)
    query <- rnorm(250)
    data[80] <- NA_real_
    query[170] <- Inf

    result <- matrixprofiler:::mpxab_na_rcpp_parallel(data, query, 20L, 1, FALSE, TRUE, FALSE)

    expect_false(any(c("profile_index", "pib") %in% names(result)))
    expect_true(all(is.na(result$matrix_profile[!result$valid_window_a])))
    expect_true(all(is.na(result$mpb[!result$valid_window_b])))
  })
})

test_that("parallel mpxab_na honors partial sample size", {
  with_mpx_test_threads({
    set.seed(2129)
    data <- rnorm(500)
    query <- rnorm(400)
    data[200] <- NA_real_
    query[300] <- Inf

    set.seed(2130)
    serial <- matrixprofiler:::mpxab_na_rcpp(data, query, 18L, 0.25, TRUE, TRUE, FALSE)
    set.seed(2130)
    parallel <- matrixprofiler:::mpxab_na_rcpp_parallel(data, query, 18L, 0.25, TRUE, TRUE, FALSE)

    expect_true(parallel$partial)
    expect_equal(parallel$matrix_profile, serial$matrix_profile, tolerance = 1e-10)
    expect_equal(parallel$mpb, serial$mpb, tolerance = 1e-10)
  })
})

test_that("mpx_na_rcpp invalidates constant windows deterministically", {
  data <- c(rep(3, 12), seq_len(30))
  window_size <- 6L
  result <- matrixprofiler:::mpx_na_rcpp(data, window_size, 0.5, 1, TRUE, TRUE, FALSE)

  expect_false(any(result$valid_window[1:7]))
  expect_true(all(is.na(result$matrix_profile[1:7])))
  expect_true(all(result$valid_window[8:length(result$valid_window)]))
})

test_that("mpx_na_rcpp normalizes windows directly for large offsets", {
  data <- 1e8 + sin(seq(0, 20, length.out = 100))
  window_size <- 12L
  oracle <- brute_force_mpx_na(data, window_size)
  result <- matrixprofiler:::mpx_na_rcpp(data, window_size, 0.5, 1, TRUE, TRUE, FALSE)

  expect_true(all(result$valid_window))
  expect_lt(max(abs(result$matrix_profile - oracle$matrix_profile)), 1e-5)
})

test_that("parallel mpx_na_rcpp also preserves large-offset precision", {
  with_mpx_test_threads({
    data <- 1e8 + sin(seq(0, 20, length.out = 100))
    window_size <- 12L
    oracle <- brute_force_mpx_na(data, window_size)
    result <- matrixprofiler:::mpx_na_rcpp_parallel(data, window_size, 0.5, 1, TRUE, TRUE, FALSE)

    expect_identical(result$valid_window, oracle$valid_window)
    expect_lt(max(abs(result$matrix_profile - oracle$matrix_profile)), 1e-5)
  })
})

test_that("mpx_na_rcpp validates its numerical controls", {
  data <- rnorm(30)

  expect_error(
    matrixprofiler:::mpx_na_rcpp(data, 5L, -0.5, 1, TRUE, TRUE, FALSE),
    "finite and non-negative", fixed = TRUE
  )
  expect_error(
    matrixprofiler:::mpx_na_rcpp(data, 5L, 0.5, 1.1, TRUE, TRUE, FALSE),
    "between 0 and 1", fixed = TRUE
  )
  expect_error(
    matrixprofiler:::mpx_na_rcpp_parallel(data, 5L, -0.5, 1, TRUE, TRUE, FALSE),
    "finite and non-negative", fixed = TRUE
  )
  expect_error(
    matrixprofiler:::mpx_na_rcpp_parallel(data, 5L, 0.5, 1.1, TRUE, TRUE, FALSE),
    "between 0 and 1", fixed = TRUE
  )
  expect_error(
    matrixprofiler:::mpxab_na_rcpp(data, rev(data), 5L, 1.1, TRUE, TRUE, FALSE),
    "between 0 and 1", fixed = TRUE
  )
  expect_error(
    matrixprofiler:::mpxab_na_rcpp_parallel(data, rev(data), 5L, 1.1, TRUE, TRUE, FALSE),
    "between 0 and 1", fixed = TRUE
  )
})

test_that("serial STAMP does not modify data containing non-finite values", {
  data <- c(1:12, NA_real_, 14:30)
  original <- data

  matrixprofiler:::stamp_rcpp(data, data, 5L, 0.5, 1, FALSE)

  expect_identical(data, original)
})

test_that("mpx_na_rcpp stress test agrees with STAMP and brute force", {
  set.seed(2103)

  for (iteration in seq_len(40)) {
    data_size <- sample(45:90, 1)
    window_size <- sample(5:min(15, floor(data_size / 3)), 1)
    data <- rnorm(data_size)
    missing_count <- sample(0:4, 1)
    if (missing_count > 0) {
      locations <- sample(seq_len(data_size), missing_count)
      data[locations] <- sample(c(NA_real_, NaN, Inf, -Inf), missing_count, replace = TRUE)
    }

    oracle <- brute_force_mpx_na(data, window_size)
    mpx_result <- matrixprofiler:::mpx_na_rcpp(data, window_size, 0.5, 1, TRUE, TRUE, FALSE)
    parallel_result <- matrixprofiler:::mpx_na_rcpp_parallel(data, window_size, 0.5, 1, TRUE, TRUE, FALSE)
    stamp_result <- matrixprofiler:::stamp_rcpp(data, data, window_size, 0.5, 1, FALSE)
    comparable <- is.finite(oracle$matrix_profile)

    expect_identical(mpx_result$valid_window, oracle$valid_window,
      info = paste("valid-window iteration", iteration)
    )
    expect_equal(mpx_result$matrix_profile[comparable], oracle$matrix_profile[comparable],
      tolerance = 1e-8, info = paste("brute-force iteration", iteration)
    )
    expect_equal(mpx_result$matrix_profile[comparable], stamp_result$matrix_profile[comparable],
      tolerance = 1e-8, info = paste("STAMP iteration", iteration)
    )
    expect_identical(parallel_result$valid_window, oracle$valid_window,
      info = paste("parallel valid-window iteration", iteration)
    )
    expect_equal(parallel_result$matrix_profile[comparable], oracle$matrix_profile[comparable],
      tolerance = 1e-8, info = paste("parallel brute-force iteration", iteration)
    )
    expect_true(all(is.na(mpx_result$matrix_profile[!oracle$valid_window])),
      info = paste("invalid-profile iteration", iteration)
    )
    expect_true(all(is.na(parallel_result$matrix_profile[!oracle$valid_window])),
      info = paste("parallel invalid-profile iteration", iteration)
    )
  }
})

test_that("mpxab_na stress test agrees with bidirectional STAMP and brute force", {
  set.seed(2130)

  for (iteration in seq_len(40)) {
    data_size <- sample(45:90, 1)
    query_size <- sample(40:85, 1)
    window_size <- sample(5:min(15, floor(min(data_size, query_size) / 3)), 1)
    data <- rnorm(data_size)
    query <- rnorm(query_size)
    missing_a <- sample(0:3, 1)
    missing_b <- sample(0:3, 1)
    if (missing_a > 0) {
      data[sample(seq_len(data_size), missing_a)] <-
        sample(c(NA_real_, NaN, Inf, -Inf), missing_a, replace = TRUE)
    }
    if (missing_b > 0) {
      query[sample(seq_len(query_size), missing_b)] <-
        sample(c(NA_real_, NaN, Inf, -Inf), missing_b, replace = TRUE)
    }

    oracle <- brute_force_mpxab_na(data, query, window_size)
    result <- matrixprofiler:::mpxab_na_rcpp(data, query, window_size, 1, TRUE, TRUE, FALSE)
    parallel <- matrixprofiler:::mpxab_na_rcpp_parallel(data, query, window_size, 1, TRUE, TRUE, FALSE)
    stamp_a <- matrixprofiler:::stamp_rcpp(data, query, window_size, 0, 1, FALSE)
    stamp_b <- matrixprofiler:::stamp_rcpp(query, data, window_size, 0, 1, FALSE)
    comparable_a <- is.finite(oracle$matrix_profile)
    comparable_b <- is.finite(oracle$mpb)

    expect_identical(result$valid_window_a, oracle$valid_window_a,
      info = paste("A mask iteration", iteration)
    )
    expect_identical(result$valid_window_b, oracle$valid_window_b,
      info = paste("B mask iteration", iteration)
    )
    expect_equal(result$matrix_profile[comparable_a], oracle$matrix_profile[comparable_a],
      tolerance = 1e-8, info = paste("A brute-force iteration", iteration)
    )
    expect_equal(result$mpb[comparable_b], oracle$mpb[comparable_b],
      tolerance = 1e-8, info = paste("B brute-force iteration", iteration)
    )
    expect_equal(result$matrix_profile[comparable_a], stamp_a$matrix_profile[comparable_a],
      tolerance = 1e-8, info = paste("A STAMP iteration", iteration)
    )
    expect_equal(result$mpb[comparable_b], stamp_b$matrix_profile[comparable_b],
      tolerance = 1e-8, info = paste("B STAMP iteration", iteration)
    )
    expect_identical(parallel$valid_window_a, oracle$valid_window_a,
      info = paste("parallel A mask iteration", iteration)
    )
    expect_identical(parallel$valid_window_b, oracle$valid_window_b,
      info = paste("parallel B mask iteration", iteration)
    )
    expect_equal(parallel$matrix_profile[comparable_a], oracle$matrix_profile[comparable_a],
      tolerance = 1e-8, info = paste("parallel A brute-force iteration", iteration)
    )
    expect_equal(parallel$mpb[comparable_b], oracle$mpb[comparable_b],
      tolerance = 1e-8, info = paste("parallel B brute-force iteration", iteration)
    )
  }
})
