brute_force_rfcp <- function(positive, negative, m, max_freq,
                             exclusion_radius = m) {
  normalized_windows <- function(series) {
    count <- length(series) - m + 1L
    valid <- logical(count)
    values <- vector("list", count)
    for (i in seq_len(count)) {
      window <- series[i:(i + m - 1L)]
      if (!all(is.finite(window))) {
        next
      }
      centered <- window - mean(window)
      centered_norm <- sqrt(sum(centered^2))
      if (!is.finite(centered_norm) || centered_norm == 0) {
        next
      }
      valid[[i]] <- TRUE
      values[[i]] <- centered / centered_norm
    }
    list(valid = valid, values = values)
  }

  select <- function(query, reference, valid_reference, query_index = NULL) {
    distance <- rep(NA_real_, length(valid_reference))
    for (j in which(valid_reference)) {
      if (!is.null(query_index) && abs(j - query_index) <= exclusion_radius) {
        next
      }
      distance[[j]] <- sqrt(m) * sqrt(sum((query - reference[[j]])^2))
    }
    order_candidates <- order(distance, seq_along(distance), na.last = NA)
    excluded <- logical(length(distance))
    selected_distance <- rep(NA_real_, max_freq)
    selected_index <- rep(NA_integer_, max_freq)
    accepted <- 0L
    for (j in order_candidates) {
      if (excluded[[j]]) {
        next
      }
      accepted <- accepted + 1L
      selected_distance[[accepted]] <- distance[[j]]
      selected_index[[accepted]] <- j
      left <- max(1L, j - exclusion_radius)
      right <- min(length(distance), j + exclusion_radius)
      excluded[left:right] <- TRUE
      if (accepted == max_freq) {
        break
      }
    }
    list(distance = selected_distance, index = selected_index)
  }

  positive_windows <- normalized_windows(positive)
  negative_windows <- normalized_windows(negative)
  positive_count <- length(positive_windows$valid)
  rfmp_aa <- matrix(NA_real_, max_freq, positive_count)
  rfmp_ab <- matrix(NA_real_, max_freq, positive_count)
  rfmpi_aa <- matrix(NA_integer_, max_freq, positive_count)
  rfmpi_ab <- matrix(NA_integer_, max_freq, positive_count)
  rfcp_profile <- matrix(NA_real_, max_freq, positive_count)
  rms <- rep(NA_real_, positive_count)
  clip <- sqrt(2 * m)

  for (i in which(positive_windows$valid)) {
    aa <- select(
      positive_windows$values[[i]], positive_windows$values,
      positive_windows$valid, i
    )
    ab <- select(
      positive_windows$values[[i]], negative_windows$values,
      negative_windows$valid
    )
    rfmp_aa[, i] <- aa$distance
    rfmp_ab[, i] <- ab$distance
    rfmpi_aa[, i] <- aa$index
    rfmpi_ab[, i] <- ab$index
    paired <- is.finite(aa$distance) & is.finite(ab$distance)
    rfcp_profile[paired, i] <- pmax(
      (pmin(ab$distance[paired], clip) -
        pmin(aa$distance[paired], clip)) / clip,
      0
    )
    squared <- rfcp_profile[, i]^2
    squared[!is.finite(squared)] <- 0
    rms[[i]] <- sqrt(sum(squared) / max_freq)
  }

  list(
    valid_positive = positive_windows$valid,
    valid_negative = negative_windows$valid,
    rfmp_aa = rfmp_aa,
    rfmp_ab = rfmp_ab,
    rfmpi_aa = rfmpi_aa,
    rfmpi_ab = rfmpi_ab,
    rfcp = rfcp_profile,
    rms = rms
  )
}

test_that("native RFCP matches an independent finite oracle", {
  set.seed(3101)
  positive <- rnorm(64)
  negative <- rnorm(57)
  expected <- brute_force_rfcp(positive, negative, 8L, 4L)
  result <- rfcp(
    positive, negative, 8L, 4L,
    return_profiles = TRUE, n_workers = 2L, progress = FALSE
  )

  expect_identical(result$valid_window_positive, expected$valid_positive)
  expect_identical(result$valid_window_negative, expected$valid_negative)
  expect_equal(result$rfmp_aa, expected$rfmp_aa, tolerance = 1e-9)
  expect_equal(result$rfmp_ab, expected$rfmp_ab, tolerance = 1e-9)
  expect_identical(result$rfmpi_aa, expected$rfmpi_aa)
  expect_identical(result$rfmpi_ab, expected$rfmpi_ab)
  expect_equal(result$rfcp, expected$rfcp, tolerance = 1e-9)
  expect_equal(result$rms_profile, expected$rms, tolerance = 1e-9)
})

test_that("native RFCP handles max_freq one and near-constant windows", {
  x <- seq(0, 9, length.out = 72)
  positive <- 1e-4 * (sin(x) + 0.2 * sin(3 * x))
  negative <- 1e-4 * (cos(x) - 0.1 * sin(2 * x))
  expected <- brute_force_rfcp(positive, negative, 9L, 1L)
  result <- rfcp(
    positive, negative, 9L, 1L,
    return_profiles = TRUE, n_workers = 2L, progress = FALSE
  )

  expect_equal(result$rfmp_aa, expected$rfmp_aa, tolerance = 1e-9)
  expect_equal(result$rfmp_ab, expected$rfmp_ab, tolerance = 1e-9)
  expect_identical(result$rfmpi_aa, expected$rfmpi_aa)
  expect_identical(result$rfmpi_ab, expected$rfmpi_ab)
  expect_equal(result$rms_profile, expected$rms, tolerance = 1e-9)
  expect_identical(result$requested_max_freq, 1L)
})

test_that("large block initialization through FFT matches direct distances", {
  set.seed(3103)
  m <- 512L
  max_freq <- 3L
  positive <- rnorm(2600)
  negative <- rnorm(2600)
  query <- positive[seq_len(m)]
  query <- (query - mean(query)) / sqrt(sum((query - mean(query))^2))

  direct_profile <- function(reference) {
    vapply(seq_len(length(reference) - m + 1L), function(index) {
      candidate <- reference[index:(index + m - 1L)]
      candidate <- (candidate - mean(candidate)) /
        sqrt(sum((candidate - mean(candidate))^2))
      sqrt(m * sum((query - candidate)^2))
    }, numeric(1))
  }
  greedy <- function(distance) {
    selected_distance <- rep(NA_real_, max_freq)
    selected_index <- rep(NA_integer_, max_freq)
    excluded <- logical(length(distance))
    accepted <- 0L
    for (index in order(distance, seq_along(distance), na.last = NA)) {
      if (excluded[[index]]) {
        next
      }
      accepted <- accepted + 1L
      selected_distance[[accepted]] <- distance[[index]]
      selected_index[[accepted]] <- index
      excluded[
        max(1L, index - m):min(length(distance), index + m)
      ] <- TRUE
      if (accepted == max_freq) {
        break
      }
    }
    list(distance = selected_distance, index = selected_index)
  }

  aa_profile <- direct_profile(positive)
  aa_profile[seq_len(m + 1L)] <- NA_real_
  expected_aa <- greedy(aa_profile)
  expected_ab <- greedy(direct_profile(negative))
  result <- rfcp(
    positive, negative, m, max_freq,
    query_begin = 1L, query_end = 1L,
    return_profiles = TRUE, n_workers = 2L, progress = FALSE
  )

  expect_equal(result$rfmp_aa[, 1L], expected_aa$distance, tolerance = 1e-9)
  expect_equal(result$rfmp_ab[, 1L], expected_ab$distance, tolerance = 1e-9)
  expect_identical(result$rfmpi_aa[, 1L], expected_aa$index)
  expect_identical(result$rfmpi_ab[, 1L], expected_ab$index)
})

test_that("native RFCP respects every barrier and normalization validity", {
  positive <- c(
    sin(seq(0, 3, length.out = 24)), NA_real_,
    rep(2, 9), NaN, cos(seq(0, 4, length.out = 27)), Inf,
    sin(seq(1, 5, length.out = 22))
  )
  negative <- c(
    cos(seq(0, 5, length.out = 31)), -Inf,
    rep(-1, 8), NA_real_, sin(seq(0, 6, length.out = 35))
  )
  expected <- brute_force_rfcp(positive, negative, 7L, 3L)
  result <- rfcp(
    positive, negative, 7L, 3L,
    return_profiles = TRUE, n_workers = 3L, progress = FALSE
  )

  expect_identical(result$valid_window_positive, expected$valid_positive)
  expect_identical(result$valid_window_negative, expected$valid_negative)
  expect_equal(result$rfmp_aa, expected$rfmp_aa, tolerance = 1e-9)
  expect_equal(result$rfmp_ab, expected$rfmp_ab, tolerance = 1e-9)
  expect_identical(result$rfmpi_aa, expected$rfmpi_aa)
  expect_identical(result$rfmpi_ab, expected$rfmpi_ab)
  expect_equal(result$rms_profile, expected$rms, tolerance = 1e-9)
  expect_true(all(is.na(
    result$rms_profile[!result$valid_window_positive]
  )))
})

test_that("RFCP serial, parallel, and adjacent query blocks compose", {
  set.seed(3102)
  positive <- c(rnorm(53), NA_real_, rnorm(48))
  negative <- c(rnorm(44), Inf, rnorm(51))
  full_serial <- rfcp(
    positive, negative, 9L, 5L,
    return_profiles = TRUE, n_workers = 1L, progress = FALSE
  )
  full_parallel <- rfcp(
    positive, negative, 9L, 5L,
    return_profiles = TRUE, n_workers = 4L, progress = FALSE
  )
  split <- 37L
  left <- rfcp(
    positive, negative, 9L, 5L, query_end = split,
    return_profiles = TRUE, n_workers = 2L, progress = FALSE
  )
  right <- rfcp(
    positive, negative, 9L, 5L, query_begin = split + 1L,
    return_profiles = TRUE, n_workers = 3L, progress = FALSE
  )

  expect_identical(full_parallel, full_serial)
  expect_identical(
    c(left$rms_profile, right$rms_profile),
    full_serial$rms_profile
  )
  expect_identical(cbind(left$rfmp_aa, right$rfmp_aa), full_serial$rfmp_aa)
  expect_identical(cbind(left$rfmp_ab, right$rfmp_ab), full_serial$rfmp_ab)
  expect_identical(cbind(left$rfmpi_aa, right$rfmpi_aa), full_serial$rfmpi_aa)
  expect_identical(cbind(left$rfmpi_ab, right$rfmpi_ab), full_serial$rfmpi_ab)
  expect_identical(cbind(left$rfcp, right$rfcp), full_serial$rfcp)
  combined_rms <- c(left$rms_profile, right$rms_profile)
  expect_identical(full_serial$plato_index, as.numeric(which.max(combined_rms)))
})

test_that("RFCP uses requested frequency denominator and deterministic ties", {
  positive <- rep(c(0, 1, 0, -1, 2), 12)
  negative <- rep(c(0, 1, 0, -1, -2), 10)
  expected <- brute_force_rfcp(positive, negative, 5L, 20L)
  result <- rfcp(
    positive, negative, 5L, 20L,
    return_profiles = TRUE, n_workers = 2L, progress = FALSE
  )

  expect_equal(result$rms_profile, expected$rms, tolerance = 1e-9)
  expect_identical(result$rfmpi_aa, expected$rfmpi_aa)
  expect_identical(result$rfmpi_ab, expected$rfmpi_ab)
  expect_identical(result$requested_max_freq, 20L)
  expect_lt(result$effective_max_freq, result$requested_max_freq)
  expect_equal(
    result$plato_rms,
    result$rms_profile[[result$plato_index]],
    tolerance = 1e-12
  )
  expect_identical(
    result$plato_twin_index,
    result$rfmpi_aa[1L, result$plato_index]
  )
  expect_identical(
    result$rfcp_at_plato,
    result$rfcp[, result$plato_index]
  )
  expect_identical(
    result$rfmp_aa_at_plato,
    result$rfmp_aa[, result$plato_index]
  )
  expect_identical(
    result$rfmp_ab_at_plato,
    result$rfmp_ab[, result$plato_index]
  )
})

test_that("RFCP validates query ranges and profile allocation mode", {
  compact <- rfcp(1:20, 21:40, 5L, 2L, progress = FALSE)
  expect_null(compact[["rfmp_aa"]])
  expect_null(compact[["rfmp_ab"]])
  expect_null(compact[["rfcp"]])

  expect_error(
    rfcp(1:20, 21:40, 5L, 2L, query_begin = 18L),
    "query range"
  )
  expect_error(
    rfcp(
      1:20, 21:40, 5L, 1000000L,
      return_profiles = TRUE, progress = FALSE
    ),
    "5,000,000"
  )
})
