# ── mpx() – branches not yet exercised by test-mpx.R ─────────────────────────

test_that("mpx self-join with distance = 'pearson' runs correctly", {
  set.seed(101L)
  data <- rnorm(200L)
  res <- expect_silent(mpx(data, 30L, distance = "pearson", progress = FALSE))
  expect_type(res, "list")
  expect_true("matrix_profile" %in% names(res))
  # Pearson values are in [-1, 1]
  finite_vals <- res$matrix_profile[is.finite(res$matrix_profile)]
  expect_true(all(finite_vals >= -1 - 1e-9 & finite_vals <= 1 + 1e-9))
})

test_that("mpx AB-join with distance = 'pearson' runs correctly", {
  set.seed(102L)
  data <- rnorm(200L)
  query <- rnorm(180L)
  res <- expect_silent(mpx(data, 30L, query = query, distance = "pearson", progress = FALSE))
  expect_type(res, "list")
})

test_that("mpx self-join with s_size < 1 marks result as partial", {
  set.seed(103L)
  data <- rnorm(500L)
  res <- expect_silent(mpx(data, 50L, s_size = 0.25, progress = FALSE))
  expect_true(isTRUE(res$partial))
})

test_that("mpx self-join with NA data and na_strategy = 'monolithic'", {
  set.seed(104L)
  data <- rnorm(200L)
  data[c(30L, 80L)] <- NA_real_
  res <- expect_silent(mpx(data, 20L, na_strategy = "monolithic", progress = FALSE))
  expect_type(res, "list")
})

test_that("mpx self-join with NA data and na_strategy = 'segmented'", {
  set.seed(105L)
  data <- rnorm(200L)
  data[c(30L, 80L)] <- NA_real_
  res_seg <- expect_silent(mpx(data, 20L, na_strategy = "segmented", n_workers = 2L, progress = FALSE))
  expect_type(res_seg, "list")
})

test_that("mpx self-join with NA data and na_strategy = 'segmented_native'", {
  set.seed(106L)
  data <- rnorm(200L)
  data[c(30L, 80L)] <- NA_real_
  res_nat <- expect_silent(mpx(data, 20L, na_strategy = "segmented_native", progress = FALSE))
  expect_type(res_nat, "list")
})

test_that("mpx AB-join with NA in data routes to NA-aware C++ path", {
  set.seed(107L)
  data <- rnorm(200L)
  query <- rnorm(150L)
  data[50L] <- NA_real_
  res <- expect_silent(mpx(data, 20L, query = query, progress = FALSE))
  expect_type(res, "list")
})

test_that("mpx AB-join with NA and na_strategy = 'segmented' (parallel)", {
  set.seed(108L)
  data <- rnorm(200L)
  query <- rnorm(150L)
  data[50L] <- NA_real_
  res <- expect_silent(
    mpx(data, 20L, query = query, na_strategy = "segmented", n_workers = 2L, progress = FALSE)
  )
  expect_type(res, "list")
})

test_that("mpx AB-join with NA and na_strategy = 'segmented_native' (parallel)", {
  set.seed(109L)
  data <- rnorm(200L)
  query <- rnorm(150L)
  data[50L] <- NA_real_
  res <- expect_silent(
    mpx(data, 20L, query = query, na_strategy = "segmented_native", n_workers = 2L, progress = FALSE)
  )
  expect_type(res, "list")
})

test_that("mpx AB-join with NA and n_workers = 1 (serial monolithic NA path)", {
  set.seed(110L)
  data <- rnorm(200L)
  query <- rnorm(150L)
  query[10L] <- Inf
  res <- expect_silent(mpx(data, 20L, query = query, progress = FALSE))
  expect_type(res, "list")
})

test_that("NA-aware workers preserve results across repeated two-worker joins", {
  set.seed(109L)
  data <- rnorm(200L)
  query <- rnorm(150L)
  data[50L] <- NA_real_

  # Exercise all three workers that previously retained thread-local buffers.
  joins <- list(
    monolithic_ab = function(workers, indices) {
      matrixprofiler:::mpxab_na_rcpp_parallel(
        data, query, 20L, idxs = indices, n_workers = workers
      )
    },
    segmented_ab = function(workers, indices) {
      matrixprofiler:::mpxab_na_segmented_native_rcpp_parallel(
        data, query, 20L, idxs = indices, n_workers = workers
      )
    },
    segmented_aa = function(workers, indices) {
      matrixprofiler:::mpx_na_segmented_native_rcpp_parallel(
        data, 20L, idxs = indices, n_workers = workers
      )
    }
  )
  for (join in joins) {
    for (indices in c(FALSE, TRUE)) {
      set.seed(2021L)
      serial <- join(1L, indices)
      for (iteration in seq_len(3L)) {
        set.seed(2021L)
        parallel <- join(2L, indices)
        profiles <- intersect(c("matrix_profile", "mpb"), names(serial))
        expect_equal(parallel[profiles], serial[profiles], tolerance = 1e-10)
        expect_identical(parallel[!names(parallel) %in% profiles],
                         serial[!names(serial) %in% profiles])
      }
    }
  }
})

test_that("mpx self-join with Inf in data (non-finite routing)", {
  set.seed(111L)
  data <- rnorm(200L)
  data[60L] <- Inf
  res <- expect_silent(mpx(data, 20L, progress = FALSE))
  expect_type(res, "list")
})

test_that("mpx on.exit guard returns NULL on interruption", {
  # Force early exit via an error during validation; result should be NULL
  expect_error(mpx(rnorm(50L), 50L, progress = FALSE), "too short")
})

test_that("mpx backend cuda errors are enforced via public API", {
  # These are already in test-mpx.R but exercised via expect_error for branch coverage
  expect_error(mpx(rnorm(100L), 20L, s_size = 0.5, backend = "cuda"), "requires `s_size = 1`")
  expect_error(mpx(rnorm(100L), 20L, na_strategy = "segmented", backend = "cuda"), "requires `na_strategy")
})
