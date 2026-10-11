set.seed(2021)
stamp_res <- stamp(
  data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
  progress = FALSE
)
set.seed(2021)
stamp_res_par <- stamp(
  data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
  n_workers = 2, progress = FALSE
)
set.seed(2021)
stomp_res <- stomp(
  data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
  progress = FALSE
)
set.seed(2021)
stomp_res_par <- stomp(
  data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
  n_workers = 2, progress = FALSE
)
set.seed(2021)
scrimp_res <- scrimp(
  data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
  progress = FALSE
)
set.seed(2021)
scrimp_res_par <- scrimp(
  data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5, pre_scrimp = 0,
  n_workers = 2, progress = FALSE
)
set.seed(2021)
mpx_res <- mpx(
  data = motifs_discords_small,
  window_size = 150, exclusion_zone = 0.5,
  idxs = TRUE, distance = "euclidean", progress = FALSE
)
set.seed(2021)
mpx_res_par <- mpx(
  data = motifs_discords_small,
  window_size = 150, exclusion_zone = 0.5,
  idxs = TRUE, distance = "euclidean", n_workers = 2, progress = FALSE
)

test_that("Algorithms are equal", {
  expect_snapshot(stamp_res)
  expect_equal(stamp_res, stomp_res)
  expect_equal(stomp_res, scrimp_res)
  expect_equal(scrimp_res, mpx_res)
  expect_equal(stamp_res_par, stomp_res_par)
  expect_equal(stomp_res_par, scrimp_res_par)
  expect_equal(scrimp_res_par, mpx_res_par)
})

test_that("has_cuda returns logical", {
  expect_type(has_cuda(), "logical")
  expect_length(has_cuda(), 1L)
})

test_that("Relaxed window_size allows AB-join with window_size > length(query) / 2", {
  set.seed(42)
  data <- rnorm(200)
  query <- rnorm(40) # half size would be 20, but window is 30
  w <- 30L

  res_mpx <- expect_silent(mpx(data, w, query = query, progress = FALSE))
  res_stomp <- expect_silent(stomp(data, w, query = query, progress = FALSE))
  res_stamp <- expect_silent(stamp(data, w, query = query, progress = FALSE))

  expect_length(res_mpx$matrix_profile, length(data) - w + 1L)
  expect_length(res_stomp$matrix_profile, length(data) - w + 1L)
  expect_length(res_stamp$matrix_profile, length(data) - w + 1L)
})

test_that("Relaxed window_size allows self-join with window_size > length(data) / 2 when exclusion zone allows", {
  set.seed(42)
  data <- rnorm(100) # half size is 50, window is 60 with exclusion_zone = 0.2 (ez = 12, leaving 40 > 12)
  w <- 60L

  res_mpx <- expect_silent(mpx(data, w, exclusion_zone = 0.2, progress = FALSE))
  res_stomp <- expect_silent(stomp(data, w, exclusion_zone = 0.2, progress = FALSE))
  res_stamp <- expect_silent(stamp(data, w, exclusion_zone = 0.2, progress = FALSE))

  expect_length(res_mpx$matrix_profile, length(data) - w + 1L)
  expect_length(res_stomp$matrix_profile, length(data) - w + 1L)
  expect_length(res_stamp$matrix_profile, length(data) - w + 1L)
})

test_that("Algorithms error informatively when window_size is truly impossible", {
  set.seed(42)
  data <- rnorm(50)

  # window_size >= length(data)
  expect_error(mpx(data, 50L, progress = FALSE), "Time series is too short")
  expect_error(stomp(data, 50L, progress = FALSE), "Time series is too short")
  expect_error(stamp(data, 50L, progress = FALSE), "Time series is too short")
  expect_error(scrimp(data, 50L, progress = FALSE), "Time series is too short")

  # window_size leaves no pairs outside exclusion zone
  expect_error(mpx(data, 45L, exclusion_zone = 0.5, progress = FALSE), "exclusion zone")
  expect_error(stomp(data, 45L, exclusion_zone = 0.5, progress = FALSE), "exclusion zone")
  expect_error(stamp(data, 45L, exclusion_zone = 0.5, progress = FALSE), "exclusion zone")
  expect_error(scrimp(data, 45L, exclusion_zone = 0.5, progress = FALSE), "exclusion zone")
})

