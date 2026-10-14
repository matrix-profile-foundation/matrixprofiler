test_that("mov_mean", {
  expect_silent(mov_mean_ogita <- mov_mean(
    data = motifs_discords_small, window_size = 150, type = "ogita"
  ))
  expect_snapshot_value(mov_mean_ogita, style = "serialize")

  expect_silent(mov_mean_normal <- mov_mean(
    data = motifs_discords_small, window_size = 150, type = "normal"
  ))
  expect_snapshot_value(mov_mean_normal, style = "serialize")

  expect_silent(mov_mean_weighted <- mov_mean(
    data = motifs_discords_small, window_size = 150, type = "weighted"
  ))
  expect_snapshot_value(round(mov_mean_weighted, 8), style = "json2")

  expect_silent(mov_mean_fading <- mov_mean(
    data = motifs_discords_small, window_size = 150, type = "fading"
  ))
  expect_snapshot_value(mov_mean_fading, style = "serialize")

  expect_equal(mov_mean_ogita, mov_mean_normal)
})

test_that("mov_var", {
  expect_silent(mov_var_ogita <- mov_var(
    data = motifs_discords_small, window_size = 150, type = "ogita"
  ))
  expect_snapshot_value(mov_var_ogita, style = "serialize")

  expect_silent(mov_var_normal <- mov_var(
    data = motifs_discords_small, window_size = 150, type = "normal"
  ))
  expect_snapshot_value(mov_var_normal, style = "serialize")

  expect_silent(mov_var_weighted <- mov_var(
    data = motifs_discords_small, window_size = 150, type = "weighted"
  ))
  expect_snapshot_value(round(mov_var_weighted, 8), style = "json2")

  expect_silent(mov_var_fading <- mov_var(
    data = motifs_discords_small, window_size = 150, type = "fading"
  ))
  expect_snapshot_value(mov_var_fading, style = "serialize")

  expect_equal(mov_var_ogita, mov_var_normal)
})

test_that("mov_sum", {
  expect_silent(mov_sum_ogita <- mov_sum(
    data = motifs_discords_small, window_size = 150, type = "ogita"
  ))
  expect_snapshot_value(mov_sum_ogita, style = "serialize")

  expect_silent(mov_sum_normal <- mov_sum(
    data = motifs_discords_small, window_size = 150, type = "normal"
  ))
  expect_snapshot_value(mov_sum_normal, style = "serialize")

  expect_silent(mov_sum_weighted <- mov_sum(
    data = motifs_discords_small, window_size = 150, type = "weighted"
  ))
  expect_snapshot_value(round(mov_sum_weighted, 8), style = "json2")

  expect_silent(mov_sum_fading <- mov_sum(
    data = motifs_discords_small, window_size = 150, type = "fading"
  ))
  expect_snapshot_value(mov_sum_fading, style = "serialize")

  expect_equal(mov_sum_ogita, mov_sum_normal)
})

test_that("mov_max", {
  expect_silent(mov_max_res <- mov_max(
    data = motifs_discords_small, window_size = 150
  ))
  expect_snapshot_value(mov_max_res, style = "serialize")
})

test_that("mov_min", {
  expect_silent(mov_min_res <- mov_min(
    data = motifs_discords_small, window_size = 150
  ))
  expect_snapshot_value(mov_min_res, style = "serialize")
})

test_that("mov_std", {
  expect_silent(mov_std_res <- mov_std(
    data = motifs_discords_small, window_size = 150
  ))
  expect_silent(mov_std_res_r <- mov_std(
    data = motifs_discords_small, window_size = 150, rcpp = FALSE
  ))
  expect_equal(mov_std_res, mov_std_res_r)
})

test_that("movmean_std", {
  expect_silent(movmean_std <- movmean_std(
    data = motifs_discords_small, window_size = 150
  ))
  expect_silent(movmean_std_r <- movmean_std(
    data = motifs_discords_small, window_size = 150, rcpp = FALSE
  ))
  expect_equal(movmean_std, movmean_std_r)
})

test_that("muinvn", {
  expect_silent(muinvn <- muinvn(
    data = motifs_discords_small, window_size = 150
  ))
  expect_silent(muinvn_par <- muinvn(
    data = motifs_discords_small, window_size = 150, n_workers = 2
  ))
  # Serial Rcpp sugar and parallel scalar arithmetic can round differently.
  expect_equal(muinvn, muinvn_par, tolerance = 1e-12)
  expect_snapshot_value(muinvn, style = "serialize")
})

test_that("na_strategy handling in window functions", {
  x <- c(1, 2, 3, NA, 5, 6, 7, 8, 9, 10)
  w <- 3L

  # mov_mean
  m_win <- mov_mean(x, w, na_strategy = "window_na")
  expect_equal(m_win[1], 2)
  expect_true(is.na(m_win[2]))
  expect_true(is.na(m_win[3]))
  expect_true(is.na(m_win[4]))
  expect_equal(m_win[5], 6)

  m_rm <- mov_mean(x, w, na_strategy = "remove")
  expect_equal(m_rm, mov_mean(x, w, na.rm = TRUE))
  expect_equal(m_rm[2], 2.5)

  expect_error(mov_mean(x, w, na_strategy = "error"))

  # mov_sum
  s_win <- mov_sum(x, w, na_strategy = "window_na")
  expect_equal(s_win[1], 6)
  expect_true(is.na(s_win[2]))
  expect_equal(s_win[5], 18)

  s_rm <- mov_sum(x, w, na.rm = TRUE)
  expect_equal(s_rm[2], 5)

  # mov_min / mov_max
  expect_equal(mov_min(x, w)[1], 1)
  expect_true(is.na(mov_min(x, w)[2]))
  expect_equal(mov_min(x, w, na.rm = TRUE)[2], 2)

  expect_equal(mov_max(x, w)[1], 3)
  expect_true(is.na(mov_max(x, w)[2]))
  expect_equal(mov_max(x, w, na.rm = TRUE)[2], 3)

  # mov_std / movmean_std
  std_win <- mov_std(x, w)
  expect_false(is.na(std_win[1]))
  expect_true(is.na(std_win[2]))
  expect_false(is.na(mov_std(x, w, na.rm = TRUE)[2]))

  mstd <- movmean_std(x, w)
  expect_false(is.na(mstd$avg[1]))
  expect_true(is.na(mstd$avg[2]))

  # muinvn
  mu <- muinvn(x, w)
  expect_false(is.na(mu$avg[1]))
  expect_true(is.na(mu$avg[2]))
  expect_false(is.na(mu$sig[1]))
  expect_true(is.na(mu$sig[2]))
})
