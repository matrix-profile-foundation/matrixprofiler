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
  expect_identical(muinvn, muinvn_par)
  expect_snapshot_value(muinvn, style = "serialize")
})
