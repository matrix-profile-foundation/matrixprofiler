test_that("Mass normalized", {
  expect_silent({
    pre_obj_norm <- mass_pre(data = motifs_discords_small, type = "normalized", window_size = 100)
  })
  expect_silent({
    res_norm <- mass(pre_obj_norm, data = motifs_discords_small, index = 100)
  })
  expect_type(res_norm, "list")
  expect_snapshot_value(round(res_norm$distance_profile, 8), style = "json2")
  expect_snapshot_value(round(res_norm$last_product, 8), style = "json2")
})

test_that("Mass Non-normalized", {
  expect_silent({
    pre_obj_non <- mass_pre(data = motifs_discords_small, type = "non_normalized", window_size = 100)
  })
  expect_silent({
    res_non <- mass(pre_obj_non, data = motifs_discords_small, index = 100)
  })
  expect_type(res_non, "list")
  expect_snapshot_value(round(res_non$distance_profile, 8), style = "json2")
  expect_snapshot_value(round(res_non$last_product, 8), style = "json2")
})

test_that("Mass Absolute", {
  expect_silent({
    pre_obj_abs <- mass_pre(data = motifs_discords_small, type = "absolute", window_size = 100)
  })
  expect_silent({
    res_abs <- mass(pre_obj_abs, data = motifs_discords_small, index = 100)
  })
  expect_type(res_abs, "list")
  expect_snapshot_value(round(res_abs$distance_profile, 8), style = "json2")
  expect_snapshot_value(round(res_abs$last_product, 8), style = "json2")
})

test_that("Mass Weighted", {
  expect_silent({
    pre_obj_weights <- mass_pre(
      data = motifs_discords_small, type = "weighted", window_size = 100,
      weights = 11:110
    )
  })
  expect_silent({
    res_weight <- mass(pre_obj_weights, data = motifs_discords_small, index = 100)
  })
  expect_type(res_weight, "list")
  expect_snapshot_value(round(res_weight$distance_profile, 8), style = "json2")
  expect_snapshot_value(round(res_weight$last_product, 8), style = "json2")
})

test_that("Mass parallel execution matches sequential across all variants", {
  # 1. Normalized v2
  pre_norm <- mass_pre(data = motifs_discords_small, type = "normalized", window_size = 100)
  res_v2_1 <- mass(pre_norm, data = motifs_discords_small, index = 100, version = "v2", n_workers = 1)
  res_v2_2 <- mass(pre_norm, data = motifs_discords_small, index = 100, version = "v2", n_workers = 2)
  expect_equal(res_v2_1$distance_profile, res_v2_2$distance_profile)
  expect_equal(res_v2_1$last_product, res_v2_2$last_product)

  # 2. Non-normalized
  pre_non <- mass_pre(data = motifs_discords_small, type = "non_normalized", window_size = 100)
  res_non_1 <- mass(pre_non, data = motifs_discords_small, index = 100, n_workers = 1)
  res_non_2 <- mass(pre_non, data = motifs_discords_small, index = 100, n_workers = 2)
  expect_equal(res_non_1$distance_profile, res_non_2$distance_profile)
  expect_equal(res_non_1$last_product, res_non_2$last_product)

  # 3. Absolute
  pre_abs <- mass_pre(data = motifs_discords_small, type = "absolute", window_size = 100)
  res_abs_1 <- mass(pre_abs, data = motifs_discords_small, index = 100, n_workers = 1)
  res_abs_2 <- mass(pre_abs, data = motifs_discords_small, index = 100, n_workers = 2)
  expect_equal(res_abs_1$distance_profile, res_abs_2$distance_profile)
  expect_equal(res_abs_1$last_product, res_abs_2$last_product)

  # 4. Weighted
  pre_w <- mass_pre(data = motifs_discords_small, type = "weighted", window_size = 100, weights = 11:110)
  res_w_1 <- mass(pre_w, data = motifs_discords_small, index = 100, n_workers = 1)
  res_w_2 <- mass(pre_w, data = motifs_discords_small, index = 100, n_workers = 2)
  expect_equal(res_w_1$distance_profile, res_w_2$distance_profile)
  expect_equal(res_w_1$last_product, res_w_2$last_product)
})
