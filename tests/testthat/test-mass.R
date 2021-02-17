if (!testthat:::on_cran()) {
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
}
