if (!testthat:::on_cran()) {
  test_that("Mass Pre normalized", {
    expect_silent({
      pre_obj_norm <- mass_pre(data = motifs_discords_small, type = "normalized", window_size = 150)
    })
    expect_silent({
      res_norm <- mass(pre_obj_norm, data = motifs_discords_small, index = 100)
    })
    expect_type(res_norm, "list")
    # expect_snapshot_value(res_norm, style = "serialize")
  })

  # TODO: check abs and weighted output

  # test_that("Mass Pre absolute", {
  #   expect_silent({
  #     pre_obj_abs <- mass_pre(data = motifs_discords_small, type = "absolute", window_size = 150)
  #   })
  #   expect_silent({
  #     res_abs <- mass(pre_obj_abs, data = motifs_discords_small, index = 100)
  #   })
  #   expect_type(res_abs, "list")
  #   expect_snapshot_value(res_abs, style = "serialize")
  # })
  #
  # test_that("Mass Pre weighted", {
  #   expect_silent({
  #     pre_obj_weights <- mass_pre(
  #       data = motifs_discords_small, type = "weighted", window_size = 150,
  #       weights = 11:160
  #     )
  #   })
  #   expect_silent({
  #     res_weight <- mass(pre_obj_weights, data = motifs_discords_small, index = 100)
  #   })
  #   expect_type(res_weight, "list")
  #   expect_snapshot_value(res_weight, style = "serialize")
  # })
}
