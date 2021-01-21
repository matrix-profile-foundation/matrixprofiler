if (!testthat:::on_cran()) {
  weights <- 11:210
  pre_obj_norm <- NULL
  pre_obj_abs <- NULL
  pre_obj_weights <- NULL

  test_that("Mass Pre normalized", {
    expect_silent(
      pre_obj_norm <<- mass_pre(data = motifs_discords_small, type = "normalized", window_size = 200)
    )
  })

  test_that("Mass Pre absolute", {
    expect_silent(
      pre_obj_abs <<- mass_pre(data = motifs_discords_small, type = "absolute", window_size = 200)
    )
  })

  test_that("Mass Pre weighted", {
    expect_silent(
      pre_obj_weights <<- mass_pre(
        data = motifs_discords_small, type = "weighted", window_size = 200,
        weights = weights
      )
    )
  })

  test_that("Mass normalized", {
    expect_silent(
      mass(pre_obj_norm, data = motifs_discords_small)
    )
  })

  test_that("Mass absolute", {
    expect_silent(
      mass(pre_obj_abs, data = motifs_discords_small)
    )
  })

  test_that("Mass weighted", {
    expect_silent(
      mass(pre_obj_weights, data = motifs_discords_small)
    )
  })
}
