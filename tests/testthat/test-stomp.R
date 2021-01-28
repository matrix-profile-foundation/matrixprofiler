if (!testthat:::on_cran()) {
  stomp_res <- NULL
  stomp_res_par <- NA

  test_that("Stomp", {
    expect_silent(stomp_res <<- stomp(
      data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
      progress = FALSE
    ))
    expect_type(stomp_res, "list")
    expect_snapshot_value(stomp_res, style = "serialize")
  })

  test_that("Stomp Parallel", {
    expect_silent(stomp_res_par <<- stomp(
      data = motifs_discords_small, window_size = 150, exclusion_zone = 0.5,
      n_workers = 2, progress = FALSE
    ))
    expect_type(stomp_res_par, "list")
    # expect_snapshot_value(stomp_res_par, style = "serialize")
  })

  test_that("Stomps are equal", {
    expect_equal(stomp_res, stomp_res_par)
  })
}
