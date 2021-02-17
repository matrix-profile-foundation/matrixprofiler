if (!testthat:::on_cran()) {
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
    expect_equal(stamp_res, stomp_res)
    expect_equal(stomp_res, scrimp_res)
    expect_equal(scrimp_res, mpx_res)
    expect_equal(stamp_res_par, stomp_res_par)
    expect_equal(stomp_res_par, scrimp_res_par)
    expect_equal(scrimp_res_par, mpx_res_par)
  })
}
