if (!testthat:::on_cran()) {

  test_that("Conversion between Euclidean Distance and Correlation.", {
    set.seed(2021)
    corr <- runif(1000)

    ed <- corr_ed(corr, 100)
    newcorr <- ed_corr(ed, 100)

    expect_equal(corr, newcorr)
  })

  test_that("General math functions.", {
    set.seed(2021)
    data_mode <- c(round(rnorm(500, mean = 10, sd = 2)), round(rnorm(500, mean = 5, sd = 2.5)))
    data_sd <- std(motifs_discords_small)
    test_mode <- mode(data_mode)
    data_znorm <- znorm(motifs_discords_small)
    data_normalized <- normalize(motifs_discords_small, 1, 5)
    data_norm_zone <- zero_one_norm(motifs_discords_small)
    zeroc <- zero_crossings(motifs_discords_small) # Counts number of zero-crossings
    data_complex <- complexity(motifs_discords_small) # Computes the complexity index of the data
    split <- binary_split(50)
    paa_obj <- paa(motifs_discords_small, 5)
    ipaa_obj <- ipaa(paa_obj, 5)
  })
}
