# test_that("Stomp", {
#   withr::local_options(digits = 22)
#   obj <- helper_objects()

#   stomp_res <- expect_silent(stomp(
#     data = obj$ref_data, window_size = obj$w, exclusion_zone = 0.5,
#     progress = FALSE
#   ))

#   stomp_res %>%
#     expect_type("list") %>%
#     expect_length(4) %>%
#     expect_named(c("matrix_profile", "profile_index", "partial", "ez"))

#   expect_snapshot(stomp_res)
# })

# test_that("Stomp Parallel", {
#   withr::local_options(digits = 22)
#   obj <- helper_objects()

#   stomp_res_par <- expect_silent(stomp(
#     data = obj$ref_data, window_size = obj$w, exclusion_zone = 0.5, n_workers = 2,
#     progress = FALSE
#   ))

#   stomp_res_par %>%
#     expect_type("list") %>%
#     expect_length(4) %>%
#     expect_named(c("matrix_profile", "profile_index", "partial", "ez"))
#   expect_snapshot(stomp_res_par)
# })

# test_that("Stomps are equal", {
#   obj <- helper_objects()
#   stomp_res <- expect_silent(stomp(
#     data = obj$ref_data, window_size = obj$w, exclusion_zone = 0.5,
#     progress = FALSE
#   ))

#   stomp_res_par <- expect_silent(stomp(
#     data = obj$ref_data, window_size = obj$w, exclusion_zone = 0.5, n_workers = 2,
#     progress = FALSE
#   ))

#   expect_equal(stomp_res, stomp_res_par)
# })

# test_that("Left Right Profiles", {
#   obj <- helper_objects()
#   stomp_res <- expect_silent(stomp(
#     data = obj$ref_data, window_size = obj$w, exclusion_zone = 0.5,
#     progress = FALSE, left_right_profile = TRUE
#   ))

#   join <- pmin(stomp_res$left_matrix_profile, stomp_res$right_matrix_profile)

#   expect_equal(stomp_res$matrix_profile, join)
# })
