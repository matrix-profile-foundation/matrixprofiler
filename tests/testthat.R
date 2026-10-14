library(testthat)
library(matrixprofiler)

if (identical(Sys.getenv("GITHUB_ACTIONS"), "true")) {
  # Print progress before completion so native crashes can be localized in CI.
  test_check("matrixprofiler", reporter = MultiReporter$new(list(
    SummaryReporter$new(),
    CheckReporter$new()
  )))
} else {
  test_check("matrixprofiler")
}

