matrixprofiler_default_options <- list(
  matrixprofiler.verbose = 2,
  matrixprofiler.exclusion_zone = 1 / 2,
  matrixprofiler.workers = 2
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(matrixprofiler_default_options) %in% names(op))
  if (any(toset)) {
    options(matrixprofiler_default_options[toset])
  }

  tryCatch(debugme::debugme(), error = identity)

  invisible()
}

.onAttach <- function(libname, pkgname) {
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
  }

  packageStartupMessage("Welcome to Matrix ProfileR")
}

.onUnload <- function(libpath) {
  unloadNamespace("debugme")
}
