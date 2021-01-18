# nolint start
.onLoad <- function(libname, pkgname) {
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
# nolint end
