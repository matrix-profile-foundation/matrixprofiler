# Run each probe in a separate R process: a native crash must not prevent the
# remaining probes from running. Reuse the package installed by R CMD check.
args <- commandArgs(trailingOnly = TRUE)
child <- length(args) == 3L && identical(args[[1L]], "--probe")

if (child) {
  probe <- args[[2L]]
  .libPaths(c(args[[3L]], .libPaths()))
  library(matrixprofiler)

  mark <- function(text) {
    cat(text, "\n", file = stderr())
    flush(stderr())
  }
  options(error = function() {
    traceback(20)
    quit(save = "no", status = 1L, runLast = FALSE)
  })
  mark(paste("START", probe))
  set.seed(2021)
  negative <- motifs_discords_small[1:400]
  positive <- motifs_discords_small[450:750]
  ab <- function(workers, indices = FALSE) {
    mark(paste("ENTER mpxab_rcpp_parallel, workers =", workers, "indices =", indices))
    result <- matrixprofiler:::mpxab_rcpp_parallel(
      positive, negative, 50L, 1, indices, TRUE, FALSE, workers
    )
    mark("RETURN mpxab_rcpp_parallel")
    result
  }
  aa <- function(workers) {
    mark(paste("ENTER mpx_rcpp_parallel, workers =", workers))
    result <- matrixprofiler:::mpx_rcpp_parallel(
      positive, 50L, 0.5, 1, TRUE, TRUE, FALSE, workers
    )
    mark("RETURN mpx_rcpp_parallel")
    result
  }
  cp <- function(workers, positive_matrix = NULL) {
    mark(paste("ENTER contrast, workers =", workers))
    result <- contrast(negative, positive, 50L, positive_matrix = positive_matrix,
                       n_workers = workers, progress = FALSE)
    mark("RETURN contrast")
    result
  }

  switch(probe,
    ab_one_worker = ab(1L),
    ab_two_workers = ab(2L),
    ab_with_indices = ab(2L, TRUE),
    aa_one_worker = aa(1L),
    aa_two_workers = aa(2L),
    contrast_serial = cp(1L),
    contrast_parallel = cp(2L),
    contrast_precomputed = cp(2L, aa(1L)),
    sequence = {
      cp(1L)
      set.seed(2021)
      cp(2L)
    },
    tests_in_order = {
      mark("ENTER algorithms and contrast tests")
      testthat::test_dir("tests/testthat", filter = "^(algorithms|contrast)$",
                         package = "matrixprofiler", load_package = "installed",
                         reporter = "summary", stop_on_failure = TRUE)
      mark("RETURN algorithms and contrast tests")
    },
    stop("Unknown probe: ", probe)
  )
  mark(paste("DONE", probe))
} else {
  lib <- normalizePath(Sys.getenv("MATRIXPROFILER_DIAGNOSTIC_LIB",
                                 "check/matrixprofiler.Rcheck"), mustWork = TRUE)
  stopifnot(dir.exists(file.path(lib, "matrixprofiler")))
  script <- normalizePath(".github/scripts/diagnose-contrast.R", mustWork = TRUE)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  # Preserve the dependency libraries when child R sessions skip .Rprofile.
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  print(sessionInfo())
  probes <- c("ab_one_worker", "ab_two_workers", "ab_with_indices",
              "aa_one_worker", "aa_two_workers", "contrast_serial",
              "contrast_parallel", "contrast_precomputed", "sequence", "tests_in_order")
  statuses <- setNames(integer(length(probes)), probes)
  for (probe in probes) {
    cat("\n::group::Contrast diagnostic: ", probe, "\n", sep = "")
    output <- tempfile(fileext = ".log")
    statuses[[probe]] <- system2(
      rscript, c("--vanilla", shQuote(script), "--probe", probe, shQuote(lib)),
      stdout = output, stderr = output, timeout = 60
    )
    cat(readLines(output, warn = FALSE), sep = "\n")
    cat("\nProcess exit status: ", statuses[[probe]], "\n::endgroup::\n", sep = "")
    unlink(output)
  }
  print(data.frame(probe = names(statuses), exit_status = unname(statuses)), row.names = FALSE)
  # The original check remains failed, even when every isolated probe passes.
}
