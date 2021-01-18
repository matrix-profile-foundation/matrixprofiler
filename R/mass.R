tryCatch(
  {
    "!DEBUG n_workers = `n_workers`"
    if (n_workers > 1) {
      p <- RcppParallel::defaultNumThreads()
      n_workers <- min(n_workers, p)
      RcppParallel::setThreadOptions(numThreads = n_workers)
      result <- scrimp_rcpp_parallel(
        data,
        data,
        window_size,
        ez,
        as.logical(progress)
      )
      RcppParallel::setThreadOptions(numThreads = p)
    } else {
      # result <- scrimp_rcpp(
      #   data,
      #   data,
      #   window_size,
      #   ez,
      #   0.25,
      #   as.logical(progress)
      # )
    }
  },
  error = print
)
