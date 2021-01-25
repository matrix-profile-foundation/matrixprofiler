
mass_pre <- function(data, window_size, query = NULL, type = c("normalized", "absolute", "weighted"),
                     weights = NULL) {
  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing Arguments"
  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (!is.null(query)) {
    query <- as.numeric(query)
    checkmate::qassert(query, c("0", "N>=4"))
  }
  type <- match.arg(type)

  query_size <- ifelse(is.null(query), length(data),
    ifelse(length(data) > length(query), length(query),
      length(data)
    )
  )

  if (type == "weighted") {
    if (is.null(weights)) {
      stop("The `weights` argument must be provided.", call. = FALSE)
    }
    if (length(weights) != window_size) {
      stop("The `weights` must be the same size as the `window_size`.", call. = FALSE)
    }
    checkmate::qassert(weights, "N+")
  }

  if (window_size > ceiling(query_size / 2)) {
    stop("Time series is too short relative to desired window size.", call. = FALSE)
  }

  # Register anytime exit point ----------------------
  result <- NULL
  "!DEBUG Register anytime exit point"
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      result$type <- type
      return(result)
    },
    TRUE
  )

  if (is.null(query)) {
    query <- data
  }

  # Computation ------------------------------------
  "!DEBUG Computation"
  tryCatch(
    {
      result <- switch(
        type,
        normalized = mass_pre_rcpp(data, query, window_size),
        absolute = mass_pre_abs_rcpp(data, query, window_size),
        weighted = mass_pre_weighted_rcpp(data, query, window_size, weights)
      )
    },
    error = print
  )
}

mass <- function(pre_obj, data, query = data, index = 1, version = c("v3", "v2"), n_workers = 1) {
  checkmate::qassert(pre_obj, "L+")
  type <- match.arg(pre_obj$type, c("normalized", "absolute", "weighted"))
  version <- match.arg(version)
  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  query <- as.numeric(query)
  checkmate::qassert(query, "N+")
  index <- as.integer(checkmate::qassert(index, "X+"))

  if (length(data) != pre_obj$data_size) {
    stop("Argument `data` is not the same as computed in `pre_obj`.", call. = FALSE)
  }

  if (length(query) != (
    ifelse(pre_obj$type == "absolute", length(pre_obj$sumy2), length(pre_obj$query_mean)) +
      pre_obj$window_size - 1)) {
    stop("Argument `query` is not the same as computed in `pre_obj`.", call. = FALSE)
  }

  # Register anytime exit point ----------------------
  result <- NULL
  "!DEBUG Register anytime exit point"
  on.exit(
    {
      RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())

      if (is.null(result)) {
        return(invisible(NULL))
      } else {
        return(result)
      }
    },
    TRUE
  )

  # Computation ------------------------------------
  "!DEBUG Computation"
  tryCatch(
    {
      query_window <- query[index:(index + pre_obj$window_size - 1)]

      result <- switch(
        type,
        normalized = if (version == "v3") {
          if (n_workers > 1) {
            n_workers <- min(n_workers, RcppParallel::defaultNumThreads())
            RcppParallel::setThreadOptions(numThreads = n_workers)
            mass3_rcpp_parallel(query_window, data, as.integer(pre_obj$data_size),
              as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd, pre_obj$query_mean[index],
              pre_obj$query_sd[index],
              k = 4096
            )
          }
          else {
            mass3_rcpp(query_window, data, as.integer(pre_obj$data_size),
              as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd, pre_obj$query_mean[index],
              pre_obj$query_sd[index],
              k = 4096
            )
          }
        } else {
          mass2_rcpp(
            pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
            as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd, pre_obj$query_mean[index],
            pre_obj$query_sd[index]
          )
        },
        absolute = mass_absolute_rcpp(
          pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
          as.integer(pre_obj$window_size), pre_obj$sumx2, pre_obj$sumy2[index]
        ),
        weighted = mass_weighted_rcpp(
          pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
          as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd,
          pre_obj$query_mean[index], pre_obj$query_sd[index], pre_obj$data_pre, pre_obj$weight
        )
      )
    },
    error = print
  )
}
