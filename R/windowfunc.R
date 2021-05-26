# Window functions ----------------------------------------------------------------------------

# Supports:
#               NA/NaN  -Inf/Inf  Edge  Rcpp
# mov_mean       No       Unk      No   Yes
# mov_var        No       Unk      No   Yes
# mov_sum        No       Unk      No   Yes
# mov_max        No       Unk      No   Yes
# mov_min        No       Unk      No   Yes
# mov_std        No       Unk      No   Yes
# movmean_std    No       Unk      No   Yes
# muinvn         No       Unk      No   Yes

#' Several moving window functions
#'
#' These functions do not handle NA values
#'
#' @param data A `vector` or a column `matrix` of `numeric`.
#' @param window_size An `integer`. The size of the rolling window.
#' @param type A `string`. Select between several algorithms. Default is `ogita` (See details).
#' @param eps A `numeric`. Used only for fading algorithms (See details), otherwise has no effect.
#'
#' @details
#' Some functions may use different algorithms to compute the results.
#' The available types are:
#' 1. **ogita**: This is the default. It uses the Ogita *et al.*, Accurate Sum, and Dot Product for precision. It is not
#' the fastest algorithm, but the time spent vs. guarantee of precision worth it.
#' 2. **normal**: This uses the `cumsum` method that is faster, but unreliable in some situations (I have to find the
#' references, but is true).
#' 3. **weighted**: This uses Rodrigues P., *et al.* algorithm that uses a weighted window for online purposes. The
#' `eps` argument controls the factor. (The function is not online yet)
#' 4. **fading**: This also uses Rodrigues P., *et al.* algorithm that in this case, uses a fading factor, also for
#' online purposes. he `eps` argument controls the factor. (The function is not online yet)
#'
#' @details
#' Another important detail is that the *standard deviation* we use for all computations is the *population* (i.e.:
#' divided by `n`), not the *sample* (i.e.: divided by `n - 1`). That is why we also provide the internally the
#' `:::std()` function that computes the *population*, differently from `stats::sd()` that is the *sample* kind. Further
#' more, `movmean_std()` shall be used when you need both results in one computation. This is faster than call
#' `mov_mean()` followed by `mov_std()`. Finally, `muinvn()` is kept like that for historical reasons, as it is the
#' function used by `mpx()`. It returns the `sig` (stable inverse centered norm) instead of `std` (`sig` is equals to
#' `1 / (std * sqrt(window_size))`).
#'
#' @return `mov_mean()` returns a `vector` with moving `avg`.
#' @export
#' @rdname windowfunc
#' @order 1
#' @examples
#' mov <- mov_mean(motifs_discords_small, 50)
mov_mean <- function(data, window_size, type = c("ogita", "normal", "weighted", "fading"), eps = 0.90) {
  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing Arguments"
  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }
  type <- match.arg(type)
  checkmate::qassert(eps, "N1")

  # Register anytime exit point ----------------------
  "!DEBUG Register anytime exit point"
  result <- NULL
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      return(result)
    },
    TRUE
  )

  # Computation ------------------------------------
  "!DEBUG Computation"
  tryCatch(
    {
      result <- switch(type,
        ogita = movsum_ogita_rcpp(data, window_size) / window_size,
        normal = movmean_rcpp(data, window_size),
        weighted = movmean_weighted_rcpp(data, window_size, eps),
        fading = movmean_fading_rcpp(data, window_size, eps)
      )
    },
    error = print
  )
}

#' Fast implementation of moving variance
#'
#' @return `mov_var()` returns a `vector` with moving `var`.
#'
#' @export
#' @rdname windowfunc
#' @order 2
#' @examples
#' mov <- mov_var(motifs_discords_small, 50)
mov_var <- function(data, window_size, type = c("ogita", "normal", "weighted", "fading"), eps = 0.90) {
  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing Arguments"
  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }
  type <- match.arg(type)
  checkmate::qassert(eps, "N1")

  # Register anytime exit point ----------------------
  "!DEBUG Register anytime exit point"
  result <- NULL
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      return(result)
    },
    TRUE
  )

  # Computation ------------------------------------
  "!DEBUG Computation"
  tryCatch(
    {
      result <- switch(type,
        ogita = movvar_rcpp(data, window_size),
        normal = movvar2_rcpp(data, window_size),
        weighted = movvar_weighted_rcpp(data, window_size, eps),
        fading = movvar_fading_rcpp(data, window_size, eps)
      )
    },
    error = print
  )
}

#' Fast implementation of moving sum
#'
#' @return `mov_sum()` returns a `vector` with moving `sum`.
#' @export
#' @rdname windowfunc
#' @order 3
#' @examples
#' mov <- mov_sum(motifs_discords_small, 50)
mov_sum <- function(data, window_size, type = c("ogita", "normal", "weighted", "fading"), eps = 0.90) {
  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing Arguments"
  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }
  type <- match.arg(type)
  checkmate::qassert(eps, "N1")

  # Register anytime exit point ----------------------
  "!DEBUG Register anytime exit point"
  result <- NULL
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      return(result)
    },
    TRUE
  )

  # Computation ------------------------------------
  "!DEBUG Computation"
  tryCatch(
    {
      result <- switch(type,
        ogita = movsum_ogita_rcpp(data, window_size),
        normal = movsum_rcpp(data, window_size),
        weighted = movsum_weighted_rcpp(data, window_size, eps),
        fading = movsum_fading_rcpp(data, window_size, eps)
      )
    },
    error = print
  )
}

#' Fast implementation of moving max
#'
#' @return `mov_max()` returns a `vector` with moving `max`.
#'
#' @export
#' @rdname windowfunc
#' @order 4
#' @examples
#' mov <- mov_max(motifs_discords_small, 50)
mov_max <- function(data, window_size) {
  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing Arguments"
  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  # Register anytime exit point ----------------------
  "!DEBUG Register anytime exit point"
  result <- NULL
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      return(result)
    },
    TRUE
  )

  # Computation ------------------------------------
  "!DEBUG Computation"
  tryCatch(
    {
      result <- movmax_rcpp(data, window_size)
    },
    error = print
  )
}

#' Fast implementation of moving min
#'
#' @return `mov_min()` returns a `vector` with moving `min`.
#' @export
#' @rdname windowfunc
#' @order 5
#' @examples
#' mov <- mov_min(motifs_discords_small, 50)
mov_min <- function(data, window_size) {
  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing Arguments"
  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  # Register anytime exit point ----------------------
  "!DEBUG Register anytime exit point"
  result <- NULL
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      return(result)
    },
    TRUE
  )

  # Computation ------------------------------------
  "!DEBUG Computation"
  tryCatch(
    {
      result <- movmin_rcpp(data, window_size)
    },
    error = print
  )
}


#' Fast implementation of moving standard deviation
#'
#' @param rcpp A `logical`. If `TRUE` will use the Rcpp implementation, otherwise will use the R implementation,
#' that may or not be slower.
#'
#' @return `mov_std()` returns a `vector` with moving `sd`.
#' @export
#' @rdname windowfunc
#' @order 6
#' @examples
#' mov <- mov_std(motifs_discords_small, 50)
mov_std <- function(data, window_size, rcpp = TRUE) {
  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing Arguments"
  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  # Computation ------------------------------------
  "!DEBUG Computation"
  tryCatch(
    {
      if (rcpp) {
        data_sd <- movstd_rcpp(data, window_size)
      } else {
        # Improve the numerical analysis by subtracting off the series mean
        # this has no effect on the standard deviation.
        data <- data - mean(data)

        data_sum <- cumsum(c(sum(data[1:window_size]), diff(data, window_size)))
        data_mean <- data_sum / window_size

        data2 <- data^2
        data2_sum <- cumsum(c(sum(data2[1:window_size]), diff(data2, window_size)))
        data_sd2 <- (data2_sum / window_size) - (data_mean^2) # variance
        data_sd <- sqrt(data_sd2)
      }
    },
    error = print
  )

  return(data_sd)
}

#' Fast implementation of moving average and moving standard deviation
#'
#' @return `movmean_std()` returns a `list` with `vectors` of the moving `avg`, `sd`, `sig`, `sum` and `sqrsum`.
#' @export
#' @rdname windowfunc
#' @order 7
#' @examples
#' mov <- movmean_std(motifs_discords_small, 50)
movmean_std <- function(data, window_size, rcpp = TRUE) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  if (rcpp) {
    result <- movmean_std_rcpp(data, window_size)
  } else {
    mov_sum <- cumsum(c(sum(data[1:window_size]), diff(data, window_size)))
    data2 <- data^2
    mov2_sum <- cumsum(c(sum(data2[1:window_size]), diff(data2, window_size)))
    mov_mean <- mov_sum / window_size
    # Improve the numerical analysis by subtracting off the series mean
    # this has no effect on the standard deviation.
    dmean <- mean(data)
    data <- data - dmean

    data_sum <- cumsum(c(sum(data[1:window_size]), diff(data, window_size)))
    data_mean <- data_sum / window_size
    data2 <- data^2
    data2_sum <- cumsum(c(sum(data2[1:window_size]), diff(data2, window_size)))
    data_sd2 <- (data2_sum / window_size) - (data_mean^2) # variance
    data_sd2[data_sd2 < 0] <- 0
    data_sd <- sqrt(data_sd2) # std deviation
    data_sig <- sqrt(1 / (data_sd2 * window_size))

    result <- list(avg = mov_mean, sd = data_sd, sig = data_sig, sum = mov_sum, sqrsum = mov2_sum)
  }

  return(result)
}

#' Fast implementation of moving average and moving sigma
#'
#' @param n_workers An `integer`. The number of threads using for computing. Defaults to `1`.
#'
#' @return `muinvn()` returns a `list` with `vectors` of moving `avg` and `sig`.
#' @export
#' @rdname windowfunc
#' @order 8
#' @examples
#' mov <- muinvn(motifs_discords_small, 50)
muinvn <- function(data, window_size, n_workers = 1) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }
  if (n_workers > 1) {
    p <- RcppParallel::defaultNumThreads()
    n_workers <- min(n_workers, p)
    RcppParallel::setThreadOptions(numThreads = n_workers)
    result <- muinvn_rcpp_parallel(data, window_size)
  } else {
    result <- muinvn_rcpp(data, window_size)
  }

  return(result)
}

#' Computes the number of times the data crossed the 'zero' line inside a rolling window
#'
#' @return `zero_crossing()` returns a `vector` of times the data crossed the 'zero' line inside a rolling window.
#' @export
#' @rdname windowfunc
#' @order 9
#' @examples
#' zero_cross <- zero_crossing(motifs_discords_small, 50)
zero_crossing <- function(data, window_size) {
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  result <- zero_crossing_rcpp(data, window_size)
  return(result)
}
