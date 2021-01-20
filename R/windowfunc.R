# Window functions ----------------------------------------------------------------------------

# Supports:
#               NA/NaN  -Inf/Inf  Edge  Rcpp
# movmin          Unk     Unk      No   Yes
# movmax          Unk     Unk      No   Yes
# fast_movsd      No      Unk      No   No
# fast_movavg     No      Unk      No   No
# fast_avg_sd     No      Unk      No   No

#' Fast implementation of moving average
#'
#' This function does not handle NA values
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window_size moving sd window size
#' @param rcpp a `logical`. Uses rcpp implementation.
#'
#' @return Returns a `vector` with the moving average
#' @export
#'
#' @examples
#'
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
      result <- switch(
        type,
        ogita = movsum_ogita_rcpp(data, window_size) / window_size,
        normal = movmean_rcpp(data, window_size),
        weighted = movmean_weighted_rcpp(data, window_size, eps),
        fading = movmean_fading_rcpp(data, window_size, eps)
      )
    },
    error = print
  )
}

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
      result <- switch(
        type,
        ogita = movvar_rcpp(data, window_size),
        normal = movvar2_rcpp(data, window_size),
        weighted = movvar_weighted_rcpp(data, window_size, eps),
        fading = movvar_fading_rcpp(data, window_size, eps)
      )
    },
    error = print
  )
}

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
      result <- switch(
        type,
        ogita = movsum_ogita_rcpp(data, window_size),
        normal = movsum_rcpp(data, window_size),
        weighted = movsum_weighted_rcpp(data, window_size, eps),
        fading = movsum_fading_rcpp(data, window_size, eps)
      )
    },
    error = print
  )
}

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
#' This function does not handle NA values
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window_size moving sd window size
#' @param rcpp a `logical`. Uses rcpp implementation.
#'
#' @return Returns a `vector` with the moving standard deviation
#' @export
#'
#' @examples
#'
mov_std <- function(data, window_size, rcpp = TRUE) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  if (rcpp) {
    return(movstd_rcpp(data, window_size))
  }

  # Improve the numerical analysis by subtracting off the series mean
  # this has no effect on the standard deviation.
  data <- data - mean(data)

  data_sum <- cumsum(c(sum(data[1:window_size]), diff(data, window_size)))
  data_mean <- data_sum / window_size

  data2 <- data^2
  data2_sum <- cumsum(c(sum(data2[1:window_size]), diff(data2, window_size)))
  data_sd2 <- (data2_sum / window_size) - (data_mean^2) # variance
  data_sd <- sqrt(data_sd2)

  return(data_sd)
}

#' Fast implementation of moving average and moving standard deviation
#'
#' This function does not handle NA values
#'
#' @inheritParams fast_movsd
#'
#' @return Returns a `list` with `avg` and `sd` `vector`s
#' @export

movmean_std <- function(data, window_size, rcpp = TRUE) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  if (rcpp) {
    return(movmean_std_rcpp(data, window_size))
  }

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

  return(list(avg = mov_mean, sd = data_sd, sig = data_sig, sum = mov_sum, sqrsum = mov2_sum))
}

# DO NOT Handles NA's
muinvn <- function(data, window_size, rcpp = TRUE, n_workers = 1) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  if (rcpp) {
    if (n_workers > 1) {
      p <- RcppParallel::defaultNumThreads()
      n_workers <- min(n_workers, p)
      RcppParallel::setThreadOptions(numThreads = n_workers)
      return(muinvn_rcpp_parallel(data, window_size))
    } else {
      return(muinvn_rcpp(data, window_size))
    }
  }

  data_sum <- mov_sum(data, window_size)
  data_mean <- data_sum / window_size
  data2 <- data^2
  data2_sum <- mov_sum(data2, window_size)
  sig <- 1 / sqrt(data2_sum - data_mean^2 * window_size)

  # std is equals to 1 / (sig * sqrt(w))
  # sig is equals to 1 / (std * sqrt(w))

  return(list(avg = data_mean, sig = sig))
}
