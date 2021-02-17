# Math functions ------------------------------------------------------------------------------
#
# Supports:
#               NA/NaN  -Inf/Inf
# std             Unk      Unk
# mode            Unk      Unk
# znorm           Unk      Unk
# paa             Unk      Unk
# ipaa            Unk      Unk


#' Converts euclidean distances into correlation values
#'
#' @param x a `vector` of `numeric`.
#' @param w the window size
#' @param rcpp A `logical`. If `TRUE` will use the Rcpp implementation, otherwise will use the R implementation,
#' that may or not be slower.
#'
#' @return Returns the converted values
#'
#' @keywords internal
#' @noRd
ed_corr <- function(x, w, rcpp = TRUE) {
  if (rcpp) {
    return(ed_corr_rcpp(x, w))
  } else {
    return((2 * w - x^2) / (2 * w))
  }
}

#' Converts correlation values into euclidean distances
#'
#' @inheritParams ed_corr
#'
#' @return Returns the converted values
#'
#' @keywords internal
#' @noRd
corr_ed <- function(x, w, rcpp = TRUE) {
  if (rcpp) {
    return(corr_ed_rcpp(x, w))
  } else {
    sqrt(2 * w * (1 - ifelse(x > 1, 1, x)))
  }
}

#' Calculates the mode of a vector
#'
#' @inheritParams ed_corr
#'
#' @return the mode
#' @keywords internal
#' @noRd

mode <- function(x, rcpp = FALSE) {
  x <- as.integer(x)

  # Rcpp is not faster
  if (rcpp) {
    return(mode_rcpp(x))
  }

  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#' Population SD, as R always calculate with n-1 (sample), here we fix it
#'
#' @param data A `vector` of `numeric`.
#' @param na.rm A logical. If `TRUE` remove the `NA` values from the computation.
#' @param rcpp A logical. If `TRUE` will use the Rcpp implementation, otherwise will use the R implementation, that may
#'   or not be slower.
#'
#' @return Returns the corrected standard deviation from sample to population
#' @keywords internal
#' @noRd
#'

std <- function(data, na.rm = FALSE, rcpp = TRUE) { # nolint

  # Rcpp is faster
  if (rcpp) {
    return(std_rcpp(data, na.rm))
  }

  sdx <- stats::sd(data, na.rm)

  if (is.na(sdx)) {
    return(NA)
  }

  return(sqrt((length(data) - 1) / length(data)) * sdx)
}

#' Normalizes data for mean Zero and Standard Deviation One
#'
#' @inheritParams std
#'
#' @return Returns the normalized data
#' @keywords internal
#' @noRd
#'
znorm <- function(data, rcpp = TRUE) {
  # Rcpp is faster
  if (rcpp) {
    return(znorm_rcpp(data))
  }

  data_mean <- mean(data)
  data_dev <- std(data)

  if (is.na(data_dev) || data_dev <= 0.01) {
    return(data - data_mean)
  }
  else {
    (data - data_mean) / (data_dev)
  }
}

#' Normalizes data to be between min and max
#'
#'
#' @param min_lim A number
#' @param max_lim A number
#' @param rcpp A logical. If `TRUE` will use the Rcpp implementation, otherwise will use the R implementation, that may
#'   or not be slower.
#' @param data A `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the normalized data
#' @keywords internal
#' @noRd
#'
normalize <- function(data, min_lim = 0, max_lim = 1, rcpp = FALSE) {
  if (rcpp) {
    na <- sort(which(is.na(data)))
    data <- data[!is.na(data)]
    data <- normalize_rcpp(data, min_lim, max_lim)
    for (n in na) {
      data <- append(data, NA, n - 1)
    }
    return(data)
  }

  min_val <- min(data, na.rm = TRUE)
  max_val <- max(data, na.rm = TRUE)

  a <- (max_lim - min_lim) / (max_val - min_val)
  b <- max_lim - a * max_val
  data <- a * data + b

  data[data < min_lim] <- min_lim
  data[data > max_lim] <- max_lim

  return(data)
}

#' Computes the complexity index of the data
#'
#' @param data a `vector` of `numeric`
#'
#' @return Returns the complexity index of the data provided (normally a subset)
#' @keywords internal
#' @noRd
#'
complexity <- function(data) {
  return(sqrt(sum(diff(data)^2)))
}

#' Binary Split algorithm
#'
#' Creates a vector with the indexes of binary split.
#'
#' @param n size of the vector
#' @param rcpp A logical. If `TRUE` will use the Rcpp implementation, otherwise will use the R implementation, that may
#'   or not be slower.
#'
#' @return Returns a `vector` with the binary split indexes
#' @keywords internal
#' @noRd

binary_split <- function(n, rcpp = TRUE) {
  if (rcpp) {
    res <- binary_split_rcpp(as.integer(n))
    return(as.integer(res))
  }

  if (n < 2) {
    return(1)
  }

  split <- function(lb, ub, m) {
    if (lb == m) {
      l <- NULL
      r <- c(m + 1, ub)
    } else if (ub == m) {
      l <- c(lb, m - 1)
      r <- NULL
    } else {
      l <- c(lb, m - 1)
      r <- c(m + 1, ub)
    }

    return(list(l = l, r = r))
  }

  idxs <- vector(mode = "numeric", length = n)
  intervals <- list()

  idxs[1] <- 1 # We always begin by explore the first integer
  intervals[[1]] <- c(2, n) # After exploring the first integer, we begin splitting the interval 2:n
  i <- 2

  while (length(intervals) > 0) {
    lb <- intervals[[1]][1]
    ub <- intervals[[1]][2]
    mid <- floor((lb + ub) / 2)
    intervals[[1]] <- NULL

    idxs[i] <- mid
    i <- i + 1

    if (lb == ub) {
      next
    } else {
      lr <- split(lb, ub, mid)
      if (!is.null(lr$l)) {
        intervals[[length(intervals) + 1]] <- lr$l
      }
      if (!is.null(lr$r)) {
        intervals[[length(intervals) + 1]] <- lr$r
      }
    }
  }
  return(as.integer(idxs))
}

#' Piecewise Aggregate Approximation of time series
#'
#' @param data a `vector` of `numeric`
#' @param p factor of PAA reduction (2 == half of size)
#'
#' @return PAA result
#' @keywords internal
#' @noRd

paa <- function(data, p) {
  # Rcpp ?
  paa_data <- as.vector(data)
  len <- length(paa_data)

  p <- round(abs(p))
  paa_size <- ceiling(len / p)

  if (len == paa_size) {
    return(data)
  } else {
    if ((len %% p) > 0) {
      pad <- p - (len %% p)
      paa_data <- c(paa_data, rep.int(NA, pad))
    }
    res <- colMeans(matrix(paa_data, nrow = p, byrow = FALSE), na.rm = TRUE)
  }

  if (is.matrix(data)) {
    return(as.matrix(res))
  } else {
    return(res)
  }
}

#' Resample data to the original size, with interpolation
#'
#' @param data a `vector` of `numeric`
#' @param p factor of PAA reduction (2 == half of size)
#'
#' @keywords internal
#' @noRd

ipaa <- function(data, p) {
  # Rcpp ?
  if (is.null(data)) {
    return(NULL)
  }

  paa_data <- as.vector(data)
  paa_size <- length(paa_data)
  size <- paa_size * p

  res <- rep.int(NA, size)

  j <- 1
  for (i in seq_len(size)) {
    if (((i - 1) %% p) == 0) {
      res[i] <- data[j]
      j <- j + 1
    } else {
      res[i] <- res[i - 1]
    }
  }

  if (is.matrix(data)) {
    return(as.matrix(res))
  } else {
    return(res)
  }
}
