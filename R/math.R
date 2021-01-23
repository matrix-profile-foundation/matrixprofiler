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
#' @param x a `vector` or a column `matrix` of `numeric`.
#' @param w the window size
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

#' Population SD, as R always calculate with n-1 (sample), here we fix it
#'
#' @inheritParams fast_movsd
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

#' Calculates the mode of a vector
#'
#' @param x
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

#' Normalizes data for mean Zero and Standard Deviation One
#'
#' @inheritParams fast_movsd
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
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param min the minimum value
#' @param max the maximum value
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

#' Normalizes data between Zero and One
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the normalized data.
#' @keywords internal
#' @noRd
#'
zero_one_norm <- function(data) {
  data <- round(data, 10)

  data <- data - min(data[!is.infinite(data) & !is.na(data)])
  data <- data / max(data[!is.infinite(data) & !is.na(data)])

  return(data)
}

#' Counts number of zero-crossings
#'
#' Count the number of zero-crossings from the supplied time-domain input vector. A simple method is
#' applied here that can be easily ported to a real-time system that would minimize the number of
#' if-else conditionals.
#'
#' @param data is the input time-domain signal (one dimensional).
#'
#' @return Returns the amount of zero-crossings in the input signal.
#' @author sparafucile17 06/27/04
#' @references <https://www.dsprelated.com/showcode/179.php>
#' @keywords internal
#' @noRd
#'
zero_crossings <- function(data) {
  # initial value
  count <- 0

  data <- as.matrix(data)

  # error checks
  if (length(data) == 1) {
    stop("Input signal must have more than one element.")
  }

  if ((ncol(data) != 1) && (nrow(data) != 1)) {
    stop("Input must be one-dimensional.")
  }

  # force signal to be a vector oriented in the same direction
  data <- as.vector(data)

  num_samples <- length(data)

  for (i in 2:num_samples) {
    # Any time you multiply to adjacent values that have a sign difference
    # the result will always be negative.  When the signs are identical,
    # the product will always be positive.
    if ((data[i] * data[i - 1]) < 0) {
      count <- count + 1
    }
  }

  return(as.integer(count))
}

#' Computes the complexity index of the data
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the complexity index of the data provided (normally a subset)
#' @keywords internal
#' @noRd
#'
complexity <- function(data) {
  return(sqrt(sum(diff(data)^2)))
}

# #' Distance between two matrices
# #'
# #' Computes the Euclidean distance between rows of two matrices.
# #'
# #' @param x a `matrix`.
# #' @param y a `matrix`.
# #'
# #' @return Returns a `matrix` of size m x n if x is of size m x k and y is of size n x k.
# #' @keywords internal
# #' @noRd

# diff2 <- function(x, y) {
#   # Rcpp ?
#   if (!is.numeric(x) || !is.numeric(y)) {
#     stop("`x` and `y` must be numeric vectors or matrices.")
#   }
#   if (is.vector(x)) {
#     dim(x) <- c(1, length(x))
#   }
#   if (is.vector(y)) {
#     dim(y) <- c(1, length(y))
#   }
#   if (ncol(x) != ncol(y)) {
#     stop("`x` and `y` must have the same number of columns.")
#   }
#   m <- nrow(x)
#   n <- nrow(y)
#   xy <- x %*% t(y)
#   xx <- matrix(rep(apply(x * x, 1, sum), n), m, n, byrow = FALSE)
#   yy <- matrix(rep(apply(y * y, 1, sum), m), m, n, byrow = TRUE)
#   sqrt(pmax(xx + yy - 2 * xy, 0))
# }

#' Binary Split algorithm
#'
#' Creates a vector with the indexes of binary split.
#'
#' @param n size of the vector
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

# #' Bubble up algorithm
# #'
# #' Bubble up algorithm.
# #'
# #' @param data a vector of values
# #' @param len size of data
# #'
# #' @return Doesnt return. Not used for now
# #' @keywords internal
# #' @noRd

# bubble_up <- function(data, len) {
#   # Rcpp ?
#   pos <- len

#   while (pos > 1 && data[pos %/% 2] < data[pos]) {
#     t <- data[pos]
#     data[pos] <- data[pos %/% 2]
#     data[pos %/% 2] <- t
#     pos <- pos %/% 2
#   }

#   return(data)
# }

#' Piecewise Aggregate Approximation of time series
#'
#' @param data time series
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
#' @param data time series
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
