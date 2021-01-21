# Math functions ------------------------------------------------------------------------------
#
# Supports:
#               NA/NaN  -Inf/Inf
# std             Unk      Unk
# mode            Unk      Unk
# znorm           Unk      Unk
# diff2           Unk      Unk
# bubble_up       Unk      Unk
# paa             Unk      Unk
# ipaa            Unk      Unk
# min_mp_idx      Unk      Unk


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

#' Distance between two matrices
#'
#' Computes the Euclidean distance between rows of two matrices.
#'
#' @param x a `matrix`.
#' @param y a `matrix`.
#'
#' @return Returns a `matrix` of size m x n if x is of size m x k and y is of size n x k.
#' @keywords internal
#' @noRd

diff2 <- function(x, y) {
  # Rcpp ?
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("`x` and `y` must be numeric vectors or matrices.")
  }
  if (is.vector(x)) {
    dim(x) <- c(1, length(x))
  }
  if (is.vector(y)) {
    dim(y) <- c(1, length(y))
  }
  if (ncol(x) != ncol(y)) {
    stop("`x` and `y` must have the same number of columns.")
  }
  m <- nrow(x)
  n <- nrow(y)
  xy <- x %*% t(y)
  xx <- matrix(rep(apply(x * x, 1, sum), n), m, n, byrow = FALSE)
  yy <- matrix(rep(apply(y * y, 1, sum), m), m, n, byrow = TRUE)
  sqrt(pmax(xx + yy - 2 * xy, 0))
}

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
    return(binary_split_rcpp(as.integer(n)))
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
  return(idxs)
}

#' Bubble up algorithm
#'
#' Bubble up algorithm.
#'
#' @param data a vector of values
#' @param len size of data
#'
#' @return Doesnt return. Not used for now
#' @keywords internal
#' @noRd

bubble_up <- function(data, len) {
  # Rcpp ?
  pos <- len

  while (pos > 1 && data[pos %/% 2] < data[pos]) {
    t <- data[pos]
    data[pos] <- data[pos %/% 2]
    data[pos %/% 2] <- t
    pos <- pos %/% 2
  }

  return(data)
}

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

#' Get index of the minimum value from a matrix profile and its nearest neighbor
#'
#' @param .mp a `MatrixProfile` object.
#' @param n_dim number of dimensions of the matrix profile
#' @param valid check for valid numbers
#'
#' @return returns a `matrix` with two columns: the minimum and the nearest neighbor
#' @export
#'
#' @examples
#' w <- 50
#' data <- mp_gait_data
#' mp <- tsmp(data, window_size = w, exclusion_zone = 1 / 4, verbose = 0)
#' min_val <- min_mp_idx(mp)
min_mp_idx <- function(.mp, n_dim = NULL, valid = TRUE) {
  if (!is.null(n_dim)) {
    .mp$mp <- .mp$mp[, n_dim, drop = FALSE]
    .mp$pi <- .mp$pi[, n_dim, drop = FALSE]
  }

  n_dim <- ncol(.mp$mp)
  mp_size <- nrow(.mp$mp)
  min <- apply(.mp$mp, 2, which.min) # support for multidimensional matrix profile

  if (any(min == 1) && any(is.infinite(.mp$mp[1, (min == 1)]))) {
    return(NA)
  }

  nn_min <- NULL

  for (i in seq_len(n_dim)) {
    nn_min <- c(nn_min, .mp$pi[min[i], i])
  }

  if (valid) {
    if (all(nn_min > 0 & nn_min <= mp_size) &&
      all(!is.infinite(diag(.mp$mp[nn_min, ], names = FALSE)))) {
      return(cbind(min, nn_min, deparse.level = 0))
    }

    for (i in seq_len(n_dim)) {
      .mp$mp[min[i], i] <- Inf
    }

    stop <- FALSE
    while (!stop) {
      min <- apply(.mp$mp, 2, which.min)

      if (any(min == 1) && any(is.infinite(.mp$mp[1, (min == 1)]))) {
        stop <- TRUE
      } else {
        nn_min <- NULL

        for (i in seq_len(n_dim)) {
          nn_min <- c(nn_min, .mp$pi[min[i], i])
        }

        if (all(nn_min > 0 & nn_min <= mp_size) &&
          all(!is.infinite(diag(.mp$mp[nn_min, ], names = FALSE)))) {
          return(cbind(min, nn_min, deparse.level = 0))
        } else {
          for (i in seq_len(n_dim)) {
            .mp$mp[min[i], i] <- Inf
          }
        }
      }
    }

    return(NA)
  } else {
    return(cbind(min, nn_min, deparse.level = 0))
  }
}

# AV Aux functions --------------------------------------------------------------------------------

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

  return(count)
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

# Misc -------------------------------------------------------------------------------------------
#' Set/changes the data included in TSMP object.
#'
#' This may be useful if you want to include the data lately or remove the included data (set as `NULL`).
#'
#' @param .mp a TSMP object.
#' @param data a `matrix` (for one series) or a `list` of matrices (for two series).
#'
#' @return Returns silently the original TSMP object with changed data.
#' @export
#'
#' @examples
#' mp <- tsmp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
#' mp <- set_data(mp, NULL)
set_data <- function(.mp, data) {
  if (!is.null(data)) {
    if (!is.list(data)) {
      data <- list(data)
    }

    data <- lapply(data, as.matrix)

    # data should be this size
    data_size <- (nrow(.mp$mp) + .mp$w - 1)

    for (i in seq_len(length(data))) {
      if (nrow(data[[i]]) != data_size) {
        warning("Warning: data size is ", nrow(data[[i]]), ", but should be ", data_size, " for this matrix profile.")
      }
    }
  }

  .mp$data <- data

  invisible(.mp)
}

#' Get the data included in a TSMP object, if any.
#'
#' @param .mp a TSMP object.
#'
#' @return Returns the data as `matrix`. If there is more than one series, returns a `list`.
#' @export
#'
#' @examples
#' mp <- tsmp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
#' get_data(mp)
get_data <- function(.mp) {
  if (length(.mp$data) == 1) {
    return(invisible(as.matrix(.mp$data[[1]])))
  } else {
    return(invisible(as.matrix(.mp$data)))
  }
}

#' Add class on front or move it to front if already exists
#'
#' @param classes result from `class()`
#' @param new_class string with the new class to add or update (put on front).
#'
#' @return Returns a vector with classes names.
#' @keywords internal
#' @noRd
#'

update_class <- function(classes, new_class) {
  classes <- classes[!(classes == new_class)]
  classes <- c(new_class, classes)

  return(classes)
}

#' Remove a `TSMP` class from an object
#'
#' @param x a `TSMP` object
#' @param class `character` string with the class name
#'
#' @return the object without the class
#' @export
#'
#' @examples
#' w <- 50
#' data <- mp_gait_data
#' mp <- tsmp(data, window_size = w, exclusion_zone = 1 / 4, verbose = 0)
#' mp <- find_chains(mp)
#' # Remove the "Chain" class information
#' mp <- remove_class(mp, "Chain")
remove_class <- function(x, class) {
  switch(class,
    "Chain" = {
      x$chain <- NULL
    },
    "Discord" = {
      x$discord <- NULL
    },
    "Motif" = {
      x$motif <- NULL
    },
    "MultiMotif" = {
      x$motif <- NULL
    },
    "AnnotationVector" = {
      x$av <- NULL
    },
    "ArcCount" = {
      x$cac <- NULL
    },
    "Fluss" = {
      x$fluss <- NULL
    },
    "Salient" = {
      x$salient <- NULL
    }
  )

  class(x) <- class(x)[class(x) != class]

  return(x)
}

# nolint start

#' Convert a TSMP object into another if possible
#'
#' The base Classes are `MatrixProfile` and `MultiMatrixProfile`, but as other functions are used,
#' classes are pushed behind, since the last output normally is the most significant. If you want,
#' for example, to plot the Matrix Profile from a `Fluss` object, you may use `as.matrixprofile()`
#' to cast it back.
#'
#' @param .mp a TSMP object.
#'
#' @return Returns the object with the new class, if possible.
#'
#' @describeIn as.matrixprofile Cast an object changed by another function back to `MatrixProfile`.
#' @export
#' @examples
#'
#' w <- 50
#' data <- mp_gait_data
#' mp <- tsmp(data, window_size = w, exclusion_zone = 1 / 4, verbose = 0)
#' mp <- find_motif(mp)
#' class(mp) # first class will be "Motif"
#'
#' plot(mp) # plots a motif plot
#'
#' plot(as.matrixprofile(mp)) # plots a matrix profile plot
as.matrixprofile <- function(.mp) {
  if (!("MatrixProfile" %in% class(.mp))) {
    stop("This object cannot be a `MatrixProfile`.")
  }

  class(.mp) <- update_class(class(.mp), "MatrixProfile")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMatrixProfile`.
#' @export
#'

as.multimatrixprofile <- function(.mp) {
  if (!("MultiMatrixProfile" %in% class(.mp))) {
    stop("This object cannot be a `MultiMatrixProfile`.")
  }

  class(.mp) <- update_class(class(.mp), "MultiMatrixProfile")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `PMP`.
#' @export
#'

as.pmp <- function(.mp) {
  if (!("PMP" %in% class(.mp))) {
    stop("This object cannot be a `PMP`.")
  }

  class(.mp) <- update_class(class(.mp), "PMP")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMatrixProfile`.
#' @export
#'

as.valmod <- function(.mp) {
  if (!("Valmod" %in% class(.mp))) {
    stop("This object cannot be a `Valmod`.")
  }

  class(.mp) <- update_class(class(.mp), "Valmod")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `Fluss`.
#' @export

as.fluss <- function(.mp) {
  if (!("Fluss" %in% class(.mp))) {
    stop("This object cannot be a `Fluss`.")
  }

  class(.mp) <- update_class(class(.mp), "Fluss")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `Chain`.
#' @export

as.chain <- function(.mp) {
  if (!("Chain" %in% class(.mp))) {
    stop("This object cannot be a `Chain`.")
  }

  class(.mp) <- update_class(class(.mp), "Chain")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `Discord`.
#' @export

as.discord <- function(.mp) {
  if (!("Discord" %in% class(.mp))) {
    stop("This object cannot be a `Discord`.")
  }

  class(.mp) <- update_class(class(.mp), "Discord")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `Motif`.
#' @export

as.motif <- function(.mp) {
  if (!("Motif" %in% class(.mp))) {
    stop("This object cannot be a `Motif`.")
  }

  class(.mp) <- update_class(class(.mp), "Motif")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMotif`.
#' @export

as.multimotif <- function(.mp) {
  if (!("MultiMotif" %in% class(.mp))) {
    stop("This object cannot be a `MultiMotif`.")
  }

  class(.mp) <- update_class(class(.mp), "MultiMotif")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `ArcCount`.
#' @export

as.arccount <- function(.mp) {
  if (!("ArcCount" %in% class(.mp))) {
    stop("This object cannot be a `ArcCount`.")
  }

  class(.mp) <- update_class(class(.mp), "ArcCount")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `Salient`.
#' @export

as.salient <- function(.mp) {
  if (!("Salient" %in% class(.mp))) {
    stop("This object cannot be a `Salient`.")
  }

  class(.mp) <- update_class(class(.mp), "Salient")

  return(.mp)
}

# nolint end
