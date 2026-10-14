#' Multivariate STOMP Algorithm
#'
#' Computes the Matrix Profile and Profile Index for Multivariate Time Series.
#'
#' @param data a `matrix` of `numeric`, where each column is a time series. Accepts `vector`, `list` and `data.frame` too.
#' @param window_size an `int` with the size of the sliding window.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is `0.5`).
#' @param must_dim an `int` or `vector` of which dimensions to forcibly include (default is `NULL`).
#' @param exc_dim an `int` or `vector` of which dimensions to exclude (default is `NULL`).
#' @param n_workers an `int`. Number of workers for parallel execution (default is `1`).
#'   For parallel execution across processes, prefer a PSOCK cluster and use
#'   `n_workers = 1` within each worker. See [matrixprofiler] for fork limitations.
#' @param progress a `logical`. If `TRUE` will show a progress bar. (Default is `TRUE`).
#'
#' @return Returns a `list` with the matrix profile `matrix_profile`, profile index `profile_index`,
#'   left and right matrix profiles `left_matrix_profile`, `right_matrix_profile` and profile indices
#'   `left_profile_index`, `right_profile_index`, window size `w`, number of dimensions `n_dim`,
#'   exclusion zone `ez`, must dimensions `must` and excluded dimensions `exc`.
#'
#' @export
#' @family matrix profile computations
#' @references * Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif
#'   Discovery.
#' @references * Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A
#'   New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1-27.
#' @references Website: <https://sites.google.com/view/mstamp/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' mp <- mstomp(motifs_discords_small, 50, progress = FALSE)
#'
mstomp <- function(data, window_size, exclusion_zone = 0.5, must_dim = NULL, exc_dim = NULL,
                   n_workers = 1L, progress = TRUE) {
  # Parse arguments ---------------------------------
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    }
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data_size <- nrow(data)
    n_dim <- ncol(data)
  } else if (is.list(data)) {
    data_size <- length(data[[1]])
    n_dim <- length(data)
    for (i in seq_len(n_dim)) {
      len <- length(data[[i]])
      if (len < data_size) {
        data[[i]] <- c(data[[i]], rep(NA_real_, data_size - len))
      }
    }
    data <- sapply(data, cbind)
  } else if (is.vector(data)) {
    data_size <- length(data)
    n_dim <- 1L
    data <- as.matrix(data)
  } else {
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list.", call. = FALSE)
  }

  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (window_size > data_size / 2) {
    stop("Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("`window_size` must be at least 4.", call. = FALSE)
  }

  checkmate::qassert(exclusion_zone, "N+")
  if (!is.null(must_dim)) {
    must_dim <- as.integer(must_dim)
    if (any(must_dim > n_dim) || any(must_dim < 1)) {
      stop("`must_dim` must be less then the total dimension.", call. = FALSE)
    }
  }
  if (!is.null(exc_dim)) {
    exc_dim <- as.integer(exc_dim)
    if (any(exc_dim > n_dim) || any(exc_dim < 1)) {
      stop("`exc_dim` must be less then the total dimension.", call. = FALSE)
    }
  }
  if (length(intersect(must_dim, exc_dim)) > 0) {
    stop("The same dimension is presented in both the exclusion dimension and must have dimension.", call. = FALSE)
  }

  n_workers <- as.integer(checkmate::qassert(n_workers, "X+"))
  if (n_workers > 1L && .Platform$OS.type != "windows") {
    is_child <- tryCatch(
      get("isChild", envir = asNamespace("parallel"))(),
      error = function(e) FALSE
    )
    if (isTRUE(is_child)) {
      warning(
        "Forked process detected: reducing n_workers to 1. ",
        "Inherited multithreaded library state may still be unsafe. ",
        "Use a PSOCK cluster for parallel execution across processes.",
        call. = FALSE
      )
      n_workers <- 1L
    }
  }

  checkmate::qassert(progress, "B+")

  result <- mstomp_rcpp(
    data_ref = data,
    window_size = window_size,
    ez = exclusion_zone,
    must_dim = must_dim,
    exc_dim = exc_dim,
    n_workers = n_workers,
    progress = progress
  )

  return(result)
}
