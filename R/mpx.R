#' Matrix Profile Computation
#'
#' MPX is by far the fastest implementation with the caveat that is not anytime as STAMP or SCRIMP.
#'
#' @param idxs (`mpx()` only) A logical. Specifies if the computation will return the Profile Index or not. Defaults to
#'   `TRUE`.
#' @param distance (`mpx()` only) A string. Currently accepts `euclidean` and `pearson`. Defaults to `euclidean`.
#'
#' @details ## mpx
#' This algorithm was developed apart from the main Matrix Profile branch that relies on Fast Fourier Transform (FFT) at
#' least in one part of the process. This algorithm doesn't use FFT at all and is several times faster. It also relies
#' on Ogita's work for better precision computing mean and standard deviation (part of the process).
#'
#' @seealso `mpxab()` for the forward and reverse join-similarity.
#'
#' @details # This document
#' Last updated on `r Sys.Date()` using R version `r getRversion()`.
#'
#' @export
#' @rdname mp_algos
#' @order 4
#' @examples
#' mp <- mpx(motifs_discords_small, 50)
mpx <- function(data, window_size, query = NULL, exclusion_zone = 0.5, idxs = TRUE,
                distance = c("euclidean", "pearson"), n_workers = 1, progress = TRUE) {


  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing Arguments"

  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (!is.null(query)) {
    query <- as.numeric(query)
    checkmate::qassert(query, c("0", "N>=4"))
  }
  checkmate::qassert(exclusion_zone, "N+")
  checkmate::qassert(idxs, "B+")
  distance <- match.arg(distance)
  if (distance == "euclidean") {
    dist <- TRUE
  } else {
    dist <- FALSE
  }
  n_workers <- as.integer(checkmate::qassert(n_workers, "X+"))
  checkmate::qassert(progress, "B+")

  ez <- exclusion_zone
  result <- NULL

  query_size <- ifelse(is.null(query), length(data),
    ifelse(length(data) > length(query), length(query),
      length(data)
    )
  )

  if (window_size > ceiling(query_size / 2)) {
    stop("Time series is too short relative to desired window size.", call. = FALSE)
  }


  # Register anytime exit point ----------------------
  "!DEBUG Register anytime exit point"
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      result$ez <- ez
      return(result)
    },
    TRUE
  )

  # Computation ------------------------------------
  "!DEBUG Computation"
  if (is.null(query)) {
    ## Self-Join ====================================
    "!DEBUG Self-Join"
    tryCatch(
      {
        "!DEBUG n_workers = `n_workers`"
        if (n_workers > 1) {
          p <- RcppParallel::defaultNumThreads()
          n_workers <- min(n_workers, p)
          RcppParallel::setThreadOptions(numThreads = n_workers)
          result <- mpx_rcpp_parallel(
            data,
            window_size,
            ez,
            as.logical(idxs),
            as.logical(dist),
            as.logical(progress)
          )
          RcppParallel::setThreadOptions(numThreads = p)
        } else {
          result <- mpx_rcpp(
            data,
            window_size,
            ez,
            as.logical(idxs),
            as.logical(dist),
            as.logical(progress)
          )
        }
      },
      error = print
    )
    "!DEBUG End Self-Join"
  } else {
    ## AB-Join ====================================
    "!DEBUG AB-Join"
    ez <- 0

    tryCatch(
      {
        "!DEBUG n_workers = `n_workers`"
        if (n_workers > 1) {
          p <- RcppParallel::defaultNumThreads()
          n_workers <- min(n_workers, p)
          RcppParallel::setThreadOptions(numThreads = n_workers)
          result <- mpxab_rcpp_parallel(
            data,
            query,
            window_size,
            as.logical(idxs),
            as.logical(dist),
            as.logical(progress)
          )
          RcppParallel::setThreadOptions(numThreads = p)
        } else {
          result <- mpxab_rcpp(
            data,
            query,
            window_size,
            as.logical(idxs),
            as.logical(dist),
            as.logical(progress)
          )
        }
      },
      error = print
    )
    "!DEBUG End AB-Join"
  }
}
