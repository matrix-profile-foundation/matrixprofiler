#' Matrix Profile Computation
#'
#' STOMP is a faster implementation with the caveat that is not anytime as STAMP or SCRIMP.
#'
#' @details ## stomp
#'  The STOMP uses a faster implementation to compute the Matrix Profile and Profile Index. It can be stopped earlier by
#'  the user, but the result is not considered anytime, just incomplete. For a anytime algorithm, use `stamp()` or
#'  `scrimp()`.
#'
#' @family matrix profile computations

#' @export
#' @rdname mp_algos
#' @order 2
#' @examples
#' mp <- stomp(motifs_discords_small, 50)
stomp <- function(data, window_size, query = NULL, exclusion_zone = 0.5, n_workers = 1, progress = TRUE) {

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
          result <- stomp_rcpp_parallel(
            data,
            data,
            window_size,
            ez,
            as.logical(progress)
          )
          RcppParallel::setThreadOptions(numThreads = p)
        } else {
          result <- stomp_rcpp(
            data,
            data,
            window_size,
            ez,
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
          result <- stomp_rcpp_parallel(
            data,
            query,
            window_size,
            ez,
            as.logical(progress)
          )
          RcppParallel::setThreadOptions(numThreads = p)
        } else {
          result <- stomp_rcpp(
            data,
            query,
            window_size,
            ez,
            as.logical(progress)
          )
        }
      },
      error = print
    )
    "!DEBUG End AB-Join"
  }
}
