#' Matrix Profile Computation
#'
#' SCRIMP is a faster implementation, like STOMP, but has the ability to return anytime results as STAMP.
#'
#' @details ## scrimp
#' The SCRIMP algorithm was the anytime solution for stomp. It is as fast as stomp but allows the user to cancel the
#' computation and get an approximation of the final result. The exact SCRIMP phase shares the NA-aware diagonal
#' recurrence used by MPX. In serial mode, pre-scrimp first seeds the profile from sampled distance profiles; the
#' multithreaded implementation skips this pre-scrimp stage.
#'
#' @family matrix profile computations

#' @export
#' @rdname mp_algos
#' @order 3
#' @examples
#' mp <- scrimp(motifs_discords_small, 50)
scrimp <- function(data, window_size, query = NULL, exclusion_zone = 0.5, s_size = 1.0, pre_scrimp = 0.25,
                   n_workers = 1, progress = TRUE) {
  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing Arguments"

  data <- as.numeric(data)
  checkmate::qassert(data, "n+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (!is.null(query)) {
    query <- as.numeric(query)
    checkmate::qassert(query, c("0", "n>=4"))
  }
  checkmate::qassert(exclusion_zone, "N+")
  checkmate::qassert(s_size, "N1")
  if (s_size < 0 || s_size > 1) {
    stop("`s_size` must be between 0 and 1.", call. = FALSE)
  }
  checkmate::qassert(pre_scrimp, "N1")
  if (pre_scrimp < 0 || pre_scrimp > 1) {
    stop("`pre_scrimp` must be between 0 and 1.", call. = FALSE)
  }
  n_workers <- as.integer(checkmate::qassert(n_workers, "X+"))
  checkmate::qassert(progress, "B+")

  ez <- exclusion_zone
  result <- NULL

  query_size <- ifelse(is.null(query), length(data), min(length(data), length(query)))
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
          on.exit(RcppParallel::setThreadOptions(numThreads = p), add = TRUE)
          n_workers <- min(n_workers, p)
          RcppParallel::setThreadOptions(numThreads = n_workers)
          result <- scrimp_rcpp_parallel(
            data,
            data,
            window_size,
            ez,
            s_size,
            as.logical(progress)
          )
        } else {
          result <- scrimp_rcpp(
            data,
            data,
            window_size,
            ez,
            s_size,
            pre_scrimp,
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
          warning("Parallel `scrimp` AB-join not implemented yet. Using 1 thread.", call. = FALSE)
          result <- scrimpab_rcpp(
            data,
            query,
            window_size,
            s_size,
            as.logical(progress)
          )
        } else {
          result <- scrimpab_rcpp(
            data,
            query,
            window_size,
            s_size,
            as.logical(progress)
          )
        }
      },
      error = print
    )
    "!DEBUG End AB-Join"
  }
}
