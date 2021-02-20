#' Matrix Profile Computation
#'
#' SCRIMP is a faster implementation, like STOMP, but has the ability to return anytime results as STAMP.
#'
#' @details ## scrimp
#' The SCRIMP algorithm was the anytime solution for stomp. It is as fast as stomp but allows the user to cancel the
#' computation and get an approximation of the final result. This implementation uses the SCRIMP++ code. This means
#' that, at first, it will compute the pre-scrimp (a very fast and good approximation), and continue improving with
#' scrimp. The exception is if you use multithreading, that skips the pre-scrimp stage.
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
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (!is.null(query)) {
    query <- as.numeric(query)
    checkmate::qassert(query, c("0", "N>=4"))

    query_size <- length(query)
    data_size <- length(data)

    if (query_size != data_size) {
      stop("Data and query must have the same size using the `scrimp` algorithm. (for now)")
    }

    if (window_size > ceiling(query_size / 2)) {
      stop("Time series is too short relative to desired window size.", call. = FALSE)
    }
  }
  checkmate::qassert(exclusion_zone, "N+")
  n_workers <- as.integer(checkmate::qassert(n_workers, "X+"))
  checkmate::qassert(progress, "B+")

  ez <- exclusion_zone
  result <- NULL

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
          result <- scrimp_rcpp_parallel(
            data,
            data,
            window_size,
            ez,
            s_size,
            as.logical(progress)
          )
          RcppParallel::setThreadOptions(numThreads = p)
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
          warnings("Parallel `scrimp` AB-join not implemented yet. Using 1 thread.")
          # p <- RcppParallel::defaultNumThreads()
          # n_workers <- min(n_workers, p)
          # RcppParallel::setThreadOptions(numThreads = n_workers)
          result <- scrimpab_rcpp(
            data,
            query,
            window_size,
            as.logical(progress)
          )
          # RcppParallel::setThreadOptions(numThreads = p)
        } else {
          result <- scrimpab_rcpp(
            data,
            query,
            window_size,
            as.logical(progress)
          )
        }
      },
      error = print
    )
    "!DEBUG End AB-Join"
  }
}
