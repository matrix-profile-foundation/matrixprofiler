#' Univariate STOMP algorithm
#'
#' Computes the Matrix Profile and Profile Index for Univariate Time Series.
#'
#' @details
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its
#' generality, versatility, simplicity and scalability. In particular it has implications for time
#' series motif discovery, time series joins, shapelet discovery (classification), density
#' estimation, semantic segmentation, visualization, rule discovery, clustering etc. `verbose`
#' changes how much information is printed by this function; `0` means nothing, `1` means text, `2`
#' adds the progress bar, `3` adds the finish sound. `exclusion_zone` is used to avoid  trivial
#' matches; if a query data is provided (join similarity), this parameter is ignored.
#'
#' @param data Required. Any 1-dimension series of numbers (`matrix`, `vector`, `ts` etc.) (See details).
#' @param window_size Required. An integer defining the rolling window size.
#' @param query Optional. Another 1-dimension series of numbers for an AB-join similarity. Default is `NULL` (See
#'   details).
#' @param exclusion_zone A numeric. Defines the size of the area around the rolling window that will be ignored to avoid
#'   trivial matches. Default is `0.5`, i.e., half of the `window_size`.
#' @param n_workers An integer. The number of threads using for computing. Defaults to `1`.
#' @param progress A logical. If `TRUE` (the default) will show a progress bar. Useful for long computations. (See
#'   details)
#'
#' @export
#'
#' @family matrix profile computations
#'
#' @describeIn stomp Single thread version.
#'
#' @references * Zhu Y, Zimmerman Z, Senobari NS, Yeh CM, Funning G. Matrix Profile II : Exploiting
#'   a Novel Algorithm and GPUs to Break the One Hundred Million Barrier for Time Series Motifs and
#'   Joins. Icdm. 2016 Jan 22;54(1):739-48.
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' \donttest{
#' mp <- stomp(runif(200), window_size = 30)
#' }
#'
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

  # Register anytime exit point
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
