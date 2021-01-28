#' Fast implementation of Matrix Profile, without FFT
#'
#' Computes the Matrix Profile and Profile Index for Univariate Time Series.
#'
#' @param data Required. Any 1-dimension series of numbers (`matrix`, `vector`, `ts` etc.) (See details).
#' @param window_size Required. An integer defining the rolling window size.
#' @param query Optional. Another 1-dimension series of numbers for an AB-join similarity. Default is `NULL` (See
#'   details).
#' @param exclusion_zone A numeric. Defines the size of the area around the rolling window that will be ignored to avoid
#'   trivial matches. Default is `0.5`, i.e., half of the `window_size`.
#' @param idxs A logical. Specifies if the computation will return the Profile Index or not. Defaults to `TRUE`.
#' @param distance A string. Currently accepts `euclidean` and `pearson`. Defaults to `euclidean`.
#' @param n_workers An integer. The number of threads using for computing. Defaults to `1`.
#' @param progress A logical. If `TRUE` (the default) will show a progress bar. Useful for long computations. (See
#'   details)
#'
#' @details This algorithm was developed apart from the main Matrix Profile branch that relies on Fast Fourier Transform
#'   (FFT) at least in one part of the process. This algorithm doesn't use FFT and is several times faster. It also
#'   relies on Ogita's work to better precision computing mean and standard deviation (part of the process). About
#'   `progress`, it is really recommended to use it as feedback for long computations. It indeed adds some (neglectable)
#'   overhead, but the benefit of knowing that your computer is still computing is much bigger than the seconds you may
#'   lose in the final benchmark. About `n_workers`, for Windows systems, this package uses TBB for multithreading, and
#'   Linux and macOS, use TinyThread++. This may or not raise some issues in the future, so we must be aware of slower
#'   processing due to different mutexes implementations or even unexpected crashes. The Windows version is usually more
#'   reliable. The `data` and `query` parameters will be internally converted to a single vector using `as.numeric()`,
#'   thus, bear in mind that a multidimensional matrix may not work as you expect, but most 1-dimensional data types
#'   will work normally. If `query` is provided, expect the same pre-procesment done for `data`; in addition,
#'   `exclusion_zone` will be ignored and set to `0`. Both `data` and `query` doesn't need to have the same size and
#'   they can be interchanged if both are provided. The difference will be in the returning object. AB-Join returns the
#'   Matrix Profile 'A' and 'B' i.e., the distance between a rolling window from query to data and from data to query.
#'
#' @return Returns a list with the Matrix Profile, Profile Index (if `idxs` is `TRUE`), and some information about the
#'   settings used to build it.
#' @export
#'
#' @family matrix profile computations
#'
#' @examples
#' \donttest{
#' mp <- mpx(runif(200), window_size = 30)
#' }
#'
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
