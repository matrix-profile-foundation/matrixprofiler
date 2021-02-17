#' Matrix Profile Computation
#'
#' STAMP Computes the best so far Matrix Profile and Profile Index for Univariate Time Series.
#'
#' @details
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its generality,
#' versatility, simplicity and scalability. In particular it has implications for time series motif discovery, time
#' series joins, shapelet discovery (classification), density estimation, semantic segmentation, visualization, rule
#' discovery, clustering etc.
#'
#' @details
#' `progress`, it is really recommended to use it as feedback for long computations. It indeed adds some
#' (neglectable) overhead, but the benefit of knowing that your computer is still computing is much bigger than the
#' seconds you may lose in the final benchmark. About `n_workers`, for Windows systems, this package uses TBB for
#' multithreading, and Linux and macOS, use TinyThread++. This may or not raise some issues in the future, so we must be
#' aware of slower processing due to different mutexes implementations or even unexpected crashes. The Windows version
#' is usually more reliable. The `data` and `query` parameters will be internally converted to a single vector using
#' `as.numeric()`, thus, bear in mind that a multidimensional matrix may not work as you expect, but most 1-dimensional
#' data types will work normally. If `query` is provided, expect the same pre-procesment done for `data`; in addition,
#' `exclusion_zone` will be ignored and set to `0`. Both `data` and `query` doesn't need to have the same size and they
#' can be interchanged if both are provided. The difference will be in the returning object. AB-Join returns the Matrix
#' Profile 'A' and 'B' i.e., the distance between a rolling window from query to data and from data to query.
#'
#' @details ## stamp
#'  The anytime STAMP computes the Matrix Profile and Profile Index in such manner that it can be stopped before its
#'  complete calculation and return the best so far results allowing ultra-fast approximate solutions.
#'
#' @param data Required. Any 1-dimension series of numbers (`matrix`, `vector`, `ts` etc.) (See details).
#' @param window_size Required. An integer defining the rolling window size.
#' @param query (not yet on `scrimp()`) Optional. Another 1-dimension series of numbers for an AB-join similarity.
#'   Default is `NULL` (See details).
#' @param exclusion_zone A numeric. Defines the size of the area around the rolling window that will be ignored to avoid
#'   trivial matches. Default is `0.5`, i.e., half of the `window_size`.
#' @param s_size A numeric. Used on anytime algorithms (stamp, scrimp, mpx) if only part of the computation is needed.
#' Default is `1.0` (means 100%).
#' @param pre_scrimp A numeric. If not zero, pre_scrimp is computed, using a fraction of the data. Default is `0.25`.
#' This parameter is ignored when using multithread or AB-join.
#' @param n_workers An integer. The number of threads using for computing. Defaults to `1`.
#' @param progress A logical. If `TRUE` (the default) will show a progress bar. Useful for long computations. (See
#'   details)
#'
#' @return Returns a `list` with the `matrix_profile`, `profile_index` (if `idxs` is `TRUE` in `mpx()`), and some
#'  information about the settings used to build it, like `ez` and `partial` when the algorithm is finished early.
#'
#' @family matrix profile computations
#'
#' @seealso `mass()` for the underlying algorithm that finds best match of a query.
#'
#' @references * Yeh CCM, Zhu Y, Ulanova L, Begum N, Ding Y, Dau HA, et al. Matrix profile I: All pairs similarity joins
#'   for time series: A unifying view that includes motifs, discords and shapelets. Proc - IEEE Int Conf Data Mining,
#'   ICDM. 2017;1317-22.
#' @references * Zhu Y, Imamura m, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A New Primitive for Time
#'   Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1-27.
#' @references * Zhu Y, Zimmerman Z, Senobari NS, Yeh CM, Funning G. Matrix Profile II : Exploiting a Novel Algorithm
#'   and GPUs to Break the One Hundred Million Barrier for Time Series Motifs and Joins. Icdm. 2016 Jan 22;54(1):739-48.
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @export
#' @rdname mp_algos
#' @order 1
#'
#' @examples
#' mp <- stamp(motifs_discords_small, 50)
stamp <- function(data, window_size, query = NULL, exclusion_zone = 0.5, s_size = 1.0, n_workers = 1, progress = TRUE) {

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
          result <- stamp_rcpp_parallel(
            data,
            data,
            window_size,
            ez,
            s_size,
            as.logical(progress)
          )
          RcppParallel::setThreadOptions(numThreads = p)
        } else {
          result <- stamp_rcpp(
            data,
            data,
            window_size,
            ez,
            s_size,
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
          result <- stamp_rcpp_parallel(
            data,
            query,
            window_size,
            ez,
            s_size,
            as.logical(progress)
          )
          RcppParallel::setThreadOptions(numThreads = p)
        } else {
          result <- stamp_rcpp(
            data,
            query,
            window_size,
            ez,
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
