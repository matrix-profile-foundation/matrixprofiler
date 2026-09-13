#' Relative Frequency Contrast Profile
#'
#' Computes an exact, NA-aware Relative Frequency Contrast Profile (RFCP)
#' between a positive series and a negative series. For every valid positive
#' window, RFCP compares ordered non-overlapping neighbor ranks from the
#' positive self-join and the positive-to-negative AB-join.
#'
#' @param positive_data Numeric positive/query time series.
#' @param negative_data Numeric negative/reference time series.
#' @param window_size Integer subsequence length, at least `2` and no larger
#'   than either input series.
#' @param max_freq Integer number of non-overlapping neighbor ranks requested.
#' @param exclusion_radius Integer exclusion radius between accepted neighbor
#'   starts. `NULL` uses `window_size`, the RFCP protocol setting.
#' @param query_begin First positive query start to process, using one-based
#'   inclusive coordinates. Defaults to the first start.
#' @param query_end Last positive query start to process, using one-based
#'   inclusive coordinates. `NULL` processes through the final start.
#' @param return_profiles Whether to return the complete RFMP, RFCP, and index
#'   matrices for the requested query block. This validation mode is limited
#'   to 5,000,000 rank-query cells. Defaults to `FALSE`.
#' @param n_workers Number of native worker threads. Defaults to `1`.
#' @param progress Logical progress flag. Reserved for native progress
#'   reporting; defaults to `TRUE`.
#'
#' @details
#' `NA`, `NaN`, `Inf`, and `-Inf` are barriers. Windows touching a barrier,
#' constant windows, and otherwise non-normalizable windows cannot be queries
#' or neighbors. Returned neighbor indices are one-based global coordinates.
#' Candidates are ordered by increasing distance and then increasing index.
#' After accepting a neighbor at `j`, the closed interval
#' `[j - exclusion_radius, j + exclusion_radius]` is excluded. The same
#' neighborhood around the query is excluded before AA selection.
#'
#' For rank `k`, distances are clipped at `sqrt(2 * window_size)` before
#' calculating `max((AB - AA) / sqrt(2 * window_size), 0)`. Missing ranks
#' contribute zero while the RMSC denominator remains the requested
#' `max_freq`. `effective_max_freq` is diagnostic: it is the largest number of
#' paired finite AA/AB ranks observed for any completed query and never changes
#' the denominator.
#'
#' The default mode materializes only the RMSC values for the requested query
#' block and the rank vectors for the selected Plato. Blocks use global query
#' coordinates and can be concatenated for checkpointed execution. With
#' `return_profiles = TRUE`, RFMP matrices contain unclipped z-normalized
#' Euclidean distances and RFCP contains the corresponding clipped contrasts.
#'
#' @return A list containing `rms_profile`, complete positive and negative
#'   validity masks, the selected Plato index and first AA twin, the RFCP/RFMP
#'   rank vectors at that Plato, query-range and completion metadata, and
#'   optionally the complete rank-by-query matrices.
#'
#' @export
#' @examples
#' positive <- c(sin(seq(0, 4 * pi, length.out = 80)), NA,
#'               sin(seq(0, 4 * pi, length.out = 80)))
#' negative <- cos(seq(0, 10 * pi, length.out = 180))
#' result <- rfcp(positive, negative, 16, max_freq = 3,
#'                n_workers = 2, progress = FALSE)
rfcp <- function(positive_data, negative_data, window_size, max_freq,
                 exclusion_radius = NULL, query_begin = 1L, query_end = NULL,
                 return_profiles = FALSE, n_workers = 1L, progress = TRUE) {
  positive_data <- as.numeric(positive_data)
  negative_data <- as.numeric(negative_data)
  if (!length(positive_data) || !length(negative_data)) {
    stop("`positive_data` and `negative_data` must be non-empty numeric vectors.",
      call. = FALSE
    )
  }
  window_size <- as.integer(checkmate::qassert(window_size, "X1[2,)"))
  max_freq <- as.integer(checkmate::qassert(max_freq, "X1[1,)"))
  if (window_size > min(length(positive_data), length(negative_data))) {
    stop("`window_size` must not exceed either input length.", call. = FALSE)
  }
  if (is.null(exclusion_radius)) {
    exclusion_radius <- window_size
  }
  exclusion_radius <- as.integer(checkmate::qassert(exclusion_radius, "X1[0,)"))
  n_positive_windows <- length(positive_data) - window_size + 1L
  query_begin <- as.integer(checkmate::qassert(query_begin, "X1[1,)"))
  if (is.null(query_end)) {
    query_end <- n_positive_windows
  }
  query_end <- as.integer(checkmate::qassert(query_end, "X1[1,)"))
  if (query_begin > query_end || query_end > n_positive_windows) {
    stop("The inclusive query range must lie within the positive window starts.",
      call. = FALSE
    )
  }
  checkmate::qassert(return_profiles, "B1")
  n_workers <- as.integer(checkmate::qassert(n_workers, "X1[1,)"))
  checkmate::qassert(progress, "B1")

  previous_threads <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = previous_threads), add = TRUE)
  n_workers <- min(n_workers, previous_threads)
  RcppParallel::setThreadOptions(numThreads = n_workers)

  rfcp_na_segmented_native_rcpp_parallel(
    positive_data,
    negative_data,
    window_size,
    max_freq,
    exclusion_radius,
    as.numeric(query_begin - 1L),
    as.numeric(query_end),
    as.logical(return_profiles),
    as.logical(progress)
  )
}
