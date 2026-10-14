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
#' @param n_workers Number of native worker threads for the CPU backend.
#'   Ignored when `backend = "cuda"`. Defaults to `1`.
#'   For parallel execution across processes, prefer a PSOCK cluster and use
#'   `n_workers = 1` within each worker. See [matrixprofiler] for fork limitations.
#' @param backend Computation backend: `"cpu"` (default) or `"cuda"`. The CUDA
#'   backend is optional and requires a local build with CUDA enabled (see
#'   `inst/cuda/build_cuda_rcpp.sh`).
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
#' Native query tiles use
#' `min(32768, max(1024, 16 * window_size))`, further capped at the ceiling of
#' the requested query count divided by the effective worker count. The first
#' tile starts at `query_begin`, so checkpoint blocks do not perform discarded
#' warm-up rows.
#' The CUDA backend executes the full join query-by-query with optional profile
#' matrices returned identically to the CPU backend.
#'
#' @return A list with RFCP outputs:
#' \itemize{
#'   \item `rms_profile`: numeric vector of root-mean-square contrasts.
#'   \item `valid_window_positive`: logical mask of valid query windows.
#'   \item `valid_window_negative`: logical mask of valid reference windows.
#'   \item `plato_index`: one-based index of the highest-contrast query.
#'   \item `plato_rms`: contrast score at `plato_index`.
#'   \item `plato_twin_index`: one-based index of the nearest self-join neighbor
#'     at `plato_index`.
#'   \item `rfcp_at_plato`: clipped contrasts across accepted neighbor ranks.
#'   \item `rfmp_aa_at_plato`: positive self-join distances at plato.
#'   \item `rfmp_ab_at_plato`: negative AB distances at plato.
#'   \item `rfmpi_aa_at_plato`: positive neighbor indices at plato.
#'   \item `rfmpi_ab_at_plato`: negative neighbor indices at plato.
#'   \item `requested_max_freq`: requested rank count.
#'   \item `effective_max_freq`: highest paired rank observed.
#'   \item `query_begin`, `query_end`: one-based query range.
#'   \item `completed_query_end`: end of completed queries.
#'   \item `partial`: whether computation stopped early.
#'   \item `rfmp_aa`, `rfmp_ab`, `rfcp`, `rfmpi_aa`, `rfmpi_ab`: optional
#'     full profile matrices when `return_profiles = TRUE`.
#' }
#'
#' @export
#' @examples
#' pos <- c(sin(seq(0, 4 * pi, length.out = 100)))
#' neg <- c(cos(seq(0, 10 * pi, length.out = 150)))
#' res <- rfcp(pos, neg, window_size = 10, max_freq = 2)
#' res$plato_index
rfcp <- function(positive_data, negative_data, window_size, max_freq,
                 exclusion_radius = NULL, query_begin = 1L, query_end = NULL,
                 return_profiles = FALSE, n_workers = 1L, progress = TRUE,
                 backend = "cpu") {
  positive_data <- as.numeric(positive_data)
  negative_data <- as.numeric(negative_data)
  if (!length(positive_data) || !length(negative_data)) {
    stop("Both positive_data and negative_data must be non-empty numeric vectors.",
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
  backend <- match.arg(backend, c("cpu", "cuda"))

  if (identical(backend, "cuda")) {
    return(rfcp_gpu_rcpp(
      positive_data,
      negative_data,
      window_size,
      max_freq,
      exclusion_radius,
      as.numeric(query_begin - 1L),
      as.numeric(query_end),
      as.logical(return_profiles)
    ))
  }

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

  rfcp_na_segmented_native_rcpp_parallel(
    positive_data,
    negative_data,
    window_size,
    max_freq,
    exclusion_radius,
    as.numeric(query_begin - 1L),
    as.numeric(query_end),
    n_workers,
    as.logical(return_profiles),
    as.logical(progress)
  )
}
