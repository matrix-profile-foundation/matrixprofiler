#' MPdist - Distance between Time Series using Matrix Profile
#'
#' MPdist is a distance measure which considers two time series to be similar
#' if they share many similar subsequences, regardless of the order of matching subsequences.
#'
#' @param ref_data a `matrix` or a `vector`. The reference data.
#' @param query_data a `matrix` or a `vector`. The query data.
#' @param window_size an `int`. Size of the sliding window. Must be at least 4.
#' @param type the type of result: `"simple"` (default) returns a scalar distance, `"vector"` returns distance between all sliding windows.
#' @param thr threshold for MPdist. (Default is `0.05`).
#' @param n_workers an `int`. Number of threads for parallel computation. (Default is `1`).
#'
#' @return Returns the distance scalar (`"simple"`) or a numeric vector (`"vector"`).
#' @export
#' @family distance measure
#' @references * Gharghabi S, Imani S, Bagnall A, Darvishzadeh A, Keogh E. Matrix Profile XII:
#' MPdist: A Novel Time Series Distance Measure to Allow Data Mining in More Challenging Scenarios.
#' In: 2018 IEEE International Conference on Data Mining (ICDM). 2018.
#' @references Website: <https://sites.google.com/site/mpdistinfo/>
#' @examples
#' ref_data <- motifs_discords_small[1:200]
#' qe_data <- motifs_discords_small[201:400]
#' w <- 30
#' deq <- mpdist(ref_data, qe_data, w)
mpdist <- function(ref_data, query_data, window_size, type = c("simple", "vector"), thr = 0.05, n_workers = 1L) {
  type <- match.arg(type)

  ref_data <- as.numeric(ref_data)
  query_data <- as.numeric(query_data)
  checkmate::qassert(ref_data, "N+")
  checkmate::qassert(query_data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (window_size < 4L) {
    stop("`window_size` must be at least 4.", call. = FALSE)
  }
  checkmate::qassert(thr, "N1(0,1)")
  n_workers <- as.integer(checkmate::qassert(n_workers, "X+"))

  if (length(ref_data) < length(query_data)) {
    temp <- ref_data
    ref_data <- query_data
    query_data <- temp
  }

  if (window_size > length(query_data)) {
    stop("`window_size` cannot be larger than query length.", call. = FALSE)
  }

  if (type == "simple" || length(ref_data) == length(query_data)) {
    dist <- mpdist_simple_rcpp(ref_data, query_data, window_size, thr, n_workers)
    return(dist)
  } else {
    dist_vect <- mpdist_vect_rcpp(ref_data, query_data, window_size, thr, n_workers)
    return(dist_vect)
  }
}
