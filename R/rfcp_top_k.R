#' Top-K Relative Frequency Contrast Profile
#'
#' Discovers the top-K contrasting patterns (Platos) between a positive time
#' series and a negative time series using an iterative Relative Frequency
#' Contrast Profile (RFCP) procedure.
#'
#' @param positive_data Numeric positive/query time series.
#' @param negative_data Numeric negative/reference time series.
#' @param window_size Integer subsequence length, at least `2` and no larger
#'   than either input series.
#' @param k_platos Integer number of distinct Platos to discover. Defaults to `2L`.
#' @param max_freq Integer scalar or integer vector specifying the number of
#'   non-overlapping neighbor ranks requested for each Plato iteration. If a
#'   scalar is provided, it is used for all iterations. If a vector is provided,
#'   its length must equal `k_platos`, allowing dynamic re-evaluation of
#'   frequency for rarer morphologies. Defaults to `5L`.
#' @param exclusion_radius Integer exclusion radius between accepted neighbor
#'   starts and for padding concatenated segments. `NULL` uses `window_size`.
#' @param min_rms Numeric threshold for root-mean-square contrast. Iteration
#'   stops early if the detected Plato has an RMS contrast below this value.
#'   Defaults to `1e-4`.
#' @param n_workers Number of native worker threads for the CPU backend.
#'   Ignored when `backend = "cuda"`. Defaults to `1L`.
#' @param backend Computation backend: `"cpu"` (default) or `"cuda"`.
#' @param progress Logical progress flag for each RFCP iteration. Defaults to `FALSE`.
#'
#' @details
#' In each iteration \eqn{i \in \{1, \dots, K\}}, the algorithm:
#' \enumerate{
#'   \item Executes \code{\link{rfcp}} using the current augmented negative series
#'     and the requested frequency \code{max_freq[i]}.
#'   \item Identifies the highest-contrast query window (\eqn{\text{Plato}_i})
#'     outside the exclusion zone of previously discovered Platos and their
#'     associated neighbors.
#'   \item If the contrast score is below \code{min_rms} or no valid candidate
#'     remains, iteration terminates early.
#'   \item Extracts the subsequences of \eqn{\text{Plato}_i} and all its accepted
#'     positive self-join neighbors (\code{rfmpi_aa_at_plato}), expanded by the
#'     exclusion zone:
#'     \deqn{[\max(1, j - \text{exclusion\_radius}), \min(N^{(+)}, j + \text{window\_size} - 1 + \text{exclusion\_radius})]}
#'   \item Merges overlapping intervals to maintain signal continuity, and appends
#'     the resulting segments to the negative series separated by \code{NA}
#'     barriers. Because \code{rfcp} is NA-aware, this suppresses the contrast of
#'     that morphology and its trivial overlaps in subsequent iterations without
#'     introducing spurious cross-boundary windows.
#' }
#'
#' @return An object of class \code{"rfcp_top_k"}, which is a list containing:
#' \itemize{
#'   \item \code{k_requested}: requested number of Platos.
#'   \item \code{k_found}: number of Platos discovered before stopping.
#'   \item \code{window_size}: subsequence length.
#'   \item \code{exclusion_radius}: exclusion radius used.
#'   \item \code{max_freq}: vector of frequencies used for each Plato.
#'   \item \code{platos}: list of discovered Plato objects, each containing:
#'     \itemize{
#'       \item \code{k}: Plato iteration rank (1, 2, ...).
#'       \item \code{plato_index}: one-based start index in \code{positive_data}.
#'       \item \code{plato_rms}: contrast score at \code{plato_index}.
#'       \item \code{plato_subsequence}: numeric vector containing the pattern.
#'       \item \code{plato_twin_index}: one-based start index of nearest self-join neighbor.
#'       \item \code{rfcp_at_plato}: clipped contrasts across accepted ranks.
#'       \item \code{rfmp_aa_at_plato}: positive self-join distances.
#'       \item \code{rfmp_ab_at_plato}: negative AB distances.
#'       \item \code{rfmpi_aa_at_plato}: positive neighbor indices.
#'       \item \code{rfmpi_ab_at_plato}: negative neighbor indices.
#'       \item \code{effective_max_freq}: highest paired rank observed.
#'     }
#'   \item \code{negative_data_final}: the negative series augmented with the
#'     suppressed Plato occurrences and exclusion zones.
#' }
#'
#' @export
#' @examples
#' set.seed(42)
#' # Create positive series with 5 sine pulses and 2 cosine pulses
#' t <- seq(0, 4 * pi, length.out = 20)
#' sine_pulse <- sin(t)
#' cos_pulse <- cos(t)
#' pos <- rnorm(300, sd = 0.05)
#' for (p in c(20, 60, 100, 140, 180)) pos[p:(p + 19)] <- pos[p:(p + 19)] + sine_pulse
#' for (p in c(220, 260)) pos[p:(p + 19)] <- pos[p:(p + 19)] + cos_pulse * 1.5
#' neg <- rnorm(300, sd = 0.05)
#'
#' top_platos <- rfcp_top_k(pos, neg, window_size = 20, k_platos = 2, max_freq = c(4, 2))
#' print(top_platos)
rfcp_top_k <- function(positive_data, negative_data, window_size,
                       k_platos = 2L, max_freq = 5L,
                       exclusion_radius = NULL, min_rms = 1e-4,
                       n_workers = 1L, backend = "cpu", progress = FALSE) {
  positive_data <- as.numeric(positive_data)
  negative_data <- as.numeric(negative_data)
  if (!length(positive_data) || !length(negative_data)) {
    stop("Both positive_data and negative_data must be non-empty numeric vectors.",
      call. = FALSE
    )
  }

  window_size <- as.integer(checkmate::qassert(window_size, "X1[2,)"))
  if (window_size > min(length(positive_data), length(negative_data))) {
    stop("`window_size` must not exceed either input length.", call. = FALSE)
  }

  k_platos <- as.integer(checkmate::qassert(k_platos, "X1[1,)"))

  if (is.null(exclusion_radius)) {
    exclusion_radius <- window_size
  }
  exclusion_radius <- as.integer(checkmate::qassert(exclusion_radius, "X1[0,)"))

  checkmate::qassert(min_rms, "N1[0,)")
  checkmate::qassert(max_freq, "X+[1,)")

  if (length(max_freq) == 1L) {
    max_freq_vec <- rep(as.integer(max_freq), k_platos)
  } else if (length(max_freq) == k_platos) {
    max_freq_vec <- as.integer(max_freq)
  } else {
    stop("`max_freq` must be either a scalar or a vector of length `k_platos`.",
      call. = FALSE
    )
  }

  backend <- match.arg(backend, c("cpu", "cuda"))
  n_workers <- as.integer(checkmate::qassert(n_workers, "X1[1,)"))
  checkmate::qassert(progress, "B1")

  n_positive_windows <- length(positive_data) - window_size + 1L
  excluded_query_mask <- logical(n_positive_windows)

  current_negative <- negative_data
  platos_found <- vector("list", k_platos)
  count_found <- 0L

  for (k_iter in seq_len(k_platos)) {
    current_freq <- max_freq_vec[[k_iter]]

    res <- rfcp(
      positive_data = positive_data,
      negative_data = current_negative,
      window_size = window_size,
      max_freq = current_freq,
      exclusion_radius = exclusion_radius,
      return_profiles = FALSE,
      n_workers = n_workers,
      progress = progress,
      backend = backend
    )

    # Determine best candidate outside the excluded query mask
    rms_prof <- res$rms_profile
    if (is.null(rms_prof) || !length(rms_prof)) {
      break
    }

    # Mask previously excluded queries
    rms_prof[excluded_query_mask] <- NA_real_

    valid_candidates <- which(is.finite(rms_prof))
    if (!length(valid_candidates)) {
      break
    }

    best_candidate_idx <- valid_candidates[[which.max(rms_prof[valid_candidates])]]
    best_candidate_rms <- rms_prof[[best_candidate_idx]]

    if (is.na(best_candidate_rms) || best_candidate_rms < min_rms) {
      break
    }

    count_found <- count_found + 1L

    # If the unexcluded best candidate matches the native Plato, we can reuse
    # res$rfmpi_aa_at_plato directly; otherwise, run a single-query lookup
    if (identical(best_candidate_idx, as.integer(res$plato_index))) {
      aa_indices <- res$rfmpi_aa_at_plato
      aa_distances <- res$rfmp_aa_at_plato
      ab_distances <- res$rfmp_ab_at_plato
      ab_indices <- res$rfmpi_ab_at_plato
      rfcp_vec <- res$rfcp_at_plato
      plato_twin <- res$plato_twin_index
      effective_freq <- res$effective_max_freq
    } else {
      # Query single window to retrieve exact rank neighbors
      single_res <- rfcp(
        positive_data = positive_data,
        negative_data = current_negative,
        window_size = window_size,
        max_freq = current_freq,
        exclusion_radius = exclusion_radius,
        query_begin = best_candidate_idx,
        query_end = best_candidate_idx,
        return_profiles = FALSE,
        n_workers = 1L,
        progress = FALSE,
        backend = backend
      )
      aa_indices <- single_res$rfmpi_aa_at_plato
      aa_distances <- single_res$rfmp_aa_at_plato
      ab_distances <- single_res$rfmp_ab_at_plato
      ab_indices <- single_res$rfmpi_ab_at_plato
      rfcp_vec <- single_res$rfcp_at_plato
      plato_twin <- single_res$plato_twin_index
      effective_freq <- single_res$effective_max_freq
    }

    subseq <- positive_data[best_candidate_idx:(best_candidate_idx + window_size - 1L)]

    platos_found[[count_found]] <- list(
      k = count_found,
      plato_index = best_candidate_idx,
      plato_rms = best_candidate_rms,
      plato_subsequence = subseq,
      plato_twin_index = plato_twin,
      rfcp_at_plato = rfcp_vec,
      rfmp_aa_at_plato = aa_distances,
      rfmp_ab_at_plato = ab_distances,
      rfmpi_aa_at_plato = aa_indices,
      rfmpi_ab_at_plato = ab_indices,
      effective_max_freq = effective_freq
    )

    # Collect all occurrences belonging to this Plato's morphology
    valid_aa_neighbors <- aa_indices[is.finite(aa_indices) & aa_indices > 0L]
    all_occurrences <- unique(c(best_candidate_idx, as.integer(valid_aa_neighbors)))

    # Mark the exclusion zone of all occurrences in the query mask
    for (occ in all_occurrences) {
      left_mask <- max(1L, occ - exclusion_radius)
      right_mask <- min(n_positive_windows, occ + exclusion_radius)
      excluded_query_mask[left_mask:right_mask] <- TRUE
    }

    # Extract concatenated segments with exclusion zone padding
    intervals <- matrix(0L, nrow = length(all_occurrences), ncol = 2L)
    for (idx in seq_along(all_occurrences)) {
      occ <- all_occurrences[[idx]]
      seg_start <- max(1L, occ - exclusion_radius)
      seg_end <- min(length(positive_data), occ + window_size - 1L + exclusion_radius)
      intervals[idx, ] <- c(seg_start, seg_end)
    }

    # Merge overlapping intervals to preserve continuous segments
    intervals <- intervals[order(intervals[, 1L]), , drop = FALSE]
    merged_intervals <- list()
    cur_start <- intervals[1L, 1L]
    cur_end <- intervals[1L, 2L]

    if (nrow(intervals) > 1L) {
      for (r in 2L:nrow(intervals)) {
        next_start <- intervals[r, 1L]
        next_end <- intervals[r, 2L]
        if (next_start <= cur_end + 1L) {
          cur_end <- max(cur_end, next_end)
        } else {
          merged_intervals[[length(merged_intervals) + 1L]] <- c(cur_start, cur_end)
          cur_start <- next_start
          cur_end <- next_end
        }
      }
    }
    merged_intervals[[length(merged_intervals) + 1L]] <- c(cur_start, cur_end)

    # Append merged segments separated by NA barriers to current_negative
    segments_to_add <- vector("list", length(merged_intervals))
    for (m_idx in seq_along(merged_intervals)) {
      s_bounds <- merged_intervals[[m_idx]]
      segments_to_add[[m_idx]] <- c(
        NA_real_,
        positive_data[s_bounds[[1L]]:s_bounds[[2L]]],
        NA_real_
      )
    }

    current_negative <- c(current_negative, unlist(segments_to_add, use.names = FALSE))
  }

  if (count_found < k_platos) {
    platos_found <- platos_found[seq_len(count_found)]
  }

  structure(
    list(
      k_requested = k_platos,
      k_found = count_found,
      window_size = window_size,
      exclusion_radius = exclusion_radius,
      max_freq = max_freq_vec[seq_len(count_found)],
      platos = platos_found,
      negative_data_final = current_negative
    ),
    class = "rfcp_top_k"
  )
}

#' @export
print.rfcp_top_k <- function(x, ...) {
  cat(sprintf("Top-%d Relative Frequency Contrast Profile (RFCP)\n", x$k_found))
  cat(sprintf(
    "Window size: %d | Exclusion radius: %d | Requested: %d\n",
    x$window_size, x$exclusion_radius, x$k_requested
  ))
  if (x$k_found == 0L) {
    cat("  No Platos found meeting the contrast threshold.\n")
    return(invisible(x))
  }
  for (i in seq_len(x$k_found)) {
    p <- x$platos[[i]]
    twin_str <- if (is.na(p$plato_twin_index) || p$plato_twin_index <= 0L) {
      "none"
    } else {
      as.character(p$plato_twin_index)
    }
    cat(sprintf(
      "  Plato %d: index = %d, RMS = %.6f, max_freq = %d (twin = %s)\n",
      p$k, p$plato_index, p$plato_rms, x$max_freq[[i]], twin_str
    ))
  }
  invisible(x)
}
