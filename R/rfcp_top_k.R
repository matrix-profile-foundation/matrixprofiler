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
#' @param adaptive_max_freq Logical. If `TRUE`, dynamically adapts `max_freq`
#'   for subsequent iterations based on the observed support (number of valid
#'   positive self-join neighbors with distance below \eqn{\sqrt{2 \cdot \text{window\_size}}})
#'   of the current Plato. Defaults to `FALSE`.
#' @param min_freq Integer minimum rank frequency floor when `adaptive_max_freq = TRUE`.
#'   Defaults to `2L`.
#' @param n_workers Number of native worker threads for the CPU backend.
#'   Ignored when `backend = "cuda"`. Defaults to `1L`.
#'   For parallel execution across processes, prefer a PSOCK cluster and use
#'   `n_workers = 1` within each worker. See [matrixprofiler] for fork limitations.
#' @param backend Computation backend: `"cpu"` (default) or `"cuda"`.
#' @param progress Logical progress flag for each RFCP iteration. Defaults to `FALSE`.
#' @param query_block_size Optional integer query block size for chunked
#'   execution. If \code{NULL} (default), queries are processed monolithically in
#'   a single call. If specified, queries are processed in contiguous blocks of
#'   size \code{query_block_size}.
#' @param progress_callback Optional callback function \code{function(event)}
#'   invoked periodically with structured progress events (\code{"round_start"},
#'   \code{"block_done"}, \code{"round_done"}).
#' @param checkpoint_callback Optional callback function \code{function(state)}
#'   invoked after each block and Plato with a serializable state list.
#' @param resume_state Optional state list from a previous execution
#'   (\code{checkpoint_callback}) used to resume execution idempotently.
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
#'   \item If \code{adaptive_max_freq = TRUE}, re-evaluates the frequency for
#'     the next iteration based on the count of finite AA neighbors with distance
#'     below \eqn{\sqrt{2 \cdot \text{window\_size}}}, bounded below by \code{min_freq}.
#' }
#'
#' @return An object of class \code{"rfcp_top_k"}, which is a list containing:
#' \itemize{
#'   \item \code{k_requested}: requested number of Platos.
#'   \item \code{k_found}: number of Platos discovered before stopping.
#'   \item \code{window_size}: subsequence length.
#'   \item \code{exclusion_radius}: exclusion radius used.
#'   \item \code{adaptive_max_freq}: whether adaptive frequency was enabled.
#'   \item \code{min_freq}: minimum frequency floor used.
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
                       adaptive_max_freq = FALSE, min_freq = 2L,
                       n_workers = 1L, backend = "cpu", progress = FALSE,
                       query_block_size = NULL,
                       progress_callback = NULL,
                       checkpoint_callback = NULL,
                       resume_state = NULL) {
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
  checkmate::assert_flag(adaptive_max_freq)
  min_freq <- as.integer(checkmate::qassert(min_freq, "X1[1,)"))
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

  if (!is.null(query_block_size)) {
    query_block_size <- as.integer(checkmate::qassert(query_block_size, "X1[1,)"))
  }
  if (!is.null(progress_callback)) {
    checkmate::assert_function(progress_callback)
  }
  if (!is.null(checkpoint_callback)) {
    checkmate::assert_function(checkpoint_callback)
  }
  if (!is.null(resume_state)) {
    checkmate::assert_list(resume_state)
  }

  initial_n_positive <- length(positive_data)
  initial_n_negative <- length(negative_data)
  n_positive_windows <- initial_n_positive - window_size + 1L

  if (!is.null(resume_state)) {
    sig <- resume_state$signature
    if (is.null(sig) || !is.list(sig) ||
      !identical(as.integer(sig$n_positive), initial_n_positive) ||
      !identical(as.integer(sig$n_negative_initial), initial_n_negative) ||
      !identical(as.integer(sig$window_size), window_size) ||
      !identical(as.integer(sig$k_platos), k_platos) ||
      !identical(as.integer(sig$exclusion_radius), exclusion_radius) ||
      !identical(as.logical(sig$adaptive_max_freq), adaptive_max_freq) ||
      !identical(as.integer(sig$min_freq), min_freq) ||
      !identical(as.character(sig$backend), backend)) {
      stop("resume_state is incompatible with current input parameters.", call. = FALSE)
    }

    count_found <- as.integer(resume_state$count_found)
    platos_found <- vector("list", k_platos)
    if (count_found > 0L && !is.null(resume_state$platos_found)) {
      platos_found[seq_len(count_found)] <- resume_state$platos_found[seq_len(count_found)]
    }
    current_negative <- resume_state$current_negative
    excluded_query_mask <- resume_state$excluded_query_mask
    max_freq_vec <- as.integer(resume_state$max_freq_vec)
    k_start <- as.integer(resume_state$k_iter)

    current_round_blocks <- resume_state$current_round_blocks
    if (is.null(current_round_blocks)) {
      current_round_blocks <- list()
    }
    completed_cursor <- as.integer(resume_state$completed_query_cursor)
    if (is.null(completed_cursor) || is.na(completed_cursor) || completed_cursor < 0L) {
      completed_cursor <- 0L
    }
    cursor <- completed_cursor + 1L
  } else {
    count_found <- 0L
    platos_found <- vector("list", k_platos)
    current_negative <- negative_data
    excluded_query_mask <- logical(n_positive_windows)
    k_start <- 1L
    current_round_blocks <- list()
    cursor <- 1L
  }

  make_state <- function(k_iter, cursor_completed, round_blocks) {
    list(
      version = 1L,
      algorithm = "rfcp_top_k",
      signature = list(
        n_positive = initial_n_positive,
        n_negative_initial = initial_n_negative,
        window_size = window_size,
        k_platos = k_platos,
        exclusion_radius = exclusion_radius,
        adaptive_max_freq = adaptive_max_freq,
        min_freq = min_freq,
        backend = backend
      ),
      k_iter = as.integer(k_iter),
      count_found = as.integer(count_found),
      platos_found = platos_found[seq_len(count_found)],
      max_freq_vec = as.integer(max_freq_vec),
      excluded_query_mask = excluded_query_mask,
      current_negative = current_negative,
      current_round_blocks = round_blocks,
      completed_query_cursor = as.integer(cursor_completed)
    )
  }

  block_step <- if (is.null(query_block_size)) n_positive_windows else query_block_size
  n_blocks <- ceiling(n_positive_windows / block_step)

  if (k_start <= k_platos) {
    for (k_iter in seq.int(k_start, k_platos)) {
      if (k_iter > k_start) {
        current_round_blocks <- list()
        cursor <- 1L
      }

      current_freq <- max_freq_vec[[k_iter]]
      t0_round <- proc.time()[["elapsed"]]

      if (!is.null(progress_callback)) {
        progress_callback(list(
          event = "round_start",
          k = k_iter,
          max_freq = current_freq,
          total_queries = n_positive_windows,
          n_blocks = n_blocks
        ))
      }

      while (cursor <= n_positive_windows) {
        begin <- cursor
        end <- min(begin + block_step - 1L, n_positive_windows)

        t0_blk <- proc.time()[["elapsed"]]

        res_block <- rfcp(
          positive_data = positive_data,
          negative_data = current_negative,
          window_size = window_size,
          max_freq = current_freq,
          exclusion_radius = exclusion_radius,
          query_begin = begin,
          query_end = end,
          return_profiles = FALSE,
          n_workers = n_workers,
          progress = if (is.null(query_block_size)) progress else FALSE,
          backend = backend
        )

        elapsed_blk <- proc.time()[["elapsed"]] - t0_blk

        block_entry <- list(
          query_begin = begin,
          query_end = end,
          rms_profile = as.numeric(res_block$rms_profile),
          plato_index = res_block$plato_index,
          plato_rms = res_block$plato_rms,
          rfmpi_aa_at_plato = res_block$rfmpi_aa_at_plato,
          rfmp_aa_at_plato = res_block$rfmp_aa_at_plato,
          rfmpi_ab_at_plato = res_block$rfmpi_ab_at_plato,
          rfmp_ab_at_plato = res_block$rfmp_ab_at_plato,
          rfcp_at_plato = res_block$rfcp_at_plato,
          plato_twin_index = res_block$plato_twin_index,
          effective_max_freq = res_block$effective_max_freq
        )
        current_round_blocks[[length(current_round_blocks) + 1L]] <- block_entry
        cursor <- end + 1L

        if (!is.null(progress_callback)) {
          progress_callback(list(
            event = "block_done",
            k = k_iter,
            block_index = length(current_round_blocks),
            n_blocks = n_blocks,
            query_begin = begin,
            query_end = end,
            elapsed_seconds = elapsed_blk,
            block_plato_rms = res_block$plato_rms
          ))
        }

        if (!is.null(checkpoint_callback)) {
          checkpoint_callback(make_state(
            k_iter = k_iter,
            cursor_completed = end,
            round_blocks = current_round_blocks
          ))
        }
      }

      # Determine best candidate outside the excluded query mask
      full_rms_profile <- unlist(lapply(current_round_blocks, `[[`, "rms_profile"), use.names = FALSE)
      if (is.null(full_rms_profile) || !length(full_rms_profile)) {
        break
      }

      # Mask previously excluded queries
      full_rms_profile[excluded_query_mask] <- NA_real_

      valid_candidates <- which(is.finite(full_rms_profile))
      if (!length(valid_candidates)) {
        break
      }

      best_candidate_idx <- valid_candidates[[which.max(full_rms_profile[valid_candidates])]]
      best_candidate_rms <- full_rms_profile[[best_candidate_idx]]

      if (is.na(best_candidate_rms) || best_candidate_rms < min_rms) {
        break
      }

      count_found <- count_found + 1L

      target_block <- Find(
        function(b) best_candidate_idx >= b$query_begin && best_candidate_idx <= b$query_end,
        current_round_blocks
      )
      if (!is.null(target_block) && identical(best_candidate_idx, as.integer(target_block$plato_index))) {
        aa_indices <- target_block$rfmpi_aa_at_plato
        aa_distances <- target_block$rfmp_aa_at_plato
        ab_distances <- target_block$rfmp_ab_at_plato
        ab_indices <- target_block$rfmpi_ab_at_plato
        rfcp_vec <- target_block$rfcp_at_plato
        plato_twin <- target_block$plato_twin_index
        effective_freq <- target_block$effective_max_freq
      } else {
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

      # Adapt max_freq dynamically for subsequent iterations if enabled
      if (isTRUE(adaptive_max_freq) && k_iter < k_platos) {
        clip_distance <- sqrt(2 * window_size)
        valid_support <- sum(
          aa_distances < clip_distance & is.finite(aa_distances) & aa_distances >= 0,
          na.rm = TRUE
        )
        next_freq <- max(
          min_freq,
          min(current_freq, as.integer(valid_support))
        )
        max_freq_vec[(k_iter + 1L):k_platos] <- next_freq
      }

      elapsed_round <- proc.time()[["elapsed"]] - t0_round

      if (!is.null(progress_callback)) {
        progress_callback(list(
          event = "round_done",
          k = k_iter,
          plato_index = best_candidate_idx,
          plato_rms = best_candidate_rms,
          effective_max_freq = effective_freq,
          next_max_freq = if (k_iter < k_platos) max_freq_vec[[k_iter + 1L]] else current_freq,
          elapsed_seconds = elapsed_round
        ))
      }

      current_round_blocks <- list()
      cursor <- 1L

      if (!is.null(checkpoint_callback)) {
        checkpoint_callback(make_state(
          k_iter = k_iter + 1L,
          cursor_completed = 0L,
          round_blocks = list()
        ))
      }
    }
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
      adaptive_max_freq = adaptive_max_freq,
      min_freq = min_freq,
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
    "Window size: %d | Exclusion radius: %d | Requested: %d | Adaptive: %s\n",
    x$window_size, x$exclusion_radius, x$k_requested,
    ifelse(isTRUE(x$adaptive_max_freq), sprintf("TRUE (min_freq=%d)", x$min_freq), "FALSE")
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
