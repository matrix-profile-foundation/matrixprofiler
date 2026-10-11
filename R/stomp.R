#' Matrix Profile Computation
#'
#' STOMP is a faster implementation with the caveat that is not anytime as STAMP or SCRIMP.
#'
#' @param left_right_profile (`stomp()` only) A boolean. If `TRUE`, the function will return the left and right profiles.
#'
#' @details ## stomp
#'  The STOMP uses a faster implementation to compute the Matrix Profile and Profile Index. It can be stopped earlier by
#'  the user, but the result is not considered anytime, just incomplete. For a anytime algorithm, use `stamp()` or
#'  `scrimp()`.
#'
#' @family matrix profile computations

#' @export
#' @rdname mp_algos
#' @order 2
#' @examples
#' mp <- stomp(motifs_discords_small, 50)
stomp <- function(data, window_size, query = NULL, exclusion_zone = 0.5, n_workers = 1, progress = TRUE, left_right_profile = FALSE) {
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
  n_workers <- as.integer(checkmate::qassert(n_workers, "X+"))
  if (n_workers > 1L && .Platform$OS.type != "windows") {
    is_child <- tryCatch(
      get("isChild", envir = asNamespace("parallel"))(),
      error = function(e) FALSE
    )
    if (isTRUE(is_child)) {
      warning("Multithreaded STOMP execution inside a forked process is unsafe and will use 1 worker to prevent crashes.", call. = FALSE)
      n_workers <- 1L
    }
  }
  checkmate::qassert(progress, "B+")
  checkmate::qassert(left_right_profile, "B1")

  ez <- exclusion_zone
  result <- NULL

  min_size <- if (is.null(query)) length(data) else min(length(data), length(query))
  if (window_size >= min_size) {
    stop("Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (is.null(query) && exclusion_zone > 0) {
    ez_points <- round(window_size * exclusion_zone)
    if ((length(data) - window_size) <= ez_points) {
      stop("Time series is too short relative to desired window size and exclusion zone.", call. = FALSE)
    }
  }

  # Register anytime exit point ----------------------
  "!DEBUG Register anytime exit point"
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      if (isTRUE(result$partial)) {
        warning(
          "STOMP computation was interrupted by the user: returning an incomplete matrix profile covering only the initial windows processed.",
          call. = FALSE
        )
      }
      result$ez <- ez
      return(result)
    },
    TRUE
  )

  # Computation ------------------------------------
  "!DEBUG Computation"
  if (is.null(query)) {
    ## Self-Join ====================================
    # Can return LMP or RMP
    "!DEBUG Self-Join"
    if (n_workers > 1) {
      result <- stomp_rcpp_parallel(
        data,
        data,
        window_size,
        ez,
        n_workers,
        as.logical(progress),
        as.logical(left_right_profile)
      )
    } else {
      result <- stomp_rcpp(
        data,
        data,
        window_size,
        ez,
        as.logical(progress),
        as.logical(left_right_profile)
      )
    }
    "!DEBUG End Self-Join"
  } else {
    ## AB-Join ====================================
    # Can't return LMP or RMP

    if (isTRUE(left_right_profile)) {
      warning("AB-Join can't return Left of Right Profiles, and this argument will be ignored.")
    }

    "!DEBUG AB-Join"
    ez <- 0

    if (n_workers > 1) {
      result <- stomp_rcpp_parallel(
        data,
        query,
        window_size,
        ez,
        n_workers,
        as.logical(progress),
        as.logical(FALSE)
      )
    } else {
      result <- stomp_rcpp(
        data,
        query,
        window_size,
        ez,
        as.logical(progress),
        as.logical(FALSE)
      )
    }
    "!DEBUG End AB-Join"
  }
}
