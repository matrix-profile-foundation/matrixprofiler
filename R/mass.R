#' Precomputes several values used on MASS
#'
#' @param window_size Required. An integer defining the rolling window size.
#' @param type This changes how the MASS algorithm will compare the rolling window and the data. (See details).
#' @param weights Optional. It is used when the `type` is `weighted`, and has to be the same size as the `window_size`.
#'
#' @details
#' There are currently four ways to compare the window with the data:
#'
#' 1. **normalized**: this normalizes the data and the query window. This is the most frequently used.
#' 2. **non_normalized**: this won't normalize the query window. The data still being normalized.
#' 3. **absolute**: this won't normalize both the data and the query window.
#' 4. **weighted**: this normalizes the data and query window, and also apply a weight vector on the query.
#'
#' @return
#' `mass_pre()` returns a `list` with several precomputations to be used on MASS later. **Attention** use this before
#' `mass()`.
#' @export
#' @rdname mass
#' @order 2
#'
#' @examples
#' pre <- mass_pre(motifs_discords_small, 50)
#' dist_profile <- mass(pre, motifs_discords_small)
mass_pre <- function(data, window_size, query = NULL, type = c("normalized", "non_normalized", "absolute", "weighted"),
                     weights = NULL) {
  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing Arguments"
  if (!is.null(query) && length(window_size) > 1 && length(query) == 1) {
    tmp <- window_size
    window_size <- query
    query <- tmp
  }
  data <- as.numeric(data)
  checkmate::qassert(data, "n+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (!is.null(query)) {
    query <- as.numeric(query)
    checkmate::qassert(query, c("0", "n>=4"))
  }
  type <- match.arg(type)

  if (type == "weighted") {
    if (is.null(weights)) {
      stop("The `weights` argument must be provided.", call. = FALSE)
    }
    if (length(weights) != window_size) {
      stop("The `weights` must be the same size as the `window_size`.", call. = FALSE)
    }
    checkmate::qassert(weights, "N+")
  }

  # Register anytime exit point ----------------------
  result <- NULL
  "!DEBUG Register anytime exit point"
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      attr(result, "type") <- type
      return(result)
    },
    TRUE
  )

  if (is.null(query)) {
    query <- data
  }

  # Computation ------------------------------------
  "!DEBUG Computation"
  result <- switch(type,
    normalized = mass_pre_rcpp(data, query, window_size),
    non_normalized = mass_pre_weighted_rcpp(data, query, window_size, rep(1, window_size)),
    absolute = mass_pre_abs_rcpp(data, query, window_size),
    weighted = mass_pre_weighted_rcpp(data, query, window_size, weights)
  )
}

#' Computes the Distance between the 'data' and the 'query'.
#'
#' This algorithm will use a rolling window, to computes the distance thorough the whole data. This means that the
#' minimum distance found is the *motif* and the maximum distance is the *discord* on that time series. **Attention**
#' you need first to create an object using `mass_pre()`. Read below.
#'
#' @param pre_obj Required. This is the object resulting from `mass_pre()`. The is no *MASS* without a *pre*.
#' @param data Required. Any 1-dimension series of numbers (`matrix`, `vector`, `ts` etc.)
#' @param query Optional. Accepts the same types as `data` and is used for join-similarity. Defaults to `data` for
#'   self-similarity. **IMPORTANT** Use the same data used on `mass_pre()`.
#' @param index An `integer`. This is the index of the rolling window that will be used. Must be between `1` and
#'   `length(data) - window_size + 1`.
#' @param version A `string`. Chooses the version of MASS what will be used. Ignored if `mass_pre()` is not the
#'   "normalized" type.
#' @param n_workers An `integer`. The number of threads used for computing. Defaults to `1`.
#'   For parallel execution across processes, prefer a PSOCK cluster and use
#'   `n_workers = 1` within each worker. See [matrixprofiler] for fork limitations.
#'
#' @return
#' `mass()` returns a `list` with the `distance_profile` and the `last_product` that is only useful for computing the
#' Matrix Profile.
#' @export
#' @rdname mass
#' @order 1
#' @examples
#' pre <- mass_pre(motifs_discords_small, 50)
#' dist_profile <- mass(pre, motifs_discords_small)
mass <- function(pre_obj, data, query = data, index = 1, version = c("v3", "v2"), n_workers = 1) {
  checkmate::qassert(pre_obj, "L+")
  pre_type <- attr(pre_obj, "type")
  if (is.null(pre_type)) {
    pre_type <- pre_obj$type
  }
  if (is.null(pre_type)) {
    pre_type <- "normalized"
  }
  type <- match.arg(pre_type, c("normalized", "non_normalized", "absolute", "weighted"))
  version <- match.arg(version)
  data <- as.numeric(data)
  checkmate::qassert(data, "n+")
  query <- as.numeric(query)
  checkmate::qassert(query, "n+")
  index <- as.integer(checkmate::qassert(index, "X+(0,)"))
  n_workers <- as.integer(checkmate::qassert(n_workers, "X+"))
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
  if (length(data) != pre_obj$data_size) {
    stop("Argument `data` is not the same as computed in `pre_obj`.", call. = FALSE)
  }

  if (length(query) != (
    ifelse(type == "absolute", length(pre_obj$sumy2), length(pre_obj$query_mean)) +
      pre_obj$window_size - 1)) {
    stop("Argument `query` is not the same as computed in `pre_obj`.", call. = FALSE)
  }

  # Register anytime exit point ----------------------
  result <- NULL
  "!DEBUG Register anytime exit point"
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      return(result)
    },
    TRUE
  )

  # Computation ------------------------------------
  "!DEBUG Computation"
  query_window <- query[index:(index + pre_obj$window_size - 1)]

  result <- switch(type,
    normalized = if (version == "v3") {
      if (n_workers > 1) {
        mass3_rcpp_parallel(
          query_window, data, as.integer(pre_obj$data_size),
          as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd, pre_obj$query_mean[index],
          pre_obj$query_sd[index],
          grain = 4096,
          n_workers = n_workers
        )
      } else {
        mass3_rcpp(
          query_window, data, as.integer(pre_obj$data_size),
          as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd, pre_obj$query_mean[index],
          pre_obj$query_sd[index],
          grain = 4096
        )
      }
    } else {
      mass2_rcpp(
        pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
        as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd, pre_obj$query_mean[index],
        pre_obj$query_sd[index],
        n_workers = n_workers
      )
    },
    non_normalized = mass_weighted_rcpp(
      pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
      as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd,
      pre_obj$query_mean[index], pre_obj$query_sd[index], pre_obj$data_pre, rep(1, pre_obj$window_size), FALSE,
      n_workers = n_workers
    ),
    absolute = mass_absolute_rcpp(
      pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
      as.integer(pre_obj$window_size), pre_obj$sumx2, pre_obj$sumy2[index],
      n_workers = n_workers
    ),
    weighted = mass_weighted_rcpp(
      pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
      as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd,
      pre_obj$query_mean[index], pre_obj$query_sd[index], pre_obj$data_pre, pre_obj$weight, TRUE,
      n_workers = n_workers
    )
  )
}
