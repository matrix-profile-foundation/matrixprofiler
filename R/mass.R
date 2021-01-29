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
  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  if (!is.null(query)) {
    query <- as.numeric(query)
    checkmate::qassert(query, c("0", "N>=4"))
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
      result$type <- type
      return(result)
    },
    TRUE
  )

  if (is.null(query)) {
    query <- data
  }

  # Computation ------------------------------------
  "!DEBUG Computation"
  tryCatch(
    {
      result <- switch(
        type,
        normalized = mass_pre_rcpp(data, query, window_size),
        non_normalized = mass_pre_weighted_rcpp(data, query, window_size, rep(1, window_size)),
        absolute = mass_pre_abs_rcpp(data, query, window_size),
        weighted = mass_pre_weighted_rcpp(data, query, window_size, weights)
      )
    },
    error = print
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
#' @param n_workers An `integer` The number of threads using for computing. Defaults to `1`.
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
  type <- match.arg(pre_obj$type, c("normalized", "non_normalized", "absolute", "weighted"))
  version <- match.arg(version)
  data <- as.numeric(data)
  checkmate::qassert(data, "N+")
  query <- as.numeric(query)
  checkmate::qassert(query, "N+")
  index <- as.integer(checkmate::qassert(index, "X+(0,)"))

  if (length(data) != pre_obj$data_size) {
    stop("Argument `data` is not the same as computed in `pre_obj`.", call. = FALSE)
  }

  if (length(query) != (
    ifelse(pre_obj$type == "absolute", length(pre_obj$sumy2), length(pre_obj$query_mean)) +
      pre_obj$window_size - 1)) {
    stop("Argument `query` is not the same as computed in `pre_obj`.", call. = FALSE)
  }

  # Register anytime exit point ----------------------
  result <- NULL
  "!DEBUG Register anytime exit point"
  on.exit(
    {
      RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())

      if (is.null(result)) {
        return(invisible(NULL))
      } else {
        return(result)
      }
    },
    TRUE
  )

  # Computation ------------------------------------
  "!DEBUG Computation"
  tryCatch(
    {
      query_window <- query[index:(index + pre_obj$window_size - 1)]

      result <- switch(
        type,
        normalized = if (version == "v3") {
          if (n_workers > 1) {
            n_workers <- min(n_workers, RcppParallel::defaultNumThreads())
            RcppParallel::setThreadOptions(numThreads = n_workers)
            mass3_rcpp_parallel(query_window, data, as.integer(pre_obj$data_size),
              as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd, pre_obj$query_mean[index],
              pre_obj$query_sd[index],
              k = 4096
            )
          }
          else {
            mass3_rcpp(query_window, data, as.integer(pre_obj$data_size),
              as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd, pre_obj$query_mean[index],
              pre_obj$query_sd[index],
              k = 4096
            )
          }
        } else {
          mass2_rcpp(
            pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
            as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd, pre_obj$query_mean[index],
            pre_obj$query_sd[index]
          )
        },
        non_normalized = mass_weighted_rcpp(
          pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
          as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd,
          pre_obj$query_mean[index], pre_obj$query_sd[index], pre_obj$data_pre, rep(1, pre_obj$window_size), FALSE
        ),
        absolute = mass_absolute_rcpp(
          pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
          as.integer(pre_obj$window_size), pre_obj$sumx2, pre_obj$sumy2[index]
        ),
        weighted = mass_weighted_rcpp(
          pre_obj$data_fft, query_window, as.integer(pre_obj$data_size),
          as.integer(pre_obj$window_size), pre_obj$data_mean, pre_obj$data_sd,
          pre_obj$query_mean[index], pre_obj$query_sd[index], pre_obj$data_pre, pre_obj$weight, TRUE
        )
      )
    },
    error = print
  )
}
