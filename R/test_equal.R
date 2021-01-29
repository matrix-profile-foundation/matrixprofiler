#' Test similarity
#'
#' Compute if two Matrix Profiles and Profile Indexes are similar. Used only for internal purposes.
#'
#' @param mpa_obj First object
#' @param mpb_obj Second object
#' @param tolerance tolerance
#'
#' @keywords internal

test_equal <- function(mpa_obj, mpb_obj, tolerance = .Machine$double.eps) {
  msg <- NULL
  can_compare <- TRUE

  obja_name <- substitute(mpa_obj)
  objb_name <- substitute(mpb_obj)

  #### Test for lengths ----
  "!DEBUG Test for lengths"

  if (is.null(mpa_obj$matrix_profile) || is.null(mpb_obj$matrix_profile)) {
    msg <- paste(msg, "One of the matrix profiles doesn't exists\n")
    can_compare <- FALSE
  } else if (!(length(mpa_obj$matrix_profile) == length(mpb_obj$matrix_profile))) {
    msg <- paste(
      msg, "Matrix Profile lengths differs,", obja_name, "=", length(mpa_obj$matrix_profile),
      objb_name, "=", length(mpb_obj$matrix_profile), "\n"
    )
    can_compare <- FALSE
  }

  if (is.null(mpa_obj$profile_index) || is.null(mpb_obj$profile_index)) {
    msg <- paste(msg, "One of the matrix profiles doesn't exists\n")
    can_compare <- FALSE
  } else if (!(length(mpa_obj$profile_index) == length(mpb_obj$profile_index))) {
    msg <- paste(
      msg, "Profile Indexes lengths differs,", obja_name, "=", length(mpa_obj$profile_index),
      objb_name, "=", length(mpb_obj$profile_index), "\n"
    )
    can_compare <- FALSE
  }

  #### Test matrix profile ----
  "!DEBUG Test matrix profile"

  if (any(is.null(mpa_obj$matrix_profile))) {
    msg <- paste(msg, obja_name, "has NULL values in Matrix Profile\n")
    can_compare <- FALSE
  }

  if (any(is.null(mpb_obj$matrix_profile))) {
    msg <- paste(msg, objb_name, "has NULL values in Matrix Profile\n")
    can_compare <- FALSE
  }

  if (anyNA(mpa_obj$matrix_profile)) {
    msg <- paste(msg, obja_name, "has NA's values in Matrix Profile\n")
    can_compare <- FALSE
  } else if (any(mpa_obj$matrix_profile < 0)) {
    msg <- paste(msg, obja_name, "has Negative values in Matrix Profile\n")
    can_compare <- FALSE
  }

  if (anyNA(mpb_obj$matrix_profile)) {
    msg <- paste(msg, objb_name, "has NA's values in Matrix Profile\n")
    can_compare <- FALSE
  } else if (any(mpb_obj$matrix_profile < 0)) {
    msg <- paste(msg, objb_name, "has Negative values in Matrix Profile\n")
    can_compare <- FALSE
  }

  if (anyNA(mpa_obj$profile_index)) {
    msg <- paste(msg, obja_name, "has NA's values in Profile Index\n")
    can_compare <- FALSE
  } else if (any(mpa_obj$profile_index < 0)) {
    msg <- paste(msg, obja_name, "has Negative values in Profile Index\n")
    can_compare <- FALSE
  }

  if (anyNA(mpb_obj$profile_index)) {
    msg <- paste(msg, objb_name, "has NA's values in Profile Index\n")
    can_compare <- FALSE
  } else if (any(mpb_obj$profile_index < 0)) {
    msg <- paste(msg, objb_name, "has Negative values in Profile Index\n")
    can_compare <- FALSE
  }

  if (can_compare) {
    #### Do Compare matrix profile ----
    "!!DEBUG Do Compare matrix profile"

    mp_diff <- abs(mpa_obj$matrix_profile - mpb_obj$matrix_profile)
    pi_diff <- abs(mpa_obj$profile_index - mpb_obj$profile_index)

    diff_idxs <- (mp_diff > tolerance)
    if (any(diff_idxs)) {
      msg <- paste(msg, "> ", objb_name, "> ***Matrix Profile differs***\n")
      if (sum(diff_idxs) > 10) {
        msg <- paste(msg, ">> Too much indexes are different\n")
      } else {
        idxs <- base::which(diff_idxs)
        difs <- mp_diff[idxs]
        msg <- paste(msg, ">>", paste(idxs, round(difs, 5), collapse = ", ", sep = ":"), "\n")
      }
    }

    #### Test profile index ----
    "!!DEBUG Test profile index"

    diff_idxs <- (pi_diff > 0.01)
    if (any(diff_idxs)) {
      msg <- paste(msg, "> ", objb_name, "> Profile Index differs\n")
      if (sum(diff_idxs) > 10) {
        msg <- paste(msg, ">> Too much indexes are different\n")
      } else {
        idxs <- base::which(diff_idxs)
        difs <- pi_diff[idxs]
        msg <- paste(msg, ">>", paste(idxs, difs, collapse = ", ", sep = ":"), "\n")
      }
    }
  }

  if (!is.null(msg)) {
    message(msg)
    return(FALSE)
  } else {
    return(TRUE)
  }
}
