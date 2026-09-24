# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# Box-Cox transform and the empirical proposal update.

# Step 7 -----------------------------------------------------------------------

#' Box-Cox transform a vector, optionally estimating lambda
#'
#' Applies the Box-Cox power transformation after shifting `x` by `delta` to
#' ensure strict positivity.  When `lambda` is `NULL` it is estimated by
#' maximising the Pearson correlation between the sorted transformed values and
#' their expected normal order statistics (normal scores), using
#' `stats::optimize()` over the interval \[−3, 3\].
#'
#' Transformation:
#' * lambda ≠ 0: `(x + delta)^lambda − 1) / lambda`
#' * lambda = 0: `log(x + delta)`
#'
#' @param x Numeric vector (length ≥ 2), no missing values.
#' @param lambda `NULL` to estimate, or a finite numeric scalar.
#' @param delta `NULL` to compute as `|min(c(x, mustAdmit))| + 1e-6`, or a
#'   non-negative scalar.  Must ensure `x + delta > 0`.
#' @param mustAdmit Optional numeric. Extra values the chosen shift must also
#'   place strictly inside the transform domain, without being transformed
#'   themselves. The next proposal centre is passed here: under `recenter` it
#'   is the best candidate, which need not be one of the retained vectors and
#'   can sit below their minimum.
#' @return Named list: `transformed` (numeric vector), `lambda` (scalar),
#'   `delta` (scalar).
#' @noRd
sirBoxCox <- function(x, lambda = NULL, delta = NULL, mustAdmit = NULL) {
  checkmate::assertNumeric(x, min.len = 2L, any.missing = FALSE, finite = TRUE)

  if (is.null(delta)) {
    # mustAdmit widens the shift and nothing else: it is not transformed here.
    # Choosing the shift from the sample alone let a recentred mean below the
    # sample minimum fall outside the domain, and .sirBcTransformMu() then
    # aborted on a run that was otherwise proceeding normally.
    checkmate::assertNumeric(
      mustAdmit,
      any.missing = FALSE, finite = TRUE, null.ok = TRUE
    )
    delta <- abs(min(c(x, mustAdmit))) + 1e-6
  }
  checkmate::assertNumber(delta, lower = 0, finite = TRUE)

  x_shifted <- x + delta
  if (any(x_shifted <= 0)) {
    cli::cli_abort(
      "{.arg delta} is too small: {.code x + delta} must be strictly positive."
    )
  }

  .bc <- function(xs, lam) {
    if (abs(lam) < 1e-10) log(xs) else (xs^lam - 1) / lam
  }

  if (is.null(lambda)) {
    if (stats::var(x) < .Machine$double.eps) {
      lambda <- 1 # degenerate: all identical, any lambda works
    } else {
      nscores <- qnorm(stats::ppoints(length(x)))
      obj <- function(lam) {
        xt <- sort(.bc(x_shifted, lam))
        v <- suppressWarnings(cor(xt, nscores))
        if (is.na(v)) 0 else -v
      }
      lambda <- stats::optimize(obj, interval = c(-3, 3))$minimum
    }
  }
  checkmate::assertNumber(lambda, finite = TRUE)

  list(transformed = .bc(x_shifted, lambda), lambda = lambda, delta = delta)
}

#' Invert a Box-Cox transformation
#'
#' Recovers the original scale given the transformed values, lambda, and delta
#' returned by `sirBoxCox()`.
#'
#' @param x_transformed Numeric vector of transformed values.
#' @param lambda Scalar lambda used in the forward transform.
#' @param delta Scalar shift used in the forward transform.
#' @return Numeric vector on the original scale.
#' @noRd
sirBoxCoxInverse <- function(x_transformed, lambda, delta) {
  checkmate::assertNumeric(x_transformed, any.missing = FALSE, finite = TRUE)
  checkmate::assertNumber(lambda, finite = TRUE)
  checkmate::assertNumber(delta, finite = TRUE)

  if (abs(lambda) < 1e-10) {
    exp(x_transformed) - delta
  } else {
    base <- lambda * x_transformed + 1
    if (any(base <= 0)) {
      cli::cli_abort(
        "Inverse Box-Cox: {.code lambda * x_transformed + 1} must be positive."
      )
    }
    base^(1 / lambda) - delta
  }
}

# Step 8 -----------------------------------------------------------------------

#' Update the proposal covariance from resampled parameter vectors
#'
#' Computes the empirical covariance of a resampled parameter matrix.  When
#' `boxcox = TRUE`, each column is first Box-Cox transformed via `sirBoxCox()`
#' to reduce skewness before the covariance is calculated; the resulting
#' covariance is in the transformed space and the per-column lambda/delta values
#' are returned so callers can back-transform sampled points if needed.
#'
#' @param resampledMat Numeric matrix with at least 2 rows (resampled parameter
#'   vectors, one per row).  Column names are preserved in the output.
#' @param boxcox Logical.  If `TRUE` (default), apply column-wise Box-Cox
#'   before computing covariance.
#' @param capCorrelation Numeric in \[0, 1\]. Caps absolute correlations in
#'   the updated proposal. Default `0.8`.
#' @return Named list:
#'   \describe{
#'     \item{`covMat`}{Symmetric positive-(semi)definite covariance matrix.}
#'     \item{`boxcoxParams`}{Data frame with columns `param`, `lambda`, `delta`
#'       (one row per column of `resampledMat`), or `NULL` when
#'       `boxcox = FALSE`.}
#'   }
#' @noRd
sirUpdateProposal <- function(
  resampledMat,
  boxcox = TRUE,
  capCorrelation = 0.8,
  centre = NULL,
  rankDeficiency = c("abort", "repair")
) {
  rankDeficiency <- match.arg(rankDeficiency)
  checkmate::assertMatrix(
    resampledMat,
    mode = "numeric",
    min.rows = 2L,
    min.cols = 1L
  )
  checkmate::assertFlag(boxcox)
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)

  # Before anything else: the retained vectors must be able to support a
  # full-rank covariance. Checked on the untransformed matrix, because that is
  # where the statistical support actually lives; the Box-Cox map below is
  # per-coordinate and monotone, so it cannot add support.
  rankInfo <- .sirCheckProposalRank(
    resampledMat,
    what = "retained",
    onDeficient = rankDeficiency
  )
  rankRepaired <- isTRUE(rankInfo$deficient > 0L)

  param_names <- colnames(resampledMat)
  n_col <- ncol(resampledMat)

  if (boxcox) {
    # The shift for each column must admit the centre that will later be
    # transformed with these same parameters, so the two cannot disagree.
    bc_list <- lapply(seq_len(n_col), function(j) {
      nm <- if (is.null(param_names)) NULL else param_names[[j]]
      admit <- if (is.null(centre) || is.null(nm) || !nm %in% names(centre)) {
        NULL
      } else {
        unname(centre[[nm]])
      }
      sirBoxCox(resampledMat[, j], mustAdmit = admit)
    })

    trans_mat <- matrix(
      unlist(lapply(bc_list, `[[`, "transformed")),
      nrow = nrow(resampledMat),
      ncol = n_col,
      dimnames = list(NULL, param_names)
    )

    bc_params <- data.frame(
      param = if (is.null(param_names)) seq_len(n_col) else param_names,
      lambda = vapply(bc_list, `[[`, numeric(1L), "lambda"),
      delta = vapply(bc_list, `[[`, numeric(1L), "delta"),
      stringsAsFactors = FALSE
    )
  } else {
    trans_mat <- resampledMat
    bc_params <- NULL
  }

  cov_mat <- cov(trans_mat)
  dimnames(cov_mat) <- list(param_names, param_names)
  cov_mat <- .sirCapCovCorrelation(cov_mat, capCorrelation = capCorrelation)
  # A deficient covariance needs a floor big enough to matter. The default
  # 1e-12 is a roundoff guard and would leave an unsupported direction at
  # whatever near-zero value it already had, which is not positive definite in
  # any useful sense. Repairing lifts it to the same relative tolerance the
  # rank verdict was made at, so the direction becomes a token rather than a
  # singularity -- and stays orders of magnitude below the supported ones.
  repaired <- if (rankRepaired) {
    .sirEnsurePosDef(cov_mat, relTol = 1e-8, report = TRUE)
  } else {
    .sirEnsurePosDef(cov_mat, report = TRUE)
  }

  list(
    covMat = repaired$covMat,
    boxcoxParams = bc_params,
    posDefAdjusted = repaired$adjusted,
    rankRepaired = rankRepaired
  )
}

# Step 9 -----------------------------------------------------------------------

# Transform a named mu vector into the sampling space using stored Box-Cox params.
.sirBcTransformMu <- function(mu, bc_state) {
  if (is.null(bc_state)) {
    return(mu)
  }
  bc_mu <- vapply(
    seq_along(mu),
    function(j) {
      nm <- names(mu)[j]
      row <- bc_state[bc_state$param == nm, ]
      if (nrow(row) != 1L) {
        cli::cli_abort(
          "Missing Box-Cox parameters for SIR parameter {.val {nm}}."
        )
      }
      xs <- mu[j] + row$delta
      if (!is.finite(xs) || xs <= 0) {
        cli::cli_abort(
          "Box-Cox shift for SIR parameter {.val {nm}} is not positive."
        )
      }
      lam <- row$lambda
      if (abs(lam) < 1e-10) log(xs) else (xs^lam - 1) / lam
    },
    numeric(1L)
  )
  setNames(bc_mu, names(mu))
}

# log|det J_T(x)| for the Box-Cox map y = T(x), evaluated per row of `mat`.
#
# Box-Cox is applied one coordinate at a time, so J_T is diagonal:
#
#   T(x)   = ((x + delta)^lambda - 1) / lambda   (log(x + delta) when lambda = 0)
#   dT/dx  = (x + delta)^(lambda - 1)            (both cases)
#
# hence log|det J_T(x)| = sum_j (lambda_j - 1) * log(x_j + delta_j).
#
# This is what turns the sampling density q_y into the density actually induced
# on the original parameter scale, q_x(x) = q_y(T(x)) * |det J_T(x)|. The
# importance weights divide by q_x, not q_y, so that the retained sample
# targets the original-scale normalized likelihood rather than a
# parameterization-dependent tilt of it. PsN omits this term; see NEWS.
#
# `mat` holds ORIGINAL-scale values (the back-transformed draws), not the
# sampled coordinates. Rows whose shifted value is not strictly positive get
# `NA`, which the weight calculation then drops along with other failures.
.sirBcLogJacobian <- function(mat, bc_state) {
  if (is.null(bc_state)) {
    return(rep(0, nrow(mat)))
  }
  out <- numeric(nrow(mat))
  for (j in seq_len(ncol(mat))) {
    nm <- colnames(mat)[j]
    row <- bc_state[bc_state$param == nm, ]
    if (nrow(row) != 1L) {
      cli::cli_abort(
        "Missing Box-Cox parameters for SIR parameter {.val {nm}}."
      )
    }
    # lambda == 1 is a pure shift: no contribution, and skipping it keeps the
    # common no-transform case exactly zero rather than zero-ish.
    if (isTRUE(all.equal(row$lambda, 1))) {
      next
    }
    shifted <- mat[, j] + row$delta
    contrib <- ifelse(
      is.finite(shifted) & shifted > 0,
      (row$lambda - 1) * log(shifted),
      NA_real_
    )
    out <- out + contrib
  }
  out
}

.sirBcInverseMatrix <- function(mat, bc_state) {
  if (is.null(bc_state)) {
    return(mat)
  }
  out <- mat
  for (j in seq_len(ncol(mat))) {
    nm <- colnames(mat)[j]
    row <- bc_state[bc_state$param == nm, ]
    if (nrow(row) != 1L) {
      cli::cli_abort(
        "Missing Box-Cox parameters for SIR parameter {.val {nm}}."
      )
    }
    out[, j] <- vapply(
      mat[, j],
      function(v) {
        tryCatch(
          sirBoxCoxInverse(v, row$lambda, row$delta),
          error = function(e) NA_real_
        )
      },
      numeric(1L)
    )
  }
  out
}
