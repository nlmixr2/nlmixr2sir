# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# Matrix helpers shared by the proposal machinery.

# Sampling Importance Resampling (SIR) ----------------------------------------
#
# Reference: Dosne et al. (2013), PAGE 22, Abstract 2907.
# Algorithm mirrors PsN's sir tool (psn_ref/sir, sir.pm, sir_userguide.tex).

# Helpers ----------------------------------------------------------------------

# Apply bounds elementwise across rows of a matrix; returns logical vector.
.sirInBounds <- function(mat, lower, upper) {
  apply(mat, 1L, function(x) all(x >= lower & x <= upper))
}

.sirCapCovCorrelation <- function(covMat, capCorrelation = 0.8) {
  checkmate::assertMatrix(covMat, mode = "numeric")
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)

  cov_mat <- (covMat + t(covMat)) / 2
  if (capCorrelation >= 1 || ncol(cov_mat) < 2L) {
    return(cov_mat)
  }

  param_names <- colnames(cov_mat)
  sd_vals <- sqrt(pmax(diag(cov_mat), 0))
  corr_mat <- suppressWarnings(cov2cor(cov_mat))
  corr_mat[!is.finite(corr_mat)] <- 0
  diag(corr_mat) <- 1

  lo <- lower.tri(corr_mat)
  corr_mat[lo] <- pmax(-capCorrelation, pmin(capCorrelation, corr_mat[lo]))
  corr_mat[upper.tri(corr_mat)] <- t(corr_mat)[upper.tri(corr_mat)]

  capped <- outer(sd_vals, sd_vals) * corr_mat
  dimnames(capped) <- list(param_names, param_names)
  capped
}

# Repair floating-point roundoff on a covariance that is ALREADY full rank.
#
# Both this and .sirCheckProposalRank() work in STANDARDIZED coordinates -- the
# correlation matrix -- rather than in the raw parameter coordinates.
#
# Why that matters: pharmacometric parameters do not share a unit. A clearance,
# a log-scale THETA, a small OMEGA element and a residual SD can differ by many
# orders of magnitude. An eigenvalue test in raw coordinates then measures the
# spread of the units as much as the spread of the information, so the same
# statistical problem gets a different verdict depending on how a parameter
# happens to be expressed. Concretely, the old code turned diag(c(1, 1e-14))
# into roughly diag(c(1, 1e-12)) -- inflating one parameter's variance a
# hundredfold because a different parameter happened to have variance one.
#
# Repairing the correlation matrix and mapping back with the original marginal
# standard deviations fixes that: the marginal variances are preserved exactly,
# and only the correlation structure is conditioned.
#
# A global rescale of the whole matrix cannot expose this, which is why the
# tests in test-sir-coordinate-invariance.R rescale each column independently.
.sirEnsurePosDef <- function(covMat, relTol = 1e-12, report = FALSE) {
  checkmate::assertMatrix(covMat, mode = "numeric")
  checkmate::assertNumber(relTol, lower = 0, finite = TRUE)
  dim_names <- dimnames(covMat)
  cov_mat <- (covMat + t(covMat)) / 2

  sds <- sqrt(pmax(diag(cov_mat), 0))
  pos <- is.finite(sds) & sds > 0

  # A zero-variance coordinate carries no uncertainty and cannot be
  # standardized. The Wishart fallback produces one for an OMEGA element
  # estimated at zero, so this is a real case rather than a defensive branch.
  # Such a parameter must stay fixed, not acquire uncertainty from the repair,
  # so its row and column are held at zero and it sits out the correlation
  # step entirely.
  if (any(!pos)) {
    cov_mat[!pos, ] <- 0
    cov_mat[, !pos] <- 0
  }

  if (!any(pos)) {
    cli::cli_abort(c(
      "Proposal covariance has no positive variance.",
      "i" = "Every parameter would be fixed at its point estimate."
    ))
  }

  sub <- cov_mat[pos, pos, drop = FALSE]
  s <- sds[pos]
  scale_mat <- outer(s, s)
  corr <- sub / scale_mat
  corr <- (corr + t(corr)) / 2

  eig <- eigen(corr, symmetric = TRUE)
  if (!all(is.finite(eig$values))) {
    cli::cli_abort("Proposal covariance has non-finite eigenvalues.")
  }
  max_eigen <- max(eig$values)
  if (max_eigen <= 0) {
    cli::cli_abort(c(
      "Proposal correlation matrix has no positive eigenvalue.",
      "i" = "Every retained vector is identical in every parameter."
    ))
  }

  floor_at <- relTol * max_eigen
  if (min(eig$values) >= floor_at) {
    dimnames(cov_mat) <- dim_names
    if (report) {
      return(list(
        covMat = cov_mat,
        adjusted = FALSE,
        floor = floor_at,
        method = "correlation-eigenvalue-floor",
        threshold = relTol,
        magnitude = 0
      ))
    }
    return(cov_mat)
  }

  values <- pmax(eig$values, floor_at)
  repaired <- eig$vectors %*% diag(values, nrow = length(values)) %*%
    t(eig$vectors)
  repaired <- (repaired + t(repaired)) / 2
  # Flooring the eigenvalues moves the diagonal off one. Renormalising to unit
  # diagonal is what guarantees the marginal variances come back untouched.
  repaired <- stats::cov2cor(repaired)
  cov_mat[pos, pos] <- repaired * scale_mat
  cov_mat <- (cov_mat + t(cov_mat)) / 2
  dimnames(cov_mat) <- dim_names

  if (report) {
    # A repaired covariance is a fact about the run that belongs in its
    # provenance, even when the repair is only roundoff.
    return(list(
      covMat = cov_mat,
      adjusted = TRUE,
      floor = floor_at,
      method = "correlation-eigenvalue-floor",
      threshold = relTol,
      # Largest absolute change to any entry, in the original coordinates.
      # The flag alone says a repair happened; this says whether it mattered.
      magnitude = max(abs(cov_mat - covMat))
    ))
  }
  cov_mat
}

# Refuse a set of parameter vectors that cannot support a full-rank covariance.
#
# The empirical covariance of m vectors in p dimensions has rank at most m - 1,
# so m > p is necessary; repeated or collinear vectors can leave it deficient
# even when m > p. Forcing such a matrix positive definite does not recover the
# missing information -- it fabricates variance in unsupported directions and
# the next iteration then proposes along them. Stopping is the honest
# behaviour: the fix is more retained vectors, not a repaired matrix.
#
# Rank is judged on the CORRELATION matrix, for the reason given above
# .sirEnsurePosDef(): in raw coordinates the verdict moves with the parameter
# units. A full-rank two-column sample whose second column was expressed in
# units 1e-6 smaller was previously rejected as rank one.
.sirCheckProposalRank <- function(mat, what = "retained", rankTol = 1e-8,
                                  onDeficient = c("abort", "repair")) {
  onDeficient <- match.arg(onDeficient)
  checkmate::assertMatrix(mat, mode = "numeric", min.rows = 1L, min.cols = 1L)
  m <- nrow(mat)
  p <- ncol(mat)

  if (m <= p) {
    cli::cli_abort(c(
      "Too few {what} vectors to estimate a full-rank proposal covariance.",
      "x" = "{m} vector{?s} for {p} parameter{?s}; the covariance has rank at most {m - 1}.",
      "i" = "More than {p} {what} vectors are needed.",
      "i" = "Increase {.arg nResample} (or supply more vectors) and rerun."
    ))
  }

  sds <- apply(mat, 2L, stats::sd)
  nm <- colnames(mat)
  if (is.null(nm)) {
    nm <- paste0("column ", seq_len(p))
  }
  # A constant column is an unsupported direction no matter what its units are,
  # and it cannot be standardized. Reported by name, because "rank deficient"
  # alone leaves the user hunting for which parameter never moved.
  constant <- !is.finite(sds) | sds <= 0
  if (any(constant)) {
    stuck <- nm[constant]
    cli::cli_abort(c(
      "Some {what} parameters do not vary, so the proposal covariance is rank deficient.",
      "x" = "Constant across every {what} vector: {.val {stuck}}.",
      "i" = "A parameter with no spread cannot be given a proposal distribution.",
      "i" = "Check the rejection diagnostics and {.arg capResampling}, and raise {.arg nResample}."
    ))
  }

  # cor() is exactly the covariance of the standardized columns, so this is the
  # rank of the centred sample with every parameter put on a common footing.
  corr <- stats::cor(mat)
  corr <- (corr + t(corr)) / 2
  eig <- eigen(corr, symmetric = TRUE, only.values = TRUE)$values
  if (!all(is.finite(eig)) || max(eig) <= 0) {
    cli::cli_abort(c(
      "The {what} vectors have rank zero: no parameter varies across them.",
      "i" = "Check the resampling diagnostics and {.arg capResampling}."
    ))
  }

  rank <- sum(eig > rankTol * max(eig))
  deficient <- p - rank
  if (rank < p) {
    # Two causes, and the remedies are opposite. The message used to offer only
    # the first, which was wrong on QR model N029 and sent the reader looking
    # in the wrong place: that fit's own covariance is numerically singular
    # (smallest eigenvalue 3.22e-10), so its proposal was degenerate before any
    # sampling and raising nResample from 200 to 500 changed nothing.
    causes <- c(
      "i" = "With {m} vector{?s} for {p} parameter{?s}, repeated or collinear draws are one cause: raise {.arg nResample}, or lower {.arg capResampling} to reduce repeats.",
      "i" = "The other is a degenerate proposal, where the fit itself does not identify every direction. Check {.code eigen(fit$cov)$values}: if the smallest is at or near zero, more vectors cannot help."
    )
    if (identical(onDeficient, "abort")) {
      cli::cli_abort(c(
        "The {what} vectors are rank deficient.",
        "x" = "Numerical rank {rank} for {p} parameters: {deficient} direction{?s} unsupported.",
        causes,
        "i" = "Forcing this positive definite would invent uncertainty the sample does not support, so this stops rather than guessing.",
        "i" = "To proceed anyway with the {deficient} unsupported direction{?s} pinned to a token variance, set {.code runSIRControl(rankDeficiency = \"repair\")}."
      ))
    }
    cli::cli_warn(c(
      "Proceeding with a rank-deficient {what} covariance because {.code rankDeficiency = \"repair\"}.",
      "x" = "Numerical rank {rank} for {p} parameters: {deficient} direction{?s} unsupported.",
      "!" = "This INVENTS uncertainty the sample does not support: {deficient} unsupported direction{?s} get a token variance, and later iterations will propose along {?it/them}.",
      "i" = "Intervals for parameters loading on those {deficient} direction{?s} are not evidence from the data.",
      causes
    ))
  }

  invisible(list(rank = rank, nVectors = m, nParams = p, deficient = deficient))
}


# PsN's math.pm round(): half away from zero, truncating toward zero first.
# R's round() is banker's rounding and differs on an exact .5 -- round(20.5)
# is 20 in R and 21 in PsN -- so the sample-count adjustments cannot use it
# and stay bit-comparable with PsN's own unit-test oracles.
.sirRound <- function(x) {
  intPart <- trunc(x)
  rem <- x - intPart
  as.integer(ifelse(
    rem >= 0,
    ifelse(rem >= 0.5, intPart + 1, intPart),
    ifelse(abs(rem) >= 0.5, intPart - 1, intPart)
  ))
}
