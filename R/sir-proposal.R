# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# Building the proposal: reading fit$cov, fallback uncertainty, assembly,
# bounds, and drawing from it.

# Step 1 -----------------------------------------------------------------------

#' Extract and inflate the proposal covariance matrix for SIR
#'
#' Pulls the parameter covariance matrix from a fitted nlmixr2 model and
#' optionally inflates variances and caps correlations to widen the proposal
#' distribution for the first SIR iteration.
#'
#' Inflation is applied per parameter type: `thetaInflation` scales the
#' variance of every THETA parameter, `omegaInflation` scales every OMEGA
#' element.  Off-diagonal elements are recomputed so that correlations are
#' preserved after inflation.  `capCorrelation` then hard-clips all absolute
#' pairwise correlations.
#'
#' @param fit An nlmixr2 fit with a successful covariance step.
#' @param thetaInflation Non-negative scalar. Variance multiplier for THETA
#'   parameters. Default `1` (no inflation).
#' @param omegaInflation Non-negative scalar. Variance multiplier for OMEGA
#'   elements. Default `1`.
#' @param sigmaInflation Not used in nlmixr2 (residual error parameters are
#'   THETAs). Retained for API compatibility with PsN. Default `1`.
#' @param capCorrelation Numeric in \[0, 1\]. All absolute pairwise
#'   correlations are capped at this value after inflation. Set to `1` to
#'   disable. Default `0.8`.
#' @param cov The covariance to build the proposal from. Defaults to
#'   `fit$cov`; a seed covariance from elsewhere (another entry of
#'   `fit$covList`, say) is named with nlmixr2est's covariance names or SIR's.
#' @return A named, symmetric covariance matrix with the same row/column order
#'   as `cov`.
#' @importFrom stats cor
#' @noRd
sirGetProposalCov <- function(
  fit,
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8,
  cov = fit$cov
) {
  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)

  cov_mat <- cov
  if (is.null(cov_mat) || nrow(cov_mat) == 0L) {
    # There is no automatic route out of this one. A missing OMEGA or SIGMA
    # block can be filled, because the Wishart-style approximation derives it
    # from the OMEGA estimates and the subject count -- real information. THETA
    # uncertainty has no such source: a fit that never ran a covariance step
    # carries none, and inventing it from an assumed RSE would fabricate
    # exactly the quantity SIR is meant to measure. So the user has to say
    # where the proposal comes from.
    cli::cli_abort(c(
      "No covariance matrix in {.arg fit}, so SIR has no proposal to start from.",
      "i" = "A failed covariance step is a normal reason to run SIR, so supply the proposal another way:",
      "*" = "{.code runSIRControl(rseTheta =)} builds a diagonal proposal from relative standard errors.",
      "*" = "{.code runSIRControl(covmatInput =)} takes a covariance matrix, a NONMEM {.file .cov} file, or {.val identity}.",
      "*" = "{.code runSIRControl(rawresInput =)} seeds the proposal from existing parameter vectors.",
      "i" = "Or refit with a covariance step, for example {.code foceiControl(covMethod = \"r\")}."
    ))
  }

  ps <- .sirParamSpace(fit)
  # fullCovName, not covName: a seed that is not fit$cov may carry parameters
  # fit$cov does not. SIR names are accepted too.
  idx <- match(rownames(cov_mat), ps$fullCovName)
  idx[is.na(idx)] <- match(rownames(cov_mat)[is.na(idx)], ps$sirName)
  if (anyNA(idx)) {
    cli::cli_abort(c(
      "The covariance carries parameter{?s} SIR cannot place: {.val {rownames(cov_mat)[is.na(idx)]}}.",
      "i" = "This means {.fn .sirParamSpace} and the covariance disagree about the parameter set."
    ))
  }

  # One inflation factor per row/col. Residual-error parameters get
  # sigmaInflation: they used to fall into the theta branch, which made
  # sigmaInflation unreachable.
  inflation <- .sirInflationVector(
    ps,
    thetaInflation = thetaInflation,
    omegaInflation = omegaInflation,
    sigmaInflation = sigmaInflation
  )[idx]

  # Preserve correlations: new SD_i = old SD_i * sqrt(inflation_i), which is
  # the same as scaling each entry by sqrt(inflation_i * inflation_j) and
  # avoids a cov2cor() that would fail on a zero variance.
  new_cov <- cov_mat * outer(sqrt(inflation), sqrt(inflation))

  # Renamed to SIR names here so that the fit$cov naming convention stays
  # confined to .sirParamSpace() and this one function.
  sir_names <- ps$sirName[idx]
  dimnames(new_cov) <- list(sir_names, sir_names)

  .sirCapCovCorrelation(new_cov, capCorrelation = capCorrelation)
}

# Fallback uncertainty for parameters `fit$cov` does not carry.
#
# nlmixr2sse has its own OMEGA-draw machinery (per-block inverse-Wishart and
# a log-Cholesky joint draw). Promoting a shared implementation to
# nlmixr2utils was considered and declined for now (PLAN.md 8.2): the cost
# would land in two packages that have users, for two callers. Revisit if a
# third caller appears.
#
# Reached only when the covariance step failed, `covMethod = ""`, or the user
# asks for `omegaFallback = "wishart"`. Under nlmixr2est 7 defaults `fit$cov`
# carries OMEGA and none of this is used.
#
# OMEGA uses the Wishart-style approximation with `df = nsub - 1` by default:
#   Var(omega_ii) = 2 * omega_ii^2 / df
#   Var(omega_ij) = (omega_ii * omega_jj + omega_ij^2) / df
# THETA and sigma use the fit's own SE when `parFixedDf` reports one, and
# `sigmaFallbackRse` percent of the estimate otherwise.
.sirFallbackSe <- function(
  fit,
  ps = .sirParamSpace(fit),
  sigmaFallbackRse = 30,
  omegaDf = NULL,
  nSub = NULL,
  parFixedSe = TRUE
) {
  checkmate::assertNumber(sigmaFallbackRse, lower = 0, finite = TRUE)
  se <- stats::setNames(rep(NA_real_, nrow(ps)), ps$sirName)

  isOmega <- ps$kind %in% c("omegaDiag", "omegaOffdiag")
  if (any(isOmega)) {
    if (is.null(nSub)) {
      nSub <- .sirNSubjects(fit)
    }
    df <- if (is.null(omegaDf)) max(nSub - 1L, 1L) else omegaDf
    checkmate::assertNumber(df, lower = 1, finite = TRUE)
    om <- fit$omega
    r <- ps$neta1[isOmega]
    cc <- ps$neta2[isOmega]
    se[isOmega] <- sqrt(ifelse(
      r == cc,
      2 * om[cbind(r, r)]^2 / df,
      (om[cbind(r, r)] * om[cbind(cc, cc)] + om[cbind(r, cc)]^2) / df
    ))
  }

  isThetaish <- ps$kind %in% c("theta", "sigma")
  if (any(isThetaish)) {
    # parFixedDf reports the installed covariance. When the proposal is seeded
    # from a different one, those SEs describe neither, so they are not used.
    pf <- if (parFixedSe) {
      tryCatch(fit$parFixedDf, error = function(e) NULL)
    } else {
      NULL
    }
    se[isThetaish] <- vapply(
      which(isThetaish),
      function(i) {
        nm <- ps$sirName[i]
        if (!is.null(pf) && nm %in% rownames(pf) && !is.na(pf[nm, "SE"])) {
          return(pf[nm, "SE"])
        }
        abs(ps$est[i]) * sigmaFallbackRse / 100
      },
      numeric(1L)
    )
  }
  se
}

# Reconstruct the full symmetric omega matrix from a named numeric vector of
# SIR parameter values. Entries absent from `vals` keep their value in
# `base_omega`.
.sirReconstructOmega <- function(ps, vals, base_omega) {
  om <- ps[ps$kind %in% c("omegaDiag", "omegaOffdiag"), , drop = FALSE]
  mat <- base_omega
  for (i in seq_len(nrow(om))) {
    nm <- om$sirName[i]
    if (nm %in% names(vals)) {
      mat[om$neta1[i], om$neta2[i]] <- vals[[nm]]
      mat[om$neta2[i], om$neta1[i]] <- vals[[nm]]
    }
  }
  mat
}

# Derive per-subject count for chi-squared omega SE approximation.
.sirNSubjects <- function(fit) {
  n <- tryCatch(as.integer(fit$nsub), error = function(e) NULL)
  if (!is.null(n) && length(n) == 1L && !is.na(n)) {
    return(n)
  }
  dat <- tryCatch(fit$origData, error = function(e) NULL)
  if (!is.null(dat) && "ID" %in% names(dat)) {
    return(length(unique(dat$ID)))
  }
  cli::cli_abort("Cannot determine subject count from {.arg fit}.")
}

# Build the initial (or continuing) SIR proposal.
#
# Two routes, reported back as `omegaRoute`:
#
#   "cov"      -- `proposalCov` already covers the whole SIR parameter vector,
#                 so it is used as-is. Under nlmixr2est 7 defaults that
#                 includes OMEGA *and its correlations with THETA*, which the
#                 old block-diagonal construction discarded even when they
#                 were available. This is also the route taken from iteration
#                 2 onwards, where the proposal is the previous empirical
#                 covariance.
#   "wishart"  -- anything `proposalCov` does not cover is filled in on the
#                 diagonal from `.sirFallbackSe()`. Needed for fits with
#                 `covFull = FALSE`, a failed covariance step, or
#                 `covMethod = ""`, and forced by `omegaFallback = "wishart"`.
.sirInitialProposal <- function(
  fit,
  mu,
  proposalCov,
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8,
  omegaFallback = c("cov", "wishart"),
  sigmaFallbackRse = 30,
  omegaDf = NULL,
  ps = .sirParamSpace(fit),
  parFixedSe = TRUE
) {
  omegaFallback <- match.arg(omegaFallback)
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, min.len = 1L)
  checkmate::assertMatrix(proposalCov, mode = "numeric")
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)

  fullNames <- ps$sirName
  nFull <- length(fullNames)

  muFull <- .sirProposalMu(fit, ps)
  shared <- intersect(names(mu), fullNames)
  muFull[shared] <- mu[shared]
  if (anyNA(muFull)) {
    cli::cli_abort(
      "No estimate available for SIR parameter{?s} {.val {fullNames[is.na(muFull)]}}."
    )
  }

  covNames <- intersect(colnames(proposalCov), fullNames)
  if (identical(omegaFallback, "wishart")) {
    # Keep only the THETA/sigma block; OMEGA is rebuilt from the fallback.
    covNames <- intersect(
      covNames,
      ps$sirName[ps$kind %in% c("theta", "sigma")]
    )
  }
  missingNames <- setdiff(fullNames, covNames)
  route <- if (length(missingNames) == 0L) "cov" else "wishart"

  covFull <- matrix(
    0,
    nrow = nFull,
    ncol = nFull,
    dimnames = list(fullNames, fullNames)
  )
  if (length(covNames) > 0L) {
    covFull[covNames, covNames] <- proposalCov[covNames, covNames, drop = FALSE]
  }

  if (length(missingNames) > 0L) {
    se <- .sirFallbackSe(
      fit,
      ps = ps,
      sigmaFallbackRse = sigmaFallbackRse,
      omegaDf = omegaDf,
      parFixedSe = parFixedSe
    )
    if (anyNA(se[missingNames])) {
      naNames <- missingNames[is.na(se[missingNames])]
      cli::cli_abort(c(
        "No uncertainty available for SIR parameter{?s} {.val {naNames}}.",
        "i" = "Missing OMEGA and SIGMA blocks can be approximated, but THETA uncertainty cannot be derived from the fit alone.",
        "i" = "Supply it with {.code runSIRControl(rseTheta =)}, {.code covmatInput}, or {.code rawresInput}, or refit with a covariance step."
      ))
    }
    diag(covFull)[match(missingNames, fullNames)] <- se[missingNames]^2
  }

  # Inflate by kind, preserving correlations. Applied to the whole assembled
  # matrix, not just the fallback blocks: under covFull the proposal arrives
  # as one covariance and inflating only the fallback would silently drop
  # thetaInflation/omegaInflation entirely.
  inflation <- .sirInflationVector(
    ps,
    thetaInflation = thetaInflation,
    omegaInflation = omegaInflation,
    sigmaInflation = sigmaInflation
  )
  if (any(inflation != 1)) {
    # cov_ij * sqrt(infl_i * infl_j) is the same rescaling as
    # corr_ij * (sd_i sqrt(infl_i)) * (sd_j sqrt(infl_j)), but goes through no
    # cov2cor() and so survives a zero variance -- which a Wishart fallback
    # produces for an OMEGA element estimated at zero.
    covFull <- covFull * outer(sqrt(inflation), sqrt(inflation))
    dimnames(covFull) <- list(fullNames, fullNames)
  }

  covFull <- .sirCapCovCorrelation(covFull, capCorrelation = capCorrelation)
  repaired <- .sirEnsurePosDef(covFull, report = TRUE)
  covFull <- repaired$covMat
  posDefAdjusted <- repaired$adjusted

  list(
    mu = muFull,
    covMat = covFull,
    paramSpace = ps,
    paramNames = fullNames,
    omegaRoute = route,
    fallbackNames = missingNames,
    posDefAdjusted = posDefAdjusted,
    # The full repair record, not just the flag. This is the covariance every
    # iteration-1 candidate is drawn from, so whether and how much it was
    # altered belongs in the run's provenance alongside the per-iteration
    # repairs of later empirical updates.
    initialRepair = repaired[c("adjusted", "method", "threshold", "magnitude")]
  )
}

# Bounds for the SIR parameter vector. OMEGA elements are not bounded here --
# they are constrained by the positive-definiteness check instead, which keeps
# the rejection diagnostics separable.
.sirParamBounds <- function(ps) {
  list(
    lower = stats::setNames(ps$lower, ps$sirName),
    upper = stats::setNames(ps$upper, ps$sirName),
    boundedNames = ps$sirName[ps$kind %in% c("theta", "sigma")]
  )
}

.sirOmegaPd <- function(ps, vals, base_omega) {
  omega_mat <- .sirReconstructOmega(ps, vals, base_omega)
  tryCatch(
    {
      chol(omega_mat)
      TRUE
    },
    error = function(e) FALSE
  )
}

.sirSampleFullProposal <- function(
  mu,
  covMat,
  n,
  lower,
  upper,
  ps,
  baseOmega,
  thetaNames,
  sigmaNames = character(0),
  boxcoxState = NULL,
  maxAttemptFactor = 10L
) {
  p <- length(mu)
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, min.len = 1L)
  checkmate::assertMatrix(covMat, mode = "numeric", nrows = p, ncols = p)
  checkmate::assertCount(n, positive = TRUE)

  n <- as.integer(n)
  max_attempts <- maxAttemptFactor * n
  param_names <- names(mu)
  sample_mat <- matrix(
    NA_real_,
    nrow = n,
    ncol = p,
    dimnames = list(NULL, param_names)
  )
  original_mat <- sample_mat

  n_filled <- 0L
  n_attempted <- 0L
  inverse_rejected <- 0L
  theta_rejected <- 0L
  omega_rejected <- 0L
  sigma_rejected <- 0L

  while (n_filled < n && n_attempted < max_attempts) {
    batch_n <- min(n - n_filled, max_attempts - n_attempted)
    draws <- mvtnorm::rmvnorm(batch_n, mean = mu, sigma = covMat)
    colnames(draws) <- param_names
    n_attempted <- n_attempted + batch_n

    original <- .sirBcInverseMatrix(draws, boxcoxState)
    # !is.finite(), not is.na(): the inverse transform can overflow to Inf
    # without raising, and is.na(Inf) is FALSE. An infinite value then
    # passed this filter and was tested against the bounds -- where
    # Inf > Inf is FALSE, so an unbounded parameter sailed through and was
    # charged to whichever check happened to reject it next, or to none.
    # A non-finite draw is an inverse failure and is classified as one here,
    # before any parameter-specific check.
    has_bad <- apply(!is.finite(original), 1L, any)
    theta_out <- rep(FALSE, nrow(original))
    if (length(thetaNames) > 0L) {
      theta_mat <- original[, thetaNames, drop = FALSE]
      theta_out <- apply(
        sweep(theta_mat, 2L, lower[thetaNames], "<") |
          sweep(theta_mat, 2L, upper[thetaNames], ">"),
        1L,
        any
      )
    }
    sigma_out <- rep(FALSE, nrow(original))
    if (length(sigmaNames) > 0L) {
      sigma_mat <- original[, sigmaNames, drop = FALSE]
      sigma_out <- apply(
        sweep(sigma_mat, 2L, lower[sigmaNames], "<") |
          sweep(sigma_mat, 2L, upper[sigmaNames], ">"),
        1L,
        any
      )
    }
    in_bounds <- !has_bad & !theta_out & !sigma_out
    omega_pd <- vapply(
      seq_len(nrow(original)),
      function(i) {
        if (has_bad[i] || !in_bounds[i]) {
          return(FALSE)
        }
        .sirOmegaPd(ps, original[i, ], baseOmega)
      },
      logical(1L)
    )

    ok <- !has_bad & in_bounds & omega_pd
    inverse_rejected <- inverse_rejected + sum(has_bad)
    theta_rejected <- theta_rejected + sum(!has_bad & theta_out)
    omega_rejected <- omega_rejected + sum(!has_bad & in_bounds & !omega_pd)
    sigma_rejected <- sigma_rejected + sum(!has_bad & sigma_out)

    n_take <- min(sum(ok), n - n_filled)
    if (n_take > 0L) {
      idx <- seq(n_filled + 1L, n_filled + n_take)
      keep <- which(ok)[seq_len(n_take)]
      sample_mat[idx, ] <- draws[keep, , drop = FALSE]
      original_mat[idx, ] <- original[keep, , drop = FALSE]
      n_filled <- n_filled + n_take
    }
  }

  if (n_filled < n) {
    cli::cli_warn(c(
      "Only {n_filled} of {n} requested full SIR samples were valid after {max_attempts} draw attempts.",
      "i" = "Consider widening the proposal or parameter bounds."
    ))
    sample_mat <- sample_mat[seq_len(n_filled), , drop = FALSE]
    original_mat <- original_mat[seq_len(n_filled), , drop = FALSE]
  }

  list(
    samples = original_mat,
    samplesForPdf = sample_mat,
    nAttempted = n_attempted,
    inverseRejected = inverse_rejected,
    thetaRejected = theta_rejected,
    omegaRejected = omega_rejected,
    sigmaRejected = sigma_rejected
  )
}

#' Run one SIR iteration
#'
#' Orchestrates sampling, OFV evaluation, weight computation, resampling, and
#' proposal update for a single iteration of the SIR algorithm.  In the first
#' iteration `boxcoxState` should be `NULL`; thereafter pass the `boxcoxState`
#' element returned by the previous call.
#'
#' @param fit nlmixr2 fit object.
#' @param mu Named numeric vector of proposal means in the **original**
#'   parameter scale.  Iteration 1 may pass theta means only; missing
#'   omega/sigma entries are initialized from the fit.
#' @param proposalCov Covariance matrix in the **sampling** scale.  Iteration 1
#'   may pass theta covariance only; later iterations pass the full empirical
#'   proposal covariance returned by the previous call.
#' @param nSamples Positive integer. Number of parameter vectors to draw.
#' @param nResample Positive integer. Target number of vectors to resample.
#' @param iterNum Positive integer. Iteration index (used for file naming).
#' @param requestedSamples Positive integer. User-requested sample count before
#'   PsN-style attempted-sample compensation. Defaults to `nSamples`.
#' @param capResampling Numeric >= 1. Max appearances per sample. Default `1`
#'   samples without replacement; values greater than `1` allow limited
#'   replacement up to the cap.
#' @param recenter Logical. If `TRUE` and any sampled dOFV < 0, shift `mu` to
#'   the best vector. Default `TRUE`.
#' @param boxcox Logical. Apply per-column Box-Cox before updating the
#'   proposal. Default `TRUE`.
#' @param directory Optional path used by higher-level callers that manage
#'   output files for the overall SIR run.
#' @param workers Passed to `.withWorkerPlan()` for parallel OFV evaluation.
#' @param rxThreads rxode2 OpenMP threads per worker; see [runSIR()].
#' @param boxcoxState `NULL` or a data frame (`param`, `lambda`, `delta`) from
#'   a previous call.  When non-`NULL` the proposal is assumed to be in
#'   Box-Cox space and samples are back-transformed before OFV evaluation.
#' @param thetaInflation,omegaInflation,sigmaInflation Non-negative variance
#'   multipliers for theta, fallback omega, and fallback sigma proposal blocks.
#' @param capCorrelation Maximum absolute proposal correlation after covariance
#'   updates.
#' @param omegaFallback How OMEGA uncertainty is obtained. `"cov"` (the
#'   default) takes it from `fit$cov`, including its correlations with
#'   THETA, and falls back automatically when `fit$cov` does not carry
#'   OMEGA. `"wishart"` always uses the Wishart-style approximation and
#'   gives a block-diagonal proposal.
#' @param sigmaFallbackRse Percent relative standard error for sigma-like
#'   fallback uncertainty when no SE is available.
#' @param omegaDf Optional degrees of freedom for the omega Wishart-style
#'   fallback; defaults to `nsub - 1`. Unused when `omegaFallback = "cov"`
#'   and `fit$cov` carries OMEGA.
#' @param isLastIteration Logical. If `TRUE`, update the proposal covariance on
#'   the original scale even when `boxcox = TRUE`.
#' @return Named list with elements `resampledMat`, `newMu`, `newCov`,
#'   `iterSummary`, `boxcoxState`, `rawResults`.
#' @noRd

# Per-parameter variance inflation factors, in SIR parameter order.
#
# Port of PsN setup_inflation() (lib/tool/sir.pm). Each of the three arguments
# is either a single value applied to the whole class, or one value per
# *diagonal* element of that class -- PsN's `inflate_only_diagonal = 1`
# semantics. An OMEGA off-diagonal is never given a factor directly: it gets
# `sqrt(infl_i) * sqrt(infl_j)` from the two diagonals it connects, which is
# what leaves the correlation unchanged when the factors are equal.
.sirInflationVector <- function(
  ps,
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1
) {
  out <- stats::setNames(rep(1, nrow(ps)), ps$sirName)

  .check <- function(given, label) {
    checkmate::assertNumeric(
      given,
      lower = 0,
      finite = TRUE,
      any.missing = FALSE,
      min.len = 1L,
      .var.name = label
    )
  }

  .fill <- function(kinds, given, label) {
    diag_i <- which(ps$kind %in% kinds[["diag"]])
    off_i <- which(ps$kind %in% kinds[["off"]])
    .check(given, label)

    if (length(diag_i) == 0L) {
      if (length(given) != 1L || !isTRUE(all.equal(given[[1L]], 1))) {
        cli::cli_abort(
          "{.arg {label}} was given, but the model has no estimated {kinds[['what']]}."
        )
      }
      return(invisible(NULL))
    }
    if (length(given) != 1L && length(given) != length(diag_i)) {
      cli::cli_abort(c(
        "{.arg {label}} must be length 1 or one value per {kinds[['what']]}.",
        "i" = "The model has {length(diag_i)} estimated {kinds[['what']]}, but {length(given)} value{?s} {?was/were} given."
      ))
    }

    if (length(given) == 1L) {
      out[c(diag_i, off_i)] <<- given
      return(invisible(NULL))
    }

    out[diag_i] <<- given
    if (length(off_i) > 0L) {
      # sqrt(infl_i) * sqrt(infl_j), keyed by eta index
      byEta <- stats::setNames(given, as.character(ps$neta1[diag_i]))
      out[off_i] <<- sqrt(byEta[as.character(ps$neta1[off_i])]) *
        sqrt(byEta[as.character(ps$neta2[off_i])])
    }
    invisible(NULL)
  }

  .fill(
    list(diag = "theta", off = character(0), what = "THETA"),
    thetaInflation,
    "thetaInflation"
  )
  .fill(
    list(
      diag = "sigma",
      off = character(0),
      what = "residual error parameters"
    ),
    sigmaInflation,
    "sigmaInflation"
  )
  .fill(
    list(diag = "omegaDiag", off = "omegaOffdiag", what = "OMEGA diagonals"),
    omegaInflation,
    "omegaInflation"
  )
  out
}
