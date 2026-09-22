# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# One SIR iteration, and the PsN-style attempted-sample adjustment.

sirRunIteration <- function(
  fit,
  mu,
  proposalCov,
  nSamples,
  nResample,
  iterNum,
  capResampling = 1,
  recenter = TRUE,
  boxcox = TRUE,
  directory = NULL,
  workers = NULL,
  rxThreads = NULL,
  boxcoxState = NULL,
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8,
  omegaFallback = c("cov", "wishart"),
  sigmaFallbackRse = 30,
  omegaDf = NULL,
  requestedSamples = nSamples,
  isLastIteration = FALSE,
  referenceOfv = NULL,
  rankDeficiency = c("abort", "repair")
) {
  rankDeficiency <- match.arg(rankDeficiency)
  omegaFallback <- match.arg(omegaFallback)

  # The OFV that dOFV is measured against. It starts at the fitted optimum and
  # moves with the centre whenever recentring finds a better one, so later
  # iterations are not scored against a superseded optimum (PsN does the same
  # in lib/tool/sir.pm). Passed explicitly so recovery and added iterations
  # resume from the reference the run actually reached.
  if (is.null(referenceOfv)) {
    referenceOfv <- fit$objf
  }
  checkmate::assertNumber(referenceOfv, finite = TRUE)

  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, min.len = 1L)
  checkmate::assertMatrix(proposalCov, mode = "numeric")
  checkmate::assertCount(nSamples, positive = TRUE)
  checkmate::assertCount(nResample, positive = TRUE)
  checkmate::assertCount(iterNum, positive = TRUE)
  checkmate::assertCount(requestedSamples, positive = TRUE)
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)
  checkmate::assertFlag(recenter)
  checkmate::assertFlag(boxcox)
  checkmate::assertFlag(isLastIteration)
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)
  checkmate::assertNumber(sigmaFallbackRse, lower = 0, finite = TRUE)
  if (!is.null(directory)) {
    checkmate::assertString(directory)
  }

  proposal <- .sirInitialProposal(
    fit = fit,
    mu = mu,
    proposalCov = proposalCov,
    thetaInflation = thetaInflation,
    omegaInflation = omegaInflation,
    sigmaInflation = sigmaInflation,
    capCorrelation = capCorrelation,
    omegaFallback = omegaFallback,
    sigmaFallbackRse = sigmaFallbackRse,
    omegaDf = omegaDf
  )

  ps <- proposal$paramSpace
  param_names <- proposal$paramNames
  theta_names <- ps$sirName[ps$kind == "theta"]
  sigma_names <- ps$sirName[ps$kind == "sigma"]
  bounds <- .sirParamBounds(ps)

  # ---- 1-3. Sample full proposal vectors and reject invalid draws ----
  bc_mu <- .sirBcTransformMu(proposal$mu, boxcoxState)
  sampled <- .sirSampleFullProposal(
    mu = bc_mu,
    covMat = proposal$covMat,
    n = nSamples,
    lower = bounds$lower,
    upper = bounds$upper,
    ps = ps,
    baseOmega = fit$omega,
    thetaNames = theta_names,
    sigmaNames = sigma_names,
    boxcoxState = boxcoxState
  )

  param_mat <- sampled$samples
  samples_for_weights <- sampled$samplesForPdf
  n_collected <- nrow(param_mat)
  if (n_collected == 0L) {
    cli::cli_abort(c(
      "No valid SIR samples were collected.",
      "i" = "Check proposal covariance, bounds, and omega positive-definiteness diagnostics."
    ))
  }

  # ---- 4. Evaluate OFV ----
  ofv_vals <- sirEvalOFV(
    fit,
    param_mat,
    workers = workers,
    rxThreads = rxThreads
  )
  dofv <- ofv_vals - referenceOfv

  # ---- 5. Handle failures ----
  n_failed <- sum(is.na(dofv))
  n_success <- n_collected - n_failed
  if (n_success == 0L) {
    eval_errors <- attr(ofv_vals, "evalErrors")
    cli::cli_abort(c(
      "All SIR OFV evaluations failed.",
      "i" = "No valid importance weights can be computed.",
      if (length(eval_errors) > 0L) {
        c("x" = "First error: {eval_errors[[1L]]}")
      }
    ))
  }
  n_resample_adj <- .sirAdjustedResamples(
    requestedResamples = nResample,
    requestedSamples = requestedSamples,
    successfulCount = n_success
  )
  if (n_resample_adj != nResample) {
    cli::cli_warn(c(
      "{n_success}/{requestedSamples} requested SIR samples had usable OFV evaluations.",
      "i" = "Adjusting nResample to {n_resample_adj}."
    ))
  }
  if (n_resample_adj < 1L) {
    cli::cli_abort(c(
      "Turnout scaled the resample count below one.",
      "i" = "{n_success}/{requestedSamples} samples usable, {nResample} requested resamples.",
      "i" = "Increase {.arg nSamples} or {.arg nResample}."
    ))
  }

  # ---- 6. Compute weights in the same full proposal scale used for sampling ----
  #
  # When Box-Cox is active the draws above are transformed coordinates, so the
  # normal density of `samples_for_weights` is not the density induced on the
  # original parameter scale. The Jacobian converts it, and is passed relative
  # to the centre so that relPDF stays 1 there.
  log_jacobian <- if (is.null(boxcoxState)) {
    NULL
  } else {
    mu_row <- matrix(
      proposal$mu[param_names],
      nrow = 1L,
      dimnames = list(NULL, param_names)
    )
    .sirBcLogJacobian(param_mat, boxcoxState) -
      .sirBcLogJacobian(mu_row, boxcoxState)
  }
  weights <- sirCalcWeights(
    samples_for_weights,
    bc_mu,
    proposal$covMat,
    dOFV = dofv,
    logJacobian = log_jacobian
  )

  # ---- 7. Resample ----
  #
  # n_resample_adj was scaled for turnout using the count of samples with a
  # usable OFV. The binding constraint is narrower: only samples with non-zero
  # resampling probability can actually be drawn, and a finite dOFV does not
  # guarantee a finite importance ratio. Clamp here rather than letting
  # sirResample() abort, so the run degrades the same graceful way it does for
  # turnout -- and so the message names the real cause.
  n_usable <- sum(is.finite(weights$prob_resample) & weights$prob_resample > 0)
  max_draws <- n_usable * as.integer(floor(capResampling))
  if (n_resample_adj > max_draws) {
    cli::cli_warn(c(
      "Only {n_usable} of {n_collected} SIR samples have non-zero resampling probability.",
      "i" = "Reducing resamples from {n_resample_adj} to {max_draws}.",
      "i" = "A wider proposal, more samples, or a higher {.arg capResampling} would avoid this."
    ))
    n_resample_adj <- max_draws
  }
  if (n_resample_adj < 1L) {
    cli::cli_abort(c(
      "No SIR samples have non-zero resampling probability.",
      "i" = "Check the OFV failures and the proposal covariance."
    ))
  }

  # A schedule can be feasible as requested and infeasible as realized:
  # evaluations fail, and the turnout and cap clamps above cut the retained
  # count further. The empirical covariance needs more than p vectors, so stop
  # here rather than inside .sirCheckProposalRank(), which sees only the final
  # count and advises raising nResample -- the wrong remedy when the real cause
  # is that most candidates never evaluated.
  n_param_cols <- ncol(param_mat)
  if (n_resample_adj <= n_param_cols) {
    cli::cli_abort(c(
      "Too few usable candidates to update the proposal.",
      "x" = "{n_usable} of {n_collected} candidate{?s} could be scored and weighted, for {n_param_cols} parameter{?s}.",
      "i" = "{n_resample_adj} vector{?s} would be retained; the covariance would have rank at most {max(n_resample_adj - 1L, 0L)}.",
      "i" = "Increase {.arg nSamples}, or check why candidates are failing to evaluate."
    ))
  }

  resampled <- sirResample(
    param_mat,
    weights,
    m = n_resample_adj,
    capResampling = capResampling
  )

  # ---- 7b. Weight degeneracy ----
  weight_diag <- .sirWeightDiagnostics(
    weights$prob_resample,
    nSuccessful = n_success
  )
  # n_success is the denominator of the efficiency ratio, and the warning uses
  # it to say when that ratio is measured on too few samples to trust.
  .sirWarnWeightDegeneracy(weight_diag, iterNum, nSamples = n_success)

  # ---- 8. Recenter ----
  new_mu <- proposal$mu
  new_reference_ofv <- referenceOfv
  if (recenter) {
    valid_dofv <- ifelse(is.na(dofv), Inf, dofv)
    best_idx <- which.min(valid_dofv)
    if (isTRUE(valid_dofv[best_idx] < 0)) {
      new_mu <- param_mat[best_idx, ]
      # Move the reference with the centre. Within this iteration a constant
      # dOFV shift cancels in the normalised weights, so this changes nothing
      # here -- it is later iterations and the chi-square diagnostic that would
      # otherwise stay pinned to the old optimum.
      new_reference_ofv <- ofv_vals[[best_idx]]
      cli::cli_inform(
        "  Iter {iterNum}: recentered mu (dOFV = {round(dofv[best_idx], 4)})."
      )
    }
  } else if (any(dofv < 0, na.rm = TRUE)) {
    cli::cli_warn(
      "At least one SIR sample had negative dOFV, but {.arg recenter} is FALSE."
    )
  }

  # ---- 9. Update full proposal in original scale ----
  updated <- sirUpdateProposal(
    resampled$samples[, param_names, drop = FALSE],
    boxcox = boxcox && !isLastIteration,
    capCorrelation = capCorrelation,
    # The centre the next iteration will transform with these very parameters.
    # Under recenter it is the best candidate, not necessarily a retained row.
    centre = new_mu,
    rankDeficiency = rankDeficiency
  )
  new_cov <- updated$covMat
  new_bc_state <- updated$boxcoxParams # NULL when boxcox = FALSE

  # ---- 10. Raw results data frame ----
  raw_df <- .sirBuildRawResults(
    paramMat = param_mat,
    weights = weights,
    dOFV = dofv,
    resampled = resampled,
    mu = proposal$mu,
    capResampling = capResampling
  )

  # ---- Summary ----
  iter_summary <- data.frame(
    iter = iterNum,
    nSamples = requestedSamples,
    nAttempted = nSamples,
    nDrawAttempts = sampled$nAttempted,
    nCollected = n_collected,
    nSuccessful = n_success,
    nFailed = n_failed,
    nResample = nResample,
    nResampled = nrow(resampled$samples),
    ess = weight_diag$ess,
    essFraction = weight_diag$essFraction,
    maxWeight = weight_diag$maxWeight,
    perplexity = weight_diag$perplexity,
    nNonNegligible = weight_diag$nNonNegligible,
    # Two different repairs, deliberately separate columns. posDefAdjusted is
    # the empirical update this iteration produced for the NEXT one;
    # proposalRepaired is the covariance this iteration actually drew from.
    posDefAdjusted = isTRUE(updated$posDefAdjusted),
    proposalRepaired = isTRUE(proposal$initialRepair$adjusted),
    proposalRepairMagnitude = proposal$initialRepair$magnitude %||% NA_real_,
    minDOFV = if (all(is.na(dofv))) NA_real_ else min(dofv, na.rm = TRUE),
    meanDOFV = if (all(is.na(dofv))) NA_real_ else mean(dofv, na.rm = TRUE),
    nNegativeDOFV = sum(dofv < 0, na.rm = TRUE),
    thetaRejected = sampled$thetaRejected,
    omegaRejected = sampled$omegaRejected,
    sigmaRejected = sampled$sigmaRejected,
    inverseRejected = sampled$inverseRejected
  )

  list(
    resampledMat = resampled$samples,
    newMu = new_mu,
    newCov = new_cov,
    # The full record for the covariance this iteration sampled from. runSIR()
    # keeps iteration 1's as the run's initial-proposal provenance.
    proposalRepair = proposal$initialRepair,
    referenceOfv = referenceOfv,
    newReferenceOfv = new_reference_ofv,
    iterSummary = iter_summary,
    boxcoxState = new_bc_state,
    rawResults = raw_df
  )
}

# Port of PsN update_attempted_samples() (lib/tool/sir.pm). Compensates the
# next iteration's sample count for samples lost to failed OFV evaluation, so
# the requested count is what actually survives.
#
# Matches PsN exactly: triggers on loss only, at turnout <= 0.95 inclusive, and
# rounds half away from zero. The previous implementation used a strict `<`,
# `ceiling()`, and a max() clamp against the requested count, which gave 112
# where PsN's own oracle says 111.
.sirAdjustedAttemptedSamples <- function(
  requestedSamples,
  previousAttempted = NULL,
  previousSuccessful = NULL
) {
  checkmate::assertCount(requestedSamples, positive = TRUE)
  if (is.null(previousAttempted) || is.null(previousSuccessful)) {
    return(as.integer(requestedSamples))
  }
  checkmate::assertCount(previousAttempted, positive = TRUE)
  checkmate::assertCount(previousSuccessful, positive = TRUE)

  previousTurnout <- previousSuccessful / previousAttempted
  if (previousTurnout > 1) {
    cli::cli_abort(c(
      "More successful samples than attempted in the previous iteration.",
      "i" = "{previousSuccessful} successful of {previousAttempted} attempted."
    ))
  }
  if (previousTurnout <= 0.95) {
    return(.sirRound(requestedSamples / previousTurnout))
  }
  as.integer(requestedSamples)
}

# Port of PsN update_actual_resamples() (lib/tool/sir.pm). Scales the resample
# count by this iteration's turnout, on a gain *or* a loss of at least 5%.
#
# `turnout` is measured against the originally requested sample count, not the
# compensated attempted count -- PsN's oracle pins this: 109 successful of a
# requested 100 gives turnout 1.09 even though 109 were attempted.
.sirAdjustedResamples <- function(
  requestedResamples,
  requestedSamples,
  successfulCount
) {
  checkmate::assertCount(requestedResamples, positive = TRUE)
  checkmate::assertCount(requestedSamples, positive = TRUE)
  checkmate::assertCount(successfulCount)

  turnout <- successfulCount / requestedSamples
  if (abs(turnout - 1) >= 0.05) {
    return(.sirRound(requestedResamples * turnout))
  }
  as.integer(requestedResamples)
}
