# Part of nlmixr2sir.
# Preflight checks that the objective SIR scores candidates with is the same
# objective that produced the input fit.

# Estimation methods whose objective sirEvalOFV() is known to reproduce.
#
# sirEvalOFV() scores every candidate by building a fresh
# nlmixr2(est = "focei", maxOuterIterations = 0) call, carrying the fit's own
# likelihood-relevant settings forward. That reproduces any method on the
# deterministic FOCEi ladder, so such a fit is on one fixed target and the
# importance weights mean what they claim.
#
# It does not reproduce a different likelihood approximation. A SAEM fit stores
# an objective computed by Gaussian quadrature: on theo_sd that is 208.512,
# against 205.820 from the FOCEi reevaluation -- a 2.69 unit gap, which is a
# different surface rather than numerical noise. Scoring candidates on it would
# make every dOFV, the chi-square diagnostic, and recentring meaningless.
#
# The deterministic conditional-estimation ladder. Every one of these runs on
# the same FOCEi engine and differs only in settings sirEvalOFV() now carries
# (fo, nAGQ, foce, interaction, muModel), so one evaluator reproduces all of
# them. Verified on theo_sd: each method's objective is recovered at its own
# fitted estimates to within the preflight tolerance.
#
# Deliberately absent:
#
#   saem            stochastic MCMC E-step, and its stored objective is a
#                   different approximation entirely -- 208.512 by Gaussian
#                   quadrature against 205.820 from FOCEi on theo_sd.
#   imp/impmap/qrpem    No way to evaluate at fixed parameters. Measured on
#                   nlmixr2est 7.0.3, theo_sd, est = "impmap":
#
#                     * impmapControl() has no EONLY analogue and no
#                       maxOuterIterations -- nothing that suppresses the
#                       M-step the way PsN's EONLY=1 does;
#                     * nIter = 0, the obvious candidate, SEGFAULTS (exit 139,
#                       reproducible with plain nlmixr2(), nlmixr2sir not
#                       loaded);
#                     * nIter >= 1 is the wrong operation regardless. nIter
#                       counts EM iterations and every one runs an M-step that
#                       UPDATES the population parameters, so it does not score
#                       the candidate -- it takes an estimation step away from
#                       it. objf goes 117.440 (nIter=1), 116.860 (nIter=2),
#                       converging back to the fit's own 116.829 as it
#                       re-estimates.
#
#                   Separately, an imp-family fit carries two objectives that
#                   disagree -- fit$objf 116.829 against fit$env$impObj 117.836
#                   -- so even with a working evaluation mode, which one the
#                   dOFVs are measured against would still need deciding.
#
#                   Supporting these needs an expectation-only mode upstream,
#                   not a workaround here.
#   npag/npb        the mixing distribution is not a normal Omega, so the
#                   whole proposal construction does not apply.
#   emvi/fbvi/vae   variational bounds, not the marginal likelihood.
#   nlm-family      no random effects, so there is no integral and no Omega.
#
# Adding a method means validating an evaluator that reproduces ITS objective,
# not adding a string here. The preflight enforces that at run time.
.sirSupportedEstimationMethods <- local({
  base <- c("focei", "foce", "focep", "laplace", "agq")
  sort(c("fo", "foi", base, paste0("m", base), paste0("i", base)))
})

.sirSupportedEstimation <- function(est) {
  supported <- .sirSupportedEstimationMethods
  if (length(est) != 1L || is.na(est)) {
    cli::cli_abort(c(
      "Cannot determine the estimation method of {.arg fit}.",
      "i" = "SIR supports {.val {supported}}."
    ))
  }
  if (!est %in% supported) {
    cli::cli_abort(c(
      "SIR does not support {.val {est}} fits.",
      "x" = "Candidate objectives are evaluated with {.code est = \"focei\"}, which does not reproduce the {.val {est}} objective.",
      "i" = "Supported: {.val {supported}}.",
      "i" = "Refit with {.code est = \"focei\"} to run SIR on this model."
    ))
  }
  invisible(est)
}

#' Verify the candidate evaluator reproduces the fit's objective
#'
#' Re-evaluates the OFV at the fitted centre and compares it with the stored
#' `fit$objf`. Importance sampling assumes every target value is an evaluation
#' of one fixed target: this is the per-run evidence that the evaluator used
#' for candidates is on the same surface as the objective the dOFVs are
#' measured against.
#'
#' @param fit An nlmixr2 fit object.
#' @param workers,rxThreads Passed to `sirEvalOFV()`.
#' @param objfTolerance Non-negative scalar. The check passes when the absolute
#'   *or* relative difference is within this tolerance.
#' @return Invisibly, a list with `stored`, `reevaluated`, `absDiff`, `relDiff`.
#' @noRd
.sirCheckObjective <- function(
  fit,
  workers = NULL,
  rxThreads = NULL,
  objfTolerance = 1e-4,
  stencil = TRUE,
  stencilTolerance = 1
) {
  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertNumber(objfTolerance, lower = 0, finite = TRUE)

  .sirSupportedEstimation(.sirFitEst(fit))

  stored <- fit$objf
  if (!checkmate::testNumber(stored, finite = TRUE)) {
    cli::cli_abort(c(
      "{.code fit$objf} is not a finite objective function value.",
      "i" = "SIR measures every dOFV against it, so it must be available."
    ))
  }

  ps <- .sirParamSpace(fit)
  mu <- .sirProposalMu(fit, ps)
  centre <- matrix(mu, nrow = 1L, dimnames = list(NULL, names(mu)))
  reevaluated <- unname(sirEvalOFV(
    fit,
    centre,
    workers = workers,
    rxThreads = rxThreads
  )[[1L]])

  if (!is.finite(reevaluated)) {
    cli::cli_abort(c(
      "The objective could not be reevaluated at the fitted estimates.",
      "i" = "Every SIR candidate is scored the same way, so none would succeed."
    ))
  }

  abs_diff <- abs(reevaluated - stored)
  rel_diff <- abs_diff / max(abs(stored), .Machine$double.eps)

  # The stencil: perturb each parameter a little either way and check that the
  # fitted centre is still a local minimum of the evaluator's surface.
  #
  # Agreement at one point does not establish that two objectives are the same
  # function -- they can cross at the centre and diverge everywhere else. If the
  # evaluator is on a different surface, the fitted estimates are generally not
  # its optimum, so some small perturbation lowers the objective. That is
  # exactly what a SAEM fit scored under FOCEi would show.
  stencil <- if (isTRUE(stencil)) {
    .sirObjectiveStencil(
      fit, ps, mu, reevaluated,
      workers = workers, rxThreads = rxThreads
    )
  } else {
    NULL
  }

  out <- list(
    stored = stored,
    reevaluated = reevaluated,
    absDiff = abs_diff,
    relDiff = rel_diff,
    stencil = stencil
  )

  # Absolute tolerance only. A relative tolerance on the raw OFV is the wrong
  # scale: the objective carries additive constants and grows with the number
  # of observations, while the weights depend on DIFFERENCES in OFV. On an
  # objective of 50,000 a relative 1e-4 would wave through five OFV units.
  if (abs_diff > objfTolerance) {
    cli::cli_abort(c(
      "SIR cannot reproduce the fit's objective at its own estimates.",
      "x" = "Stored {.code fit$objf}: {format(stored, digits = 10)}",
      "x" = "Reevaluated at the same estimates: {format(reevaluated, digits = 10)}",
      "i" = "Absolute difference {format(abs_diff, digits = 4)}; tolerance {objfTolerance} (absolute).",
      "i" = "Candidates would be scored on a different surface from the dOFV reference, so the importance weights would not be meaningful.",
      "i" = "Raise {.code runSIRControl(objfTolerance =)} only if this difference is understood and acceptable."
    ))
  }

  # A probe that fails outright means the evaluator cannot score points near the
  # centre, so it will not score candidates either.
  if (!is.null(stencil) && stencil$nFailed > 0L) {
    cli::cli_abort(c(
      "The objective could not be evaluated near the fitted estimates.",
      "x" = "{stencil$nFailed} of {stencil$nProbes} probe{?s} failed.",
      "i" = "Candidates are scored the same way, so most would fail too."
    ))
  }

  if (!is.null(stencil) && is.finite(stencil$minDOFV)) {
    worst <- which.min(stencil$dOFV)
    if (stencil$minDOFV < -stencilTolerance) {
      cli::cli_abort(c(
        "The fitted estimates are not a local optimum of the objective SIR evaluates.",
        "x" = "Perturbing {.val {stencil$param[[worst]]}} by a fraction of its estimate lowered the objective by {format(-stencil$minDOFV, digits = 4)}.",
        "i" = "Agreement at the centre alone does not make two objectives the same function; this probes the surface around it.",
        "i" = "The usual cause is a fit whose likelihood settings cannot be reproduced by the evaluator.",
        "i" = "Tolerance is {.code runSIRControl(objfStencilTolerance =)}, currently {stencilTolerance}."
      ))
    }
    if (stencil$minDOFV < -objfTolerance) {
      cli::cli_warn(c(
        "The fitted estimates are not quite a local optimum of the objective.",
        "i" = "Best probe improved the objective by {format(-stencil$minDOFV, digits = 4)} ({.val {stencil$param[[worst]]}}).",
        "i" = "Usually means the fit stopped a little short of convergence; {.code recenter = TRUE} will move to a better point if one is sampled."
      ))
    }
  }

  invisible(out)
}

# Evaluate the objective either side of the fitted centre along every
# parameter, and report the largest improvement found.
#
# Steps are a small fraction of each estimate's own magnitude, so the stencil
# is unit-free, and are clipped into the parameter's bounds so no probe leaves
# the admissible region.
.sirObjectiveStencil <- function(fit, ps, mu, centreOfv, relStep = 1e-3,
                                 workers = NULL, rxThreads = NULL) {
  nm <- names(mu)
  bounds <- .sirParamBounds(ps)
  steps <- pmax(abs(as.numeric(mu)) * relStep, 1e-6)

  rows <- list()
  labels <- character()
  for (j in seq_along(mu)) {
    for (sgn in c(-1, 1)) {
      v <- mu
      cand <- v[[j]] + sgn * steps[[j]]
      lo <- bounds$lower[[nm[[j]]]]
      hi <- bounds$upper[[nm[[j]]]]
      if (is.finite(lo) && cand <= lo) next
      if (is.finite(hi) && cand >= hi) next
      v[[j]] <- cand
      rows[[length(rows) + 1L]] <- v
      labels <- c(labels, nm[[j]])
    }
  }
  if (length(rows) == 0L) {
    return(NULL)
  }

  mat <- do.call(rbind, rows)
  colnames(mat) <- nm
  ofv <- sirEvalOFV(fit, mat, workers = workers, rxThreads = rxThreads)
  dofv <- as.numeric(ofv) - centreOfv

  list(
    nProbes = nrow(mat),
    param = labels,
    dOFV = dofv,
    nFailed = sum(!is.finite(dofv)),
    minDOFV = if (all(!is.finite(dofv))) NA_real_ else min(dofv, na.rm = TRUE)
  )
}
