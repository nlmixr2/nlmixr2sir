# Part of nlmixr2sir.
# Preflight checks that the objective SIR scores candidates with is the same
# objective that produced the input fit.

# Estimation methods whose objective sirEvalOFV() is known to reproduce.
#
# sirEvalOFV() scores every candidate by re-evaluating the fit's OWN method
# with maxOuterIterations = 0, carrying the fit's own control object forward
# with only the evaluation fields overridden. That reproduces any method on the
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
#   npag/npb        the mixing distribution is not a normal Omega, so the
#                   whole proposal construction does not apply.
#   emvi/fbvi/vae   variational bounds, not the marginal likelihood.
#   nlm-family      no random effects, so there is no integral and no Omega.
#
# Adding a method means validating an evaluator that reproduces ITS objective,
# not adding a string here. The preflight enforces that at run time.
# The importance-sampling family is admitted on a DIFFERENT basis from the
# ladder above, and the difference matters if this list is ever revisited.
#
# The ladder is admitted because sirEvalOFV() reproduces each method's own
# objective. imp/impmap/qrpem are admitted because their objective is not their
# own: nlmixr2est recomputes it as a nested FOCEi evaluation at the converged
# estimates for every such fit (.impmapRecomputeObjf(), nlmixr2est
# R/impmap.R:1106), so fit$objf is already a FOCEi number. SIR scores their
# candidates as FOCEi to match -- see .sirImpEvalControl() in R/sir-eval.R.
#
# That is not a shortcut. Measured on theo_sd at nlmixr2est 7.1.0, one eta and
# three etas, all agreeing bit-for-bit with fit$objf:
#
#   one eta    193.6046289649, against $impObj 193.9889503885
#   three etas 116.8319956005, against $impObj 117.8403322041
#
# The one-eta case is the demanding one: it is where the eta-Hessian defect
# that forced the upstream recompute bites, and the raw C++ objective is ~19.96
# units low there. Reproducing 193.60 rather than 173.63 is the evidence that
# SIR is on the recomputed surface.
#
# No minimum nlmixr2est version is asserted for this. SIR never takes the
# nIter = 0 evaluation path that 7.1.0 added -- it scores as FOCEi -- so what it
# relies on is .impmapRecomputeObjf() running unconditionally, which is older
# than 7.1.0 and is verified empirically by the per-run preflight anyway. A
# version that published the raw importance-sampling objective instead would
# miss by ~19.96 units on a one-eta model and be refused.
.sirSupportedEstimationMethods <- local({
  base <- c("focei", "foce", "focep", "laplace", "agq")
  sort(c(
    "fo", "foi", base, paste0("m", base), paste0("i", base),
    "imp", "impmap", "qrpem"
  ))
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
      "x" = "Candidates are scored by re-evaluating the fit's own method at fixed parameters, and {.val {est}} cannot be evaluated that way.",
      "i" = "Supported: {.val {supported}}.",
      "i" = "Refit with a deterministic method such as {.code est = \"focei\"} to run SIR on this model."
    ))
  }
  # Announced, not silent. Every other supported method scores candidates with
  # itself; this family does not, and a user who chose impmap deliberately is
  # owed the fact that the numbers come from somewhere else. Once per run --
  # this is the preflight, which runs once, not the evaluator, which runs per
  # candidate.
  if (est %in% .sirImpFamilyMethods) {
    cli::cli_warn(c(
      "SIR will score {.val {est}} candidates with FOCEi.",
      "i" = "{.code fit$objf} on an {.val {est}} fit is already a FOCEi re-evaluation at the converged estimates, not the importance-sampling objective. nlmixr2est recomputes it that way for every fit of this family.",
      "i" = "The importance-sampling objective is {.code fit$env$impObj}, and SIR does not use it.",
      "i" = "Scoring as FOCEi reproduces {.code fit$objf} exactly and skips an E-step whose result would be discarded. This is expected, not a problem with your fit."
    ))
  }
  invisible(est)
}

# How the objective preflight decides what is worth refusing.
#
# The check compares fit$objf with a fresh re-evaluation at the same estimates.
# That difference is a CONSTANT across candidates -- every candidate is scored
# by the same evaluator, so it shifts every dOFV equally, and a constant shift
# in dOFV multiplies every weight by the same factor and divides out of the
# normalised weights. R/sir-iterate.R already says exactly this where recentring
# moves the reference.
#
# Since runSIR() now measures dOFV against the RE-EVALUATED centre rather than
# fit$objf, that constant is zero by construction and the difference no longer
# reaches the weights at all. What remains of this check is its real job:
# noticing that the evaluator is on a different SURFACE from the fit.
#
# Measured on Rik Schoemaker's QR models (P9-PROGRESS.md), the difference was
# 10x to 116x larger than the part that does not cancel, so aborting on it
# refused sound runs over the term that provably vanishes.
#
# The threshold is RELATIVE because that is the scale on which the two things
# this check must tell apart actually separate:
#
#   convergence slack, 20 QR models   <= 1.4e-5 of the objective
#   a SAEM fit scored under FOCEi      2.3e-2 of the objective (2.69 on theo_sd)
#
# 1600x apart relatively; only 10x apart absolutely (0.27 against 2.69), which
# is why the absolute default that preceded this could not separate them.
# 1e-3 sits ~70x above the worst observed slack and ~23x below the SAEM
# signature.
.sirObjfAbsFloor <- 1e-2

# Never tighter than the floor: a relative tolerance alone would be absurd on a
# small objective, and 1e-2 is comfortably achievable at theo_sd scale.
.sirObjfAbortThreshold <- function(stored, objfTolerance) {
  # An explicit zero is honoured rather than floored. It is the only way to ask
  # for exact agreement, the suite uses it to exercise the mismatch path, and
  # silently overriding a user's zero with 1e-2 would be the kind of quiet
  # substitution this package refuses elsewhere.
  if (objfTolerance <= 0) {
    return(0)
  }
  max(objfTolerance * abs(stored), .sirObjfAbsFloor)
}

# Purely informational, and absolute on purpose: it exists to tell the user
# that reported dOFVs are measured against a centre that differs from the
# published fit$objf by this much. Below ~0.1 OFV units nobody needs telling.
.sirObjfWarnTolerance <- 0.1

# A fixed pseudo-random unit direction that does NOT consume the session's RNG
# stream. Drawing from the global stream here would shift every subsequent
# sample the run draws, so the whole run's answer would depend on whether the
# noise floor was measured.
# Run `expr` with the session's RNG state restored afterwards.
#
# Wrapped around the WHOLE noise measurement, not merely the direction draw:
# the evaluator goes through the worker plan, which consumes from the stream
# itself. Without this, whether the noise floor was measured would change every
# sample the run subsequently draws.
.sirWithPreservedRng <- function(expr) {
  hasSeed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (hasSeed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(
    {
      if (is.null(old)) {
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        }
      } else {
        assign(".Random.seed", old, envir = globalenv())
      }
    },
    add = TRUE
  )
  force(expr)
}

.sirFixedDirection <- function(n, seed = 20260921L) {
  set.seed(seed)
  v <- stats::rnorm(n)
  v / sqrt(sum(v^2))
}

# Measure the evaluator's noise floor: how much the scored objective wanders
# for reasons that are NOT the objective changing.
#
# This is the term the old centre check never measured --
# e(candidate) - e(centre) -- and it is the only one that reaches the weights.
#
# Method: walk a transect through parameter space along a fixed direction,
# scaled by each parameter's own proposal SD, so the separations are the ones
# candidates actually have. The true objective is smooth along a straight line,
# so whatever a low-order polynomial in the step cannot absorb is evaluator
# slack. Order 4 over 15 points leaves 10 residual degrees of freedom and was
# measured to sit within ~20% of the order-6 plateau.
#
# Measured floors ranged 1e-4 to 2.8e-2 across the QR models -- 280x, and
# model-specific, which is exactly why it has to be measured per run rather
# than assumed.
.sirObjectiveNoise <- function(fit,
                               ps = .sirParamSpace(fit),
                               mu = NULL,
                               nPoints = 15L,
                               span = 0.25,
                               workers = NULL,
                               rxThreads = NULL) {
  .sirWithPreservedRng(.sirObjectiveNoiseImpl(
    fit, ps, mu, nPoints, span, workers, rxThreads
  ))
}

.sirObjectiveNoiseImpl <- function(fit, ps, mu, nPoints, span, workers,
                                   rxThreads) {
  none <- list(noise = NA_real_, nOk = 0L, range = NA_real_)
  if (is.null(mu)) {
    mu <- tryCatch(.sirProposalMu(fit, ps), error = function(e) NULL)
  }
  if (is.null(mu)) {
    return(none)
  }
  cov <- tryCatch(
    suppressMessages(sirGetProposalCov(fit)),
    error = function(e) NULL
  )
  if (!is.matrix(cov) || is.null(rownames(cov))) {
    return(none)
  }
  nm <- intersect(names(mu), rownames(cov))
  if (length(nm) == 0L) {
    return(none)
  }
  sdv <- sqrt(diag(cov)[nm])
  if (!all(is.finite(sdv)) || all(sdv == 0)) {
    return(none)
  }
  dir <- .sirFixedDirection(length(nm))
  steps <- seq(-span, span, length.out = nPoints)
  mat <- t(vapply(
    steps,
    function(s) {
      v <- mu
      v[nm] <- mu[nm] + s * dir * sdv
      v
    },
    numeric(length(mu))
  ))
  colnames(mat) <- names(mu)

  ofv <- tryCatch(
    sirEvalOFV(fit, mat, workers = workers, rxThreads = rxThreads),
    error = function(e) rep(NA_real_, nPoints)
  )
  ok <- is.finite(ofv)
  # Order 4 needs 5 coefficients; insist on a few residual degrees of freedom
  # rather than reporting a number fitted to nothing.
  if (sum(ok) < 9L) {
    return(none)
  }
  # Order 6, not 4. The residual is only "noise" once the polynomial has
  # absorbed the shape, and a quartic demonstrably cannot: on QR model N021 a
  # +/-1 SD transect swings 1208 OFV units, where the quartic residual was
  # 73.94 and still falling steeply with order -- pure unfitted curvature, and
  # it produced a false refusal before this was corrected.
  fitPoly <- tryCatch(
    stats::lm(ofv[ok] ~ stats::poly(steps[ok], 6L)),
    error = function(e) NULL
  )
  if (is.null(fitPoly)) {
    return(none)
  }
  noise <- stats::sd(stats::residuals(fitPoly))
  span_range <- diff(range(ofv[ok]))

  # Refuse to report a number the fit does not support. Genuine slack is orders
  # of magnitude below the transect's own variation (9e-5 of it on N001,
  # 3.5e-4 on N021 at this span); anything approaching a few percent means the
  # polynomial is still fitting the objective, not the slack in it.
  if (is.finite(span_range) && span_range > 0 && noise > 0.02 * span_range) {
    return(none)
  }
  list(noise = noise, nOk = sum(ok), range = span_range)
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
  objfTolerance = 1e-3,
  noise = TRUE,
  noiseTolerance = 1,
  warnTolerance = .sirObjfWarnTolerance,
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
    stencil = stencil,
    noise = NA_real_,
    abortThreshold = NA_real_
  )

  # Absolute tolerance only. A relative tolerance on the raw OFV is the wrong
  # scale: the objective carries additive constants and grows with the number
  # of observations, while the weights depend on DIFFERENCES in OFV. On an
  # objective of 50,000 a relative 1e-4 would wave through five OFV units.
  abortThreshold <- .sirObjfAbortThreshold(stored, objfTolerance)
  # Local alias: cli treats a `{}` expression starting with a dot as a style
  # name, so `{.sirObjfAbsFloor}` errors rather than interpolating.
  absFloor <- .sirObjfAbsFloor
  out$abortThreshold <- abortThreshold
  if (abs_diff > abortThreshold) {
    cli::cli_abort(c(
      "SIR cannot reproduce the fit's objective at its own estimates.",
      "x" = "Stored {.code fit$objf}: {format(stored, digits = 10)}",
      "x" = "Reevaluated at the same estimates: {format(reevaluated, digits = 10)}",
      "i" = "Absolute difference {format(abs_diff, digits = 4)}; threshold {format(abortThreshold, digits = 4)} ({objfTolerance} of the objective, floored at {absFloor}).",
      "i" = "A difference this large relative to the objective is the signature of a different likelihood surface, not of convergence slack, and candidates would not be scored on the surface the dOFVs are measured against.",
      "i" = "If the fit and the evaluator are on the same surface, refit with a higher {.code sigdig}: each additional digit has been measured to shrink the gap about 3-4 fold."
    ))
  }

  # Informational only. Since runSIR() measures dOFV against the re-evaluated
  # centre, this difference no longer reaches the weights -- but it does mean
  # reported dOFVs are not measured against the published fit$objf, and at this
  # magnitude that is worth saying once.
  if (abs_diff > warnTolerance) {
    cli::cli_warn(c(
      "The evaluator's objective at the fitted estimates differs from {.code fit$objf} by {format(abs_diff, digits = 4)}.",
      "i" = "Stored {.code fit$objf}: {format(stored, digits = 10)}; reevaluated: {format(reevaluated, digits = 10)}.",
      "i" = "dOFVs are measured against the reevaluated centre, so this constant does not reach the importance weights -- but reported dOFVs will not line up with {.code fit$objf}.",
      "i" = "This is inner-solve convergence slack. A higher {.code sigdig} shrinks it about 3-4 fold per digit."
    ))
  }

  # The part that does NOT cancel, and the only part that reaches the weights.
  if (isTRUE(noise)) {
    nz <- .sirObjectiveNoise(fit, ps, mu, workers = workers, rxThreads = rxThreads)
    out$noise <- nz$noise
    # Reported, never refused. An earlier draft aborted here and immediately
    # produced a false refusal on a sound model, which is precisely what this
    # preflight must not do: a diagnostic whose own estimator can be wrong has
    # no business stopping a run. If the number is untrustworthy
    # .sirObjectiveNoise() returns NA and nothing is said at all.
    if (is.finite(nz$noise) && nz$noise > noiseTolerance) {
      cli::cli_warn(c(
        "The objective carries a noise floor of {format(nz$noise, digits = 4)} OFV units.",
        "i" = "This is the part of the evaluation error that does NOT cancel between a candidate and the centre, so it enters each weight as {.code exp(-noise/2)}: about {format(100 * (1 - exp(-nz$noise / 2)), digits = 2)}%.",
        "i" = "The run continues. A higher {.code sigdig} shrinks it about 3-4 fold per digit."
      ))
    }
  }

  # Inside tolerance but worth saying out loud. A gap of this size does not
  # threaten the result -- it moves a weight by well under a percent -- but it
  # is the signal that the fit is closer to the threshold than most, and the
  # remedy is cheap.
  if (abs_diff > .sirObjfWarnTolerance) {
    cli::cli_warn(c(
      "The fit reproduces its own objective to {format(abs_diff, digits = 4)}, which is larger than usual.",
      "i" = "Stored {.code fit$objf}: {format(stored, digits = 10)}; reevaluated: {format(reevaluated, digits = 10)}.",
      "i" = "Within {.code objfTolerance} ({objfTolerance}), so the run continues, and a gap this size moves an importance weight by well under one percent.",
      "i" = "This is inner-solve convergence slack. Refitting with a higher {.code sigdig} shrinks it about 3-4 fold per digit if you want it smaller."
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
    # Absolute, and deliberately NOT objfTolerance: that is now a fraction of
    # the objective, so using it here would fire this warning whenever a probe
    # improved by 0.001 -- well inside the measured noise floor on some models.
    if (stencil$minDOFV < -warnTolerance) {
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
