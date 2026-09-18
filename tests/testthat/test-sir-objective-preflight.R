# B2: the candidate objective must be the fit's objective ---------------------
#
# sirEvalOFV() scores every candidate with a freshly built FOCEi call at
# maxOuterIterations = 0. That is the workflow analogue of NONMEM MAXEVAL=0,
# but it is not by itself proof that the resulting objective is the same
# surface that produced fit$objf. Importance sampling assumes one fixed target:
# a candidate-dependent difference between evaluators changes the shape of the
# target and so the retained distribution.
#
# The preflight settles it empirically, per run, by re-evaluating the fitted
# centre and comparing with the stored objective before any sampling happens.

# A stand-in evaluator that scores every point `delta` off the stored objective,
# held ETAs or not: a genuine mismatch on every platform.
.sirOffsetEval <- function(delta) {
  function(fit, paramSamples, workers = NULL, rxThreads = NULL, fixEtas = NULL) {
    rep(fit$objf + delta, nrow(paramSamples))
  }
}

test_that("the preflight accepts a fit whose centre reproduces its objective", {
  skip_on_cran()
  fit <- theoFit()
  expect_no_error(.sirCheckObjective(fit, workers = 1L))
})

test_that("the preflight returns the stored and reevaluated objectives", {
  skip_on_cran()
  fit <- theoFit()
  res <- .sirCheckObjective(fit, workers = 1L)
  expect_named(
    res,
    c("stored", "reevaluated", "absDiff", "relDiff", "candidateCentre",
      "innerNoise", "stencil")
  )
  expect_equal(res$stored, fit$objf, tolerance = 1e-12)
  expect_lt(res$absDiff, 1e-3)
  expect_equal(res$innerNoise, res$candidateCentre - res$stored)
})

test_that("the identity check holds the fit's ETAs, so inner-optimization noise does not fail it", {
  skip_on_cran()
  # A FOCE fit whose ETAs re-optimize to a slightly different objective at the
  # same THETA: 1.5e-4 away on theo_sd, above the 1e-4 tolerance. At the fit's
  # own ETAs the objective is the stored one to rounding.
  fit <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    theoOneCmt, nlmixr2data::theo_sd, est = "foce",
    control = list(print = 0L, covMethod = "", calcTables = FALSE)
  )))
  r <- .sirCheckObjective(fit, workers = 1L, stencil = FALSE)
  expect_lt(r$absDiff, 1e-8)
  expect_true(is.finite(r$innerNoise))
})

test_that("a candidate-style centre far from fit$objf warns", {
  skip_on_cran()
  # The identity holds at the fit's own ETAs, but re-optimized ETAs land 5 OFV
  # units away: every candidate dOFV would carry that offset, so say so.
  fit <- theoFit()
  local_mocked_bindings(
    sirEvalOFV = function(fit, paramSamples, workers = NULL, rxThreads = NULL,
                          fixEtas = NULL) {
      if (is.null(fixEtas)) fit$objf + 5 else fit$objf
    }
  )
  expect_warning(
    r <- .sirCheckObjective(fit, workers = 1L, stencil = FALSE),
    "away from"
  )
  expect_equal(r$innerNoise, 5)
  expect_lt(r$absDiff, 1e-12)

  # Within the stencil tolerance it is noise, and passes quietly.
  local_mocked_bindings(
    sirEvalOFV = function(fit, paramSamples, workers = NULL, rxThreads = NULL,
                          fixEtas = NULL) {
      if (is.null(fixEtas)) fit$objf + 1e-3 else fit$objf
    }
  )
  expect_no_warning(.sirCheckObjective(fit, workers = 1L, stencil = FALSE))
})

test_that("the abort says so when the fit's ETAs could not be held", {
  skip_on_cran()
  fit <- theoFit()
  local_mocked_bindings(
    .sirFitEtaMat = function(fit) NULL,
    sirEvalOFV = .sirOffsetEval(1e-3)
  )
  err <- tryCatch(
    .sirCheckObjective(fit, workers = 1L, stencil = FALSE),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "could not be held", fixed = TRUE)
  expect_false(grepl("estimates and ETAs", err, fixed = TRUE))
})

test_that("holding the ETAs still refuses a different surface", {
  skip_on_cran()
  # FO scored as FOCEi at the FO fit's own ETAs is a different objective, and
  # the check must say so -- holding the ETAs must not make it pass trivially.
  fitFo <- theoFitFo()
  local_mocked_bindings(.sirEvalMethod = function(fit) "focei")
  expect_error(
    .sirCheckObjective(fitFo, workers = 1L, stencil = FALSE),
    "cannot reproduce"
  )
})

test_that("the preflight aborts when the centre does not reproduce the objective", {
  skip_on_cran()
  fit <- theoFit()
  # An evaluator 1e-3 off the stored objective. This is the mismatch path: the
  # message must name both values so the user can judge the gap. (This used to
  # rely on objfTolerance = 0 and a nonzero numerical difference, but at the
  # fit's own ETAs the objective reproduces exactly on some platforms --
  # Windows among them -- so there was no difference to find.)
  local_mocked_bindings(sirEvalOFV = .sirOffsetEval(1e-3))
  expect_error(
    .sirCheckObjective(fit, workers = 1L),
    "objective"
  )
  err <- tryCatch(
    .sirCheckObjective(fit, workers = 1L),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, format(fit$objf, digits = 10), fixed = TRUE)
})

test_that("unsupported estimation methods are rejected before sampling", {
  skip_on_cran()
  # The check is on the recorded method, so it fires without needing a real fit
  # of that kind. saem is the case that matters: its objective comes from
  # Gaussian quadrature, and reevaluating it under FOCEi shifts the OFV by
  # ~2.7 units on theo_sd -- a different surface, not numerical noise.
  expect_error(.sirSupportedEstimation("saem"), "saem")
  expect_error(.sirSupportedEstimation("nlme"), "nlme")
  expect_error(.sirSupportedEstimation(NA_character_), "determine")
  expect_silent(.sirSupportedEstimation("focei"))
})

test_that("runSIR runs the objective preflight before sampling", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()
  local_mocked_bindings(sirEvalOFV = .sirOffsetEval(1e-3))
  expect_error(
    .sirQuiet(runSIR(
      fit,
      nSamples = 16L,
      nResample = 8L,
      directory = dir,
      control = runSIRControl(recover = FALSE, workers = 1L)
    )),
    "objective"
  )
  # It aborted before doing any work, so there are no iteration artifacts.
  expect_false(file.exists(file.path(dir, "raw_results.csv")))
})

# B2: the stencil, and carrying the fit's likelihood settings ----------------

test_that("the evaluator carries the fit's likelihood-relevant controls", {
  skip_on_cran()
  fit <- theoFit()
  ec <- .sirEvalControl(fit)
  expect_s3_class(ec, "foceiControl")
  # Evaluation-only overrides are always ours.
  expect_equal(ec$maxOuterIterations, 0L)
  # foceiControl() normalises covMethod = "" to integer 0, meaning no step.
  expect_equal(as.integer(ec$covMethod), 0L)
  # Likelihood-defining settings come from the fit, not from the defaults.
  expect_identical(ec$interaction, fit$control$interaction)
  expect_identical(ec$addProp, fit$control$addProp)
})

test_that("the preflight probes the surface around the centre", {
  skip_on_cran()
  fit <- theoFit()
  res <- .sirCheckObjective(fit, workers = 1L)
  expect_false(is.null(res$stencil))
  # Two probes per parameter, less any clipped by a bound.
  expect_lte(res$stencil$nProbes, 2L * nrow(.sirParamSpace(fit)))
  expect_gt(res$stencil$nProbes, 0L)
  expect_equal(res$stencil$nFailed, 0L)
  # The fitted estimates are a local optimum, so no probe improves much on it.
  expect_gt(res$stencil$minDOFV, -1)
})

test_that("the stencil can be switched off", {
  skip_on_cran()
  res <- .sirCheckObjective(theoFit(), workers = 1L, stencil = FALSE)
  expect_null(res$stencil)
})

test_that("the preflight tolerance is absolute, not relative", {
  skip_on_cran()
  fit <- theoFit()
  # A relative rule would wave through a large absolute gap on a large
  # objective. The weights depend on differences in OFV, so only the absolute
  # scale is meaningful.
  local_mocked_bindings(sirEvalOFV = .sirOffsetEval(1e-3))
  expect_error(
    .sirCheckObjective(fit, workers = 1L, stencil = FALSE),
    "absolute"
  )
})

# Fit `m` with `est` and check the evaluator reproduces ITS objective, both at
# the fitted centre and at a candidate far from it.
#
# The off-centre half is the part that matters. Agreement at the centre is
# necessary but weak: an AGQ fit whose integration bounds were dropped agreed
# there to 1e-08 and was 6490 OFV units out at tka = -20, because bounds only
# bite far from the mode. A test that only reproduces the centre certifies
# nothing about the surface candidates are actually drawn from.
.sirExpectReproduces <- function(fit, est, offParam = "tka", offValue = -20) {
  testthat::expect_equal(nlmixr2sir:::.sirFitEst(fit), est)

  r <- nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE)
  testthat::expect_lt(r$absDiff, 1e-4)

  ps <- nlmixr2sir:::.sirParamSpace(fit)
  mu <- nlmixr2sir:::.sirProposalMu(fit, ps)
  off <- mu
  off[[offParam]] <- offValue
  mat <- matrix(off, nrow = 1L, dimnames = list(NULL, names(off)))
  got <- unname(nlmixr2sir:::sirEvalOFV(fit, mat, workers = 1L)[[1L]])

  # Reference: the fit's own control, with only the evaluation fields changed.
  ctl <- fit$control
  ctl$maxOuterIterations <- 0L
  ctl$calcTables <- FALSE
  ctl$covMethod <- ""
  ctl$compress <- FALSE
  ctl$print <- 0L
  want <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    rxode2::ini(fit$ui, off), nlmixr2est::getData(fit),
    est = est, control = ctl
  )))$objf

  testthat::expect_true(is.finite(got))
  testthat::expect_equal(got, unname(want), tolerance = 1e-6)
  invisible(got)
}

test_that("unsupported methods are excluded from the allowlist", {
  supported <- nlmixr2sir:::.sirSupportedEstimationMethods
  # Stochastic and non-FOCEi-family methods stay out: their objectives are not
  # reproduced by this evaluator, which is the whole point of the allowlist.
  for (e in c("saem", "imp", "impmap", "qrpem", "npag", "npb", "vae", "emvi")) {
    expect_false(e %in% supported, info = e)
  }
})

test_that("the evaluator reproduces each deterministic method's objective", {
  skip_on_cran()
  # Each rung of the ladder, fitted and then put through the real evaluator --
  # not merely asserted to be present in a character vector.
  for (e in c("focei", "foce", "fo", "foi", "focep", "laplace", "agq")) {
    expect_true(e %in% nlmixr2sir:::.sirSupportedEstimationMethods, info = e)
    fit <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
      theoOneCmt, nlmixr2data::theo_sd, est = e,
      control = list(print = 0L, covMethod = "", calcTables = FALSE)
    )))
    .sirExpectReproduces(fit, e)
    # At the fit's own ETAs these reproduce to rounding, far inside 1e-4, so a
    # drift towards the tolerance shows up here first.
    r <- nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE)
    expect_lt(r$absDiff, 1e-8)
  }
})

test_that("the evaluator reproduces the muModel variants", {
  skip_on_cran()
  for (e in c("mfocei", "ifocei")) {
    expect_true(e %in% nlmixr2sir:::.sirSupportedEstimationMethods, info = e)
    fit <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
      theoOneCmt, nlmixr2data::theo_sd, est = e,
      control = list(print = 0L, covMethod = "", calcTables = FALSE)
    )))
    .sirExpectReproduces(fit, e)
  }
})

test_that("non-default AGQ integration bounds reach the evaluator", {
  skip_on_cran()
  # The regression guard for the dropped-bounds defect. agqLow only changes the
  # objective away from the mode, so this asserts both that the setting arrives
  # AND that it makes a difference -- otherwise the test would pass against an
  # evaluator that ignored it.
  fit <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    theoOneCmt, nlmixr2data::theo_sd, est = "agq",
    control = nlmixr2est::agqControl(
      nAGQ = 2, agqLow = -100, agqHi = Inf,
      print = 0L, covMethod = "", calcTables = FALSE
    )
  )))
  expect_equal(fit$control$agqLow, -100)

  ctl <- nlmixr2sir:::.sirEvalControl(fit)
  expect_equal(ctl$agqLow, -100)
  expect_equal(ctl$agqHi, Inf)
  expect_equal(ctl$maxOuterIterations, 0L)

  .sirExpectReproduces(fit, "agq")

  # Discriminating: the default bounds give a materially different answer at
  # the same off-centre point, so carrying them is not a no-op.
  ps <- nlmixr2sir:::.sirParamSpace(fit)
  mu <- nlmixr2sir:::.sirProposalMu(fit, ps)
  off <- mu
  off[["tka"]] <- -20
  wrong <- ctl
  wrong$agqLow <- -Inf
  got <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(
    rxode2::ini(fit$ui, off), nlmixr2data::theo_sd, est = "agq", control = wrong
  )))$objf
  mat <- matrix(off, nrow = 1L, dimnames = list(NULL, names(off)))
  right <- unname(nlmixr2sir:::sirEvalOFV(fit, mat, workers = 1L)[[1L]])
  expect_gt(abs(got - right), 100)
})

test_that("the evaluator changes only the evaluation fields of the control", {
  skip_on_cran()
  # The invariant that replaced the hand-picked allowlist: everything except
  # the documented overrides is carried through untouched, so a setting nobody
  # anticipated cannot be silently dropped.
  fit <- theoFit()
  ctl <- nlmixr2sir:::.sirEvalControl(fit)
  overrides <- names(nlmixr2sir:::.sirEvalOverrides)
  shared <- setdiff(intersect(names(fit$control), names(ctl)), overrides)
  expect_gt(length(shared), 100)
  for (nm in shared) {
    expect_equal(ctl[[nm]], fit$control[[nm]], info = nm)
  }
  # The overrides are asserted by intent rather than by literal value, because
  # they are stored in the form foceiControl() normalises them to, not the form
  # they are written in: covMethod "" becomes integer 0 and print 0 becomes
  # NULL. Comparing against the raw list would pin the wrong thing.
  expect_equal(ctl$maxOuterIterations, 0L)   # no estimation
  expect_false(isTRUE(ctl$calcTables))       # no output tables
  expect_false(isTRUE(ctl$compress))
  expect_equal(as.integer(ctl$covMethod), 0L) # no covariance step
  expect_true(is.null(ctl$print) || isTRUE(ctl$print == 0)) # quiet
  # Every override key is still present -- assigning NULL must not delete one.
  for (nm in overrides) {
    expect_true(nm %in% names(ctl), info = nm)
  }
})

test_that("the evaluator scores a fit with its own method, not always focei", {
  skip_on_cran()
  fitFo <- theoFitFo()
  expect_equal(.sirFitEst(fitFo), "fo")

  # The evaluator must reproduce the FO objective at the FO estimates.
  r <- nlmixr2sir:::.sirCheckObjective(fitFo, workers = 1L, stencil = FALSE)
  expect_lt(r$absDiff, 1e-4)

  # Discriminating half. The premise is that FO and FOCEi are genuinely
  # different surfaces on this model, so "reproduces its own objective" is not
  # a claim that holds trivially. The gap is about 1.5 OFV units here (it is
  # larger on models with more etas), which is far above the 1e-4 tolerance the
  # reproduction is asserted at.
  gap <- abs(fitFo$objf - theoFit()$objf)
  expect_gt(gap, 0.5)
  expect_lt(r$absDiff, gap / 100)

  # Before the evaluator dispatched on the fit's own method, a hardcoded
  # est = "focei" scored an FO fit 24 units low on a three-eta model -- and
  # below FOCEi's own minimum, because the etas were being estimated rather
  # than held at zero, which is the one thing FO must not do.
  expect_equal(unname(r$reevaluated), unname(fitFo$objf), tolerance = 1e-6)
})

test_that("the evaluator control is built with the fit's own constructor", {
  skip_on_cran()
  expect_equal(nlmixr2sir:::.sirEvalMethod(theoFitFo()), "fo")
  expect_equal(nlmixr2sir:::.sirEvalMethod(theoFit()), "focei")
  ctl <- nlmixr2sir:::.sirEvalControl(theoFitFo())
  expect_equal(ctl$maxOuterIterations, 0L)
  expect_true(isTRUE(as.logical(ctl$fo)))
})
