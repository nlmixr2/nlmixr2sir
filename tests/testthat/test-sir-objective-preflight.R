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

test_that("the preflight accepts a fit whose centre reproduces its objective", {
  skip_on_cran()
  fit <- theoFit()
  expect_no_error(.sirCheckObjective(fit, workers = 1L))
})

test_that("the preflight returns the stored and reevaluated objectives", {
  skip_on_cran()
  fit <- theoFit()
  res <- .sirCheckObjective(fit, workers = 1L)
  expect_named(res, c("stored", "reevaluated", "absDiff", "relDiff", "stencil"))
  expect_equal(res$stored, fit$objf, tolerance = 1e-12)
  expect_lt(res$absDiff, 1e-3)
})

test_that("the preflight aborts when the centre does not reproduce the objective", {
  skip_on_cran()
  fit <- theoFit()
  # A tolerance tight enough that even the genuine numerical difference between
  # the stored and reevaluated objective fails it. This is the mismatch path:
  # the message must name both values so the user can judge the gap.
  expect_error(
    .sirCheckObjective(fit, workers = 1L, objfTolerance = 0),
    "objective"
  )
  err <- tryCatch(
    .sirCheckObjective(fit, workers = 1L, objfTolerance = 0),
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
  expect_error(
    suppressMessages(runSIR(
      fit,
      nSamples = 16L,
      nResample = 8L,
      directory = dir,
      control = runSIRControl(
        recover = FALSE,
        workers = 1L,
        objfTolerance = 0
      )
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
  expect_error(
    .sirCheckObjective(fit, workers = 1L, objfTolerance = 0, stencil = FALSE),
    "absolute"
  )
})

test_that("the deterministic ladder is accepted, not just focei", {
  skip_on_cran()
  supported <- nlmixr2sir:::.sirSupportedEstimationMethods
  for (e in c("focei", "foce", "fo", "foi", "focep", "laplace", "agq")) {
    expect_true(e %in% supported, info = e)
  }
  # Stochastic and non-FOCEi-family methods stay out: their objectives are not
  # reproduced by this evaluator, which is the whole point of the allowlist.
  for (e in c("saem", "imp", "impmap", "qrpem", "npag", "npb", "vae", "emvi")) {
    expect_false(e %in% supported, info = e)
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
