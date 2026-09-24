# P7: the importance-sampling family scores on FOCEi, by documented exception.
#
# imp/impmap/qrpem were excluded from .sirSupportedEstimationMethods because
# nlmixr2est 7.0.3 had no way to evaluate at fixed parameters. 7.1.0 added one
# (nIter = 0, the EONLY=1 analogue), which re-opened the question of HOW such a
# candidate should be scored.
#
# The answer is not the obvious one. nlmixr2est recomputes the objective of
# EVERY imp-family fit as a nested FOCEi evaluation at the converged estimates
# (.impmapRecomputeObjf(), nlmixr2est R/impmap.R:1106), because the in-C++
# finalize leaves the eta-Hessian without its data term. So fit$objf on an imp
# fit has never been the importance-sampling objective -- it is a FOCEi number,
# and $impObj is where the importance-sampling one lives.
#
# SIR therefore scores imp candidates directly as FOCEi. That is a deliberate
# exception to the rule that est selects the objective, and these tests pin
# both halves of it: that the exception produces the right number, and that it
# is genuinely an exception rather than the ordinary carried-control path.

# Muffle ONLY the scoring notice, so that the tests below are about what they
# say they are about. Deliberately not suppressWarnings(): any other warning
# still surfaces and still fails a suite that is expected to be quiet. The
# notice itself is asserted directly, further down.
withoutScoringNotice <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("score", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

test_that("the imp family is admitted to the allowlist", {
  supported <- nlmixr2sir:::.sirSupportedEstimationMethods
  for (e in c("imp", "impmap", "qrpem")) {
    expect_true(e %in% supported, info = e)
    # Admitted, but not silently: see the scoring-notice tests below.
    expect_warning(nlmixr2sir:::.sirSupportedEstimation(e), "FOCEi", info = e)
  }
})

test_that("an imp-family fit is scored as focei, not as its own method", {
  skip_on_cran()
  fit <- impmapFit()
  # The recorded method is untouched: the allowlist and every provenance record
  # still see impmap. Only the evaluator redirects.
  expect_equal(nlmixr2sir:::.sirFitEst(fit), "impmap")
  expect_equal(nlmixr2sir:::.sirEvalMethod(fit), "focei")
})

test_that("the imp-family evaluation control is a bare foceiControl carrying sigdig", {
  skip_on_cran()
  fit <- impmapFit()
  ctl <- nlmixr2sir:::.sirEvalControl(fit)
  expect_s3_class(ctl, "foceiControl")
  expect_equal(ctl$maxOuterIterations, 0L)
  expect_equal(ctl$sigdig, fit$foceiControl$sigdig)
  expect_false(isTRUE(ctl$calcTables))
  expect_false(isTRUE(ctl$compress))
  expect_equal(as.integer(ctl$covMethod), 0L)

  # Discriminating: isample and nIter exist only on impmapControl. Their
  # absence is what proves the fit's own control was NOT carried forward, which
  # is the whole exception. An evaluator that took the ordinary path would
  # carry them.
  expect_true("isample" %in% names(fit$control))
  expect_true("nIter" %in% names(fit$control))
  expect_false("isample" %in% names(ctl))
  expect_false("nIter" %in% names(ctl))
})

test_that("the evaluator reproduces a one-eta imp fit's objective exactly", {
  skip_on_cran()
  # One eta is the case that matters: the eta-Hessian defect that forced the
  # upstream recompute is specific to it, and the published C++ objective is
  # ~19.96 units low there (173.63 against 193.60 on this model). Reproducing
  # fit$objf therefore proves SIR is on the recomputed FOCEi surface and not on
  # the raw importance-sampling one.
  fit <- impmapFit()
  expect_equal(nrow(fit$omega), 1L)
  r <- withoutScoringNotice(
    nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE)
  )
  # Measured bit-identical, so this is asserted far tighter than the 1e-4 the
  # preflight itself allows.
  expect_equal(unname(r$reevaluated), unname(fit$objf), tolerance = 1e-12)

  # Discriminating half: $impObj is a genuinely different number on the same
  # fit, so "reproduces its objective" is not a claim that holds for any
  # objective this fit carries.
  impObj <- fit$env$impObj
  expect_gt(abs(impObj - fit$objf), 0.1)
  expect_lt(r$absDiff, abs(impObj - fit$objf) / 100)
})

test_that("the evaluator reproduces a three-eta imp fit's objective exactly", {
  skip_on_cran()
  # The control case. Models with 2+ random effects are unaffected by the
  # Hessian defect, so this checks the exception is not one-eta-specific.
  fit <- impmapFitThreeEta()
  expect_gt(nrow(fit$omega), 1L)
  r <- withoutScoringNotice(
    nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE)
  )
  expect_equal(unname(r$reevaluated), unname(fit$objf), tolerance = 1e-12)
})

test_that("an imp fit survives the objective stencil", {
  skip_on_cran()
  # The stencil probes either side of every parameter. It is the check that the
  # evaluator tracks the fit's surface AWAY from the centre, not just at it --
  # which is what the agqLow defect showed can differ.
  #
  # The near-optimum warning is expected and is NOT imp-specific: theoFit(), a
  # plain focei fixture, raises the same one. These fixtures are deliberately
  # small and stop a little short of convergence. Asserted rather than muffled
  # so that a genuinely new warning here would still fail the suite.
  # The threshold is injected: since P9 the near-optimum warning fires at 0.1
  # OFV units rather than 1e-2, and this fixture's best probe improves by about
  # 0.027 -- real, but below the level worth reporting now that the measured
  # noise floor can itself reach 0.028.
  fit <- impmapFit()
  expect_warning(
    withoutScoringNotice(nlmixr2sir:::.sirCheckObjective(
      fit, workers = 1L, warnTolerance = 1e-3
    )),
    "local optimum"
  )
})

test_that("the imp evaluator reproduces the objective more exactly than focei's own path", {
  skip_on_cran()
  # Not a vanity comparison: it is the evidence that scoring as FOCEi is the
  # RIGHT quantity rather than merely a close one. The imp route reproduces
  # fit$objf to exactly zero, because fit$objf was produced by this very
  # calculation; the ordinary carried-control focei route lands within its
  # tolerance but not on the nose.
  impDiff <- withoutScoringNotice(nlmixr2sir:::.sirCheckObjective(
    impmapFit(), workers = 1L, stencil = FALSE
  ))$absDiff
  foceiDiff <- nlmixr2sir:::.sirCheckObjective(
    theoFit(), workers = 1L, stencil = FALSE
  )$absDiff
  expect_equal(impDiff, 0)
  expect_gt(foceiDiff, 0)
})

test_that("runSIR completes on an imp fit", {
  skip_on_cran()
  fit <- impmapFit()
  dir <- withr::local_tempdir()
  res <- withoutScoringNotice(.sirQuiet(runSIR(
    fit,
    nSamples = 16L,
    nResample = 8L,
    directory = dir,
    control = runSIRControl(
      recover = FALSE, workers = 1L, objfStencil = FALSE
    )
  )))
  expect_s3_class(res, "nlmixr2SIR")
  expect_equal(nrow(res), length(fit$theta) + 1L)
})

test_that("the stochastic methods stay excluded", {
  # The re-admission is specific to the imp family and rests on the recompute.
  # Nothing else acquired one.
  supported <- nlmixr2sir:::.sirSupportedEstimationMethods
  for (e in c("saem", "npag", "npb", "vae", "emvi", "fbvi")) {
    expect_false(e %in% supported, info = e)
    expect_error(nlmixr2sir:::.sirSupportedEstimation(e), e)
  }
})


# The scoring notice -----------------------------------------------------------
#
# Admitting the family silently would leave a user to infer from nothing that
# their impmap fit was scored by a method they did not choose. The redirect is
# defensible but it is not obvious, so it is announced once per run.

test_that("an imp-family fit warns that candidates are scored as FOCEi", {
  for (e in c("imp", "impmap", "qrpem")) {
    w <- tryCatch(
      {
        nlmixr2sir:::.sirSupportedEstimation(e)
        NULL
      },
      warning = function(x) conditionMessage(x)
    )
    expect_false(is.null(w), info = e)
    # Names the method it was given, and the method it will use instead.
    expect_match(w, e, fixed = TRUE, info = e)
    expect_match(w, "FOCEi", info = e)
  }
})

test_that("the scoring notice does not fire for the deterministic ladder", {
  # Discriminating: a notice that fired for everything would carry no
  # information. Only the redirected family is announced.
  for (e in c("focei", "foce", "fo", "foi", "focep", "laplace", "agq")) {
    expect_silent(nlmixr2sir:::.sirSupportedEstimation(e))
  }
})

test_that("the scoring notice fires once per run, not once per candidate", {
  skip_on_cran()
  # It lives on the preflight, which runs once. A notice attached to the
  # evaluator would fire nSamples times and be worse than useless.
  fit <- impmapFit()
  dir <- withr::local_tempdir()
  w <- character(0)
  withCallingHandlers(
    .sirQuiet(runSIR(
      fit, nSamples = 16L, nResample = 8L, directory = dir,
      control = runSIRControl(
        recover = FALSE, workers = 1L, objfStencil = FALSE
      )
    )),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(grep("FOCEi", w), 1L)
})

# Provenance -------------------------------------------------------------------

test_that("the fingerprint records the method candidates were scored with", {
  skip_on_cran()
  # estMethod alone cannot describe an imp run: it says impmap, while the
  # numbers were produced by FOCEi. Both are recorded.
  fit <- impmapFit()
  ps <- nlmixr2sir:::.sirParamSpace(fit)
  fp <- nlmixr2sir:::.sirRunFingerprint(
    fit, ps, nlmixr2sir:::.sirSchedule(16L, 8L),
    runSIRControl(objfStencil = FALSE, workers = 1L)
  )
  expect_equal(fp$estMethod, "impmap")
  expect_equal(fp$evalMethod, "focei")
})

test_that("a change in scoring method invalidates a saved run", {
  skip_on_cran()
  # evalMethod is derived from estMethod today, so it adds no discriminating
  # power right now. It is in the identity fields so that if the mapping ever
  # changes, a state file written under the old mapping is refused rather than
  # silently reused.
  expect_true("evalMethod" %in% nlmixr2sir:::.sirIdentityFields)

  fit <- impmapFit()
  ps <- nlmixr2sir:::.sirParamSpace(fit)
  current <- nlmixr2sir:::.sirRunFingerprint(
    fit, ps, nlmixr2sir:::.sirSchedule(16L, 8L),
    runSIRControl(objfStencil = FALSE, workers = 1L)
  )
  saved <- current
  saved$evalMethod <- "laplace"
  expect_true("evalMethod" %in% nlmixr2sir:::.sirCompareFingerprints(saved, current))
})
