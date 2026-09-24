# The objective preflight after P9.
#
# Four changes, all resting on one decomposition of what a candidate's dOFV
# actually contains:
#
#   dOFV_i = [ f(x_i) - f(xhat) ]            the signal
#          + [ e(x_i) - e(xhat) ]            real noise
#          + [ e(xhat) - e_fit(xhat) ]       a CONSTANT across candidates
#
# The preflight used to abort on the third term. That term is identical for
# every candidate, and a constant shift in dOFV multiplies every weight by the
# same factor, which divides out of the normalised weights -- as R/sir-iterate.R
# already says where recentring moves the reference. Measured on Rik
# Schoemaker's QR models, that term was 10x to 116x larger than the second one,
# so the check was refusing runs over the part that provably cancels while
# saying nothing about the part that does not.
#
# PsN does not check this at all (lib/tool/sir.pm): its reference is the
# original .lst OFV, candidates are scored through MAXEVAL=0, the two are never
# compared, and the centre's dOFV is hardcoded to zero. So:
#
#  1. the initial reference is now the RE-EVALUATED centre, which makes the
#     third term zero by construction rather than by assertion -- stricter than
#     PsN, which sets it to zero by fiat;
#  2. the abort is relative to the objective, because slack and a genuinely
#     different surface separate by 1600x on that scale and only 10x on an
#     absolute one;
#  3. the noise floor -- the second term -- is measured and gated instead;
#  4. the old warning's weight claim, which was simply wrong, is gone.

test_that("the abort tolerance is relative, with an absolute floor", {
  expect_equal(runSIRControl()$objfTolerance, 1e-3)
  # A relative tolerance alone would be absurdly tight on a small objective, so
  # it never goes below the floor.
  expect_equal(nlmixr2sir:::.sirObjfAbsFloor, 1e-2)
  thr <- nlmixr2sir:::.sirObjfAbortThreshold
  expect_equal(thr(117, 1e-3), 0.117)      # theo_sd scale
  expect_equal(thr(19367, 1e-3), 19.367)   # QR scale
  expect_equal(thr(0.5, 1e-3), 1e-2)       # floor wins
})

test_that("the relative abort still separates slack from a different surface", {
  thr <- nlmixr2sir:::.sirObjfAbortThreshold(117, 1e-3)
  # A SAEM fit scored under FOCEi is 2.69 OFV units out on theo_sd; that is the
  # case this check exists for and it must still abort.
  expect_gt(2.69, thr)
  # The worst convergence slack measured on the QR models, scaled to theo_sd's
  # objective, must not.
  expect_lt(1.2e-3, thr)
})

test_that("the preflight aborts only on egregious disagreement", {
  skip_on_cran()
  fit <- theoFit()
  # Ordinary slack passes silently at the default.
  expect_no_error(nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE))
  # A surface-scale disagreement still aborts.
  expect_error(
    nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE,
                                    objfTolerance = 0),
    "reproduce"
  )
})

test_that("the abort message no longer claims a weight impact it cannot know", {
  skip_on_cran()
  err <- tryCatch(
    nlmixr2sir:::.sirCheckObjective(theoFit(), workers = 1L, stencil = FALSE,
                                    objfTolerance = 0),
    error = function(e) conditionMessage(e)
  )
  # The old text asserted "well under one percent", which is false at a 0.27
  # gap (exp(-0.135) is 13%). It was harmless only because the gap cancels --
  # a different reason from the one it gave.
  expect_false(grepl("under one percent", err, fixed = TRUE))
  expect_match(err, "sigdig")
})

test_that("runSIR measures dOFV against the re-evaluated centre", {
  skip_on_cran()
  # The substance of change 1. The initial reference must be the evaluator's
  # own centre, not fit$objf, so that iteration 1 is on the same footing as
  # every later iteration.
  fit <- theoFit()
  pre <- nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE)
  dir <- withr::local_tempdir()
  res <- .sirQuiet(runSIR(
    fit, nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(recover = FALSE, workers = 1L, objfStencil = FALSE,
                            objfNoise = FALSE)
  ))
  expect_equal(attr(res, "initialReferenceOfv"), unname(pre$reevaluated),
               tolerance = 1e-12)
  # Discriminating: it is NOT fit$objf, and the two genuinely differ here.
  expect_gt(abs(unname(pre$reevaluated) - fit$objf), 0)
})

test_that("the evaluator noise floor is measured and is deterministic", {
  skip_on_cran()
  fit <- theoFit()
  ps <- nlmixr2sir:::.sirParamSpace(fit)
  a <- nlmixr2sir:::.sirObjectiveNoise(fit, ps, workers = 1L)
  b <- nlmixr2sir:::.sirObjectiveNoise(fit, ps, workers = 1L)
  expect_true(is.finite(a$noise))
  expect_gte(a$noise, 0)
  # Same answer twice: the direction is drawn from a fixed local seed, so the
  # measurement cannot wander between runs -- and must not disturb the run's
  # own RNG stream either.
  expect_equal(a$noise, b$noise)
})

test_that("measuring the noise floor does not disturb the sampling stream", {
  skip_on_cran()
  # If it consumed from the global stream, every run's samples would shift.
  fit <- theoFit()
  ps <- nlmixr2sir:::.sirParamSpace(fit)
  set.seed(99); before <- runif(3)
  set.seed(99); invisible(nlmixr2sir:::.sirObjectiveNoise(fit, ps, workers = 1L))
  after <- runif(3)
  expect_equal(before, after)
})

test_that("the noise floor is reported, never refused", {
  skip_on_cran()
  # An earlier draft aborted above the tolerance. It then refused QR model
  # N021 over a floor of 73.94 OFV units that was not noise at all: a +/-1 SD
  # transect on that model swings 1208 units, and a quartic could not absorb
  # the shape, so the residual was unfitted curvature. A diagnostic whose own
  # estimator can be wrong must not stop a run.
  fit <- theoFit()
  expect_warning(
    nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE,
                                    noise = TRUE, noiseTolerance = 0),
    "noise floor"
  )
  expect_no_error(suppressWarnings(
    nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE,
                                    noise = TRUE, noiseTolerance = 0)
  ))
  # An ordinary floor is not even mentioned.
  expect_no_warning(
    nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE,
                                    noise = TRUE, noiseTolerance = 1)
  )
})

test_that("an untrustworthy noise estimate is withheld rather than reported", {
  skip_on_cran()
  # Wide span: the polynomial cannot absorb the objective's shape, so the
  # residual is curvature. The estimator must say nothing rather than report a
  # number that would be read as noise.
  fit <- theoFit()
  ps <- nlmixr2sir:::.sirParamSpace(fit)
  wide <- nlmixr2sir:::.sirObjectiveNoise(fit, ps, span = 50, workers = 1L)
  narrow <- nlmixr2sir:::.sirObjectiveNoise(fit, ps, workers = 1L)
  expect_true(is.finite(narrow$noise))
  expect_true(is.na(wide$noise))
})
