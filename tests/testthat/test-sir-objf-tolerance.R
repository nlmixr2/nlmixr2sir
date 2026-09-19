# The objective-reproduction tolerance.
#
# The default was 1e-4, which was not a property of any fit. Measured across
# models, the reproduction error -- the gap between fit$objf and a fresh
# re-evaluation at the same estimates -- spans about 1e-6 to 1.2e-3, and 1e-4
# sat in the middle of that range. So it misfired erratically: theoFit()
# reproduces to 8.2e-5 and passed, while a near-identical one-eta fit
# reproduced to 1.107e-4 and aborted. The package's own vignette and its own
# runSIR() example both aborted, which is to say runSIR() did not work on
# ordinary models.
#
# What causes the gap is inner-solve convergence slack: fit$objf comes from the
# final outer iteration's inner solve with warm-started etas, and the
# re-evaluation solves the inner problem fresh. Only `sigdig` predicts its size
# (~3.4x per digit). Eta count does not -- flat from 1 to 12 etas -- and
# neither does design collinearity. Both were tested and falsified.
#
# The threshold is therefore set from what the error AFFECTS rather than from
# any model of its size. A dOFV error of d moves an importance weight by
# exp(-d/2): at 1e-2 that is 0.5%, against dOFV of order 1 to 10. It still
# catches the two mismatches the preflight exists for by more than two orders
# of magnitude -- a SAEM fit scored under FOCEi is 2.69 units out on theo_sd,
# and the dropped-agqLow defect was 6490.

test_that("the default tolerance is 1e-2", {
  expect_equal(runSIRControl()$objfTolerance, 1e-2)
  expect_equal(formals(nlmixr2sir:::.sirCheckObjective)$objfTolerance, 1e-2)
})

test_that("the warning band sits an order of magnitude below the abort", {
  # One order of magnitude: wide enough that an ordinary fit does not trip it,
  # narrow enough that a fit approaching the abort is announced first.
  expect_equal(nlmixr2sir:::.sirObjfWarnTolerance, 1e-3)
  expect_lt(nlmixr2sir:::.sirObjfWarnTolerance, runSIRControl()$objfTolerance)
})

test_that("an ordinary fit reproduces its objective silently", {
  skip_on_cran()
  # theoFit() reproduces to ~8e-5, comfortably inside the band.
  # expect_no_warning, not expect_silent: the evaluator emits worker-plan
  # messages, which are not the subject here.
  expect_no_warning(
    nlmixr2sir:::.sirCheckObjective(theoFit(), workers = 1L, stencil = FALSE)
  )
})

test_that("a fit in the warning band is announced but not refused", {
  skip_on_cran()
  # threeEtaFit() reproduces to ~1.15e-3 -- past the warning threshold, well
  # inside the abort. This is the case the old default got wrong: it is a
  # perfectly ordinary three-eta FOCEi fit at the default sigdig.
  fit <- threeEtaFit()
  # Assign inside, rather than relying on what expect_warning() returns.
  r <- NULL
  expect_warning(
    r <- nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE),
    "reproduc"
  )
  expect_gt(r$absDiff, nlmixr2sir:::.sirObjfWarnTolerance)
  expect_lt(r$absDiff, runSIRControl()$objfTolerance)
})

test_that("the warning names sigdig, the one lever known to work", {
  skip_on_cran()
  w <- tryCatch(
    {
      nlmixr2sir:::.sirCheckObjective(threeEtaFit(), workers = 1L, stencil = FALSE)
      NULL
    },
    warning = function(x) conditionMessage(x)
  )
  expect_false(is.null(w))
  expect_match(w, "sigdig")
})

test_that("the abort message points at sigdig rather than at the tolerance alone", {
  skip_on_cran()
  # Raising objfTolerance hides the problem; refitting at a higher sigdig is
  # the only lever measured to reduce the gap, so the message has to say so.
  err <- tryCatch(
    nlmixr2sir:::.sirCheckObjective(
      theoFit(), workers = 1L, objfTolerance = 0, stencil = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "sigdig")
})

test_that("the vignette's own model now passes the preflight", {
  skip_on_cran()
  # The regression this whole change exists for. This is the fit the package
  # vignette runs SIR on; under the old default it aborted, which failed
  # R CMD build.
  fit <- threeEtaFit()
  expect_no_error(suppressWarnings(
    nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE)
  ))
})

test_that("a genuinely different surface is still refused", {
  skip_on_cran()
  # The discriminating half. Raising the default only helps if it still catches
  # what the check exists for. An offset of 2.69 -- the measured SAEM-under-
  # FOCEi gap on theo_sd -- must still abort, with 269x margin over the new
  # default.
  fit <- theoFit()
  r <- suppressWarnings(
    nlmixr2sir:::.sirCheckObjective(fit, workers = 1L, stencil = FALSE)
  )
  expect_lt(r$absDiff, 2.69 / 100)
  expect_error(
    nlmixr2sir:::.sirCheckObjective(
      fit, workers = 1L, objfTolerance = 1e-8, stencil = FALSE
    ),
    "reproduce"
  )
})
