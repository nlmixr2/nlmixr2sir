# Registering SIR as a covariance method on the fit.
#
# nlmixr2est 7.1.0 made setCov() an S3 generic dispatched on the covariance
# method, and names SIR as the worked example. Most of what that needs already
# existed here: runSIR() registers its empirical covariance in fit$env$covList
# under "sir" (.sirRegisterCov()), and setCov() installs a cached covariance
# BEFORE it dispatches, so setCov(fit, "sir") already worked after a run.
#
# Two things did not. "sir" was absent from setCovAllMethods(), so nobody could
# discover it; and on a fit that had never been through runSIR(), setCov()
# dispatched to setCov.default() and said "covariance method 'sir' not
# supported" -- which is wrong. It is supported. There is just nothing to
# install yet, and the remedy is to run runSIR().
#
# setCov.sir() deliberately does NOT run SIR. nlmixr2est's contract says a
# method that cannot compute the covariance should stop(), and computing this
# one means a full SIR run: minutes to hours, needing a schedule and a
# directory, and throwing away the diagnostics that are the point of the
# runSIR() result. A setter should not silently do that.

test_that("sir is discoverable as a covariance method", {
  expect_true("sir" %in% nlmixr2est::setCovAllMethods())
})

test_that("setCov() on a fit that has not run SIR says so", {
  skip_on_cran()
  fit <- blockFit()
  expect_false("sir" %in% names(fit$env$covList))
  err <- tryCatch(
    nlmixr2est::setCov(fit, "sir"),
    error = function(e) conditionMessage(e)
  )
  # The old message claimed the method was unsupported. The useful message
  # names the thing to do instead.
  expect_match(err, "runSIR", fixed = TRUE)
  expect_false(grepl("not supported", err, fixed = TRUE))
})

test_that("setCov() installs the SIR covariance after a run", {
  skip_on_cran()
  # The regression guard: registering a method must not break the cache path
  # that already worked.
  fit <- theoFitForSetCov()
  dir <- withr::local_tempdir()
  .sirQuiet(runSIR(
    fit, nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(recover = FALSE, workers = 1L, objfStencil = FALSE)
  ))
  expect_true("sir" %in% names(fit$env$covList))

  before <- fit$covMethod
  expect_false(identical(before, "sir"))
  nlmixr2est::setCov(fit, "sir")
  expect_equal(fit$covMethod, "sir")

  # It is the SIR matrix that got installed, not merely a relabelled old one.
  expect_true(is.matrix(fit$cov))
  expect_true(all(is.finite(fit$cov)))
})

test_that("the SIR covariance is on fit$cov's own names and order", {
  skip_on_cran()
  # .sirCovAsFitCov() maps SIR's parameter names onto fit$cov's via
  # .sirParamSpace(), rather than reconstructing them by pattern the way
  # nlmixr2boot has to.
  fit <- theoFit()
  ps <- nlmixr2sir:::.sirParamSpace(fit)
  covSir <- fit$cov
  rownames(covSir) <- colnames(covSir) <-
    ps$sirName[match(rownames(fit$cov), ps$covName)]
  out <- nlmixr2sir:::.sirCovAsFitCov(fit, covSir, ps)
  expect_equal(rownames(out), rownames(fit$cov))
  expect_equal(colnames(out), colnames(fit$cov))
})

test_that("a mismatched SIR covariance is refused rather than installed", {
  skip_on_cran()
  # Fails closed: a covariance whose parameters do not correspond exactly to
  # fit$cov would silently misreport uncertainty if it were reshaped anyway.
  fit <- theoFit()
  ps <- nlmixr2sir:::.sirParamSpace(fit)
  bogus <- matrix(1, 2L, 2L, dimnames = list(c("nope1", "nope2"),
                                             c("nope1", "nope2")))
  expect_message(
    out <- nlmixr2sir:::.sirCovAsFitCov(fit, bogus, ps),
    "do not match"
  )
  expect_null(out)
})
