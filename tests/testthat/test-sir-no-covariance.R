# R2: what happens when fit$cov is incomplete, and when it is absent ----------
#
# These are two different situations and the documentation used to conflate
# them:
#
#   * fit$cov present but missing blocks (covFull = FALSE, or a partial
#     covariance). The missing OMEGA/SIGMA uncertainty is derivable -- the
#     Wishart-style approximation uses the OMEGA estimates and the subject
#     count, which are real information -- so the run proceeds and fills them.
#
#   * fit$cov absent entirely (covMethod = "", or a failed covariance step).
#     THETA uncertainty cannot be derived from a fit that never computed any.
#     Manufacturing it from an assumed RSE would invent precisely the quantity
#     SIR exists to measure, so the run stops and asks for one of the
#     alternative proposal sources instead.

test_that("a theta-only covariance is completed by the OMEGA fallback", {
  skip_on_cran()
  fit <- theoFitThetaCov()
  expect_false(is.null(fit$cov))
  expect_false(any(grepl("^om", rownames(fit$cov))))

  ps <- .sirParamSpace(fit)
  prop <- .sirInitialProposal(
    fit,
    mu = .sirProposalMu(fit, ps),
    proposalCov = sirGetProposalCov(fit)
  )

  # Every estimated parameter gets a proposal, including the OMEGA element
  # that fit$cov never carried.
  expect_setequal(prop$paramNames, ps$sirName)
  expect_true("eta.ka" %in% prop$paramNames)
  expect_gt(prop$covMat["eta.ka", "eta.ka"], 0)
  # And the run says which parameters came from the fallback.
  expect_true("eta.ka" %in% prop$fallbackNames)
})

test_that("runSIR completes end-to-end on a theta-only covariance", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  set.seed(21)
  res <- .sirQuiet(runSIR(
    theoFitThetaCov(),
    nSamples = 16L,
    nResample = 8L,
    directory = dir,
    control = runSIRControl(recover = FALSE, workers = 1L)
  ))
  expect_s3_class(res, "nlmixr2SIR")
  expect_true("eta.ka" %in% res$param)
})

test_that("an absent covariance aborts naming the alternative proposal sources", {
  skip_on_cran()
  fit <- theoFitNoCov()
  expect_null(fit$cov)

  err <- tryCatch(sirGetProposalCov(fit), error = function(e) conditionMessage(e))

  # The old message said "re-run the model with a successful covariance step",
  # which is no help at all for the models SIR exists to serve. It must point
  # at the three routes that do work.
  expect_match(err, "rseTheta")
  expect_match(err, "covmatInput")
  expect_match(err, "rawresInput")
})

test_that("the default runSIR call on a no-covariance fit is actionable", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  err <- tryCatch(
    .sirQuiet(runSIR(
      theoFitNoCov(),
      nSamples = 16L,
      nResample = 8L,
      directory = dir,
      control = runSIRControl(recover = FALSE, workers = 1L)
    )),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "rseTheta")
  expect_match(err, "covmatInput|rawresInput")
})
