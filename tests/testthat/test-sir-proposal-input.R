# Alternative proposal sources: RSE percentages and a directly supplied
# covariance. These exist so SIR can run on a fit whose covariance step failed.

test_that(".sirRseVariance uses (rse * estimate / 100)^2 on the diagonal", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  v <- .sirRseVariance(ps, rseTheta = 30)
  expected <- (30 * ps$est / 100)^2
  expect_equal(
    unname(v[ps$kind %in% c("theta", "sigma", "omegaDiag")]),
    expected[ps$kind %in% c("theta", "sigma", "omegaDiag")],
    tolerance = 1e-12
  )
})

test_that("a scalar rseTheta fills in rseOmega and rseSigma", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  expect_equal(
    .sirRseVariance(ps, rseTheta = 30),
    .sirRseVariance(ps, rseTheta = 30, rseOmega = 30, rseSigma = 30)
  )
})

test_that("a vector rseTheta does not fill in the others", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  expect_error(
    .sirRseVariance(ps, rseTheta = c(10, 20, 30)),
    "required when building the proposal from RSE"
  )
  expect_silent(
    .sirRseVariance(ps, rseTheta = c(10, 20, 30), rseOmega = 25, rseSigma = 15)
  )
})

test_that("rseOmega without rseTheta is an error", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  expect_error(.sirRseVariance(ps, rseOmega = 30), "without .*rseTheta")
})

test_that("rse vectors must be length 1 or one per estimated element", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  expect_error(
    .sirRseVariance(ps, rseTheta = c(10, 20)),
    "one value per"
  )
})

test_that(".sirRseVariance derives OMEGA off-diagonals by PsN's formula", {
  skip_on_cran()
  ps <- .sirParamSpace(blockFit())
  v <- .sirRseVariance(ps, rseTheta = 20, rseOmega = c(30, 40))
  off <- which(ps$kind == "omegaOffdiag")
  # neta1 = 2 (eta.cl), neta2 = 1 (eta.ka)
  varI <- ps$est[ps$kind == "omegaDiag" & ps$neta1 == 2L]
  varJ <- ps$est[ps$kind == "omegaDiag" & ps$neta1 == 1L]
  n <- (100 / 40)^2 + (100 / 30)^2 + 1
  expect_equal(
    unname(v[off]),
    (ps$est[off]^2 + varI * varJ) / n,
    tolerance = 1e-12
  )
})

test_that("covmatInput = 'identity' gives the identity in SIR names", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  m <- .sirProposalFromCovmatInput(ps, "identity")
  expect_equal(unname(m), diag(1, nrow(ps)))
  expect_equal(rownames(m), ps$sirName)
})

test_that("an unnamed covmatInput must cover the whole parameter vector", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  good <- diag(0.5, nrow(ps))
  expect_equal(unname(.sirProposalFromCovmatInput(ps, good)), good)
  expect_error(
    .sirProposalFromCovmatInput(ps, diag(0.5, nrow(ps) - 1L)),
    "must be"
  )
})

test_that("a named covmatInput may be a subset and is reordered", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  sub <- matrix(
    c(4, 1, 1, 9),
    2L,
    2L,
    dimnames = list(c("tcl", "tka"), c("tcl", "tka"))
  )
  m <- .sirProposalFromCovmatInput(ps, sub)
  expect_equal(m["tcl", "tcl"], 4)
  expect_equal(m["tka", "tka"], 9)
  expect_equal(m["tka", "tcl"], 1)
  expect_equal(m["tv", "tv"], 0)
})

test_that("a named covmatInput also accepts fit$cov names", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  sub <- matrix(1, 1L, 1L, dimnames = list("om.eta.ka", "om.eta.ka"))
  m <- .sirProposalFromCovmatInput(ps, sub)
  expect_equal(m["eta.ka", "eta.ka"], 1)
})

test_that("covmatInput rejects names that match no parameter", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  bad <- matrix(1, 1L, 1L, dimnames = list("nope", "nope"))
  expect_error(.sirProposalFromCovmatInput(ps, bad), "do not match")
})

test_that("runSIRControl rejects incompatible proposal sources", {
  skip_on_cran()
  expect_error(
    runSIRControl(objfStencil = FALSE, rseTheta = 30, covmatInput = "identity"),
    "alternative proposal sources"
  )
  expect_error(
    runSIRControl(objfStencil = FALSE, rseTheta = 30, thetaInflation = 2),
    "Inflation cannot be combined"
  )
})

test_that(".sirResolveInitialProposal follows PsN's dispatch order", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  expect_equal(
    .sirResolveInitialProposal(fit, ps, runSIRControl(objfStencil = FALSE))$source,
    "cov"
  )
  expect_equal(
    .sirResolveInitialProposal(fit, ps, runSIRControl(objfStencil = FALSE, rseTheta = 30))$source,
    "rse"
  )
  expect_equal(
    .sirResolveInitialProposal(
      fit,
      ps,
      runSIRControl(objfStencil = FALSE, covmatInput = "identity")
    )$source,
    "covmatInput"
  )
})

# The gate for this phase: the covariance step is no longer required.
test_that("runSIR runs on a fit with no covariance step, via rseTheta", {
  skip_on_cran()
  fit <- theoFitNoCov()
  expect_null(fit$cov)
  expect_error(sirGetProposalCov(fit), "covariance matrix")

  tmp <- tempfile("sir_nocov_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  set.seed(20260913)
  res <- .sirQuiet(runSIR(
    fit,
    nSamples = c(16L, 16L),
    nResample = c(8L, 8L),
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L, rseTheta = 30)
  ))
  expect_s3_class(res, "nlmixr2SIR")
  expect_true(all(c("tka", "add.sd", "eta.ka") %in% res$param))
})

test_that("runSIR runs on a fit with no covariance step, via covmatInput", {
  skip_on_cran()
  tmp <- tempfile("sir_ident_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  set.seed(20260913)
  res <- .sirQuiet(runSIR(
    theoFitNoCov(),
    nSamples = c(16L, 16L),
    nResample = c(8L, 8L),
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, 
      recover = FALSE,
      workers = 1L,
      covmatInput = "identity",
      thetaInflation = 0.01
    )
  ))
  expect_s3_class(res, "nlmixr2SIR")
})

# rawresInput: PsN's iteration 0. The empirical mean and covariance of the
# supplied parameter vectors become the iteration-1 proposal.

test_that(".sirRawResultsMatrix returns vectors in SIR parameter order", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  m <- .sirRawResultsMatrix(theoFit(), ps, sirRawResultsPath())
  expect_equal(colnames(m), ps$sirName)
  expect_false(anyNA(m))
  expect_gt(nrow(m), 1L)
})

test_that(".sirRawResultsMatrix reads OMEGA from the matrix, not by name", {
  skip_on_cran()
  # parseRawResultsParams() hands OMEGA back as a matrix keyed by eta name,
  # not as a column named omega(eta.ka,eta.ka), so it is read by index.
  ps <- .sirParamSpace(theoFit())
  m <- .sirRawResultsMatrix(theoFit(), ps, sirRawResultsPath())
  expect_true("eta.ka" %in% colnames(m))
  expect_true(all(m[, "eta.ka"] > 0))
})

test_that("offsetRawres drops low-numbered samples", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  n0 <- nrow(.sirRawResultsMatrix(
    theoFit(),
    ps,
    sirRawResultsPath(),
    offsetRawres = 0L
  ))
  n1 <- nrow(.sirRawResultsMatrix(
    theoFit(),
    ps,
    sirRawResultsPath(),
    offsetRawres = 1L
  ))
  # sample 0 is the reference row, dropped by the default offset of 1
  expect_equal(n0 - n1, 1L)
})

test_that("inFilter narrows the raw-results rows used", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  all_rows <- .sirRawResultsMatrix(theoFit(), ps, sirRawResultsPath())
  few <- .sirRawResultsMatrix(
    theoFit(),
    ps,
    sirRawResultsPath(),
    inFilter = function(d) d$sample <= 4L
  )
  expect_lt(nrow(few), nrow(all_rows))
})

test_that("an over-narrow filter is an error, not a silent empty proposal", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  expect_error(
    .sirRawResultsMatrix(
      theoFit(),
      ps,
      sirRawResultsPath(),
      inFilter = function(d) d$sample > 10000L
    ),
    "No raw-results rows survived"
  )
})

test_that(".sirProposalFromRawResults centres on the supplied vectors", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  got <- .sirProposalFromRawResults(fit, ps, sirRawResultsPath())
  m <- .sirRawResultsMatrix(fit, ps, sirRawResultsPath())
  expect_equal(got$mu, colMeans(m), tolerance = 1e-12)
  expect_equal(dim(got$covMat), c(nrow(ps), nrow(ps)))
  expect_equal(got$nVectors, nrow(m))
})

test_that("runSIR can be seeded from a raw-results file", {
  skip_on_cran()
  tmp <- tempfile("sir_rr_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  set.seed(2)
  res <- .sirQuiet(runSIR(
    theoFit(),
    nSamples = 16L,
    nResample = 8L,
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, 
      recover = FALSE,
      workers = 1L,
      rawresInput = sirRawResultsPath()
    )
  ))
  expect_s3_class(res, "nlmixr2SIR")
})

test_that("runSIRControl rejects more than one proposal source", {
  skip_on_cran()
  expect_error(
    runSIRControl(objfStencil = FALSE, rawresInput = "x.csv", rseTheta = 30),
    "alternative proposal sources"
  )
  expect_error(
    runSIRControl(objfStencil = FALSE, rawresInput = "x.csv", covmatInput = "identity"),
    "alternative proposal sources"
  )
  expect_warning(runSIRControl(objfStencil = FALSE, inFilter = function(d) TRUE), "no effect")
})
