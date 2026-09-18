# setCov(fit, "sir"): SIR as an nlmixr2est covariance method, cached by its
# options and by the covariance it is seeded from.

# A fresh fit per test: setCov() mutates fit$env, which would leak between
# tests through the shared fixtures.
.freshTheoFit <- function(covMethod = "r") {
  suppressMessages(suppressWarnings(nlmixr2utils::nlmixr2(
    theoOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L, covMethod = covMethod)
  )))
}

# Five parameters, so the smallest schedule that reaches full rank. The
# stencil preflight is off: this tiny fixture stops just short of the optimum
# and would warn in every test, which is covered elsewhere.
.smallSir <- function(...) {
  sirControl(
    nSamples = c(16L, 16L),
    nResample = c(8L, 8L),
    workers = 1L,
    objfStencil = FALSE,
    ...
  )
}

# Count the SIR runs setCov() makes, while still running them.
.countRuns <- function(env = parent.frame()) {
  counter <- new.env()
  counter$n <- 0L
  counter$args <- list()
  real <- .sirRunCore
  local_mocked_bindings(
    .sirRunCore = function(...) {
      counter$n <- counter$n + 1L
      counter$args[[counter$n]] <- list(...)
      real(...)
    },
    .env = env
  )
  counter
}

test_that("sirControl() has covariance-step defaults and validates", {
  ctl <- sirControl()
  expect_s3_class(ctl, "sirControl")
  expect_equal(ctl$nSamples, c(1000L, 1000L, 1000L, 2000L, 2000L))
  expect_equal(ctl$nResample, c(200L, 400L, 500L, 1000L, 1000L))
  expect_null(ctl$seedCov)
  expect_equal(ctl$rseTheta, 30)
  expect_equal(ctl$seed, 42L)
  # 1L and 1 are the same option
  expect_identical(sirControl(thetaInflation = 2L), sirControl(thetaInflation = 2))

  expect_error(sirControl(nSamples = 10, nResample = c(5, 5)))
  expect_error(sirControl(nSamples = 10, nResample = 20), "capResampling")
  expect_error(sirControl(capCorrelation = 2))
  expect_error(sirControl(seedCov = "sir"), "own result")
  expect_error(sirControl(seedCov = matrix(1, 2, 3)))
  expect_error(sirControl(rseTheta = NULL, rseOmega = 20), "rseTheta")
  expect_error(sirControl(seed = 1.5))
  expect_message(print(ctl), "nlmixr2sir covariance control")
})

test_that("rxUiDeparse round-trips sirControl()", {
  skip_on_cran()
  ctl <- sirControl(thetaInflation = 2, nSamples = c(50, 50), nResample = c(10, 20))
  x <- eval(rxode2::rxUiDeparse(ctl, "x"))
  expect_identical(x, ctl)
  expect_equal(
    deparse(rxode2::rxUiDeparse(sirControl(), "y")),
    "y <- sirControl()"
  )
})

test_that("setCov(fit, 'sir') installs SIR and keeps the result on the fit", {
  skip_on_cran()
  fit <- .freshTheoFit()
  method0 <- fit$covMethod
  cov0 <- fit$cov
  se0 <- fit$parFixedDf$SE
  expect_true("sir" %in% nlmixr2est::setCovAllMethods())

  runs <- .countRuns()
  .sirQuiet(suppressMessages(nlmixr2est::setCov(fit, "sir", control = .smallSir())))
  expect_equal(runs$n, 1L)
  expect_identical(fit$covMethod, "sir")
  expect_s3_class(fit$sir, "nlmixr2SIR")
  expect_true(method0 %in% names(fit$covList))
  expect_setequal(rownames(fit$cov), rownames(cov0))
  expect_false(inherits(try(chol(fit$cov), silent = TRUE), "try-error"))
  expect_false(isTRUE(all.equal(fit$parFixedDf$SE, se0)))
  # the run is given the resolved seed, not left to read fit$cov itself
  expect_identical(runs$args[[1]]$seedCov, cov0)
  expect_false(runs$args[[1]]$parFixedSe)
  opts <- fit$env$covOptions$sir
  expect_identical(opts$seedMethod, method0)
  expect_equal(opts$seedCov, cov0)
  expect_equal(opts$nSamples, c(16L, 16L))
  expect_null(opts$workers)

  # the same options and seed: already installed, then served from the cache
  expect_error(
    nlmixr2est::setCov(fit, "sir", control = .smallSir()),
    "no need to switch"
  )
  suppressMessages(nlmixr2est::setCov(fit, method0))
  expect_equal(fit$cov, cov0)
  suppressMessages(nlmixr2est::setCov(fit, "sir", control = .smallSir()))
  expect_equal(runs$n, 1L)
  expect_identical(fit$covMethod, "sir")

  # workers do not change the answer, so they do not recompute
  expect_error(
    nlmixr2est::setCov(fit, "sir", control = .smallSir(rxThreads = 1L)),
    "no need to switch"
  )

  # different options recompute, seeded from the recorded seed, not from SIR
  .sirQuiet(suppressMessages(
    nlmixr2est::setCov(fit, "sir", control = .smallSir(seed = 7L))
  ))
  expect_equal(runs$n, 2L)
  expect_equal(runs$args[[2]]$seedCov, cov0)
  expect_false(runs$args[[2]]$parFixedSe)
  expect_identical(fit$env$covOptions$sir$seedMethod, method0)
})

test_that("a seed that is symmetric only to rounding still reuses the cache", {
  skip_on_cran()
  # Estimation-time covariances can be off symmetric by ~1e-17; nlmixr2est
  # symmetrizes one when it reinstalls it, which must not look like a new seed.
  fit <- .freshTheoFit()
  method0 <- fit$covMethod
  cov <- fit$env$cov
  cov[1L, 2L] <- cov[1L, 2L] * (1 + 4 * .Machine$double.eps)
  assign("cov", cov, envir = fit$env)
  expect_false(isSymmetric(unname(fit$cov), tol = 0))

  runs <- .countRuns()
  .sirQuiet(suppressMessages(nlmixr2est::setCov(fit, "sir", control = .smallSir())))
  suppressMessages(nlmixr2est::setCov(fit, method0))
  expect_true(isSymmetric(unname(fit$cov), tol = 0))
  suppressMessages(nlmixr2est::setCov(fit, "sir", control = .smallSir()))
  expect_equal(runs$n, 1L)
  expect_identical(fit$covMethod, "sir")
})

test_that("a different seed covariance recomputes SIR", {
  skip_on_cran()
  fit <- .freshTheoFit()
  runs <- .countRuns()
  .sirQuiet(suppressMessages(nlmixr2est::setCov(fit, "sir", control = .smallSir())))
  expect_equal(runs$n, 1L)

  # A THETA-only covariance becomes the installed one, and so the seed
  suppressMessages(suppressWarnings(nlmixr2est::setCov(fit, "r")))
  expect_identical(fit$covMethod, "r")
  .sirQuiet(suppressMessages(nlmixr2est::setCov(fit, "sir", control = .smallSir())))
  expect_equal(runs$n, 2L)
  expect_identical(fit$env$covOptions$sir$seedMethod, "r")
  # the SIR covariance is full-shape even from a THETA-only seed
  expect_true("om.eta.ka" %in% rownames(fit$cov))
})

test_that("the same control and seed give the same covariance, and restore the RNG", {
  skip_on_cran()
  fit1 <- .freshTheoFit()
  fit2 <- .freshTheoFit()
  set.seed(99)
  before <- .Random.seed
  .sirQuiet(suppressMessages(nlmixr2est::setCov(fit1, "sir", control = .smallSir())))
  expect_identical(.Random.seed, before)
  .sirQuiet(suppressMessages(nlmixr2est::setCov(fit2, "sir", control = .smallSir())))
  expect_equal(fit1$cov, fit2$cov)
})

test_that("seedCov picks the seed by name or matrix", {
  skip_on_cran()
  fit <- .freshTheoFit()
  seen <- NULL
  local_mocked_bindings(.sirRunCore = function(...) {
    seen <<- list(...)
    stop("captured")
  })
  expect_error(
    nlmixr2est::setCov(fit, "sir", control = .smallSir(seedCov = "nonesuch")),
    "no covariance"
  )
  m <- fit$cov * 2
  expect_error(
    nlmixr2est::setCov(fit, "sir", control = .smallSir(seedCov = m)),
    "captured"
  )
  expect_identical(seen$seedCov, m)
  expect_false(seen$parFixedSe)
  # the installed covariance by name is the same seed as the default
  expect_error(
    nlmixr2est::setCov(fit, "sir", control = .smallSir(seedCov = fit$covMethod)),
    "captured"
  )
  expect_identical(seen$seedCov, fit$cov)
})

test_that("a fit without a covariance is seeded from the RSE default", {
  skip_on_cran()
  fit <- .freshTheoFit(covMethod = "")
  expect_null(fit$cov)
  # .sirQuiet() muffles messages, so they are collected inside it
  msgs <- character()
  .sirQuiet(withCallingHandlers(
    nlmixr2est::setCov(fit, "sir", control = .smallSir()),
    message = function(m) msgs <<- c(msgs, conditionMessage(m))
  ))
  expect_true(any(grepl("no covariance to seed", msgs)))
  expect_identical(fit$covMethod, "sir")
  expect_identical(fit$env$covOptions$sir$seedMethod, "rse")
  expect_equal(fit$env$covOptions$sir$rseTheta, 30)
  expect_identical(attr(fit$sir, "proposalSource"), "rse")

  fit2 <- .freshTheoFit(covMethod = "")
  expect_error(
    nlmixr2est::setCov(fit2, "sir", control = .smallSir(thetaInflation = 2)),
    "inflation"
  )
})

test_that("setCov(fit) <- runSIR() installs a finished run on its own fit", {
  skip_on_cran()
  fit <- .freshTheoFit()
  method0 <- fit$covMethod
  set.seed(3)
  sir <- .sirQuiet(suppressMessages(runSIR(
    fit,
    nSamples = 16L,
    nResample = 8L,
    control = runSIRControl(saveFiles = FALSE, workers = 1L, objfStencil = FALSE)
  )))
  # runSIR() registered it, with its options and the result itself
  expect_true("sir" %in% names(fit$covList))
  expect_identical(fit$env$covOptions$sir$source, "runSIR")
  expect_identical(fit$env$covOptions$sir$seedMethod, method0)
  expect_s3_class(fit$sir, "nlmixr2SIR")

  nlmixr2est::setCov(fit) <- sir
  expect_identical(fit$covMethod, "sir")
  expect_identical(fit$sir, sir)
  expect_identical(fit$env$covOptions$sir$source, "runSIR")
  expect_true(method0 %in% names(fit$covList))

  # a default setCov(fit, "sir") is a different computation, so it recomputes
  runs <- .countRuns()
  .sirQuiet(suppressMessages(nlmixr2est::setCov(fit, "sir", control = .smallSir())))
  expect_equal(runs$n, 1L)
  # ... seeded from what the runSIR() result started from
  expect_identical(runs$args[[1]]$seedCov, fit$env$covOptions$sir$seedCov)

  # runSIR() while "sir" is installed keeps the original seed, not SIR's own
  seed0 <- fit$env$covOptions$sir$seedCov
  set.seed(4)
  .sirQuiet(suppressMessages(runSIR(
    fit,
    nSamples = 16L,
    nResample = 8L,
    control = runSIRControl(saveFiles = FALSE, workers = 1L, objfStencil = FALSE)
  )))
  expect_identical(fit$covMethod, "sir")
  expect_identical(fit$env$covOptions$sir$seedMethod, method0)
  expect_identical(fit$env$covOptions$sir$seedCov, seed0)

  other <- threeEtaFit()
  expect_error(nlmixr2est::setCov(other) <- sir, "not run on")
})
