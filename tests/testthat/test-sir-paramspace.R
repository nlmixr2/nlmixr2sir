# .sirParamSpace is the single source of truth for the SIR parameter vector.
# These tests pin the three naming conventions it bridges, and the fit$cov
# shape assumption that broke when nlmixr2est switched to covFull = TRUE.

test_that(".sirParamSpace returns the documented columns", {
  skip_on_cran()
  expect_named(
    .sirParamSpace(theoFit()),
    c(
      "sirName",
      "covName",
      "rawName",
      "kind",
      "ntheta",
      "neta1",
      "neta2",
      "est",
      "lower",
      "upper",
      "fullCovName"
    )
  )
})

test_that(".sirParamSpace classifies residual error as sigma, not theta", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFit())
  expect_equal(ps$kind[ps$sirName == "add.sd"], "sigma")
  expect_equal(ps$kind[ps$sirName == "tka"], "theta")
})

test_that(".sirParamSpace covers every row of fit$cov, in the same order", {
  skip_on_cran()
  for (fit in list(theoFit(), threeEtaFit(), blockFit())) {
    expect_identical(.sirParamSpace(fit)$covName, rownames(fit$cov))
  }
})

test_that(".sirParamSpace covers every raw-results schema parameter column", {
  skip_on_cran()
  for (fit in list(theoFit(), threeEtaFit(), blockFit())) {
    schema <- nlmixr2utils::rawResultsSchema(fit)
    expect_setequal(
      .sirParamSpace(fit)$rawName,
      c(schema$thetaCols, schema$omegaCols, schema$sigmaCols)
    )
  }
})

test_that(".sirParamSpace names OMEGA off-diagonals in neta1, neta2 order", {
  skip_on_cran()
  ps <- .sirParamSpace(blockFit())
  off <- ps[ps$kind == "omegaOffdiag", ]
  expect_equal(nrow(off), 1L)
  expect_equal(off$sirName, "eta.cl:eta.ka")
  expect_equal(off$covName, "cov.eta.cl.eta.ka")
  expect_equal(off$rawName, "omega(eta.cl,eta.ka)")
  expect_equal(off$neta1, 2L)
  expect_equal(off$neta2, 1L)
})

test_that(".sirParamSpace bounds variances below at zero but not covariances", {
  skip_on_cran()
  ps <- .sirParamSpace(blockFit())
  expect_equal(ps$lower[ps$kind == "omegaDiag"], c(0, 0))
  expect_equal(ps$lower[ps$kind == "omegaOffdiag"], -Inf)
  expect_equal(ps$lower[ps$sirName == "add.sd"], 0)
})

test_that(".sirParamSpace OMEGA estimates agree with fit$omega", {
  skip_on_cran()
  ps <- .sirParamSpace(blockFit())
  om <- ps[ps$kind %in% c("omegaDiag", "omegaOffdiag"), ]
  expect_equal(
    om$est,
    blockFit()$omega[cbind(om$neta1, om$neta2)],
    tolerance = 1e-12
  )
})

test_that(".sirParamSpace marks covName NA when there is no covariance step", {
  skip_on_cran()
  ps <- .sirParamSpace(theoFitNoCov())
  expect_true(all(is.na(ps$covName)))
  expect_false(.sirCovHasOmega(ps))
})

test_that(".sirCovHasOmega is TRUE for a covFull fit", {
  skip_on_cran()
  expect_true(.sirCovHasOmega(.sirParamSpace(theoFit())))
  expect_true(.sirCovHasOmega(.sirParamSpace(blockFit())))
})

# Regression test for the P0 breakage: fit$theta[rownames(fit$cov)] used to
# build the proposal mean, and silently produced NA for every OMEGA row once
# nlmixr2est started reporting OMEGA in fit$cov. If the core changes shape
# again, this is where it should fail.
test_that(".sirProposalMu is complete and finite for every fit shape", {
  skip_on_cran()
  for (fit in list(theoFit(), threeEtaFit(), blockFit(), theoFitNoCov())) {
    mu <- .sirProposalMu(fit)
    expect_false(anyNA(mu))
    expect_true(all(is.finite(mu)))
    expect_named(mu, .sirParamSpace(fit)$sirName)
  }
})

test_that("fit$cov carries the OMEGA block under nlmixr2est defaults", {
  skip_on_cran()
  expect_true(any(grepl("^om\\.", rownames(theoFit()$cov))))
})
