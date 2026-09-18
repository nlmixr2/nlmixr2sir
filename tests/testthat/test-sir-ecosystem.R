# Integration with the rest of the nlmixr2 ecosystem: registering the SIR
# covariance for setCov(), and deparsing the control object back into source.

test_that("rxUiDeparse round-trips the control object", {
  skip_on_cran()
  ctl <- runSIRControl(objfStencil = FALSE, thetaInflation = 2, workers = 4)
  code <- rxode2::rxUiDeparse(ctl, "x")
  expect_type(code, "language")
  x <- eval(code)
  expect_identical(x, ctl)
})

test_that("rxUiDeparse emits only the arguments that differ from defaults", {
  skip_on_cran()
  txt <- deparse(rxode2::rxUiDeparse(runSIRControl(objfStencil = FALSE, workers = 2), "x"))
  expect_match(txt, "workers")
  expect_false(grepl("capCorrelation", txt))
  # all defaults deparse to a bare constructor call
  expect_equal(
    deparse(rxode2::rxUiDeparse(runSIRControl(objfStencil = FALSE), "y")),
    "y <- runSIRControl(objfStencil = FALSE)"
  )
})

test_that("the deparsed constructor name is one that actually exists", {
  skip_on_cran()
  # .deparseFinal() uses class(object) as the function name, so the class must
  # be the constructor's name. This is what the rename to runSIRControl was
  # for; a package-prefixed class would emit a call to a missing function.
  expect_equal(class(runSIRControl(objfStencil = FALSE)), "runSIRControl")
  expect_true(is.function(get("runSIRControl")))
})

test_that(".sirCovAsFitCov maps SIR names onto fit$cov names and order", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  # A well-conditioned covariance in SIR names, deliberately in a different
  # order from fit$cov so the reordering is exercised.
  nms <- rev(ps$sirName)
  covSir <- diag(seq_along(nms) / 100, nrow = length(nms))
  dimnames(covSir) <- list(nms, nms)

  out <- .sirCovAsFitCov(fit, covSir, ps)
  expect_false(is.null(out))
  expect_identical(rownames(out), rownames(fit$cov))
  expect_identical(colnames(out), colnames(fit$cov))
  # om.eta.ka in fit$cov is eta.ka in SIR names; the value must follow the name
  expect_equal(
    out["om.eta.ka", "om.eta.ka"],
    covSir["eta.ka", "eta.ka"]
  )
})

test_that(".sirCovAsFitCov declines a rank-deficient covariance", {
  skip_on_cran()
  # Registering a singular covariance would give setCov() an unusable matrix.
  #
  # This used to lean on iter1() being rank deficient -- it resampled 4 vectors
  # for 5 parameters. That is now rejected at source, so the singular matrix is
  # built here explicitly rather than borrowed from a broken fixture.
  fit <- theoFit()
  s <- sirSummary(iter1()$resampledMat, fit)
  nm <- colnames(attr(s, "covMatrix"))
  # Genuinely singular and symmetric: the covariance of fewer vectors than
  # there are parameters has rank at most nrow - 1, so chol() must fail.
  set.seed(5)
  cm <- stats::cov(matrix(
    stats::rnorm(3L * length(nm)),
    nrow = 3L,
    dimnames = list(NULL, nm)
  ))
  expect_true(inherits(try(chol(cm), silent = TRUE), "try-error"))

  expect_message(
    out <- .sirCovAsFitCov(fit, cm, .sirParamSpace(fit)),
    "not positive definite"
  )
  expect_null(out)
})

test_that(".sirCovAsFitCov names the covariance when there is no fit$cov", {
  skip_on_cran()
  # setCov(fit, "sir") on a fit without a covariance step installs the SIR
  # covariance under nlmixr2est's full-shape names.
  s <- sirSummary(iter1()$resampledMat, theoFit())
  out <- .sirCovAsFitCov(theoFitNoCov(), attr(s, "covMatrix"))
  expect_identical(rownames(out), c("tka", "tcl", "tv", "add.sd", "om.eta.ka"))
  expect_equal(out["om.eta.ka", "om.eta.ka"], attr(s, "covMatrix")["eta.ka", "eta.ka"])
})

test_that(".sirCovAsFitCov keeps fit$cov's order, then the parameters it lacks", {
  skip_on_cran()
  s <- sirSummary(iter1()$resampledMat, theoFit())
  fit <- theoFitThetaCov()
  out <- .sirCovAsFitCov(fit, attr(s, "covMatrix"))
  expect_identical(rownames(out)[seq_len(nrow(fit$cov))], rownames(fit$cov))
  expect_true("om.eta.ka" %in% rownames(out))
})

test_that(".sirCovAsFitCov declines on a parameter mismatch", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  partial <- attr(sirSummary(iter1()$resampledMat, fit), "covMatrix")
  partial <- partial[1:2, 1:2, drop = FALSE]
  expect_message(
    out <- .sirCovAsFitCov(fit, partial, ps),
    "do not match"
  )
  expect_null(out)
})

test_that(".sirRegisterCovList merges rather than replaces", {
  skip_on_cran()
  env <- new.env()
  fake <- list(env = env)
  m <- matrix(1, 1L, 1L, dimnames = list("a", "a"))
  .sirRegisterCovList(fake, "one", m)
  .sirRegisterCovList(fake, "two", m)
  expect_setequal(names(get("covList", envir = env)), c("one", "two"))
})

test_that(".sirRegisterCovList ignores a NULL covariance", {
  env <- new.env()
  expect_false(.sirRegisterCovList(list(env = env), "sir", NULL))
  expect_false(exists("covList", envir = env, inherits = FALSE))
})

test_that("runSIR registers a covariance that setCov() can select", {
  skip_on_cran()
  # A fit built here rather than the shared fixture, since registration
  # mutates fit$env and would leak into other tests.
  fit <- suppressMessages(nlmixr2utils::nlmixr2(
    theoOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L, covMethod = "r")
  ))
  tmp <- tempfile("sir_setcov_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  set.seed(3)
  .sirQuiet(runSIR(
    fit,
    nSamples = 16L,
    nResample = 8L,
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))

  covList <- get("covList", envir = fit$env)
  expect_true("sir" %in% names(covList))
  expect_identical(rownames(covList$sir), rownames(fit$cov))
  expect_false(inherits(try(chol(covList$sir), silent = TRUE), "try-error"))
  # the result and the options it was computed with are kept with it
  expect_s3_class(fit$sir, "nlmixr2SIR")
  expect_identical(fit$env$covOptions$sir$source, "runSIR")
})
