skip_if_not(exists("rxEventEmit", envir = asNamespace("rxode2"), inherits = FALSE),
            "rxode2 has no event bus")

test_that("runSIR emits one fitResult with a small summary and nothing else", {
  skip_on_cran()
  fit <- rxode2::rxEventScope(theoFit())
  tmp <- tempfile("sir_ev_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  rec <- new.env()
  rec$ev <- list()
  rxode2::rxEventListen("nlmixr2sir-test", function(event, ...) {
    rec$ev[[length(rec$ev) + 1L]] <- list(event = event, p = list(...))
  })
  on.exit(rxode2::rxEventUnlisten("nlmixr2sir-test"), add = TRUE)
  set.seed(20261009)
  res <- .sirQuiet(runSIR(
    fit,
    nSamples = c(16L, 16L),
    nResample = c(8L, 8L),
    directory = tmp,
    fitName = "demo",
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))
  expect_identical(vapply(rec$ev, function(e) e$event, ""), "fitResult")
  p <- rec$ev[[1]]$p
  expect_identical(p$kind, "sir")
  expect_identical(p$fit, fit)
  expect_s3_class(p$result, "nlmixr2SIR")
  expect_null(attr(p$result, "resampledMat"))
  expect_null(attr(p$result, "rawResults"))
  expect_identical(attr(p$result, "outputDir"), attr(res, "outputDir"))
  expect_true(object.size(p$result) < 1e6)
  ## the user's result is untouched
  expect_false(is.null(attr(res, "resampledMat")))
  expect_identical(rxode2::rxEventDepth(), 0L)
})

test_that("a runSIR error emits nothing and restores the depth", {
  rec <- new.env()
  rec$n <- 0L
  rxode2::rxEventListen("nlmixr2sir-test", function(event, ...) rec$n <- rec$n + 1L)
  on.exit(rxode2::rxEventUnlisten("nlmixr2sir-test"), add = TRUE)
  expect_error(runSIR(structure(list(), class = c("nlmixr2FitCore", "list")), bogus = 1))
  expect_identical(rec$n, 0L)
  expect_identical(rxode2::rxEventDepth(), 0L)
})
