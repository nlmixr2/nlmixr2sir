# B4: the fingerprint must identify the effective proposal, and fail closed
#
# The fingerprint covered the model, data, parameter names, estimates,
# objective, method, schedule and controls. It did not cover the numbers the
# run actually proposes from. For a path-based input the control digest stored
# the *path*, so replacing the file at that path left the fingerprint
# unchanged and a stale result could be returned for different scientific
# input.
#
# It also skipped any field it could not digest, which is a fail-open policy:
# recovery proceeded precisely when identity could not be established.

test_that("the fingerprint covers the resolved proposal contents", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  ctl <- runSIRControl(objfStencil = FALSE, workers = 1L)
  initial <- .sirResolveInitialProposal(fit, ps, ctl)

  fp <- .sirRunFingerprint(fit, ps, .sirSchedule(16L, 8L), ctl, initial = initial)
  expect_true("proposal" %in% names(fp))
  expect_false(is.na(fp$proposal))

  # A different proposal covariance is a different run, even with everything
  # else identical.
  moved <- initial
  moved$covMat <- moved$covMat * 2
  fp2 <- .sirRunFingerprint(fit, ps, .sirSchedule(16L, 8L), ctl, initial = moved)
  expect_false(identical(fp$proposal, fp2$proposal))
  expect_true("proposal" %in% .sirCompareFingerprints(fp, fp2))
})

test_that("the fingerprint covers the parameter schema, not just names", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  ctl <- runSIRControl(objfStencil = FALSE, workers = 1L)
  fp <- .sirRunFingerprint(fit, ps, .sirSchedule(16L, 8L), ctl)
  expect_true("paramSchema" %in% names(fp))

  # Same names, different bounds: a different estimation problem. Use a
  # parameter whose bound is finite -- most default to -Inf, where subtracting
  # one changes nothing at all.
  ps2 <- ps
  ps2$lower[1L] <- -5
  ps2$upper[1L] <- 5
  expect_false(identical(ps$lower, ps2$lower))
  fp2 <- .sirRunFingerprint(fit, ps2, .sirSchedule(16L, 8L), ctl)
  expect_true("paramSchema" %in% .sirCompareFingerprints(fp, fp2))
})

test_that("changing a covariance file's contents changes the fingerprint", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  p <- nrow(ps)
  path <- withr::local_tempfile(fileext = ".csv")

  # .sirReadCovFile() reads a whitespace-delimited square table with a NAME
  # column, which is the shape a NONMEM .cov file has.
  writeMat <- function(m) {
    df <- data.frame(NAME = ps$sirName, stringsAsFactors = FALSE)
    block <- as.data.frame(m)
    colnames(block) <- ps$sirName
    utils::write.table(
      cbind(df, block), path,
      row.names = FALSE, quote = FALSE, sep = "	"
    )
  }

  writeMat(diag(p) * 0.01)
  ctl <- runSIRControl(objfStencil = FALSE, covmatInput = path, workers = 1L)
  fpA <- .sirRunFingerprint(
    fit, ps, .sirSchedule(16L, 8L), ctl,
    initial = .sirResolveInitialProposal(fit, ps, ctl)
  )

  # Same path, different numbers.
  writeMat(diag(p) * 0.05)
  fpB <- .sirRunFingerprint(
    fit, ps, .sirSchedule(16L, 8L), ctl,
    initial = .sirResolveInitialProposal(fit, ps, ctl)
  )

  expect_true("proposal" %in% .sirCompareFingerprints(fpA, fpB))
})

test_that("recovery fails closed when an identity field cannot be verified", {
  a <- list(stateVersion = 2L, algoVersion = 1L, model = NA_character_, params = "x")
  b <- list(stateVersion = 2L, algoVersion = 1L, model = "abc", params = "x")

  # A fresh run may proceed without being able to digest everything.
  open <- .sirCompareFingerprints(a, b)
  expect_false("model" %in% open)

  # Recovery may not: it is exactly the case where identity matters.
  closed <- .sirCompareFingerprints(a, b, failClosed = TRUE)
  expect_true("model" %in% closed)
  expect_true("model" %in% attr(closed, "unverifiable"))
})

test_that("the algorithm version is an identity field", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  fp <- .sirRunFingerprint(fit, ps, .sirSchedule(16L, 8L), runSIRControl(objfStencil = FALSE, workers = 1L))
  expect_true("algoVersion" %in% names(fp))

  stale <- fp
  stale$algoVersion <- fp$algoVersion - 1L
  expect_true("algoVersion" %in% .sirCompareFingerprints(stale, fp))
})

test_that("recovery refuses a mutated fit covariance", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()
  set.seed(81)
  .sirQuiet(runSIR(
    fit, nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))

  # The proposal the run would now build is different, so the saved result is
  # not an answer to this question.
  st <- nlmixr2utils::readRunState(dir, .sirStateSchema())
  st$fingerprint$proposal <- "a-different-digest"
  nlmixr2utils::writeRunState(dir, st, .sirStateSchema())

  err <- tryCatch(
    .sirQuiet(runSIR(
      fit, nSamples = 16L, nResample = 8L, directory = dir,
      control = runSIRControl(objfStencil = FALSE, recover = TRUE, workers = 1L)
    )),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "proposal")
})

test_that("recovery refuses state whose identity cannot be established", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()
  set.seed(82)
  .sirQuiet(runSIR(
    fit, nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))

  st <- nlmixr2utils::readRunState(dir, .sirStateSchema())
  st$fingerprint$model <- NA_character_
  nlmixr2utils::writeRunState(dir, st, .sirStateSchema())

  err <- tryCatch(
    .sirQuiet(runSIR(
      fit, nSamples = 16L, nResample = 8L, directory = dir,
      control = runSIRControl(objfStencil = FALSE, recover = TRUE, workers = 1L)
    )),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "verif|model")
})

test_that("a dependency version change warns but does not block recovery", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()
  set.seed(83)
  .sirQuiet(runSIR(
    fit, nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))

  st <- nlmixr2utils::readRunState(dir, .sirStateSchema())
  st$fingerprint$pkgVersions <- "0.0/0.0/0.0/0.0"
  nlmixr2utils::writeRunState(dir, st, .sirStateSchema())

  expect_warning(
    .sirQuiet(runSIR(
      fit, nSamples = 16L, nResample = 8L, directory = dir,
      control = runSIRControl(objfStencil = FALSE, recover = TRUE, workers = 1L)
    )),
    "version"
  )
})
