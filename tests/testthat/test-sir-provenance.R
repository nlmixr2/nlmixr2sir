# B4: run identity, directory ownership, and optional file saving -------------
#
# `recover = TRUE` is the default, and a saved run used to be returned solely
# because a state file existed in the directory and recorded enough completed
# iterations. Nothing checked that the state belonged to the fit in hand, so
# pointing a different model, dataset, or schedule at the same directory
# returned a stale result labelled as belonging to the new run. That is a
# scientific-integrity failure, not merely a caching bug.
#
# Separately, an explicitly supplied directory in overwrite mode was removed
# with unlink(recursive = TRUE) with no check that nlmixr2sir had created it.

test_that("the run fingerprint is deterministic for identical inputs", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  ctl <- runSIRControl(objfStencil = FALSE, workers = 1L)
  a <- .sirRunFingerprint(fit, ps, .sirSchedule(16L, 8L), ctl)
  b <- .sirRunFingerprint(fit, ps, .sirSchedule(16L, 8L), ctl)
  expect_identical(a, b)
  expect_true(is.integer(a$stateVersion) || is.numeric(a$stateVersion))
})

test_that("the fingerprint changes when the schedule or controls change", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  base <- .sirRunFingerprint(fit, ps, .sirSchedule(16L, 8L), runSIRControl(objfStencil = FALSE, workers = 1L))

  sched <- .sirRunFingerprint(fit, ps, .sirSchedule(24L, 8L), runSIRControl(objfStencil = FALSE, workers = 1L))
  expect_false(identical(base$schedule, sched$schedule))

  ctl <- .sirRunFingerprint(
    fit, ps, .sirSchedule(16L, 8L), runSIRControl(objfStencil = FALSE, workers = 1L, boxcox = FALSE)
  )
  expect_false(identical(base$controls, ctl$controls))

  # Parallelism does not change the answer, so it must not invalidate a run.
  par <- .sirRunFingerprint(
    fit, ps, .sirSchedule(16L, 8L), runSIRControl(objfStencil = FALSE, workers = 1L, rxThreads = 2L)
  )
  expect_identical(base$controls, par$controls)
})

test_that("the fingerprint distinguishes different models and data", {
  skip_on_cran()
  a <- .sirRunFingerprint(
    theoFit(), .sirParamSpace(theoFit()), .sirSchedule(16L, 8L), runSIRControl(objfStencil = FALSE, workers = 1L)
  )
  b <- .sirRunFingerprint(
    blockFit(), .sirParamSpace(blockFit()), .sirSchedule(16L, 8L), runSIRControl(objfStencil = FALSE, workers = 1L)
  )
  expect_false(identical(a$model, b$model))
  expect_false(identical(a$params, b$params))
})

test_that(".sirCompareFingerprints reports mismatched fields by name", {
  skip_on_cran()
  fit <- theoFit()
  ps <- .sirParamSpace(fit)
  a <- .sirRunFingerprint(fit, ps, .sirSchedule(16L, 8L), runSIRControl(objfStencil = FALSE, workers = 1L))
  b <- .sirRunFingerprint(fit, ps, .sirSchedule(24L, 8L), runSIRControl(objfStencil = FALSE, workers = 1L, boxcox = FALSE))

  bad <- .sirCompareFingerprints(a, b)
  expect_true("schedule" %in% bad)
  expect_true("controls" %in% bad)
  expect_false("model" %in% bad)

  expect_length(.sirCompareFingerprints(a, a), 0L)
  # Fields can be exempted, which is how addIterations tolerates a new schedule.
  expect_false("schedule" %in% .sirCompareFingerprints(a, b, ignore = "schedule"))
})

test_that("recovery refuses a state file from a different run", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()
  set.seed(31)
  .sirQuiet(runSIR(
    fit, nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))

  # Same directory, different schedule: the stored result does not belong to
  # this request and must not be handed back.
  err <- tryCatch(
    .sirQuiet(runSIR(
      fit, nSamples = c(16L, 16L), nResample = c(8L, 8L), directory = dir,
      control = runSIRControl(objfStencil = FALSE, recover = TRUE, workers = 1L)
    )),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "schedule")
})

test_that("recovery refuses a state file from a different model", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  set.seed(32)
  .sirQuiet(runSIR(
    theoFit(), nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))

  err <- tryCatch(
    .sirQuiet(runSIR(
      blockFit(), nSamples = 16L, nResample = 8L, directory = dir,
      control = runSIRControl(objfStencil = FALSE, recover = TRUE, workers = 1L)
    )),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "model|params|does not match")
})

test_that("an unowned non-empty directory is never deleted", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  sentinel <- file.path(dir, "important.txt")
  writeLines("do not delete me", sentinel)

  err <- tryCatch(
    .sirQuiet(runSIR(
      theoFit(), nSamples = 16L, nResample = 8L, directory = dir,
      control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
    )),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "not appear to be|nlmixr2sir")
  expect_true(file.exists(sentinel))
  expect_equal(readLines(sentinel), "do not delete me")
})

test_that("a directory nlmixr2sir created carries an ownership marker", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  set.seed(33)
  .sirQuiet(runSIR(
    theoFit(), nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))
  expect_true(.sirDirIsOwned(dir))

  # A fresh run over its own directory is fine.
  expect_no_error(.sirQuiet(runSIR(
    theoFit(), nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  )))
})

test_that("saveFiles = FALSE writes nothing at all", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  before <- list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  expect_length(before, 0L)

  set.seed(34)
  res <- .sirQuiet(runSIR(
    theoFit(), nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(objfStencil = FALSE, saveFiles = FALSE, workers = 1L)
  ))

  expect_s3_class(res, "nlmixr2SIR")
  expect_null(attr(res, "outputDir"))
  after <- list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  expect_length(after, 0L)
})

test_that("saveFiles = FALSE still produces a usable result and diagnostics", {
  skip_on_cran()
  set.seed(35)
  res <- .sirQuiet(runSIR(
    theoFit(), nSamples = 16L, nResample = 8L,
    control = runSIRControl(objfStencil = FALSE, saveFiles = FALSE, workers = 1L)
  ))
  expect_true(all(c("param", "estimate", "mean", "sd") %in% names(res)))
  expect_false(is.null(attr(res, "iterations")))
  expect_no_error(print(res))
})

test_that("saveFiles = FALSE cannot be combined with addIterations", {
  skip_on_cran()
  expect_error(
    .sirQuiet(runSIR(
      theoFit(), nSamples = 16L, nResample = 8L,
      control = runSIRControl(objfStencil = FALSE, saveFiles = FALSE, addIterations = TRUE, workers = 1L)
    )),
    "saveFiles"
  )
})

test_that("a state file from an older schema version is refused", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()
  set.seed(36)
  .sirQuiet(runSIR(
    fit, nSamples = 16L, nResample = 8L, directory = dir,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))

  st <- nlmixr2utils::readRunState(dir, .sirStateSchema())
  st$fingerprint$stateVersion <- 0L
  nlmixr2utils::writeRunState(dir, st, .sirStateSchema())

  err <- tryCatch(
    .sirQuiet(runSIR(
      fit, nSamples = 16L, nResample = 8L, directory = dir,
      control = runSIRControl(objfStencil = FALSE, recover = TRUE, workers = 1L)
    )),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "stateVersion|version")
})

test_that("initial proposal repair is recorded separately from iteration repair", {
  skip_on_cran()
  fit <- theoFit()
  res <- suppressMessages(suppressWarnings(runSIR(
    fit,
    nSamples = 16L,
    nResample = 8L,
    control = runSIRControl(objfStencil = FALSE, workers = 1L, saveFiles = FALSE)
  )))

  ip <- attr(res, "initialProposalRepair")
  expect_type(ip, "list")
  expect_named(
    ip,
    c("adjusted", "method", "threshold", "magnitude"),
    ignore.order = TRUE
  )
  expect_false(is.null(ip$adjusted))
  expect_true(is.logical(ip$adjusted))
  expect_true(is.finite(ip$magnitude))
  # Distinct from the per-iteration record, which covers later empirical
  # updates rather than the covariance iteration 1 is drawn from.
  expect_true("posDefAdjusted" %in% names(attr(res, "iterationSummary")))
})
