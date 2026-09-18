# The schedule identity of an extended run.
#
# addIterations retains the previous iterations and their summaries, so the
# result describes every iteration that ran. Its identity has to describe the
# same thing. Saving only the extension schedule means a later recover() call
# presenting that short schedule matches, and receives the long result.

test_that("extending a run yields a cumulative schedule everywhere", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()

  set.seed(42)
  first <- .sirQuiet(runSIR(
    fit,
    nSamples = c(16L, 16L),
    nResample = c(8L, 8L),
    directory = dir,
    control = runSIRControl(objfStencil = FALSE, workers = 1L)
  ))
  expect_equal(nrow(attr(first, "schedule")), 2L)

  set.seed(42)
  second <- .sirQuiet(runSIR(
    fit,
    nSamples = 16L,
    nResample = 8L,
    directory = dir,
    control = runSIRControl(objfStencil = FALSE, workers = 1L, addIterations = TRUE)
  ))

  sched <- attr(second, "schedule")
  expect_equal(nrow(sched), 3L)
  expect_equal(sched$iter, 1:3)
  expect_equal(as.integer(sched$nSamples), c(16L, 16L, 16L))
  expect_equal(as.integer(sched$nResample), c(8L, 8L, 8L))

  # The identity the next run will be compared against must agree with the
  # result it describes.
  expect_equal(attr(second, "fingerprint")$schedule, "16,16,16/8,8,8")

  state <- nlmixr2utils::readRunState(dir, nlmixr2sir:::.sirStateSchema())
  expect_equal(state$completedIterations, 3L)
  expect_equal(state$fingerprint$schedule, "16,16,16/8,8,8")
  expect_equal(nrow(state$iterationSummary), 3L)
})

test_that("a recover request with the extension schedule alone is refused", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()

  set.seed(42)
  .sirQuiet(runSIR(
    fit,
    nSamples = c(16L, 16L),
    nResample = c(8L, 8L),
    directory = dir,
    control = runSIRControl(objfStencil = FALSE, workers = 1L)
  ))
  set.seed(42)
  .sirQuiet(runSIR(
    fit,
    nSamples = 16L,
    nResample = 8L,
    directory = dir,
    control = runSIRControl(objfStencil = FALSE, workers = 1L, addIterations = TRUE)
  ))

  # Before the fix the saved identity described only the one-iteration
  # extension, so this matched and returned the three-iteration result as
  # though it were the answer to a one-iteration request.
  expect_error(
    .sirQuiet(runSIR(
      fit,
      nSamples = 16L,
      nResample = 8L,
      directory = dir,
      control = runSIRControl(objfStencil = FALSE, workers = 1L)
    )),
    "does not match this recovery request"
  )
})
