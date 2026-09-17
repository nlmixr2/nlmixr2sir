# Release-readiness: a real runSIR() call, the complete artifact set it should
# leave on disk, and the recover / addIterations round trip. These are the
# paths a user actually exercises, and the ones a unit test of an internal
# helper cannot vouch for.

test_that("runSIR writes the complete artifact set", {
  skip_on_cran()
  tmp <- tempfile("sir_e2e_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  set.seed(20260913)
  res <- .sirQuiet(runSIR(
    theoFit(),
    nSamples = c(16L, 16L),
    nResample = c(8L, 8L),
    directory = tmp,
    fitName = "demo",
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))

  expect_setequal(
    basename(list.files(tmp)),
    c(
      "demo_sir.cov",
      "demo_sir.sdcorr",
      "raw_results.csv",
      "raw_results.rds",
      "raw_results_header.json",
      "sample_rejection_summary.txt",
      "sir_manifest.dcf",
      "sir_results.csv",
      "sir_seed.rds",
      "sir_state.rds",
      "summary_iterations.csv"
    )
  )
  expect_s3_class(res, "nlmixr2SIR")
  expect_equal(nrow(attr(res, "iterationSummary")), 2L)
})

test_that("the written summary round-trips back through read.csv", {
  skip_on_cran()
  tmp <- tempfile("sir_e2e_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  set.seed(1)
  res <- .sirQuiet(runSIR(
    theoFit(),
    nSamples = 16L,
    nResample = 8L,
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))
  back <- utils::read.csv(
    file.path(tmp, "sir_results.csv"),
    check.names = FALSE
  )
  expect_equal(back$param, res$param)
  expect_equal(back$estimate, res$estimate, tolerance = 1e-8)
})

test_that("recover returns the completed run instead of repeating it", {
  skip_on_cran()
  tmp <- tempfile("sir_recover_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  set.seed(42)
  first <- .sirQuiet(runSIR(
    theoFit(),
    nSamples = c(16L, 16L),
    nResample = c(8L, 8L),
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))
  expect_true(file.exists(file.path(tmp, "sir_state.rds")))

  # Same schedule, same directory, recover = TRUE: the stored result comes
  # back rather than the model being evaluated again.
  set.seed(999)
  again <- .sirQuiet(runSIR(
    theoFit(),
    nSamples = c(16L, 16L),
    nResample = c(8L, 8L),
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, recover = TRUE, workers = 1L)
  ))
  expect_equal(again$estimate, first$estimate, tolerance = 1e-12)
  expect_equal(again$sd, first$sd, tolerance = 1e-12)
  expect_equal(nrow(attr(again, "iterationSummary")), 2L)
})

test_that("addIterations extends a completed run", {
  skip_on_cran()
  tmp <- tempfile("sir_additer_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  set.seed(7)
  first <- .sirQuiet(runSIR(
    theoFit(),
    nSamples = c(16L, 16L),
    nResample = c(8L, 8L),
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))
  expect_equal(nrow(attr(first, "iterationSummary")), 2L)

  extended <- .sirQuiet(runSIR(
    theoFit(),
    nSamples = 16L,
    nResample = 8L,
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, addIterations = TRUE, workers = 1L)
  ))
  summary_df <- attr(extended, "iterationSummary")
  expect_equal(nrow(summary_df), 3L)
  expect_equal(summary_df$iter, 1:3)
  # The first two iterations are carried over unchanged, not recomputed.
  expect_equal(
    summary_df$nSuccessful[1:2],
    attr(first, "iterationSummary")$nSuccessful
  )
})

test_that("every diagnostic plot builds from a real run", {
  skip_on_cran()
  tmp <- tempfile("sir_plots_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  set.seed(5)
  res <- .sirQuiet(runSIR(
    theoFit(),
    nSamples = c(20L, 20L),
    nResample = c(10L, 10L),
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))
  for (ty in c(
    "parameters",
    "dofv",
    "resampling",
    "convergence",
    "intervals",
    "rsecor"
  )) {
    p <- suppressWarnings(plot(res, type = ty, nReplicate = 25L))
    expect_s3_class(p, "ggplot")
  }
})

test_that("print() works on a real run", {
  skip_on_cran()
  tmp <- tempfile("sir_print_")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  set.seed(6)
  res <- .sirQuiet(runSIR(
    theoFit(),
    nSamples = 16L,
    nResample = 8L,
    directory = tmp,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))
  # cli headers go to the message stream; capture.output() sees the tables.
  out <- utils::capture.output(ret <- print(res))
  expect_identical(ret, res)
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "param")
  expect_match(txt, "estimate")
  expect_match(txt, "nAttempted")
})
