# The dOFV-versus-chi-square convergence diagnostic, and PsN's summary-
# statistics parity. Ported from PsN R-scripts/sir_default.R and
# lib/tool/sir.pm empirical_statistics().

test_that(".sirDofvQuantiles stops short of 1 so the reference stays finite", {
  q <- .sirDofvQuantiles(20L)
  expect_length(q, 19L)
  expect_equal(min(q), 0)
  expect_equal(max(q), 19 / 20)
  expect_true(all(is.finite(stats::qchisq(q, df = 5))))
})

test_that(".sirDofvQuantiles refuses a grid too small to be a curve", {
  expect_error(.sirDofvQuantiles(2L), "At least 3 resamples")
})

test_that(".sirDofvCurves gives proposal and SIR curves per iteration", {
  skip_on_cran()
  cv <- .sirDofvCurves(sirObj())
  expect_setequal(cv$type, c("reference", "proposal", "SIR"))
  # one reference curve, shared across iterations
  expect_equal(sum(cv$type == "reference"), sum(cv$type == "SIR"))
  expect_true(all(is.na(cv$iteration[cv$type == "reference"])))
  expect_false(anyNA(cv$iteration[cv$type != "reference"]))
})

test_that(".sirDofvCurves reference is a chi-square on the parameter count", {
  skip_on_cran()
  cv <- .sirDofvCurves(sirObj())
  ref <- cv[cv$type == "reference", ]
  expect_equal(
    ref$dOFV,
    stats::qchisq(ref$quantile, df = length(unique(sirObj()$param))),
    tolerance = 1e-10
  )
  expect_false(is.unsorted(ref$dOFV))
})

test_that(".sirProposalTooNarrow applies PsN's 25% rule", {
  # Hand-built curves: the proposal sits below the reference at 3 of 4
  # quantiles, which is above the quarter threshold.
  curves <- rbind(
    data.frame(
      iteration = NA_integer_,
      label = "reference",
      type = "reference",
      quantile = c(0.1, 0.2, 0.3, 0.4),
      dOFV = c(1, 2, 3, 4),
      stringsAsFactors = FALSE
    ),
    data.frame(
      iteration = 1L,
      label = "proposal 1",
      type = "proposal",
      quantile = c(0.1, 0.2, 0.3, 0.4),
      dOFV = c(0, 0, 0, 5),
      stringsAsFactors = FALSE
    )
  )
  chk <- .sirProposalTooNarrow(curves)
  expect_equal(chk$fraction, 0.75)
  expect_true(chk$warn)
})

test_that(".sirProposalTooNarrow stays quiet for a wide enough proposal", {
  curves <- rbind(
    data.frame(
      iteration = NA_integer_,
      label = "reference",
      type = "reference",
      quantile = c(0.1, 0.2, 0.3, 0.4),
      dOFV = c(1, 2, 3, 4),
      stringsAsFactors = FALSE
    ),
    data.frame(
      iteration = 1L,
      label = "proposal 1",
      type = "proposal",
      quantile = c(0.1, 0.2, 0.3, 0.4),
      dOFV = c(2, 3, 4, 5),
      stringsAsFactors = FALSE
    )
  )
  chk <- .sirProposalTooNarrow(curves)
  expect_equal(chk$fraction, 0)
  expect_false(chk$warn)
})

test_that("the too-narrow warning names inflation as the remedy", {
  curves <- rbind(
    data.frame(
      iteration = NA_integer_,
      label = "reference",
      type = "reference",
      quantile = c(0.1, 0.2),
      dOFV = c(1, 2),
      stringsAsFactors = FALSE
    ),
    data.frame(
      iteration = 1L,
      label = "proposal 1",
      type = "proposal",
      quantile = c(0.1, 0.2),
      dOFV = c(0, 0),
      stringsAsFactors = FALSE
    )
  )
  expect_warning(
    .sirWarnProposalTooNarrow(.sirProposalTooNarrow(curves)),
    "inflated proposal"
  )
})

test_that(".sirDofvNoise brackets the curve it is built from", {
  skip_on_cran()
  cv <- .sirDofvCurves(sirObj())
  quant <- cv$quantile[cv$type == "reference"]
  set.seed(1)
  nb <- .sirDofvNoise(sirObj(), 1L, quant, nReplicate = 50L)
  expect_equal(nrow(nb), length(quant))
  expect_true(all(nb$low <= nb$high))
})

test_that("plot(type = 'convergence') builds a faceted ggplot", {
  skip_on_cran()
  p <- suppressWarnings(
    plot(sirObj(), type = "convergence", nReplicate = 25L)
  )
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$y, "dOFV")
  expect_equal(p$labels$x, "Quantile")
})

test_that("plot(type = 'convergence') can omit the noise band", {
  skip_on_cran()
  p <- suppressWarnings(plot(sirObj(), type = "convergence", noise = FALSE))
  expect_s3_class(p, "ggplot")
})

# Summary-statistics parity.

test_that(".sirPercentileLabels reproduces PsN's percentile set", {
  # From prediction intervals 0, 40, 80, 90, 95.
  expect_equal(
    .sirPercentileLabels(),
    c(2.5, 5, 10, 30, 50, 70, 90, 95, 97.5)
  )
})

test_that("sirSummary reports both mean and median", {
  skip_on_cran()
  s <- sirSummary(iter1()$resampledMat, theoFit())
  expect_true(all(c("mean", "p50") %in% names(s)))
  expect_equal(
    s$mean,
    unname(colMeans(iter1()$resampledMat)),
    tolerance = 1e-12
  )
})

test_that("rse_sd_scale halves the RSE of variance parameters only", {
  skip_on_cran()
  s <- sirSummary(iter1()$resampledMat, theoFit())
  ps <- .sirParamSpace(theoFit())
  kind <- ps$kind[match(s$param, ps$sirName)]

  isOmega <- kind %in% c("omegaDiag", "omegaOffdiag")
  expect_equal(s$rse_sd_scale[isOmega], s$rse[isOmega] / 2, tolerance = 1e-12)
  # THETA is not a variance, and nlmixr2's residual error is already on the SD
  # scale, so neither is rescaled.
  expect_true(all(is.na(s$rse_sd_scale[!isOmega])))
})

test_that("sirSummary attaches an sd/correlation matrix", {
  skip_on_cran()
  s <- sirSummary(iter1()$resampledMat, theoFit())
  sdcor <- attr(s, "sdCorMatrix")
  cm <- stats::cov(iter1()$resampledMat)
  expect_equal(diag(sdcor), sqrt(diag(cm)), tolerance = 1e-12)
  expect_equal(
    sdcor[lower.tri(sdcor)],
    stats::cov2cor(cm)[lower.tri(cm)],
    tolerance = 1e-12
  )
})

test_that("sirSummary records that rse is a percentage", {
  skip_on_cran()
  # PsN reports the same quantity as a fraction; the units are recorded so the
  # difference is not silent.
  expect_equal(
    attr(sirSummary(iter1()$resampledMat, theoFit()), "rseUnits"),
    "percent"
  )
})

# P2.2 intervals by iteration, P2.3 RSE/correlation, P2.5 on-disk parity.

test_that(".sirIterationIntervals covers both distributions per iteration", {
  skip_on_cran()
  iv <- .sirIterationIntervals(sirObj())
  expect_setequal(iv$type, c("proposal", "SIR"))
  expect_setequal(iv$param, colnames(attr(sirObj(), "resampledMat")))
  expect_true(all(iv$low <= iv$median))
  expect_true(all(iv$median <= iv$high))
})

test_that(".sirIterationIntervals honours the requested interval width", {
  skip_on_cran()
  narrow <- .sirIterationIntervals(sirObj(), ci = 50)
  wide <- .sirIterationIntervals(sirObj(), ci = 95)
  expect_true(all(
    (wide$high - wide$low) >= (narrow$high - narrow$low) - 1e-12
  ))
})

test_that("asymmetry is the ratio of the two half-widths", {
  skip_on_cran()
  iv <- .sirIterationIntervals(sirObj())
  expect_equal(
    iv$asymmetry,
    (iv$high - iv$median) / (iv$median - iv$low),
    tolerance = 1e-12
  )
})

test_that(".sirRseCorData keeps one triangle with RSE on the diagonal", {
  skip_on_cran()
  d <- .sirRseCorData(sirObj())
  n <- length(unique(sirObj()$param))
  expect_equal(nrow(d), n * (n + 1) / 2)
  expect_equal(sum(d$isDiagonal), n)
  # off-diagonal cells are correlations
  expect_true(all(abs(d$value[!d$isDiagonal]) <= 1 + 1e-12))
  # diagonal cells carry an asymmetry band, off-diagonal ones do not
  expect_false(anyNA(d$asymmetryBand[d$isDiagonal]))
  expect_true(all(is.na(d$asymmetryBand[!d$isDiagonal])))
})

test_that(".sirRseCorData can show the proposal instead of the posterior", {
  skip_on_cran()
  expect_equal(
    unique(.sirRseCorData(sirObj(), which = "proposal")$which),
    "proposal"
  )
  expect_equal(unique(.sirRseCorData(sirObj(), which = "SIR")$which), "SIR")
})

test_that("asymmetry bands follow PsN's breakpoints", {
  skip_on_cran()
  d <- .sirRseCorData(sirObj())
  diagRows <- d[d$isDiagonal, ]
  expected <- cut(
    diagRows$asymmetry,
    breaks = c(-Inf, 0.5, 1, 1.25, 2, Inf),
    labels = c("<0.5", "0.5-1", "1-1.25", "1.25-2", ">2"),
    right = FALSE
  )
  expect_equal(diagRows$asymmetryBand, expected)
})

test_that("plot(type = 'intervals') and plot(type = 'rsecor') build", {
  skip_on_cran()
  expect_s3_class(plot(sirObj(), type = "intervals"), "ggplot")
  expect_s3_class(plot(sirObj(), type = "rsecor"), "ggplot")
  expect_s3_class(
    plot(sirObj(), type = "rsecor", which = "proposal"),
    "ggplot"
  )
})

test_that("summary_iterations.csv leads with PsN's column names", {
  skip_on_cran()
  tmp <- tempfile("sir_cols_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  .sirWriteIterationSummary(attr(sirObj(), "iterationSummary"), tmp)
  got <- names(utils::read.csv(
    file.path(tmp, "summary_iterations.csv"),
    check.names = FALSE
  ))
  expect_identical(
    got[seq_len(10)],
    c(
      "iteration",
      "commandline.samples",
      "attempted.samples",
      "successful.samples",
      "commandline.resamples",
      "actual.resamples",
      "requested.ratio",
      "actual.ratio",
      "negative.dOFV",
      "minimum.sample.ofv"
    )
  )
})

test_that(".sirWriteCovMatrices exports the covariance and sd/correlation", {
  skip_on_cran()
  tmp <- tempfile("sir_cov_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  s <- sirSummary(iter1()$resampledMat, theoFit())
  .sirWriteCovMatrices(s, tmp, fitName = "demo")

  expect_true(file.exists(file.path(tmp, "demo_sir.cov")))
  expect_true(file.exists(file.path(tmp, "demo_sir.sdcorr")))
  back <- utils::read.delim(file.path(tmp, "demo_sir.cov"), check.names = FALSE)
  expect_equal(back$NAME, rownames(attr(s, "covMatrix")))
  expect_equal(
    as.matrix(back[, -1L]),
    unname(attr(s, "covMatrix")),
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

# R4: the synthetic reference row must not enter empirical proposal summaries --
#
# .sirBuildRawResults() prepends a row holding the proposal centre itself
# (sample_id 0, dOFV 0, no resamples). It belongs in the written raw-results
# file, because PsN writes it too, but it is not a draw from the proposal and
# must not be counted as one. PsN drops it before summarising for exactly this
# reason (R-scripts/sir_default.R).

test_that(".sirBuildRawResults marks the centre row as a reference row", {
  skip_on_cran()
  raw <- iter1()$rawResults
  expect_true("role" %in% names(raw))
  expect_identical(raw$role[raw$sample_id == 0L], "reference")
  expect_true(all(raw$role[raw$sample_id != 0L] == "sample"))
})

test_that("proposal dOFV quantiles exclude the reference row", {
  skip_on_cran()
  it <- iter1()
  raw <- it$rawResults
  curves <- .sirDofvCurves(sirObj(), quant = seq(0.01, 0.99, by = 0.01))
  got <- curves$dOFV[curves$type == "proposal"]

  # The proposal sample is the distinct candidates, with the centre removed.
  drawn <- raw[raw$role == "sample" & !duplicated(raw$sample_id), , drop = FALSE]
  expected <- unname(stats::quantile(
    drawn$deltaofv[!is.na(drawn$deltaofv)],
    probs = seq(0.01, 0.99, by = 0.01),
    na.rm = TRUE
  ))

  expect_equal(got, expected, tolerance = 1e-12)
  # The centre contributes an exact zero; keeping it shifts the low quantiles.
  expect_false(0 %in% drawn$deltaofv)
})

test_that("proposal intervals and RSEs exclude the reference row", {
  skip_on_cran()
  obj <- sirObj()
  raw <- iter1()$rawResults
  params <- colnames(iter1()$resampledMat)

  intervals <- .sirIterationIntervals(obj, ci = 95)
  prop <- intervals[intervals$type == "proposal", , drop = FALSE]
  drawn <- raw[raw$role == "sample" & !duplicated(raw$sample_id), , drop = FALSE]

  p1 <- params[[1L]]
  expect_equal(
    prop$median[prop$param == p1],
    unname(stats::quantile(drawn[[p1]], probs = 0.5, na.rm = TRUE)),
    tolerance = 1e-12
  )

  # .sirRseCorData() puts RSEs on the diagonal of its grid; those come from the
  # proposal covariance, so they move if the centre row is counted.
  rc <- .sirRseCorData(obj, which = "proposal", ci = 95)
  cm <- stats::cov(as.matrix(drawn[, params, drop = FALSE]))
  estimate <- stats::setNames(obj$estimate, obj$param)[params]
  expected_rse <- 100 * abs(sqrt(diag(cm)) / estimate)

  got <- rc$value[rc$isDiagonal][match(params, as.character(rc$row[rc$isDiagonal]))]
  expect_equal(got, unname(expected_rse), tolerance = 1e-12)
})

# R3: recentring must move the dOFV reference, not just the proposal centre ----
#
# When a candidate beats the fitted optimum, PsN moves both the centre and
# reference_ofv (lib/tool/sir.pm). Moving only the centre leaves every later
# dOFV, negative-dOFV count, and the chi-square convergence curve measured
# against an optimum the run has already superseded.

test_that("sirRunIteration accepts and returns a reference OFV", {
  skip_on_cran()
  fit <- theoFit()
  set.seed(42)
  it <- .sirQuiet(sirRunIteration(
    fit,
    mu = .sirProposalMu(fit),
    proposalCov = sirGetProposalCov(fit),
    nSamples = 16L,
    nResample = 8L,
    iterNum = 1L,
    recenter = TRUE,
    boxcox = FALSE,
    directory = NULL
  ))
  expect_true("newReferenceOfv" %in% names(it))
  expect_true(is.finite(it$newReferenceOfv))
})

test_that("recentring moves the reference OFV to the better optimum", {
  skip_on_cran()
  fit <- theoFit()
  set.seed(42)
  it <- .sirQuiet(sirRunIteration(
    fit,
    mu = .sirProposalMu(fit),
    proposalCov = sirGetProposalCov(fit),
    nSamples = 16L,
    nResample = 8L,
    iterNum = 1L,
    recenter = TRUE,
    boxcox = FALSE,
    directory = NULL
  ))
  best <- it$iterSummary$minDOFV

  if (is.finite(best) && best < 0) {
    # A better optimum was found: the reference must drop by exactly that much.
    expect_equal(it$newReferenceOfv, fit$objf + best, tolerance = 1e-8)
  } else {
    # Nothing beat the fit, so the reference must not move.
    expect_equal(it$newReferenceOfv, fit$objf, tolerance = 1e-8)
  }
})

test_that("a supplied referenceOfv is what dOFV is measured against", {
  skip_on_cran()
  fit <- theoFit()
  shifted <- fit$objf - 10

  set.seed(42)
  base <- .sirQuiet(sirRunIteration(
    fit,
    mu = .sirProposalMu(fit),
    proposalCov = sirGetProposalCov(fit),
    nSamples = 16L, nResample = 8L, iterNum = 1L,
    recenter = FALSE, boxcox = FALSE, directory = NULL
  ))
  set.seed(42)
  moved <- suppressMessages(suppressWarnings(sirRunIteration(
    fit,
    mu = .sirProposalMu(fit),
    proposalCov = sirGetProposalCov(fit),
    nSamples = 16L, nResample = 8L, iterNum = 1L,
    recenter = FALSE, boxcox = FALSE, directory = NULL,
    referenceOfv = shifted
  )))

  # Same draws, reference lowered by 10, so every dOFV rises by 10.
  expect_equal(
    moved$rawResults$deltaofv[moved$rawResults$role == "sample"],
    base$rawResults$deltaofv[base$rawResults$role == "sample"] + 10,
    tolerance = 1e-8
  )
})

test_that("runSIR persists the reference OFV across added iterations", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()

  set.seed(11)
  .sirQuiet(runSIR(
    fit,
    nSamples = 16L,
    nResample = 8L,
    directory = dir,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))

  state <- nlmixr2utils::readRunState(dir, .sirStateSchema())
  expect_true("nextReferenceOfv" %in% names(state))
  expect_true(is.finite(state$nextReferenceOfv))

  # The reference can only ever improve on where it started, never worsen.
  # Measured against the INITIAL reference, not fit$objf: since P9 the run
  # starts from the evaluator's own value at the centre, which can sit
  # marginally above the published fit$objf (8.2e-05 on this fixture) without
  # anything being wrong.
  initialRef <- nlmixr2sir:::.sirCheckObjective(
    fit, workers = 1L, stencil = FALSE
  )$reevaluated
  expect_lte(state$nextReferenceOfv, unname(initialRef) + 1e-8)

  # It must equal what the last iteration reported, not fit$objf by default.
  last <- state$iterations[[length(state$iterations)]]
  expect_equal(state$nextReferenceOfv, last$newReferenceOfv, tolerance = 1e-10)
})

test_that("added iterations resume from the saved reference OFV", {
  skip_on_cran()
  fit <- theoFit()
  dir <- withr::local_tempdir()

  set.seed(7)
  .sirQuiet(runSIR(
    fit,
    nSamples = 16L,
    nResample = 8L,
    directory = dir,
    control = runSIRControl(objfStencil = FALSE, recover = FALSE, workers = 1L)
  ))
  first <- nlmixr2utils::readRunState(dir, .sirStateSchema())$nextReferenceOfv

  set.seed(8)
  .sirQuiet(runSIR(
    fit,
    nSamples = 16L,
    nResample = 8L,
    directory = dir,
    control = runSIRControl(objfStencil = FALSE, addIterations = TRUE, workers = 1L)
  ))
  state <- nlmixr2utils::readRunState(dir, .sirStateSchema())

  # The extension starts from where the first run left off, not from fit$objf,
  # and the reference can only improve.
  expect_lte(state$nextReferenceOfv, first + 1e-8)
  added <- state$iterations[[length(state$iterations)]]
  expect_lte(added$referenceOfv, first + 1e-8)
})

test_that("the noise band keeps a candidate whose raw ratio overflowed", {
  # sirCalcWeights() normalizes in log space, so prob_resample stays finite even
  # when importance_ratio = exp(log_ir) overflows to Inf. The real resampler
  # uses prob_resample and picks the dominant candidate; the diagnostic filtered
  # on importance_ratio and threw that candidate away, so its band described a
  # different weighted population from the one that produced the result.
  raw <- data.frame(
    sample_id = 1:3,
    role = c("sample", "sample", "sample"),
    deltaofv = c(0, -10, -1500),
    importance_ratio = c(1, 1, Inf),
    prob_resample = c(2.225074e-308, 2.225074e-308, 1),
    resamples = c(0L, 1L, 2L)
  )
  x <- structure(
    data.frame(iter = 1L),
    class = c("nlmixr2SIR", "data.frame"),
    iterations = list(list(rawResults = raw))
  )

  # Several quantiles, as the plot path always passes.
  band <- .sirDofvNoise(x, 1L, quant = c(0.25, 0.5, 0.75), capResampling = 2)

  expect_false(is.null(band))
  # Dropping the overflowed row leaves only dOFVs of 0 and -10, so the band
  # cannot reach -1500. Keeping it, the band sits on the dominant candidate.
  expect_lt(min(band$low), -100)
})

test_that("the noise band works for a single quantile", {
  # vapply() returns a bare vector rather than a 1 x nReplicate matrix when
  # length(quant) == 1, and apply(curves, 1L, ...) then fails with
  # "dim(X) must have a positive length". The plot path always passes several
  # quantiles, so this never surfaced in normal use.
  raw <- data.frame(
    sample_id = 1:4,
    role = rep("sample", 4L),
    deltaofv = c(0, -2, -4, -6),
    importance_ratio = c(1, 2, 3, 4),
    prob_resample = c(0.1, 0.2, 0.3, 0.4),
    resamples = c(1L, 1L, 1L, 1L)
  )
  x <- structure(
    data.frame(iter = 1L),
    class = c("nlmixr2SIR", "data.frame"),
    iterations = list(list(rawResults = raw))
  )

  one <- .sirDofvNoise(x, 1L, quant = 0.5, capResampling = 1)
  expect_s3_class(one, "data.frame")
  expect_equal(nrow(one), 1L)
  expect_equal(one$quantile, 0.5)
  expect_true(is.finite(one$low) && is.finite(one$high))
  expect_lte(one$low, one$high)

  # The many-quantile path must keep working, and agree in shape.
  many <- .sirDofvNoise(x, 1L, quant = c(0.25, 0.5, 0.75), capResampling = 1)
  expect_equal(nrow(many), 3L)
  expect_named(one, names(many))
})
