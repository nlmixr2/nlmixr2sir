# Shared fixtures for the nlmixr2sir test suite.
#
# Fixtures are lazy and memoised: each is built on first use and the outcome --
# success or failure -- is cached. Lazy means nothing expensive is built on
# CRAN, because every test calls skip_on_cran() before touching a fixture.
# Memoised failure means a broken fixture is reported by the tests that need
# it, once each, instead of halting the rest of the file.

.sirLazy <- function(expr) {
  expr <- substitute(expr)
  env <- parent.frame()
  cached <- NULL
  function() {
    if (is.null(cached)) {
      cached <<- tryCatch(
        list(ok = TRUE, value = eval(expr, env)),
        error = function(e) list(ok = FALSE, value = e)
      )
    }
    if (!cached$ok) {
      stop(cached$value)
    }
    cached$value
  }
}

# One-compartment model on theo_sd, single eta, with a covariance step.
theoOneCmt <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- log(31.5)
    eta.ka ~ 0.6
    add.sd <- 0.7 # nolint: object_usage_linter.
  })
  model({
    ka <- exp(tka + eta.ka) # nolint: object_usage_linter.
    cl <- exp(tcl) # nolint: object_usage_linter.
    v <- exp(tv) # nolint: object_usage_linter.
    cp <- linCmt() # nolint: object_usage_linter.
    cp ~ add(add.sd)
  })
}

theoFit <- .sirLazy(suppressMessages(
  nlmixr2utils::nlmixr2(
    theoOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L, covMethod = "r")
  )
))

# Same fit with the covariance step suppressed; exercises the fallback paths.
theoFitNoCov <- .sirLazy(suppressMessages(
  nlmixr2utils::nlmixr2(
    theoOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L, covMethod = "")
  )
))

# Theta-only covariance: covFull = FALSE gives a fit$cov with no OMEGA rows.
# This is the case the automatic OMEGA fallback actually serves -- a covariance
# that is present but incomplete, as opposed to one that is absent entirely.
theoFitThetaCov <- .sirLazy(suppressMessages(suppressWarnings(
  nlmixr2utils::nlmixr2(
    theoOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L, covMethod = "r", covFull = FALSE)
  )
)))

# First-order fit. Its objective is a genuinely different surface from FOCEi's
# (127.98 against 116.80 on theo_sd), which is what makes it the fixture that
# discriminates: an evaluator that quietly scores everything as FOCEi passes
# every other method's preflight and fails this one.
theoFitFo <- .sirLazy(suppressMessages(suppressWarnings(
  nlmixr2utils::nlmixr2(
    theoOneCmt,
    nlmixr2data::theo_sd,
    est = "fo",
    control = list(print = 0L, covMethod = "")
  )
)))

# Three-eta variant, used by the tests that need more than one omega element.
threeEtaOneCmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.00
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7 # nolint: object_usage_linter.
  })
  model({
    ka <- exp(tka + eta.ka) # nolint: object_usage_linter.
    cl <- exp(tcl + eta.cl) # nolint: object_usage_linter.
    v <- exp(tv + eta.v) # nolint: object_usage_linter.
    linCmt() ~ add(add.sd)
  })
}

threeEtaFit <- .sirLazy(suppressMessages(suppressWarnings(
  nlmixr2utils::nlmixr2(
    threeEtaOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L),
    table = list(npde = TRUE, cwres = TRUE)
  )
)))

# Correlated-eta variant. The only fixture with an OMEGA off-diagonal, so it
# is what pins the cov.<eta1>.<eta2> naming bridge.
blockOneCmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.00
    tv <- 3.45
    eta.ka + eta.cl ~ c(0.6, 0.01, 0.3)
    add.sd <- 0.7 # nolint: object_usage_linter.
  })
  model({
    ka <- exp(tka + eta.ka) # nolint: object_usage_linter.
    cl <- exp(tcl + eta.cl) # nolint: object_usage_linter.
    v <- exp(tv) # nolint: object_usage_linter.
    linCmt() ~ add(add.sd)
  })
}

blockFit <- .sirLazy(suppressMessages(suppressWarnings(
  nlmixr2utils::nlmixr2(
    blockOneCmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L)
  )
)))

# A single SIR iteration on theoFit(); tiny schedule for speed.
iter1 <- .sirLazy(local({
  fit <- theoFit()
  set.seed(42)
  .sirQuiet(
    sirRunIteration(
      fit,
      mu = .sirProposalMu(fit),
      proposalCov = sirGetProposalCov(fit),
      nSamples = 16L,
      nResample = 8L,
      iterNum = 1L,
      recenter = TRUE,
      boxcox = TRUE,
      directory = NULL
    )
  )
}))

# A minimal nlmixr2SIR object for the S3 method tests.
sirObj <- .sirLazy(local({
  it <- iter1()
  out <- sirSummary(it$resampledMat, theoFit())
  class(out) <- c("nlmixr2SIR", "data.frame")
  attr(out, "iterationSummary") <- it$iterSummary
  attr(out, "iterations") <- list(it)
  attr(out, "resampledMat") <- it$resampledMat
  attr(out, "outputDir") <- tempdir()
  # runSIR() records the effective controls on its result; this stand-in does
  # the same, so diagnostics read the same provenance they would in a real run.
  attr(out, "control") <- runSIRControl(workers = 1L)
  out
}))

# A canonical raw-results file, produced by a tiny SIR run on theoFit(). Used
# to exercise the rawresInput route, including the round trip from SIR's own
# output back into a new run's proposal.
sirRawResultsPath <- .sirLazy(local({
  dir <- file.path(tempdir(), "sir_rawres_fixture")
  set.seed(1)
  # .sirQuiet(), not suppressMessages(): this fixture is lazy, so its
  # degeneracy warning would otherwise be charged to whichever test happens to
  # touch it first.
  .sirQuiet(runSIR(
    theoFit(),
    nSamples = 12L,
    nResample = 8L,
    directory = dir,
    control = runSIRControl(recover = FALSE, workers = 1L)
  ))
  file.path(dir, "raw_results.csv")
}))

# Run an expression quietly, muffling ONLY the weight-degeneracy warning.
#
# Every fit-based fixture here is deliberately tiny, and small SIR runs on a
# 5-parameter model concentrate their weights -- which is exactly what the S1
# diagnostics exist to report. In a test about schedules, ownership or
# fingerprints that warning is noise, and noise is what hides a genuinely new
# warning. It is asserted directly where it is the subject, in
# test-sir-weight-diagnostics.R.
#
# Note this is NOT suppressWarnings(): any other warning still surfaces and
# still fails a suite that is expected to be quiet. Resizing the fixtures is
# not an alternative -- essFraction is ESS/n and falls as n rises on these
# models (0.54 at n = 16 against 0.02 at n = 500), so a bigger fixture warns
# harder, not less.
.sirQuiet <- function(expr) {
  withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      if (grepl("importance weights are concentrated", conditionMessage(w),
                fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
