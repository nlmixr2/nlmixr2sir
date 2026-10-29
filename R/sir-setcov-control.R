#' Options for the SIR covariance in `setCov()`
#'
#' `sirControl()` holds the options of `nlmixr2est::setCov(fit, "sir")`, which
#' runs sampling importance resampling at the fit's estimates and installs the
#' resampled covariance as the fit's covariance. It holds only what shapes the
#' SIR result, plus how the model evaluations are spread over workers. The
#' covariance step runs in memory: nothing is written to disk, and a run cannot
#' be resumed or extended. Use [runSIR()] when you want those, or the full set
#' of diagnostics; its result can still be installed with
#' `setCov(fit) <- sirResult`.
#'
#' @section The seed covariance:
#'
#' SIR starts from a proposal built from a *seed* covariance. By default
#' (`seedCov = NULL`) this is the covariance installed on the fit, such as the
#' `"r,s"` or `"analytic"` covariance the estimation computed. When `"sir"` is
#' installed itself, the seed that SIR covariance was computed from is used
#' again, so SIR is never seeded from its own result.
#'
#' If the fit has no covariance at all, the seed is a diagonal proposal built
#' from relative standard errors (`rseTheta`, by default 30%), as in
#' `runSIRControl(rseTheta = 30)`. This is an assumed uncertainty, not an
#' estimated one, and `setCov()` says so when it uses it.
#'
#' A SIR covariance is cached on the fit like any other `setCov()`
#' covariance. It is reused only when both these options and the seed
#' covariance are the same as when it was computed. Changing the installed
#' covariance (for example with `setCov(fit, "analytic")`) changes the default
#' seed, so a later `setCov(fit, "sir")` recomputes.
#'
#' @param nSamples,nResample Samples and resamples per iteration, as in
#'   [runSIR()]. The default is the PsN schedule.
#' @inheritParams runSIRControl
#' @param seedCov The seed covariance. `NULL` (default) uses the installed
#'   covariance, as described above; a covariance-method name uses that entry of
#'   `fit$covList` (or `fit$cov` when it is installed); a named matrix is used
#'   directly, keyed by either `fit$cov` or SIR parameter names.
#' @param rseTheta,rseOmega,rseSigma Relative standard errors, as percentages,
#'   for the seed used **only when the fit has no covariance**; see
#'   [runSIRControl()]. They cannot be combined with inflation.
#' @param seed Random seed for the SIR sampling. A fixed seed makes the
#'   covariance reproducible, which is what lets it be cached. The global
#'   random-number stream is restored afterwards.
#'
#' @return An object of class `sirControl`.
#' @seealso [runSIR()], [runSIRControl()], `nlmixr2est::setCov()`
#' @examples
#' sirControl()
#' sirControl(nSamples = c(500, 500), nResample = c(100, 200), workers = 4)
#' \dontrun{
#' setCov(fit, "sir")
#' setCov(fit, "sir", control = sirControl(thetaInflation = 1.5))
#' fit$sir # the SIR result behind the installed covariance
#' }
#' @export
sirControl <- function(
  nSamples = c(1000, 1000, 1000, 2000, 2000),
  nResample = c(200, 400, 500, 1000, 1000),
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8,
  capResampling = 1,
  recenter = TRUE,
  boxcox = TRUE,
  omegaFallback = c("cov", "wishart"),
  sigmaFallbackRse = 30,
  omegaDf = NULL,
  seedCov = NULL,
  rseTheta = 30,
  rseOmega = NULL,
  rseSigma = NULL,
  seed = 42L,
  workers = NULL,
  rxThreads = NULL,
  objfTolerance = 1e-4,
  objfStencil = TRUE,
  objfStencilTolerance = 1
) {
  omegaFallback <- match.arg(omegaFallback)
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)
  .sirCheckSchedule(nSamples, nResample, capResampling)
  .sirAssertRse(rseTheta, rseOmega, rseSigma)
  if (is.null(rseTheta) && (!is.null(rseOmega) || !is.null(rseSigma))) {
    cli::cli_abort(c(
      "{.arg rseOmega}/{.arg rseSigma} were given without {.arg rseTheta}.",
      "i" = "The RSE seed is driven by {.arg rseTheta}; set it too."
    ))
  }
  checkmate::assertIntegerish(seed, len = 1L, any.missing = FALSE)
  .sirAssertSeedCov(seedCov)

  # Everything else is validated by the control the run itself uses, so the
  # two cannot disagree about what is allowed.
  runSIRControl(
    thetaInflation = thetaInflation,
    omegaInflation = omegaInflation,
    sigmaInflation = sigmaInflation,
    capCorrelation = capCorrelation,
    capResampling = capResampling,
    recenter = recenter,
    boxcox = boxcox,
    workers = workers,
    rxThreads = rxThreads,
    recover = FALSE,
    saveFiles = FALSE,
    objfTolerance = objfTolerance,
    objfStencil = objfStencil,
    objfStencilTolerance = objfStencilTolerance,
    omegaFallback = omegaFallback,
    sigmaFallbackRse = sigmaFallbackRse,
    omegaDf = omegaDf
  )

  structure(
    list(
      nSamples = as.integer(nSamples),
      nResample = as.integer(nResample),
      # Stored as doubles, so 1L and 1 key the setCov() cache the same way.
      thetaInflation = as.double(thetaInflation),
      omegaInflation = as.double(omegaInflation),
      sigmaInflation = as.double(sigmaInflation),
      capCorrelation = as.double(capCorrelation),
      capResampling = as.double(capResampling),
      recenter = recenter,
      boxcox = boxcox,
      omegaFallback = omegaFallback,
      sigmaFallbackRse = as.double(sigmaFallbackRse),
      omegaDf = if (!is.null(omegaDf)) as.double(omegaDf),
      seedCov = seedCov,
      rseTheta = if (!is.null(rseTheta)) as.double(rseTheta),
      rseOmega = if (!is.null(rseOmega)) as.double(rseOmega),
      rseSigma = if (!is.null(rseSigma)) as.double(rseSigma),
      seed = as.integer(seed),
      workers = workers,
      rxThreads = rxThreads,
      objfTolerance = objfTolerance,
      objfStencil = objfStencil,
      objfStencilTolerance = objfStencilTolerance
    ),
    class = "sirControl"
  )
}

.sirAssertSeedCov <- function(seedCov) {
  if (is.null(seedCov)) {
    return(invisible(TRUE))
  }
  if (is.character(seedCov)) {
    checkmate::assertString(seedCov, min.chars = 1L, .var.name = "seedCov")
    if (identical(seedCov, "sir")) {
      cli::cli_abort(c(
        "{.arg seedCov} cannot be {.val sir}.",
        "i" = "SIR is seeded from another covariance, never from its own result."
      ))
    }
    return(invisible(TRUE))
  }
  checkmate::assertMatrix(
    seedCov,
    mode = "numeric",
    any.missing = FALSE,
    min.rows = 1L,
    .var.name = "seedCov"
  )
  if (
    nrow(seedCov) != ncol(seedCov) ||
      is.null(rownames(seedCov)) ||
      !identical(rownames(seedCov), colnames(seedCov))
  ) {
    cli::cli_abort(
      "A {.arg seedCov} matrix must be square with matching row and column names."
    )
  }
  invisible(TRUE)
}

# The runSIRControl() a sirControl() runs under, for a resolved seed. The RSE
# options only apply when there is no seed covariance.
.sirControlToRun <- function(ctl, seed) {
  useRse <- identical(seed$method, "rse")
  if (useRse) {
    inflated <- !all(
      c(ctl$thetaInflation, ctl$omegaInflation, ctl$sigmaInflation) == 1
    )
    if (inflated) {
      cli::cli_abort(c(
        "The fit has no covariance, so the SIR seed comes from {.arg rseTheta}, which cannot be combined with inflation.",
        "i" = "Widen the RSE instead, or give a seed with {.code sirControl(seedCov =)}."
      ))
    }
  }
  runSIRControl(
    thetaInflation = ctl$thetaInflation,
    omegaInflation = ctl$omegaInflation,
    sigmaInflation = ctl$sigmaInflation,
    capCorrelation = ctl$capCorrelation,
    capResampling = ctl$capResampling,
    recenter = ctl$recenter,
    boxcox = ctl$boxcox,
    workers = ctl$workers,
    rxThreads = ctl$rxThreads,
    recover = FALSE,
    addIterations = FALSE,
    saveFiles = FALSE,
    objfTolerance = ctl$objfTolerance,
    objfStencil = ctl$objfStencil,
    objfStencilTolerance = ctl$objfStencilTolerance,
    omegaFallback = ctl$omegaFallback,
    sigmaFallbackRse = ctl$sigmaFallbackRse,
    omegaDf = ctl$omegaDf,
    rseTheta = if (useRse) ctl$rseTheta,
    rseOmega = if (useRse) ctl$rseOmega,
    rseSigma = if (useRse) ctl$rseSigma
  )
}

#' @export
print.sirControl <- function(x, ...) {
  cli::cli_h2("nlmixr2sir covariance control")
  seedTxt <- if (is.null(x$seedCov)) {
    "installed covariance"
  } else if (is.character(x$seedCov)) {
    x$seedCov
  } else {
    paste0(nrow(x$seedCov), "x", ncol(x$seedCov), " matrix")
  }
  cli::cli_dl(c(
    schedule = "nSamples {x$nSamples}; nResample {x$nResample}",
    seed = "{seedTxt} (no covariance: rseTheta {format(x$rseTheta %||% 'none')}); RNG seed {x$seed}",
    inflation = "theta {x$thetaInflation}, omega {x$omegaInflation}, sigma {x$sigmaInflation}",
    caps = "correlation {x$capCorrelation}, resampling {x$capResampling}",
    proposal = "recenter {x$recenter}, boxcox {x$boxcox}, omegaFallback {.val {x$omegaFallback}}",
    parallel = "workers {format(x$workers %||% 'current plan')}, rxThreads {format(x$rxThreads %||% 'rxode2 default')}"
  ))
  invisible(x)
}

#' @exportS3Method rxode2::rxUiDeparse
rxUiDeparse.sirControl <- function(object, var) {
  .default <- sirControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}
