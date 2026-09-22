#' Configure an `nlmixr2sir` run
#'
#' `runSIRControl()` constructs and validates the control object used by
#' [runSIR()]. It holds everything that tunes *how* SIR runs, leaving
#' [runSIR()] itself to take only the fit, the sampling schedule, and where
#' output goes.
#'
#' @param thetaInflation,omegaInflation,sigmaInflation Non-negative variance
#'   multipliers applied to the THETA, OMEGA, and residual-error blocks of the
#'   **initial** proposal. From the second iteration onwards the proposal is
#'   the previous iteration's empirical covariance and is not re-inflated.
#' @param capCorrelation Numeric in `[0, 1]`. Maximum absolute proposal
#'   correlation after covariance construction and updates.
#' @param capResampling Numeric greater than or equal to one. `1` resamples
#'   without replacement; larger values allow limited replacement.
#' @param recenter Logical. If `TRUE`, recenter the next proposal on the best
#'   sampled vector when any sample has negative dOFV.
#' @param boxcox Logical. If `TRUE`, use Box-Cox transformed empirical
#'   covariances for non-final iterations.
#' @param workers `NULL`, `"auto"`, `1`, or a positive integer. Controls
#'   parallel OFV evaluation through `future`. `NULL` leaves the current
#'   `future::plan()` unchanged, `1` forces sequential execution, a positive
#'   integer temporarily uses a multisession plan, and `"auto"` uses
#'   `future::availableCores(omit = 1L)`. SIR iterations are always sequential,
#'   because each iteration builds the next iteration's proposal.
#' @param rxThreads Integer, `"auto"`, or `NULL`; rxode2 OpenMP threads per
#'   worker. `NULL` (the default) uses the current `rxode2::getRxThreads()`
#'   value for every worker; `"auto"` divides the core count evenly across
#'   workers. Whenever `workers > 1`, `workers * rxThreads` must not exceed the
#'   machine's core count, since each worker is a separate process running its
#'   own rxode2 thread pool.
#' @param saveFiles Logical. If `TRUE` (default) the run writes its raw
#'   results, summaries, covariance files, seed and resumable state to
#'   `directory`. Set to `FALSE` to run entirely in memory, writing nothing and
#'   creating no directory; the result is returned as usual. Recovery,
#'   `addIterations`, and per-iteration seed reproduction all need the saved
#'   state, so they are unavailable when this is `FALSE` -- seed a run with
#'   `set.seed()` beforehand to reproduce it.
#' @param objfStencil Logical. If `TRUE` (default) the objective preflight also
#'   evaluates small perturbations either side of every parameter and requires
#'   the fitted estimates to remain a local optimum. Agreement at the centre
#'   alone does not establish that two objectives are the same function.
#' @param objfStencilTolerance Non-negative number. How much a probe may lower
#'   the objective before the run is refused. Smaller decreases warn instead,
#'   since a fit that stopped just short of convergence is common and harmless.
#' @param objfTolerance Non-negative number, default `1e-3`, interpreted as a
#'   **fraction of the objective** (floored at an absolute `1e-2`). Before
#'   sampling, SIR re-evaluates the objective at the fit's own estimates and
#'   compares it with `fit$objf`; the run aborts only if the absolute
#'   difference exceeds that threshold.
#'
#'   It is relative because that is the scale on which the two things this
#'   check must separate actually differ: convergence slack measured across 20
#'   population PK models stayed within `1.4e-5` of the objective, while a SAEM
#'   fit scored under FOCEi is `2.3e-2` of it. Those are 1600x apart relatively
#'   and only 10x apart absolutely.
#'
#'   The difference itself is a constant across candidates, and since dOFVs are
#'   measured against the re-evaluated centre it does not reach the importance
#'   weights at all. Refitting at a higher `sigdig` shrinks it about 3-4 fold
#'   per digit.
#' @param objfNoise Logical, default `TRUE`. Measure the evaluator's noise
#'   floor before sampling, by walking a short transect through parameter space
#'   at the proposal's own scale and taking the residual from a smooth
#'   polynomial. Unlike the reproduction difference above, this part does *not*
#'   cancel between a candidate and the centre, so it enters every weight as
#'   `exp(-noise/2)`. Costs 15 model evaluations.
#' @param objfNoiseTolerance Non-negative number, default `1` OFV unit. The run
#'   *warns* above this and never refuses: at 1 OFV unit a single weight is
#'   perturbed by about 39%, which is worth knowing, but a diagnostic whose own
#'   estimator can be wrong has no business stopping a run. When the estimate
#'   is not trustworthy -- the transect's polynomial has not absorbed the
#'   objective's shape -- nothing is reported at all.
#' @param rankDeficiency How to handle retained vectors that cannot support a
#'   full-rank covariance. `"abort"` (the default) refuses, because flooring a
#'   deficient direction's eigenvalue does not recover missing information --
#'   it fabricates variance the sample never supported, and later iterations
#'   then propose along it. `"repair"` proceeds with those directions pinned to
#'   a token variance, with a warning saying so.
#'
#'   `"repair"` exists for fits whose own covariance is degenerate, where no
#'   amount of resampling helps: a Michaelis-Menten model on single-dose data,
#'   for instance, may not identify every direction, and `eigen(fit$cov)$values`
#'   will show it. Intervals for parameters loading on a repaired direction are
#'   not evidence from the data, so this is opt-in rather than automatic.
#' @param recover Logical. If `TRUE` and the output directory holds
#'   `sir_state.rds`, resume from the last completed iteration when possible.
#' @param addIterations Logical. If `TRUE`, append the supplied schedule after
#'   an existing completed state in the output directory.
#' @param omegaFallback How OMEGA uncertainty is obtained. `"cov"` (the
#'   default) takes it from `fit$cov`, including its correlations with THETA,
#'   and degrades to `"wishart"` automatically when `fit$cov` does not carry
#'   OMEGA. `"wishart"` always uses the Wishart-style approximation, which
#'   gives a block-diagonal proposal.
#' @param sigmaFallbackRse Percent relative standard error used for
#'   residual-error uncertainty when neither `fit$cov` nor `fit$parFixedDf`
#'   reports a standard error.
#' @param omegaDf Optional degrees of freedom for the Wishart-style OMEGA
#'   fallback; defaults to `nsub - 1`. Unused when OMEGA comes from `fit$cov`.
#' @param rseTheta,rseOmega,rseSigma Relative standard errors, as percentages,
#'   used to build a diagonal proposal without a covariance step. Each is
#'   either a single value for the whole class or one value per estimated
#'   (non-fixed) element of it -- THETA, OMEGA **diagonals**, residual error.
#'   Diagonal variances come out as `(rse * estimate / 100)^2`; OMEGA
#'   off-diagonals are derived from the two diagonals they connect. Following
#'   PsN, a scalar `rseTheta` fills in an unset `rseOmega` and `rseSigma`, a
#'   vector `rseTheta` does not, and setting `rseOmega` without `rseTheta` is
#'   an error. Cannot be combined with `covmatInput` or with inflation.
#' @param covmatInput A proposal covariance supplied directly: a numeric
#'   matrix, a path to a NONMEM-style `.cov` file, or the string `"identity"`.
#'   An unnamed matrix must cover the whole SIR parameter vector in order; a
#'   named one may be a subset, keyed by either SIR or `fit$cov` names.
#'   `"identity"` together with inflation is the cheap "any diagonal proposal"
#'   route. Cannot be combined with `rseTheta`.
#' @param rawresInput A canonical raw-results file path or data frame whose
#'   parameter vectors seed the first proposal, as PsN's iteration 0 does: the
#'   empirical mean and covariance of the supplied vectors become the
#'   iteration-1 proposal. Any canonical raw-results file works, including one
#'   written by `nlmixr2boot`. Cannot be combined with `rseTheta` or
#'   `covmatInput`.
#' @param offsetRawres Integer. Skip raw-results samples numbered below this.
#'   Defaults to `1`, which drops the reference row.
#' @param inFilter Optional filter applied to the raw-results rows before they
#'   are used, in any form accepted by
#'   [nlmixr2utils::setupRawResultsFilter()]. Only meaningful with
#'   `rawresInput`.
#'
#' @return An object of class `runSIRControl`.
#' @examples
#' runSIRControl(thetaInflation = 2, workers = 4, rxThreads = 2)
#' @export
runSIRControl <- function(
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8,
  capResampling = 1,
  recenter = TRUE,
  boxcox = TRUE,
  workers = NULL,
  rxThreads = NULL,
  rankDeficiency = c("abort", "repair"),
  recover = TRUE,
  addIterations = FALSE,
  saveFiles = TRUE,
  objfTolerance = 1e-3,
  objfNoise = TRUE,
  objfNoiseTolerance = 1,
  objfStencil = TRUE,
  objfStencilTolerance = 1,
  omegaFallback = c("cov", "wishart"),
  sigmaFallbackRse = 30,
  omegaDf = NULL,
  rseTheta = NULL,
  rseOmega = NULL,
  rseSigma = NULL,
  covmatInput = NULL,
  rawresInput = NULL,
  offsetRawres = 1L,
  inFilter = NULL
) {
  omegaFallback <- match.arg(omegaFallback)

  checkmate::assertNumeric(
    thetaInflation,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    min.len = 1L
  )
  checkmate::assertNumeric(
    omegaInflation,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    min.len = 1L
  )
  checkmate::assertNumeric(
    sigmaInflation,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    min.len = 1L
  )
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)
  checkmate::assertFlag(recenter)
  checkmate::assertFlag(boxcox)
  rankDeficiency <- match.arg(rankDeficiency)
  checkmate::assertFlag(recover)
  checkmate::assertFlag(addIterations)
  checkmate::assertFlag(saveFiles)
  if (!saveFiles && isTRUE(addIterations)) {
    cli::cli_abort(c(
      "{.arg addIterations} needs {.arg saveFiles} to be {.code TRUE}.",
      "i" = "Extending a run reads the saved state of the run it extends."
    ))
  }
  checkmate::assertNumber(objfTolerance, lower = 0, finite = TRUE)
  checkmate::assertFlag(objfNoise)
  checkmate::assertNumber(objfNoiseTolerance, lower = 0, finite = TRUE)
  checkmate::assertFlag(objfStencil)
  checkmate::assertNumber(objfStencilTolerance, lower = 0, finite = TRUE)
  checkmate::assertNumber(sigmaFallbackRse, lower = 0, finite = TRUE)
  if (!is.null(omegaDf)) {
    checkmate::assertNumber(omegaDf, lower = 1, finite = TRUE)
  }
  nlmixr2utils::.validateWorkers(workers)

  for (nm in c("rseTheta", "rseOmega", "rseSigma")) {
    v <- get(nm)
    if (!is.null(v)) {
      checkmate::assertNumeric(
        v,
        lower = .Machine$double.eps,
        finite = TRUE,
        any.missing = FALSE,
        min.len = 1L,
        .var.name = nm
      )
    }
  }
  sources <- c(
    covmatInput = !is.null(covmatInput),
    rseTheta = !is.null(rseTheta),
    rawresInput = !is.null(rawresInput)
  )
  if (sum(sources) > 1L) {
    cli::cli_abort(c(
      "{.arg {names(sources)[sources]}} are alternative proposal sources; give only one.",
      "i" = "Dispatch order when several are set would be ambiguous."
    ))
  }
  checkmate::assertCount(offsetRawres)
  if (!is.null(inFilter) && is.null(rawresInput)) {
    cli::cli_warn(
      "{.arg inFilter} has no effect without {.arg rawresInput}."
    )
  }
  # PsN forbids rse_* together with inflation: the RSE already states the
  # width, so inflating it on top makes the stated RSE a fiction.
  inflated <- !all(
    c(thetaInflation, omegaInflation, sigmaInflation) == 1
  )
  if (!is.null(rseTheta) && inflated) {
    cli::cli_abort(c(
      "Inflation cannot be combined with {.arg rseTheta}.",
      "i" = "The RSE already specifies the proposal width; widen the RSE instead."
    ))
  }

  structure(
    list(
      thetaInflation = thetaInflation,
      omegaInflation = omegaInflation,
      sigmaInflation = sigmaInflation,
      capCorrelation = capCorrelation,
      capResampling = capResampling,
      recenter = recenter,
      boxcox = boxcox,
      workers = workers,
      rxThreads = rxThreads,
      rankDeficiency = rankDeficiency,
      recover = recover,
      addIterations = addIterations,
      saveFiles = saveFiles,
      objfTolerance = objfTolerance,
      objfNoise = objfNoise,
      objfNoiseTolerance = objfNoiseTolerance,
      objfStencil = objfStencil,
      objfStencilTolerance = objfStencilTolerance,
      omegaFallback = omegaFallback,
      sigmaFallbackRse = sigmaFallbackRse,
      omegaDf = omegaDf,
      rseTheta = rseTheta,
      rseOmega = rseOmega,
      rseSigma = rseSigma,
      covmatInput = covmatInput,
      rawresInput = rawresInput,
      offsetRawres = as.integer(offsetRawres),
      inFilter = inFilter
    ),
    class = "runSIRControl"
  )
}

#' @export
print.runSIRControl <- function(x, ...) {
  cli::cli_h2("nlmixr2sir control")
  cli::cli_dl(c(
    inflation = "theta {x$thetaInflation}, omega {x$omegaInflation}, sigma {x$sigmaInflation}",
    caps = "correlation {x$capCorrelation}, resampling {x$capResampling}",
    proposal = "recenter {x$recenter}, boxcox {x$boxcox}, omegaFallback {.val {x$omegaFallback}}",
    parallel = "workers {format(x$workers %||% 'current plan')}, rxThreads {format(x$rxThreads %||% 'rxode2 default')}",
    resume = "recover {x$recover}, addIterations {x$addIterations}"
  ))
  invisible(x)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# rxode2::rxUiDeparse() turns an object stored in a model UI's `meta`
# environment back into reproducible source, emitting only the arguments that
# differ from the defaults. .deparseFinal() builds the call as
# paste0(var, " <- ", class(object), "(...)"), taking the *class name* as the
# constructor name -- which is why this object's class is `runSIRControl` and
# not a package-prefixed variant.
#' @exportS3Method rxode2::rxUiDeparse
rxUiDeparse.runSIRControl <- function(object, var) {
  .default <- runSIRControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}
