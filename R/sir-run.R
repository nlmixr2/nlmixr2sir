# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# runSIR(): the top-level driver.

#' Run sampling importance resampling for an nlmixr2 fit
#'
#' `runSIR()` runs one or more sampling importance resampling iterations using
#' an nlmixr2 fit as the reference model. It follows the PsN SIR workflow where
#' practical for nlmixr2: sample parameter vectors from a proposal covariance,
#' evaluate them with population parameters fixed, compute importance weights,
#' resample, and use the empirical resampled covariance as the next proposal.
#'
#' @param fit An nlmixr2 fit object with a covariance matrix.
#' @param nSamples Integer vector. Requested number of samples per iteration.
#' @param nResample Integer vector. Requested number of resamples per
#'   iteration. Must have the same length as `nSamples`.
#' @param directory Output directory. Run artifacts, saved state and
#'   diagnostics are written here. When `NULL` (default), **nothing is written
#'   to disk**: supplying a directory is what grants permission to write, so a
#'   run that is not given one leaves the filesystem untouched. Recovery and
#'   `addIterations` read the saved state, so both need a directory.
#' @param fitName Optional fit label used in canonical raw-results metadata and
#'   in the names of the files written under `directory`. When `NULL` (default), the
#'   label is derived from the expression supplied to `fit`.
#' @param control A [runSIRControl()] object holding everything that tunes
#'   how the run behaves: inflation, caps, recentering, Box-Cox, parallelism,
#'   resume behaviour, and the OMEGA fallback.
#' @param ... Reserved for future PsN-compatible inputs. Passing a run
#'   setting here is an error; put it in `control` instead.
#' @return A data frame of final SIR summary statistics with class
#'   `c("nlmixr2SIR", "data.frame")`. Attributes include
#'   `iterationSummary`, `iterations`, `resampledMat`, `covMatrix`,
#'   `corMatrix`, and `outputDir`.
#' @examples
#' \dontrun{
#' sir <- runSIR(
#'   fit,
#'   nSamples = c(1000, 1000, 1000),
#'   nResample = c(200, 400, 500),
#'   control = runSIRControl(workers = 4, rxThreads = 2)
#' )
#' }
#' @seealso [runSIRControl()] for the run settings.
#' @export
runSIR <- function(
  fit,
  nSamples = c(1000, 1000, 1000, 2000, 2000),
  nResample = c(200, 400, 500, 1000, 1000),
  directory = NULL,
  fitName = NULL,
  control = runSIRControl(),
  ...
) {
  dots <- list(...)
  if (length(dots) > 0L) {
    cli::cli_abort(c(
      "Unsupported SIR argument(s): {.arg {names(dots)}}.",
      "i" = "Run settings now live in {.fn runSIRControl}."
    ))
  }
  checkmate::assertClass(control, "runSIRControl")

  thetaInflation <- control$thetaInflation
  omegaInflation <- control$omegaInflation
  sigmaInflation <- control$sigmaInflation
  capCorrelation <- control$capCorrelation
  capResampling <- control$capResampling
  recenter <- control$recenter
  boxcox <- control$boxcox
  workers <- control$workers
  rxThreads <- control$rxThreads
  recover <- control$recover
  addIterations <- control$addIterations
  saveFiles <- control$saveFiles
  omegaFallback <- control$omegaFallback
  sigmaFallbackRse <- control$sigmaFallbackRse
  omegaDf <- control$omegaDf

  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertIntegerish(
    nSamples,
    lower = 1,
    any.missing = FALSE,
    min.len = 1L
  )
  checkmate::assertIntegerish(
    nResample,
    lower = 1,
    any.missing = FALSE,
    len = length(nSamples)
  )
  if (is.null(fitName)) {
    fitName <- nlmixr2utils::deriveFitName(substitute(fit))
  }

  nSamples <- as.integer(nSamples)
  nResample <- as.integer(nResample)
  ps <- .sirParamSpace(fit)

  # Free check, so it goes first. Each iteration rebuilds the proposal from the
  # empirical covariance of its retained vectors, which has rank at most
  # nResample - 1; nResample <= nParams therefore cannot produce a usable
  # proposal, and there is no point evaluating a single model to find out.
  n_params <- nrow(ps)
  if (min(nResample) <= n_params) {
    bad <- nResample[nResample <= n_params]
    cli::cli_abort(c(
      "{.arg nResample} is too small for the number of estimated parameters.",
      "x" = "Requested {.val {bad}} for {n_params} parameter{?s}.",
      "i" = "The empirical proposal covariance would have rank at most {min(bad) - 1}.",
      "i" = "Use {.arg nResample} greater than {n_params}; PsN's working ratio is about 5 samples per resample."
    ))
  }

  # A candidate must be drawn before it can be retained, so nSamples bounds the
  # retained set too. Without this the run reached .sirCheckProposalRank() after
  # a full round of model evaluations and then blamed nResample, which was fine.
  if (min(nSamples) <= n_params) {
    bad <- nSamples[nSamples <= n_params]
    cli::cli_abort(c(
      "{.arg nSamples} is too small for the number of estimated parameters.",
      "x" = "Requested {.val {bad}} for {n_params} parameter{?s}.",
      "i" = "At most {min(bad)} distinct vector{?s} can be drawn, so the retained sample cannot reach full rank.",
      "i" = "Increase {.arg nSamples}; PsN's working ratio is about 5 samples per resample."
    ))
  }

  # Limited replacement: each candidate fills at most capResampling slots, so
  # nSamples draws can supply at most nSamples * cap retained vectors. Asking
  # for more silently produced a clamped run rather than saying so.
  cap_int <- max(1L, as.integer(floor(capResampling)))
  if (any(nResample > nSamples * cap_int)) {
    j <- which(nResample > nSamples * cap_int)[[1L]]
    cli::cli_abort(c(
      "{.arg nResample} cannot be met under the current {.arg capResampling}.",
      "x" = "Iteration {j} asks for {nResample[[j]]} from {nSamples[[j]]} candidate{?s} at cap {cap_int}.",
      "i" = "The cap allows at most {nSamples[[j]] * cap_int} retained vector{?s}.",
      "i" = "Lower {.arg nResample}, raise {.arg nSamples}, or raise {.code runSIRControl(capResampling =)}."
    ))
  }

  # Before anything expensive or destructive: prove that the evaluator used for
  # candidates reproduces this fit's own objective. Every dOFV is measured
  # against fit$objf, so if the two are on different surfaces the importance
  # weights are not weights for the advertised target. Runs first so an
  # unsupported fit cannot create or overwrite a run directory.
  .sirCheckObjective(
    fit,
    workers = workers,
    rxThreads = rxThreads,
    objfTolerance = control$objfTolerance,
    stencil = control$objfStencil,
    stencilTolerance = control$objfStencilTolerance
  )

  # Resolve the initial proposal up front, whatever route supplies it, so its
  # resolved numbers can go into the fingerprint. This also validates file-based
  # inputs before any model evaluation. A recovery run pays for this too, and
  # that is the point: re-reading the file is what detects a changed file at an
  # unchanged path.
  initial <- .sirResolveInitialProposal(fit, ps, control)

  # The schedule this call asks for. For a fresh run or a plain recovery it is
  # also the whole schedule; addIterations replaces it below with the prior
  # schedule plus this extension, because the result will contain both.
  request_schedule <- .sirSchedule(nSamples, nResample)
  cumulative_schedule <- request_schedule

  fingerprint <- .sirRunFingerprint(
    fit, ps, request_schedule, control,
    initial = initial
  )

  # Files are written only when the user has named somewhere to write them.
  #
  # CRAN policy is that a package must not write into the user's filespace
  # without explicit consent. runSIR() used to create a numbered
  # <fitName>_sir_<N> directory in the working directory whenever `directory`
  # was left at its NULL default, which is exactly that. Supplying `directory`
  # IS the consent, so that is now what enables persistence.
  #
  # Nothing is written without it: no directory, no state, no per-iteration
  # seed file. A single set.seed() before the call is what makes such a run
  # reproducible, and recovery and addIterations are unavailable because there
  # is nothing on disk to resume from.
  write_files <- saveFiles && !is.null(directory)

  if (!write_files) {
    if (!saveFiles && !is.null(directory)) {
      cli::cli_inform(
        "{.arg saveFiles} is {.code FALSE}; {.arg directory} is ignored and nothing is written."
      )
    } else if (saveFiles) {
      # The common case: defaults. Said once, plainly, because a user expecting
      # artifacts on disk needs to know why there are none.
      cli::cli_inform(c(
        "No {.arg directory} given, so nothing is written to disk.",
        "i" = "Pass {.arg directory} to save run artifacts, state and diagnostics.",
        "i" = "Recovery and {.code addIterations} need those files, so they are unavailable here."
      ))
    }
    if (isTRUE(addIterations)) {
      cli::cli_abort(c(
        "{.code addIterations = TRUE} needs the saved state of a previous run.",
        "x" = if (is.null(directory)) {
          "No {.arg directory} was given, so there is nothing to extend."
        } else {
          "{.arg saveFiles} is {.code FALSE}, so no state was ever written."
        },
        "i" = "Pass the {.arg directory} of the run to extend."
      ))
    }
    output_dir <- NULL
    master_seed <- NULL
    saved_state <- NULL
    recover <- FALSE
    addIterations <- FALSE
  } else {
    run_dir <- nlmixr2utils::resolveRunDir(
      "sir",
      fitName,
      restart = !(recover || addIterations),
      outputDir = directory
    )
    output_dir <- run_dir$path

    # Establish the right to write here before writing anything at all --
    # including the seed file and the manifest. resolveRunDir() returns resume
    # mode for any existing directory under the default recover = TRUE, so
    # guarding only the overwrite path let an unrelated non-empty directory
    # acquire a manifest and become deletable by the next run.
    .sirAssertClaimable(output_dir)

    if (identical(run_dir$mode, "overwrite") && dir.exists(output_dir)) {
      # Only ever clear a directory this package owns. An explicitly supplied
      # path can hold anything.
      .sirAssertSafeToClear(output_dir)
      unlink(output_dir, recursive = TRUE, force = TRUE)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    master_seed <- nlmixr2utils::withRunSeed(output_dir, prefix = "sir")
    saved_state <- if (recover || addIterations) {
      nlmixr2utils::readRunState(output_dir, .sirStateSchema())
    } else {
      NULL
    }

    # An extended run contains the prior iterations and this extension, so its
    # identity has to describe both. Recorded before the manifest is written,
    # because the manifest carries the schedule too.
    #
    # This does not weaken the identity check below: addIterations exempts the
    # schedule field from comparison, and on every other path the cumulative
    # schedule is the requested one.
    if (!is.null(saved_state) && isTRUE(addIterations)) {
      prior_schedule <- saved_state$schedule
      if (is.null(prior_schedule)) {
        # States written before the schedule was persisted. The per-iteration
        # summary has always carried the requested counts, so rebuild from it
        # rather than refusing to extend an otherwise valid run.
        prior_schedule <- data.frame(
          iter = saved_state$iterationSummary$iter,
          nSamples = as.integer(saved_state$iterationSummary$nSamples),
          nResample = as.integer(saved_state$iterationSummary$nResample)
        )
      }
      extension <- .sirSchedule(nSamples, nResample)
      extension$iter <- nrow(prior_schedule) + extension$iter
      cumulative_schedule <- rbind(prior_schedule, extension)
      fingerprint <- .sirRunFingerprint(
        fit, ps, cumulative_schedule, control,
        initial = initial
      )
    }

    # Whatever is on disk has to belong to the run being asked for. Extending a
    # run deliberately changes the schedule, so that field is exempt there.
    if (!is.null(saved_state)) {
      # failClosed: reusing saved state is exactly the situation where an
      # unverifiable identity field must block rather than be waved through.
      bad <- .sirCompareFingerprints(
        saved_state$fingerprint,
        fingerprint,
        ignore = if (isTRUE(addIterations)) "schedule" else character(),
        failClosed = TRUE
      )
      if (length(bad) > 0L) {
        .sirAbortFingerprint(
          bad,
          output_dir,
          if (isTRUE(addIterations)) "addIterations" else "recovery"
        )
      }
      # Dependency versions are recorded but not enforced: a bump does not by
      # itself invalidate a result. It can change proposal construction,
      # Box-Cox estimation or RNG behaviour, though, so it is never silent.
      savedPkgs <- saved_state$fingerprint$pkgVersions
      if (
        !is.null(savedPkgs) &&
          !is.na(savedPkgs) &&
          !identical(savedPkgs, fingerprint$pkgVersions)
      ) {
        cli::cli_warn(c(
          "The saved SIR run was produced under different package versions.",
          "i" = "Saved: {savedPkgs}",
          "i" = "Now:   {fingerprint$pkgVersions}",
          "i" = "Order is nlmixr2sir/nlmixr2est/nlmixr2utils/rxode2. Results are being reused; rerun with {.code recover = FALSE} if that is not wanted."
        ))
      }
    }

    .sirWriteManifest(output_dir, fingerprint, fitName, seed = master_seed)
  }

  if (!is.null(saved_state) && isTRUE(addIterations)) {
    iter_offset <- saved_state$completedIterations
    iter_index <- seq_along(nSamples)
    iter_numbers <- iter_offset + iter_index
    mu <- saved_state$nextMu
    proposal_cov <- saved_state$nextCov
    boxcox_state <- saved_state$nextBoxcoxState
    reference_ofv <- saved_state$nextReferenceOfv %||% fit$objf
    proposal_source <- saved_state$proposalSource %||% NA_character_
    iter_results <- saved_state$iterations
    iter_summary <- saved_state$iterationSummary
    prev_attempted <- tail(iter_summary$nAttempted, 1L)
    prev_successful <- tail(iter_summary$nSuccessful, 1L)
  } else if (!is.null(saved_state) && isTRUE(recover)) {
    completed <- saved_state$completedIterations
    if (completed >= length(nSamples) && !is.null(saved_state$result)) {
      cli::cli_inform("Recovered completed SIR run from {.path {output_dir}}.")
      return(saved_state$result)
    }
    iter_numbers <- seq.int(completed + 1L, length(nSamples))
    iter_index <- iter_numbers
    mu <- saved_state$nextMu
    proposal_cov <- saved_state$nextCov
    boxcox_state <- saved_state$nextBoxcoxState
    reference_ofv <- saved_state$nextReferenceOfv %||% fit$objf
    proposal_source <- saved_state$proposalSource %||% NA_character_
    iter_results <- saved_state$iterations
    iter_summary <- saved_state$iterationSummary
    prev_attempted <- tail(iter_summary$nAttempted, 1L)
    prev_successful <- tail(iter_summary$nSuccessful, 1L)
  } else {
    iter_numbers <- seq_along(nSamples)
    iter_index <- iter_numbers
    proposal_cov <- initial$covMat
    # The raw-results route derives its own centre and Box-Cox state from the
    # supplied vectors, the way PsN's iteration 0 does; every other route
    # centres on the fit's own estimates.
    mu <- initial$mu %||% .sirProposalMu(fit, ps)
    boxcox_state <- initial$boxcoxState
    proposal_source <- initial$source
    if (!identical(initial$source, "cov")) {
      cli::cli_inform(
        "Initial SIR proposal built from {.arg {initial$source}}, not {.code fit$cov}."
      )
    }
    reference_ofv <- fit$objf
    iter_results <- list()
    iter_summary <- data.frame()
    prev_attempted <- NULL
    prev_successful <- NULL
  }

  # A resumed run must not lose the record: its iteration 1 is not re-run.
  initial_repair <- if (is.null(saved_state)) {
    NULL
  } else {
    saved_state$initialProposalRepair
  }

  reference_ofv_history <- if (is.null(saved_state)) {
    numeric(0)
  } else {
    saved_state$referenceOfvHistory %||% numeric(0)
  }

  for (j in seq_along(iter_numbers)) {
    iter_num <- iter_numbers[[j]]
    schedule_idx <- iter_index[[j]]
    requested_samples <- nSamples[[schedule_idx]]
    attempted_samples <- .sirAdjustedAttemptedSamples(
      requested_samples,
      previousAttempted = prev_attempted,
      previousSuccessful = prev_successful
    )
    is_last <- j == length(iter_numbers)
    # Inflation widens the *initial* proposal only, as in PsN. From the second
    # iteration on the proposal is the previous empirical covariance, which
    # must not be re-inflated each time.
    inflate_now <- j == 1L && is.null(saved_state)

    cli::cli_inform(
      "Running SIR iteration {iter_num}: {attempted_samples} attempted samples, {nResample[[schedule_idx]]} requested resamples." # nolint: line_length_linter.
    )
    runIteration <- function() {
      sirRunIteration(
        fit = fit,
        mu = mu,
        proposalCov = proposal_cov,
        nSamples = attempted_samples,
        requestedSamples = requested_samples,
        nResample = nResample[[schedule_idx]],
        iterNum = iter_num,
        capResampling = capResampling,
        recenter = recenter,
        boxcox = boxcox,
        directory = output_dir,
        workers = workers,
        rxThreads = rxThreads,
        boxcoxState = boxcox_state,
        thetaInflation = if (inflate_now) thetaInflation else 1,
        omegaInflation = if (inflate_now) omegaInflation else 1,
        sigmaInflation = if (inflate_now) sigmaInflation else 1,
        capCorrelation = capCorrelation,
        omegaFallback = omegaFallback,
        sigmaFallbackRse = sigmaFallbackRse,
        omegaDf = omegaDf,
        isLastIteration = is_last,
        referenceOfv = reference_ofv
      )
    }
    # Per-iteration seeding exists so a resumed run reproduces the stream it
    # would have had. With nothing persisted there is nothing to resume, so the
    # ambient RNG drives the run and a single set.seed() reproduces it.
    iter_res <- if (write_files) {
      nlmixr2utils::withRunSeed(
        output_dir,
        key = paste0("sir-iteration-", iter_num),
        prefix = "sir",
        expr = runIteration()
      )
    } else {
      runIteration()
    }

    iter_results[[as.character(iter_num)]] <- iter_res
    iter_summary <- rbind(iter_summary, iter_res$iterSummary)
    # Iteration 1 builds its proposal from the resolved initial covariance, so
    # its repair record is the run's. Later iterations repair an empirical
    # update instead, which iter_summary already tracks per iteration.
    if (is.null(initial_repair)) {
      initial_repair <- iter_res$proposalRepair
    }
    if (write_files) {
      .sirWriteIterationSummary(iter_summary, output_dir)
      .sirWriteRejectionSummary(iter_summary, output_dir)
    }

    mu <- iter_res$newMu
    proposal_cov <- iter_res$newCov
    boxcox_state <- iter_res$boxcoxState
    reference_ofv <- iter_res$newReferenceOfv
    reference_ofv_history <- c(
      reference_ofv_history,
      stats::setNames(iter_res$newReferenceOfv, as.character(iter_num))
    )
    prev_attempted <- iter_res$iterSummary$nAttempted
    prev_successful <- iter_res$iterSummary$nSuccessful

    if (write_files) {
      nlmixr2utils::writeRunState(
      output_dir,
      list(
        fingerprint = fingerprint,
        schedule = cumulative_schedule,
        initialProposalRepair = initial_repair,
        proposalSource = proposal_source,
        referenceOfvHistory = reference_ofv_history,
        completedIterations = iter_num,
        nextMu = mu,
        nextCov = proposal_cov,
        nextBoxcoxState = boxcox_state,
        nextReferenceOfv = reference_ofv,
        iterations = iter_results,
        iterationSummary = iter_summary,
        result = NULL
      ),
      .sirStateSchema()
      )
    }
  }

  final_iter <- iter_results[[length(iter_results)]]
  summary_df <- sirSummary(final_iter$resampledMat, fit)
  cov_mat <- attr(summary_df, "covMatrix")
  cor_mat <- attr(summary_df, "corMatrix")
  sdcor_mat <- attr(summary_df, "sdCorMatrix")
  # The canonical raw results are part of the returned object either way; only
  # writing them to disk is optional.
  raw_results <- .sirCanonicalRawResults(fit, fitName, final_iter$resampledMat)
  if (write_files) {
    utils::write.csv(
      summary_df,
      file.path(output_dir, "sir_results.csv"),
      row.names = FALSE
    )
    .sirWriteCovMatrices(summary_df, output_dir, fitName = fitName)
    nlmixr2utils::writeRawResults(raw_results, output_dir)
  }
  # Make the SIR uncertainty selectable with nlmixr2est::setCov(fit, "sir").
  # This lives in the fit, not on disk, so it happens either way.
  .sirRegisterCov(fit, summary_df, ps)

  class(summary_df) <- c("nlmixr2SIR", "data.frame")
  attr(summary_df, "iterationSummary") <- iter_summary
  attr(summary_df, "iterations") <- iter_results
  attr(summary_df, "resampledMat") <- final_iter$resampledMat
  attr(summary_df, "covMatrix") <- cov_mat
  attr(summary_df, "corMatrix") <- cor_mat
  attr(summary_df, "sdCorMatrix") <- sdcor_mat
  attr(summary_df, "outputDir") <- output_dir
  attr(summary_df, "fitName") <- fitName
  attr(summary_df, "rawResults") <- raw_results
  attr(summary_df, "seed") <- master_seed
  attr(summary_df, "call") <- match.call()
  # Complete provenance on the object itself, so a result read back from an
  # .rds can still say what produced it and the diagnostics can describe the
  # algorithm that actually ran rather than assuming defaults.
  attr(summary_df, "control") <- control
  attr(summary_df, "schedule") <- cumulative_schedule
  attr(summary_df, "initialProposalRepair") <- initial_repair
  attr(summary_df, "proposalSource") <- proposal_source
  attr(summary_df, "referenceOfvHistory") <- reference_ofv_history
  attr(summary_df, "fingerprint") <- fingerprint
  # Records whether files were actually written, which is saveFiles AND a
  # directory having been supplied -- not the control flag alone.
  attr(summary_df, "saveFiles") <- write_files

  if (write_files) {
    nlmixr2utils::writeRunState(
      output_dir,
      list(
        fingerprint = fingerprint,
        schedule = cumulative_schedule,
        initialProposalRepair = initial_repair,
        proposalSource = proposal_source,
        referenceOfvHistory = reference_ofv_history,
        completedIterations = tail(iter_summary$iter, 1L),
        nextMu = mu,
        nextCov = proposal_cov,
        nextBoxcoxState = boxcox_state,
        nextReferenceOfv = reference_ofv,
        iterations = iter_results,
        iterationSummary = iter_summary,
        result = summary_df
      ),
      .sirStateSchema()
    )
  }

  summary_df
}
