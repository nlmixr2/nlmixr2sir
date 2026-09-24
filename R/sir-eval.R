# Part of nlmixr2sir. Split out of the original single-file R/sir.R.
# OFV evaluation of sampled parameter vectors.

# Step 4 -----------------------------------------------------------------------

#' Evaluate OFV for a matrix of sampled parameter vectors
#'
#' For each row of `paramSamples`, seeds the model's ini block with the
#' supplied values (theta and/or lower-triangle omega elements) and calls
#' `nlmixr2(est = <the fit's own method>)` with `maxOuterIterations = 0`, making it the
#' nlmixr2 equivalent of NONMEM `MAXEVAL=0`.  THETA column names must match
#' parameter names in `fit$iniDf`; omega columns use the SIR lower-triangle
#' proposal names.
#'
#' Parameters NOT present as columns in `paramSamples` retain their estimated
#' values from `fit`.
#'
#' @param fit An nlmixr2 fit object (carries the data internally).
#' @param paramSamples Numeric matrix, one row per sample.  Column names may
#'   include THETA names from `fit$iniDf` and SIR omega lower-triangle names.
#' @param workers Passed to `.withWorkerPlan()`: `NULL` (keep current plan),
#'   `1` (force sequential), a positive integer, or `"auto"`.
#' @param rxThreads Integer, `"auto"`, or `NULL`; rxode2 OpenMP threads per
#'   worker. Required by `.withWorkerPlan()` whenever `workers > 1`, since
#'   each worker runs its own thread pool.
#' @return Named numeric vector of length `nrow(paramSamples)`.  Entries are
#'   `NA_real_` for rows that produced an error during evaluation.
#' @noRd
# Evaluation-only overrides: what the call DOES, as opposed to what the
# objective MEANS. These are the only fields the evaluator is entitled to
# change, and everything else is taken from the fit untouched.
#
# maxOuterIterations = 0   fix the population parameters; this is the whole
#                          point of the evaluator
# calcTables, compress     output shaping, not likelihood
# covMethod = ""           a covariance step per candidate would be enormous
#                          waste and is never read
# print = 0                quiet
.sirEvalOverrides <- list(
  maxOuterIterations = 0L,
  calcTables = FALSE,
  covMethod = "",
  compress = FALSE,
  print = 0L
)

.sirFitEst <- function(fit) {
  # fit$env$est, not fit$est. A fitted object is data-frame-like, so `$est`
  # can resolve to an output-table COLUMN instead of the scalar method name:
  # an `fo` fit with calcTables = TRUE returns a 132-long vector on theo_sd,
  # one element per row of the data. fit$env$est is the scalar in every case
  # measured (focei and fo, tables on and off). fit$est is kept only as a
  # fallback for an object with no env.
  est <- tryCatch(as.character(fit$env$est), error = function(e) character(0))
  if (length(est) != 1L || is.na(est) || !nzchar(est)) {
    est <- tryCatch(as.character(fit$est), error = function(e) character(0))
  }
  if (length(est) != 1L || is.na(est) || !nzchar(est)) {
    return(NA_character_)
  }
  est
}

# The importance-sampling family. These are the one documented exception to
# "est selects the objective, so dispatch on the fit's own est".
#
# nlmixr2est recomputes the objective of EVERY imp/impmap/qrpem fit as a nested
# FOCEi evaluation at the converged estimates (.impmapRecomputeObjf(),
# nlmixr2est R/impmap.R:1106), because the in-C++ finalize leaves d(pred)/d(eta)
# at allocation residue and the eta-Hessian collapses to Omega^-1. So fit$objf
# on an imp fit is NOT the importance-sampling objective -- that is $impObj,
# which is a different number and is never what SIR scores against.
#
# Scoring such a candidate as its own method would run a MAP pass and an E-step
# over isample draws per subject and THEN the same nested FOCEi call, to arrive
# at the number the FOCEi call alone produces. Measured on theo_sd, both one-eta
# and three-eta: the two routes agree BIT-FOR-BIT, and the FOCEi route is 6-9x
# faster per candidate. The E-step is pure cost, and it is per candidate.
.sirImpFamilyMethods <- c("imp", "impmap", "qrpem")

.sirEvalMethod <- function(fit) {
  est <- .sirFitEst(fit)
  if (is.na(est)) {
    return("focei")
  }
  # Deliberate: the evaluator redirects, the RECORDED method does not. The
  # allowlist and every provenance record still see the fit's own est.
  if (est %in% .sirImpFamilyMethods) {
    return("focei")
  }
  est
}

.sirEvalControlFun <- function(est) {
  tryCatch(
    getExportedValue("nlmixr2est", paste0(est, "Control")),
    error = function(e) nlmixr2est::foceiControl
  )
}

# The evaluation control for an imp-family fit.
#
# This is the one place the fail-open-is-worse argument below INVERTS. Carrying
# the fit's whole control is right everywhere else, because an unrecognised
# setting is then preserved rather than dropped. Here it is wrong twice over:
#
#   * the reference number is not produced with the fit's control.
#     .impmapRecomputeObjf() builds a DEFAULT foceiControl() and carries only
#     sigdig across (nlmixr2est R/impmap.R:1162-1165). Reproducing fit$objf
#     means reproducing that recipe, not the impmap control the fit was run
#     with -- which is not a foceiControl at all.
#   * carrying it is measurably worse. An impmap evaluation with the fit's own
#     control reports "diag(V) had non-positive or NA entries; the non-finite
#     result may be dubious" and "NaNs produced" on $runInfo, which the same
#     evaluation with a fresh control does not, and takes 3x as long.
#
# So: a bare foceiControl() with the evaluation overrides and the fit's sigdig.
.sirImpEvalControl <- function(fit) {
  sigdig <- tryCatch(fit$foceiControl$sigdig, error = function(e) NULL)
  # 4 is nlmixr2est's own fallback in .impmapRecomputeObjf(); matching it keeps
  # the two recipes identical rather than merely similar.
  if (is.null(sigdig) || length(sigdig) != 1L || is.na(sigdig)) {
    sigdig <- 4
  }
  ctl <- tryCatch(
    do.call(
      nlmixr2est::foceiControl,
      c(.sirEvalOverrides, list(sigdig = sigdig))
    ),
    error = function(e) NULL
  )
  # Fail closed, for the same reason the general path does: a candidate scored
  # on the wrong surface is not detectable from the result.
  if (is.null(ctl)) {
    cli::cli_abort(c(
      "Could not build the FOCEi evaluation control for this importance-sampling fit.",
      "i" = "SIR scores {.val {.sirImpFamilyMethods}} candidates as FOCEi, because {.code fit$objf} is itself a FOCEi re-evaluation.",
      "i" = "This is refused rather than approximated."
    ))
  }
  ctl
}

.sirEvalControl <- function(fit) {
  # Branch on the RECORDED method, not the evaluation method: .sirEvalMethod()
  # already maps the imp family onto focei, so asking it here would lose the
  # very distinction this branch needs.
  recorded <- .sirFitEst(fit)
  if (!is.na(recorded) && recorded %in% .sirImpFamilyMethods) {
    return(.sirImpEvalControl(fit))
  }
  est <- .sirEvalMethod(fit)
  ctlFun <- .sirEvalControlFun(est)
  ctl <- fit$control
  if (!is.list(ctl) || length(ctl) == 0L) {
    ctl <- NULL
  }
  # Everything else comes from the fit's OWN control object, unchanged.
  #
  # This used to rebuild the control from a hand-picked allowlist of
  # "likelihood-relevant" fields. foceiControl() has 150 arguments, so that
  # list failed OPEN: a field nobody thought of was silently dropped and the
  # candidate scored on a different surface. agqLow/agqHi were dropped exactly
  # that way -- an AGQ fit with agqLow = -100 was reevaluated with the default
  # -Inf, which agrees at the centre to 1e-08 (so the preflight passed) and is
  # 6490 OFV units out at tka = -20. That enters the weight as
  # exp(-dOFV/2). The stencil cannot catch it either: it perturbs by a
  # thousandth of each estimate, and integration bounds only bite far away.
  #
  # Copying the fit's control and overriding only the evaluation fields
  # inverts that: an unrecognised setting is PRESERVED rather than lost, so
  # the failure mode of being wrong about this list is a slower evaluation,
  # not a wrong answer.
  # The override VALUES have to be normalised the way the constructor would
  # normalise them, not assigned raw: foceiControl() turns covMethod = "" into
  # integer 0 and print = 0 into NULL, so raw assignment would leave the
  # control in a shape the constructor never produces. Build one reference
  # control and copy its versions of exactly those fields.
  ref <- tryCatch(do.call(ctlFun, .sirEvalOverrides), error = function(e) NULL)

  # Fail closed on both counts. Either fallback -- raw unnormalised overrides,
  # or the method's bare defaults -- can put the evaluator on a different
  # surface from the one fit$objf came from, which is precisely the failure the
  # whole-control carry exists to prevent. Refusing to run is the only safe
  # answer: a wrong dOFV is invisible, a refused run is not.
  if (is.null(ref)) {
    cli::cli_abort(c(
      "Could not build the evaluation control for {.val {est}}.",
      "i" = "SIR must score candidates on the same likelihood that produced {.code fit$objf}, and cannot verify that it would.",
      "i" = "This is refused rather than approximated, because a candidate scored on the wrong surface is not detectable from the result."
    ))
  }
  if (is.null(ctl)) {
    cli::cli_abort(c(
      "{.arg fit} carries no usable control object.",
      "i" = "Its estimation settings are what the evaluator reproduces, so they cannot be guessed at.",
      "i" = "Refit, or supply a fit whose {.code fit$control} is intact."
    ))
  }

  for (nm in names(.sirEvalOverrides)) {
    # ctl[nm] <- list(v), not ctl[[nm]] <- v: the latter DELETES the element
    # when v is NULL, and the normalised `print` is NULL.
    ctl[nm] <- list(ref[[nm]])
  }
  ctl
}

sirEvalOFV <- function(fit, paramSamples, workers = NULL, rxThreads = NULL) {
  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertMatrix(
    paramSamples,
    mode = "numeric",
    min.rows = 1L,
    min.cols = 1L
  )

  ps <- .sirParamSpace(fit)
  col_names <- colnames(paramSamples)
  theta_cols <- intersect(
    col_names,
    ps$sirName[ps$kind %in% c("theta", "sigma")]
  )
  omega_cols <- intersect(
    col_names,
    ps$sirName[ps$kind %in% c("omegaDiag", "omegaOffdiag")]
  )
  has_omega <- length(omega_cols) > 0L

  if (length(theta_cols) == 0L && !has_omega) {
    cli::cli_abort(
      "No column names in {.arg paramSamples} match any parameter in {.arg fit}."
    )
  }

  base_omega <- fit$omega
  base_theta <- fit$theta
  # Built once: it is the same objective for every candidate, and that is the
  # whole point.
  evalControl <- .sirEvalControl(fit)
  evalEst <- .sirEvalMethod(fit)

  eval_one <- function(i) {
    row <- paramSamples[i, ]

    theta_vals <- base_theta
    if (length(theta_cols) > 0L) {
      theta_vals[theta_cols] <- row[theta_cols]
    }

    ini_args <- as.list(theta_vals)

    if (has_omega) {
      omega_mat <- .sirReconstructOmega(ps, row, base_omega)
      eta_names <- rownames(omega_mat)
      lt_vals <- unlist(
        lapply(seq_len(nrow(omega_mat)), function(r) {
          omega_mat[r, seq_len(r)]
        }),
        use.names = FALSE
      )
      lhs <- paste(eta_names, collapse = " + ")
      lt_txt <- format(lt_vals, scientific = TRUE, digits = 17, trim = TRUE)
      rhs <- if (length(lt_vals) == 1L) {
        lt_txt
      } else {
        paste0("c(", paste(lt_txt, collapse = ", "), ")")
      }
      omega_expr <- str2lang(paste(lhs, "~", rhs))
      ini_args <- c(ini_args, list(omega_expr))
    }

    tryCatch(
      {
        model_new <- suppressMessages(
          do.call(rxode2::ini, c(list(x = fit), ini_args))
        )
        f <- suppressMessages(
          nlmixr2est::nlmixr2(
            model_new,
            est = evalEst,
            control = evalControl
          )
        )
        list(objf = f$objf, error = NA_character_)
      },
      # The message is kept rather than discarded: a configuration problem
      # (a missing import, an unloadable model) fails every sample
      # identically, and reporting only "all evaluations failed" hides why.
      error = function(e) {
        list(objf = NA_real_, error = conditionMessage(e))
      }
    )
  }

  results <- nlmixr2utils::.withWorkerPlan(
    workers,
    rxThreads = nlmixr2utils::resolveRxThreads(workers, rxThreads),
    {
      # nolint: object_usage_linter.
      nlmixr2utils::.plap(
        # nolint: object_usage_linter.
        seq_len(nrow(paramSamples)),
        eval_one,
        .label = function(i) sprintf("sample %d", i)
      )
    }
  )

  objf <- vapply(results, function(r) r$objf, numeric(1L))
  errs <- unique(stats::na.omit(vapply(
    results,
    function(r) r$error,
    character(1L)
  )))
  if (length(errs) > 0L) {
    attr(objf, "evalErrors") <- as.character(errs)
  }
  objf
}
