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
# Fields of the fitted control that define or materially affect the value of
# the objective, as opposed to how it was optimised. Optimisation-only settings
# need not be carried when maxOuterIterations = 0, but these do: they change
# the number the evaluator returns.
#
#   interaction       FOCEi versus FOCE
#   fo                first-order: expand at eta = 0 rather than at each
#                     subject's conditional mode. This one is easy to miss --
#                     an `fo` fit shows the same interaction/nAGQ/foce as a
#                     `foce` fit and differs only here, yet scores 127.98
#                     against 116.80 on theo_sd. Omitting it would have scored
#                     every candidate on the wrong surface silently.
#   nAGQ              quadrature nodes: 0 is FOCEi, 1 is Laplace, >= 2 is AGQ
#   foce              residual-variance convention ("nonmem" vs "foce+")
#   muModel           mu-referencing regression variant (the m.../i... methods)
#   addProp           how additive and proportional error combine
#   adjLik            likelihood constant adjustment
#   badSolveObjfAdj   the penalty applied to a failed solve -- candidate
#                     dependent, so a mismatch changes the shape of the target
#   rxControl         ODE solver method and tolerances
#   sumProd, optExpression, literalFix, sigdig
#                     expression handling and derived tolerances
#
# Carrying these is what lets one FOCEi-family evaluator reproduce the whole
# deterministic ladder: every one of fo/foi/foce/focei/focep/laplace/agq runs
# on this engine and differs only in these settings. The preflight is still the
# gate -- it re-evaluates at the fitted centre and aborts on a mismatch -- so a
# field missed here shows up as a refused run rather than a wrong answer.
.sirLikelihoodControlFields <- c(
  "interaction", "fo", "nAGQ", "foce", "muModel",
  "addProp", "adjLik", "badSolveObjfAdj", "rxControl",
  "sumProd", "optExpression", "literalFix", "sigdig"
)

# Build the fixed-parameter evaluator's control by carrying the fitted model's
# likelihood-relevant settings forward, rather than accepting foceiControl()'s
# defaults for all of them. A fresh default control is a different objective
# whenever the fit used anything but the defaults.
# The method a candidate is scored with, and the constructor for its control.
#
# Carrying the control fields is not on its own enough: the `est` string itself
# selects the objective. An `fo` fit evaluated as `est = "focei"` scored
# 103.870 against its own 127.982 -- 24 units out, and *below* FOCEi's own
# minimum, because the etas were being estimated rather than held at zero.
# Evaluating it as `est = "fo"` reproduces 127.98223 exactly.
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

.sirEvalMethod <- function(fit) {
  est <- .sirFitEst(fit)
  if (is.na(est)) "focei" else est
}

.sirEvalControlFun <- function(est) {
  tryCatch(
    getExportedValue("nlmixr2est", paste0(est, "Control")),
    error = function(e) nlmixr2est::foceiControl
  )
}

.sirEvalControl <- function(fit) {
  est <- .sirEvalMethod(fit)
  ctlFun <- .sirEvalControlFun(est)
  base <- fit$control
  args <- if (is.list(base)) {
    keep <- intersect(.sirLikelihoodControlFields, names(base))
    as.list(base)[keep]
  } else {
    list()
  }
  # A thin wrapper such as foceControl() or laplaceControl() forces the args
  # that define its rung and may not accept every FOCEi field by name, so drop
  # anything it cannot take unless it forwards through `...`.
  fml <- names(formals(ctlFun))
  if (!("..." %in% fml)) {
    args <- args[intersect(names(args), fml)]
  }
  # Evaluation-only overrides. These control what the call does, not what the
  # objective means, so they are always ours.
  args$calcTables <- FALSE
  args$covMethod <- ""
  args$compress <- FALSE
  args$maxOuterIterations <- 0L
  args$print <- 0L
  tryCatch(
    do.call(ctlFun, args),
    error = function(e) {
      cli::cli_abort(c(
        "Could not reconstruct the fit's objective settings for evaluation.",
        "x" = conditionMessage(e),
        "i" = "Estimation method: {.val {est}}.",
        "i" = "SIR must score candidates on the same likelihood that produced {.code fit$objf}."
      ))
    }
  )
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
