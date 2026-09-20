# Registering the SIR covariance with the fit, so that
# nlmixr2est::setCov(fit, "sir") switches the fit's reported uncertainty to the
# SIR result. Mirrors nlmixr2boot's .registerBootCovList() /
# .bootstrapCovAsFitCov().
#
# nlmixr2boot has to guess how its parameter names correspond to fit$cov
# rownames, reconstructing "cov.<eta1>.<eta2>" and "om.<eta>" by pattern. That
# guesswork is unnecessary here: .sirParamSpace() already carries the exact
# sirName <-> covName correspondence, so the mapping is a lookup.

# Reshape the SIR covariance into fit$cov's names and order, or return NULL
# with an explanation if it cannot be done safely.
.sirCovAsFitCov <- function(fit, covSir, ps = .sirParamSpace(fit)) {
  fitCov <- fit$cov
  if (!is.matrix(fitCov) || is.null(rownames(fitCov))) {
    cli::cli_inform(c(
      "i" = "No {.code fit$cov} to match; skipping {.fn setCov} registration."
    ))
    return(NULL)
  }
  if (!is.matrix(covSir) || is.null(rownames(covSir))) {
    return(NULL)
  }

  fitNames <- rownames(fitCov)
  mapped <- ps$covName[match(rownames(covSir), ps$sirName)]

  if (anyNA(mapped) || !setequal(mapped, fitNames)) {
    cli::cli_inform(c(
      "i" = "SIR parameters do not match {.code fit$cov} exactly; skipping {.fn setCov} registration.",
      "i" = "This is expected when SIR ran on a subset, or on a fit whose covariance step failed."
    ))
    return(NULL)
  }

  out <- covSir
  dimnames(out) <- list(mapped, mapped)
  out <- out[fitNames, fitNames, drop = FALSE]

  if (inherits(try(chol(out), silent = TRUE), "try-error")) {
    cli::cli_inform(c(
      "i" = "The SIR covariance is not positive definite; skipping {.fn setCov} registration."
    ))
    return(NULL)
  }
  out
}

# Merge a covariance into fit$env$covList under `label`, leaving any other
# registered matrices in place.
.sirRegisterCovList <- function(fit, label, covMat) {
  if (is.null(covMat)) {
    return(invisible(FALSE))
  }
  env <- tryCatch(fit$env, error = function(e) NULL)
  if (!is.environment(env)) {
    return(invisible(FALSE))
  }
  existing <- if (exists("covList", envir = env, inherits = FALSE)) {
    get("covList", envir = env)
  } else {
    list()
  }
  existing[[label]] <- covMat
  assign("covList", existing, envir = env)
  invisible(TRUE)
}

# Called at the end of a run: register the empirical SIR covariance as "sir".
.sirRegisterCov <- function(
  fit,
  summary,
  ps = .sirParamSpace(fit),
  label = "sir"
) {
  covMat <- attr(summary, "covMatrix", exact = TRUE)
  if (is.null(covMat)) {
    return(invisible(FALSE))
  }
  .sirRegisterCovList(fit, label, .sirCovAsFitCov(fit, covMat, ps))
}

#' Use a SIR covariance as a fit's covariance
#'
#' Registers SIR as a covariance method for [nlmixr2est::setCov()], so that
#' `setCov(fit, "sir")` switches a fit's reported uncertainty -- its standard
#' errors, RSEs and print output -- from the asymptotic covariance to the
#' empirical one SIR produced. This is the point of running SIR in the first
#' place: the asymptotic covariance is the thing SIR exists to improve on.
#'
#' [runSIR()] registers its covariance on the fit automatically when it
#' finishes, so the usual sequence is
#'
#' ```
#' sir <- runSIR(fit, ...)
#' nlmixr2est::setCov(fit, "sir")
#' ```
#'
#' and `nlmixr2est::setCov(fit, fit$covMethod)` puts the original back, since
#' `setCov()` keeps the previous covariance in `fit$covList`.
#'
#' This method does **not** run SIR. Computing a SIR covariance means a full
#' SIR run -- minutes to hours, needing a sampling schedule and somewhere to
#' write -- and it produces convergence and weight diagnostics that a
#' covariance setter would throw away. A setter should not silently start that,
#' so a fit with no SIR covariance is refused with a pointer to [runSIR()],
#' which is what nlmixr2est asks a method that cannot compute its covariance to
#' do.
#'
#' Requires nlmixr2est >= 7.1.0, which is where `setCov()` became a generic.
#' On earlier versions `setCov(fit, "sir")` still installs a covariance
#' [runSIR()] has already registered, because that path does not dispatch, but
#' the method is neither listed by `nlmixr2est::setCovAllMethods()` nor able to
#' produce this message.
#'
#' @param fit An nlmixr2 fit that has been through [runSIR()].
#' @param method Covariance method, supplied by [nlmixr2est::setCov()].
#' @param ... Unused; present for consistency with the generic.
#' @return The SIR covariance, named and ordered like `fit$cov`, for
#'   `setCov()` to install.
#' @exportS3Method nlmixr2est::setCov
setCov.sir <- function(fit, method, ...) {
  # Normally unreachable when a covariance IS registered: setCov() installs a
  # cached one before it dispatches. Looked up here anyway so the method is
  # correct on its own terms rather than only in combination with that path.
  env <- tryCatch(fit$env, error = function(e) NULL)
  covMat <- if (is.environment(env) &&
                exists("covList", envir = env, inherits = FALSE)) {
    get("covList", envir = env)[["sir"]]
  } else {
    NULL
  }

  if (is.null(covMat)) {
    cli::cli_abort(c(
      "This fit has no SIR covariance to install.",
      "i" = "Run {.run runSIR(fit)} first; it registers its covariance on the fit when it finishes.",
      "i" = "{.fn setCov} does not run SIR itself: a SIR run needs a sampling schedule and produces convergence and weight diagnostics that installing a covariance would discard."
    ), call. = NULL)
  }
  covMat
}
