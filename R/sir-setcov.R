# SIR as an nlmixr2est covariance method.
#
# nlmixr2est::setCov() is an S3 generic on the covariance-method name, and
# caches every covariance it computes together with the options it was computed
# with. This file plugs SIR into that:
#
#   setCov(fit, "sir", control = sirControl())   computes (setCov.sir)
#   setCov(fit) <- runSIR(fit, ...)              installs a finished run
#                                                (setCovValue.nlmixr2SIR)
#   runSIR()                                     registers its result in
#                                                fit$covList (.sirRegisterCov)
#
# SIR starts from a seed covariance, so the options it is cached under are the
# sirControl() options that shape the result *plus* the seed
# (setCovOptions.sirControl). A changed seed recomputes rather than reinstalling
# a SIR covariance built on a different starting point. Every route keeps the
# SIR result itself on the fit, as fit$sir.

# The fit environment, from a fit or the environment itself.
.sirFitEnv <- function(fit) {
  if (is.environment(fit)) {
    return(fit)
  }
  tryCatch(fit$env, error = function(e) NULL)
}

# Reshape a SIR covariance (SIR names) into nlmixr2est's covariance names, or
# return NULL with an explanation if it cannot be done safely.
#
# The names are the full-shape ones (om.<eta>, cov.<eta1>.<eta2>) whether or not
# fit$cov carries OMEGA, so a fit with a THETA-only covariance, or none, can
# still install the SIR covariance. Parameters fit$cov has keep its order.
.sirCovAsFitCov <- function(fit, covSir, ps = .sirParamSpace(fit)) {
  if (!is.matrix(covSir) || is.null(rownames(covSir))) {
    return(NULL)
  }
  if (!setequal(rownames(covSir), ps$sirName)) {
    cli::cli_inform(c(
      "i" = "SIR parameters do not match the fit's estimated parameters; skipping {.fn setCov} registration."
    ))
    return(NULL)
  }

  mapped <- ps$fullCovName[match(rownames(covSir), ps$sirName)]
  out <- covSir
  dimnames(out) <- list(mapped, mapped)
  fitNames <- rownames(fit$cov)
  ord <- c(intersect(fitNames, mapped), setdiff(ps$fullCovName, fitNames))
  out <- out[ord, ord, drop = FALSE]

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
  env <- .sirFitEnv(fit)
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

# Forget the options recorded for `label`, leaving it "computed with unknown
# options" -- which nlmixr2est reuses for a request that names no control.
.sirDropCovOptions <- function(env, label) {
  rec <- if (exists("covOptions", envir = env, inherits = FALSE)) {
    get("covOptions", envir = env)
  } else {
    NULL
  }
  if (is.null(rec) || is.null(rec[[label]])) {
    return(invisible(FALSE))
  }
  rec[[label]] <- NULL
  assign("covOptions", rec, envir = env)
  invisible(TRUE)
}

# Record the options a covariance on the fit was computed with, as
# nlmixr2est's setCov() does for its own.
.sirSetCovOptions <- function(env, label, options) {
  rec <- if (exists("covOptions", envir = env, inherits = FALSE)) {
    get("covOptions", envir = env)
  } else {
    list()
  }
  rec[[label]] <- options
  assign("covOptions", rec, envir = env)
  invisible(TRUE)
}

# Seed ----------------------------------------------------------------------

# The covariance SIR starts from.
#
# Returns list(cov, method, installed): the seed matrix (NULL for the RSE
# seed), a label for it, and whether it is the covariance installed as
# fit$cov.
.sirResolveSeed <- function(fit, control) {
  env <- .sirFitEnv(fit)
  installedName <- NULL
  installedCov <- NULL
  if (
    is.environment(env) &&
      exists("cov", envir = env, inherits = FALSE) &&
      is.matrix(get("cov", envir = env))
  ) {
    installedCov <- get("cov", envir = env)
    if (exists("covMethod", envir = env, inherits = FALSE)) {
      installedName <- as.character(get("covMethod", envir = env))
    }
  }
  covList <- if (
    is.environment(env) && exists("covList", envir = env, inherits = FALSE)
  ) {
    get("covList", envir = env)
  } else {
    list()
  }

  seedCov <- control$seedCov
  if (is.matrix(seedCov)) {
    return(list(cov = seedCov, method = "matrix", installed = FALSE))
  }
  if (is.character(seedCov)) {
    if (identical(seedCov, installedName)) {
      return(list(cov = installedCov, method = seedCov, installed = TRUE))
    }
    if (!is.null(covList[[seedCov]])) {
      return(list(cov = covList[[seedCov]], method = seedCov, installed = FALSE))
    }
    avail <- setdiff(c(installedName, names(covList)), "sir")
    cli::cli_abort(c(
      "The fit has no covariance {.val {seedCov}} to seed SIR from.",
      "i" = if (length(avail) > 0L) {
        "Available: {.val {avail}}."
      } else {
        "The fit has no other covariance."
      }
    ))
  }

  if (!is.null(installedCov) && !identical(installedName, "sir")) {
    return(list(cov = installedCov, method = installedName, installed = TRUE))
  }
  if (identical(installedName, "sir")) {
    # Reuse the seed the installed SIR covariance was computed from: SIR is
    # never seeded from its own result.
    rec <- env$covOptions[["sir"]]
    seedMethod <- rec$seedMethod
    seedCov <- rec$seedCov
    if (is.null(seedMethod)) {
      # A covariance registered by runSIR() carries no options; the result
      # itself records what the run was seeded from.
      stored <- tryCatch(get("sir", envir = env, inherits = FALSE),
                         error = function(e) NULL)
      seedMethod <- attr(stored, "seedMethod", exact = TRUE)
      seedCov <- attr(stored, "seedCov", exact = TRUE)
    }
    if (!is.null(seedMethod) && !identical(seedMethod, "sir")) {
      return(list(cov = seedCov, method = seedMethod, installed = FALSE))
    }
    cli::cli_abort(c(
      "The installed {.val sir} covariance does not record the covariance it was seeded from.",
      "i" = "Choose a seed with {.code sirControl(seedCov =)}, or install another covariance first."
    ))
  }
  list(cov = NULL, method = "rse", installed = FALSE)
}

# The cache key of a SIR covariance: every option that shapes the result, and
# the seed. Parallelism and the objective preflight do not change the answer,
# so they are left out; the RSE options only matter for the RSE seed.
.sirCovKeyFields <- c(
  "nSamples", "nResample",
  "thetaInflation", "omegaInflation", "sigmaInflation",
  "capCorrelation", "capResampling", "recenter", "boxcox",
  "omegaFallback", "sigmaFallbackRse", "omegaDf", "seed"
)

# The seed as it goes into a cache key. nlmixr2est symmetrizes a covariance
# whenever it installs one, so an estimation-time covariance that is off
# symmetric by rounding (1e-17, say) comes back from fit$covList a few bits
# different. Symmetrizing here too, which is idempotent, keeps the key stable
# across any number of swaps.
.sirSeedKeyCov <- function(cov) {
  if (!is.matrix(cov)) {
    return(cov)
  }
  out <- 0.5 * (cov + t(cov))
  dimnames(out) <- dimnames(cov)
  out
}

.sirCovKey <- function(control, seed) {
  key <- unclass(control)[.sirCovKeyFields]
  rse <- if (identical(seed$method, "rse")) {
    unclass(control)[c("rseTheta", "rseOmega", "rseSigma")]
  } else {
    list(rseTheta = NULL, rseOmega = NULL, rseSigma = NULL)
  }
  c(
    key,
    rse,
    list(seedMethod = seed$method, seedCov = .sirSeedKeyCov(seed$cov))
  )
}

#' @exportS3Method nlmixr2est::setCovOptions
setCovOptions.sirControl <- function(control, fit, ...) {
  .sirCovKey(control, .sirResolveSeed(fit, control))
}

# Run `expr` under a fixed random seed, restoring the caller's stream after.
.sirWithRngSeed <- function(seed, expr) {
  genv <- globalenv()
  hadSeed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  oldSeed <- if (hadSeed) get(".Random.seed", envir = genv) else NULL
  on.exit(
    if (hadSeed) {
      assign(".Random.seed", oldSeed, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    },
    add = TRUE
  )
  set.seed(seed)
  force(expr)
}

# setCov(fit, "sir") ---------------------------------------------------------

#' SIR covariance for `setCov()`
#'
#' `setCov(fit, "sir")` runs sampling importance resampling at the fit's
#' estimates, seeded from the fit's covariance, and installs the resampled
#' covariance as the fit's covariance: the standard errors change, and the
#' previous covariance stays in `fit$covList` so `setCov()` can swap back. The
#' SIR result behind it is kept as `fit$sir`.
#'
#' `setCov(fit) <- sir` installs a result already returned by [runSIR()] the
#' same way.
#'
#' @param fit An nlmixr2 fit.
#' @param method The covariance method, `"sir"`.
#' @param control A [sirControl()] object. Only needed to change the defaults.
#' @param value A [runSIR()] result from this fit.
#' @param ... Ignored.
#' @return `setCov.sir()` installs the covariance and returns `NULL`, as
#'   `nlmixr2est::setCov()` expects of a method that installs itself.
#'   `setCovValue.nlmixr2SIR()` returns what `nlmixr2est::setCov<-` installs.
#' @seealso [sirControl()] for the options and how the result is cached.
#' @examples
#' \dontrun{
#' setCov(fit, "sir")
#' fit$parFixedDf
#' setCov(fit, "r,s") # back to the asymptotic covariance, from the cache
#' setCov(fit, "sir") # and back to SIR, also from the cache
#'
#' sir <- runSIR(fit)
#' setCov(fit) <- sir
#' }
#' @name setCov.sir
#' @exportS3Method nlmixr2est::setCov
setCov.sir <- function(fit, method, control = sirControl(), ...) {
  method <- unclass(method)
  if (!inherits(control, "sirControl")) {
    cli::cli_abort(
      "Covariance method {.val {method}} needs {.arg control} from {.fn sirControl}."
    )
  }
  seed <- .sirResolveSeed(fit, control)
  if (identical(seed$method, "rse")) {
    cli::cli_inform(c(
      "!" = "The fit has no covariance to seed SIR from.",
      "i" = "Seeding it from an assumed {.arg rseTheta} of {.val {control$rseTheta}}%, which is not an estimate of the uncertainty."
    ))
  }
  res <- .sirWithRngSeed(
    control$seed,
    .sirRunCore(
      fit,
      nSamples = control$nSamples,
      nResample = control$nResample,
      directory = NULL,
      fitName = "setCov",
      control = .sirControlToRun(control, seed),
      call = sys.call(),
      # Always the resolved seed, never fit$cov or fit$parFixedDf directly: the
      # same seed is reached both installed and from the record while "sir" is
      # installed, and one cache key must mean one computation.
      seedCov = seed$cov,
      parFixedSe = FALSE,
      register = FALSE
    )
  )
  attr(res, "sirControl") <- control
  attr(res, "covOptions") <- .sirCovKey(control, seed)
  nlmixr2est::`setCov<-`(fit, method, value = res)
  invisible(NULL)
}

# setCov(fit) <- runSIR(...) -------------------------------------------------

# A SIR result can only describe the fit it was run on.
.sirCheckOwnership <- function(fit, result) {
  fp <- attr(result, "fingerprint")
  if (!is.list(fp)) {
    cli::cli_abort(c(
      "This SIR result carries no fingerprint, so it cannot be matched to {.arg fit}.",
      "i" = "Rerun {.fn runSIR} on this fit."
    ))
  }
  ps <- .sirParamSpace(fit)
  current <- .sirRunFingerprint(
    fit,
    ps,
    attr(result, "schedule"),
    attr(result, "control")
  )
  fields <- c("model", "data", "params", "estimates", "objf")
  bad <- fields[!vapply(
    fields,
    function(f) isTRUE(all.equal(fp[[f]], current[[f]])),
    logical(1L)
  )]
  if (length(bad) > 0L) {
    cli::cli_abort(c(
      "This SIR result was not run on {.arg fit}.",
      "x" = "It differs in: {.field {bad}}."
    ))
  }
  invisible(ps)
}

#' @rdname setCov.sir
#' @exportS3Method nlmixr2est::setCovValue
setCovValue.nlmixr2SIR <- function(value, fit, method = NULL, ...) {
  if (is.null(method)) {
    method <- "sir"
  }
  ps <- .sirCheckOwnership(fit, value)
  covMat <- .sirCovAsFitCov(fit, attr(value, "covMatrix", exact = TRUE), ps)
  if (is.null(covMat)) {
    cli::cli_abort("The SIR covariance cannot be installed on {.arg fit}.")
  }
  key <- attr(value, "covOptions", exact = TRUE)
  extra <- stats::setNames(list(value), method)
  if (is.null(key)) {
    # A runSIR() result is not what any sirControl() would compute, so it gets
    # no key: unrecorded means a later plain setCov(fit, method) reinstalls it,
    # while one naming a control recomputes. setCov<- records the options it is
    # given, so the entry is dropped again afterwards through `extra`, which is
    # assigned last.
    extra[["covOptions"]] <- .sirCovOptionsWithout(fit, method)
  }
  list(
    cov = covMat,
    method = method,
    options = key,
    extra = extra
  )
}

# The fit's covOptions with `label` removed, for `extra` to reinstate.
.sirCovOptionsWithout <- function(fit, label) {
  env <- .sirFitEnv(fit)
  rec <- if (is.environment(env) &&
               exists("covOptions", envir = env, inherits = FALSE)) {
    get("covOptions", envir = env)
  } else {
    list()
  }
  rec[[label]] <- NULL
  rec
}

# runSIR() registration -----------------------------------------------------

# Called at the end of a run: make the SIR covariance selectable with
# setCov(fit, "sir"), and keep the result on the fit. When "sir" is the
# installed covariance, it is replaced in place so the fit never holds two.
.sirRegisterCov <- function(
  fit,
  result,
  ps = .sirParamSpace(fit),
  label = "sir"
) {
  env <- .sirFitEnv(fit)
  if (!is.environment(env)) {
    return(invisible(FALSE))
  }
  covMat <- .sirCovAsFitCov(
    fit,
    attr(result, "covMatrix", exact = TRUE),
    ps
  )
  if (is.null(covMat)) {
    return(invisible(FALSE))
  }
  if (identical(env$covMethod, label)) {
    ok <- tryCatch(
      {
        nlmixr2est::`setCov<-`(fit, label, value = result)
        TRUE
      },
      error = function(e) {
        cli::cli_inform(c(
          "i" = "Could not replace the installed {.val {label}} covariance: {conditionMessage(e)}"
        ))
        FALSE
      }
    )
    return(invisible(ok))
  }
  .sirRegisterCovList(fit, label, covMat)
  # Deliberately no recorded options. nlmixr2est reads an unrecorded covariance
  # as "computed with unknown options", which it reinstalls for a plain
  # setCov(fit, "sir") and recomputes for one that names a control. So a fit
  # that has been through runSIR() installs that covariance rather than paying
  # for a fresh run, and asking for particular sirControl() options still gets
  # them.
  .sirDropCovOptions(env, label)
  assign(label, result, envir = env)
  invisible(TRUE)
}
