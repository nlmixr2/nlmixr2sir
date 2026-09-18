# Alternative sources for the initial proposal covariance.
#
# `runSIR()` used to require a successful covariance step. PsN offers four ways
# in, and the three it has that we lacked exist precisely for the models where
# the covariance step fails -- which is where SIR is most wanted. Dispatch
# order follows PsN sir.pm: covmatInput, then rawresInput, then rse*, then
# fit$cov.

# Variance of each parameter's uncertainty, derived from RSE percentages.
#
# Port of PsN setup_variancevec_from_rse() (lib/tool/sir.pm:251).
#
#   diagonal / THETA : variance = (rse * estimate / 100)^2
#   OMEGA off-diag   : N        = (100/rse_i)^2 + (100/rse_j)^2 + 1
#                      variance = (cov_ij^2 + var_i * var_j) / N
#
# Note that PsN's *documentation* describes the off-diagonal rule as choosing
# the variance "so that the correlation from the final estimate is unchanged".
# Its implementation is the Wishart-style expression above, which is not the
# same thing. This follows the implementation, since that is what PsN actually
# runs.
.sirRseVariance <- function(
  ps,
  rseTheta = NULL,
  rseOmega = NULL,
  rseSigma = NULL
) {
  rse <- .sirResolveRse(rseTheta, rseOmega, rseSigma)

  out <- stats::setNames(rep(NA_real_, nrow(ps)), ps$sirName)
  perParam <- stats::setNames(rep(NA_real_, nrow(ps)), ps$sirName)

  .fillClass <- function(kind, given, label, what, explicit) {
    idx <- which(ps$kind == kind)
    if (length(idx) == 0L) {
      # Complain only about a value the caller actually supplied. A scalar
      # rseTheta fills in rseOmega/rseSigma automatically, and it must not be
      # an error for that fill-in to land on a class the model does not have.
      if (!is.null(given) && isTRUE(explicit)) {
        cli::cli_abort(
          "{.arg {label}} was given, but the model has no estimated {what}."
        )
      }
      return(invisible(NULL))
    }
    if (is.null(given)) {
      cli::cli_abort(c(
        "{.arg {label}} is required when building the proposal from RSE.",
        "i" = "The model has {length(idx)} estimated {what}."
      ))
    }
    checkmate::assertNumeric(
      given,
      lower = .Machine$double.eps,
      finite = TRUE,
      any.missing = FALSE,
      min.len = 1L,
      .var.name = label
    )
    if (length(given) != 1L && length(given) != length(idx)) {
      cli::cli_abort(c(
        "{.arg {label}} must be length 1 or one value per {what}.",
        "i" = "The model has {length(idx)} estimated {what}, but {length(given)} value{?s} {?was/were} given."
      ))
    }
    perParam[idx] <<- given
    out[idx] <<- (given * ps$est[idx] / 100)^2
    invisible(NULL)
  }

  .fillClass("theta", rse$theta, "rseTheta", "THETA", TRUE)
  .fillClass(
    "sigma",
    rse$sigma,
    "rseSigma",
    "residual error parameters",
    rse$explicit[["sigma"]]
  )
  .fillClass(
    "omegaDiag",
    rse$omega,
    "rseOmega",
    "OMEGA diagonals",
    rse$explicit[["omega"]]
  )

  off <- which(ps$kind == "omegaOffdiag")
  if (length(off) > 0L) {
    diagIdx <- which(ps$kind == "omegaDiag")
    rseByEta <- stats::setNames(
      perParam[diagIdx],
      as.character(ps$neta1[diagIdx])
    )
    varByEta <- stats::setNames(
      ps$est[diagIdx],
      as.character(ps$neta1[diagIdx])
    )
    for (i in off) {
      ei <- as.character(ps$neta1[i])
      ej <- as.character(ps$neta2[i])
      varI <- varByEta[[ei]]
      varJ <- varByEta[[ej]]
      if (is.na(varI) || is.na(varJ) || varI <= 0 || varJ <= 0) {
        cli::cli_abort(c(
          "Cannot derive an RSE variance for {.val {ps$sirName[i]}}.",
          "i" = "Both connected OMEGA diagonals must be estimated and positive."
        ))
      }
      n <- (100 / rseByEta[[ei]])^2 + (100 / rseByEta[[ej]])^2 + 1
      out[i] <- (ps$est[i]^2 + varI * varJ) / n
    }
  }

  if (anyNA(out)) {
    cli::cli_abort(
      "No RSE supplied for SIR parameter{?s} {.val {ps$sirName[is.na(out)]}}."
    )
  }
  out
}

# PsN's defaulting rules for the rse_* trio (bin/sir -rse_omega help text):
# a scalar rseTheta fills in an unset rseOmega/rseSigma; a vector rseTheta does
# not, and leaving the others unset is then an error; and setting rseOmega
# without rseTheta is an error.
.sirResolveRse <- function(rseTheta, rseOmega, rseSigma) {
  if (is.null(rseTheta) && (!is.null(rseOmega) || !is.null(rseSigma))) {
    cli::cli_abort(c(
      "{.arg rseOmega}/{.arg rseSigma} were given without {.arg rseTheta}.",
      "i" = "The RSE proposal is driven by {.arg rseTheta}; set it too."
    ))
  }
  if (is.null(rseTheta)) {
    return(list(
      theta = NULL,
      omega = NULL,
      sigma = NULL,
      explicit = c(omega = FALSE, sigma = FALSE)
    ))
  }
  explicit <- c(omega = !is.null(rseOmega), sigma = !is.null(rseSigma))
  if (length(rseTheta) == 1L) {
    if (is.null(rseOmega)) {
      rseOmega <- rseTheta
    }
    if (is.null(rseSigma)) {
      rseSigma <- rseTheta
    }
  }
  list(
    theta = rseTheta,
    omega = rseOmega,
    sigma = rseSigma,
    explicit = explicit
  )
}

# Diagonal proposal covariance from RSE percentages.
.sirProposalFromRse <- function(
  ps,
  rseTheta = NULL,
  rseOmega = NULL,
  rseSigma = NULL
) {
  v <- .sirRseVariance(
    ps,
    rseTheta = rseTheta,
    rseOmega = rseOmega,
    rseSigma = rseSigma
  )
  out <- diag(unname(v), nrow = length(v))
  dimnames(out) <- list(ps$sirName, ps$sirName)
  out
}

# Proposal covariance supplied directly: a matrix, a NONMEM-style .cov file, or
# the literal "identity". PsN's -covmat_input. "identity" plus inflation is the
# cheap "any diagonal proposal" route.
.sirProposalFromCovmatInput <- function(ps, covmatInput) {
  if (identical(covmatInput, "identity")) {
    out <- diag(1, nrow = nrow(ps))
    dimnames(out) <- list(ps$sirName, ps$sirName)
    return(out)
  }
  if (is.character(covmatInput)) {
    checkmate::assertFileExists(covmatInput, .var.name = "covmatInput")
    covmatInput <- .sirReadCovFile(covmatInput)
  }
  checkmate::assertMatrix(
    covmatInput,
    mode = "numeric",
    any.missing = FALSE,
    .var.name = "covmatInput"
  )
  if (nrow(covmatInput) != ncol(covmatInput)) {
    cli::cli_abort("{.arg covmatInput} must be square.")
  }
  nms <- rownames(covmatInput)
  if (is.null(nms)) {
    if (nrow(covmatInput) != nrow(ps)) {
      cli::cli_abort(c(
        "Unnamed {.arg covmatInput} must be {nrow(ps)}x{nrow(ps)}.",
        "i" = "It is {nrow(covmatInput)}x{ncol(covmatInput)}.",
        "i" = "Name its rows and columns to supply a subset instead."
      ))
    }
    dimnames(covmatInput) <- list(ps$sirName, ps$sirName)
    return(covmatInput)
  }
  # Named: accept either SIR names or fit$cov names, and reorder.
  idx <- match(nms, ps$sirName)
  if (anyNA(idx)) {
    idx <- match(nms, ps$covName)
  }
  if (anyNA(idx)) {
    cli::cli_abort(c(
      "{.arg covmatInput} names do not match the SIR parameters.",
      "x" = "Unmatched: {.val {nms[is.na(idx)]}}.",
      "i" = "Expected names from: {.val {ps$sirName}}."
    ))
  }
  out <- matrix(
    0,
    nrow = nrow(ps),
    ncol = nrow(ps),
    dimnames = list(ps$sirName, ps$sirName)
  )
  out[idx, idx] <- covmatInput
  out
}

# Minimal reader for a NONMEM-style .cov table: whitespace-delimited, a leading
# NAME column, one row per parameter.
.sirReadCovFile <- function(path) {
  raw <- utils::read.table(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    comment.char = "",
    check.names = FALSE
  )
  nameCol <- which(tolower(names(raw)) %in% c("name", "names"))
  if (length(nameCol) == 1L) {
    nms <- raw[[nameCol]]
    raw <- raw[, -nameCol, drop = FALSE]
  } else {
    nms <- names(raw)
  }
  out <- as.matrix(raw)
  if (nrow(out) != ncol(out)) {
    cli::cli_abort(c(
      "{.file {path}} is not a square covariance table.",
      "i" = "Read {nrow(out)} rows and {ncol(out)} numeric columns."
    ))
  }
  dimnames(out) <- list(as.character(nms), as.character(nms))
  storage.mode(out) <- "double"
  out
}

# Parameter vectors read from a canonical raw-results file, as a matrix in SIR
# parameter order. PsN's -rawres_input / -offset_rawres / -in_filter.
#
# All three pieces already exist in nlmixr2utils, and the rawName column of
# .sirParamSpace() is exactly the naming convention they speak, so this is
# wiring rather than new numerics. Any canonical raw-results file works,
# including one written by nlmixr2boot.
.sirRawResultsMatrix <- function(
  fit,
  ps,
  rawresInput,
  offsetRawres = 1L,
  inFilter = NULL
) {
  parsed <- nlmixr2utils::parseRawResultsParams(
    rawresInput,
    fit,
    offset = offsetRawres,
    filter = inFilter
  )
  if (length(parsed) == 0L) {
    cli::cli_abort(c(
      "No raw-results rows survived {.arg offsetRawres} and {.arg inFilter}.",
      "i" = "Nothing is left to build a proposal from."
    ))
  }

  out <- matrix(
    NA_real_,
    nrow = length(parsed),
    ncol = nrow(ps),
    dimnames = list(NULL, ps$sirName)
  )
  # THETA and residual error arrive as a named vector keyed by parameter name;
  # OMEGA arrives as a matrix, so it is read by (neta1, neta2) rather than by
  # name.
  flatIdx <- which(ps$kind %in% c("theta", "sigma"))
  omegaIdx <- which(ps$kind %in% c("omegaDiag", "omegaOffdiag"))

  for (i in seq_along(parsed)) {
    vals <- rep(NA_real_, nrow(ps))

    if (length(flatIdx) > 0L) {
      flat <- c(parsed[[i]]$theta, parsed[[i]]$sigma)
      j <- match(ps$rawName[flatIdx], names(flat))
      if (anyNA(j)) {
        cli::cli_abort(c(
          "Raw-results sample {parsed[[i]]$sample} is missing parameter{?s} {.val {ps$sirName[flatIdx][is.na(j)]}}.",
          "i" = "Expected raw-results column{?s} {.val {ps$rawName[flatIdx][is.na(j)]}}."
        ))
      }
      vals[flatIdx] <- unname(flat[j])
    }

    if (length(omegaIdx) > 0L) {
      om <- parsed[[i]]$omega
      need <- max(ps$neta1[omegaIdx], ps$neta2[omegaIdx])
      if (!is.matrix(om) || nrow(om) < need) {
        cli::cli_abort(c(
          "Raw-results sample {parsed[[i]]$sample} has no usable OMEGA matrix.",
          "i" = "Need at least {need} eta{?s}, got {if (is.matrix(om)) nrow(om) else 0}."
        ))
      }
      vals[omegaIdx] <- om[cbind(ps$neta1[omegaIdx], ps$neta2[omegaIdx])]
    }

    out[i, ] <- vals
  }
  if (anyNA(out)) {
    cli::cli_abort("Raw-results input produced missing parameter values.")
  }
  out
}

# Derive the iteration-1 proposal from raw-results vectors, the way PsN's
# iteration 0 does: take the empirical mean and covariance of the supplied
# vectors, Box-Cox transformed if requested.
.sirProposalFromRawResults <- function(
  fit,
  ps,
  rawresInput,
  offsetRawres = 1L,
  inFilter = NULL,
  boxcox = TRUE,
  capCorrelation = 0.8
) {
  mat <- .sirRawResultsMatrix(
    fit,
    ps,
    rawresInput = rawresInput,
    offsetRawres = offsetRawres,
    inFilter = inFilter
  )
  if (nrow(mat) < 2L) {
    cli::cli_abort(c(
      "At least two raw-results vectors are needed to form a proposal covariance.",
      "i" = "{nrow(mat)} row{?s} survived {.arg offsetRawres} and {.arg inFilter}."
    ))
  }
  # PsN validates raw-results rank before using it, and so do we. Naming the
  # raw-results source here is more use than the generic message
  # sirUpdateProposal() would raise a moment later.
  .sirCheckProposalRank(mat, what = "raw-results")
  updated <- sirUpdateProposal(
    mat,
    boxcox = boxcox,
    capCorrelation = capCorrelation
  )
  list(
    mu = stats::setNames(colMeans(mat), colnames(mat)),
    covMat = updated$covMat,
    boxcoxState = updated$boxcoxParams,
    nVectors = nrow(mat)
  )
}


# Pick the initial proposal covariance, following PsN's dispatch order. Returns
# the raw source matrix; inflation, correlation capping and the
# positive-definiteness fix are applied downstream by .sirInitialProposal().
.sirResolveInitialProposal <- function(fit, ps, control, seedCov = NULL) {
  if (!is.null(control$covmatInput)) {
    return(list(
      covMat = .sirProposalFromCovmatInput(ps, control$covmatInput),
      source = "covmatInput"
    ))
  }
  if (!is.null(control$rawresInput)) {
    got <- .sirProposalFromRawResults(
      fit,
      ps,
      rawresInput = control$rawresInput,
      offsetRawres = control$offsetRawres,
      inFilter = control$inFilter,
      boxcox = control$boxcox,
      capCorrelation = control$capCorrelation
    )
    return(list(
      covMat = got$covMat,
      mu = got$mu,
      boxcoxState = got$boxcoxState,
      source = "rawresInput"
    ))
  }
  if (!is.null(control$rseTheta)) {
    return(list(
      covMat = .sirProposalFromRse(
        ps,
        rseTheta = control$rseTheta,
        rseOmega = control$rseOmega,
        rseSigma = control$rseSigma
      ),
      source = "rse"
    ))
  }
  if (!is.null(seedCov)) {
    # A seed chosen by setCov(fit, "sir") that is not the installed fit$cov.
    return(list(
      covMat = sirGetProposalCov(
        fit,
        capCorrelation = control$capCorrelation,
        cov = seedCov
      ),
      source = "seedCov"
    ))
  }
  list(
    covMat = sirGetProposalCov(fit, capCorrelation = control$capCorrelation),
    source = "cov"
  )
}
