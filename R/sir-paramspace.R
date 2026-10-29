# Single source of truth for the SIR parameter vector.
#
# Three naming conventions meet in this package, and deriving them
# independently in each consumer is what let the code paths drift apart:
#
#   kind          sirName        covName (nlmixr2est)  rawName (nlmixr2utils)
#   theta         tka            tka                   tka
#   sigma         add.sd         add.sd                add.sd
#   omegaDiag     eta.ka         om.eta.ka             omega(eta.ka,eta.ka)
#   omegaOffdiag  eta.ka:eta.cl  cov.eta.ka.eta.cl     omega(eta.ka,eta.cl)
#
# `.sirParamSpace()` returns one row per estimated, non-fixed parameter
# carrying all three names plus everything needed to place the parameter in a
# proposal, bound it, and rebuild an OMEGA matrix from it.

# Names nlmixr2est gives OMEGA entries in fit$cov, built the same way
# nlmixr2est:::.foceiOmegaCovNames() builds them. Derived from the eta names
# rather than parsed back out of the rownames, because eta names contain dots
# and "cov.eta.ka.eta.cl" cannot be split unambiguously.
.sirOmegaCovName <- function(n1, n2) {
  ifelse(n1 == n2, paste0("om.", n1), paste0("cov.", n1, ".", n2))
}

.sirOmegaSirName <- function(n1, n2) {
  ifelse(n1 == n2, n1, paste(n1, n2, sep = ":"))
}

.sirOmegaRawName <- function(n1, n2) {
  paste0("omega(", n1, ",", n2, ")")
}

#' Describe every estimated parameter SIR works with
#'
#' @param fit An nlmixr2 fit.
#' @return A data frame with one row per estimated, non-fixed parameter and
#'   columns `sirName`, `covName`, `rawName`, `kind`, `ntheta`, `neta1`,
#'   `neta2`, `est`, `lower`, `upper` and `fullCovName`. `covName` is `NA`
#'   when `fit$cov` does not carry that parameter, which is what drives the fallback choice
#'   in `.sirInitialProposal()`. Rows are ordered THETA, sigma, then OMEGA
#'   lower-triangle by column then row.
#' @noRd
.sirParamSpace <- function(fit) {
  # Deliberately no class assertion here: the public entry points already
  # check, and nlmixr2FitCore's custom `$` method makes a lightweight list
  # stand-in impossible to use, which the fixed-parameter tests need.
  ini_df <- fit$iniDf
  fixed <- !is.na(ini_df$fix) & ini_df$fix
  cov_names <- rownames(fit$cov)
  if (is.null(cov_names)) {
    cov_names <- character(0)
  }

  # THETA-side rows. "sigma" is the residual-error subset, identified
  # structurally by iniDf$err rather than by absence from fit$cov -- the
  # latter silently changed meaning when nlmixr2est switched to covFull.
  theta_rows <- ini_df[!is.na(ini_df$ntheta) & !fixed, , drop = FALSE]
  theta_rows <- theta_rows[order(theta_rows$ntheta), , drop = FALSE]
  is_sigma <- !is.na(theta_rows$err)
  theta_part <- data.frame(
    sirName = theta_rows$name,
    covName = theta_rows$name,
    rawName = theta_rows$name,
    kind = ifelse(is_sigma, "sigma", "theta"),
    ntheta = theta_rows$ntheta,
    neta1 = NA_integer_,
    neta2 = NA_integer_,
    est = theta_rows$est,
    lower = ifelse(is.na(theta_rows$lower), -Inf, theta_rows$lower),
    upper = ifelse(is.na(theta_rows$upper), Inf, theta_rows$upper),
    stringsAsFactors = FALSE
  )
  theta_part <- theta_part[
    order(theta_part$kind != "theta", theta_part$ntheta),
    ,
    drop = FALSE
  ]

  # OMEGA-side rows: the free lower triangle, ordered by column then row.
  omega_rows <- ini_df[!is.na(ini_df$neta1) & !fixed, , drop = FALSE]
  omega_rows <- omega_rows[omega_rows$neta1 >= omega_rows$neta2, , drop = FALSE]
  omega_rows <- omega_rows[
    order(omega_rows$neta2, omega_rows$neta1),
    ,
    drop = FALSE
  ]

  if (nrow(omega_rows) == 0L) {
    out <- theta_part
  } else {
    diag_rows <- ini_df[
      !is.na(ini_df$neta1) & ini_df$neta1 == ini_df$neta2,
      ,
      drop = FALSE
    ]
    etaName <- stats::setNames(diag_rows$name, as.character(diag_rows$neta1))
    n1 <- unname(etaName[as.character(omega_rows$neta1)])
    n2 <- unname(etaName[as.character(omega_rows$neta2)])

    omega_part <- data.frame(
      sirName = .sirOmegaSirName(n1, n2),
      covName = .sirOmegaCovName(n1, n2),
      rawName = .sirOmegaRawName(n1, n2),
      kind = ifelse(
        omega_rows$neta1 == omega_rows$neta2,
        "omegaDiag",
        "omegaOffdiag"
      ),
      ntheta = NA_integer_,
      neta1 = omega_rows$neta1,
      neta2 = omega_rows$neta2,
      est = unname(fit$omega[cbind(omega_rows$neta1, omega_rows$neta2)]),
      lower = -Inf,
      upper = Inf,
      stringsAsFactors = FALSE
    )
    # A variance cannot be negative; a covariance can.
    omega_part$lower[omega_part$kind == "omegaDiag"] <- 0
    out <- rbind(theta_part, omega_part)
  }

  # fullCovName is what nlmixr2est calls the parameter in a full-shape
  # covariance, whether or not fit$cov carries it. It names a seed covariance
  # that is not fit$cov, and the SIR covariance installed with setCov().
  out$fullCovName <- out$covName
  # covName is only a claim about fit$cov if fit$cov actually carries it.
  out$covName[!(out$covName %in% cov_names)] <- NA_character_
  rownames(out) <- NULL
  out
}

# Proposal mean, in SIR parameter-vector order. Replaces the
# `fit$theta[rownames(fit$cov)]` idiom, which returned NA for every OMEGA row
# once nlmixr2est started reporting OMEGA in fit$cov.
.sirProposalMu <- function(fit, ps = .sirParamSpace(fit)) {
  stats::setNames(ps$est, ps$sirName)
}

# TRUE when fit$cov carries the OMEGA block, i.e. the Wishart fallback is not
# needed. nlmixr2est does this whenever foceiControl(covFull = TRUE), which
# has been the default since nlmixr2est 7.
.sirCovHasOmega <- function(ps) {
  omega <- ps$kind %in% c("omegaDiag", "omegaOffdiag")
  any(omega) && !anyNA(ps$covName[omega])
}
