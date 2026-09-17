# Part of nlmixr2sir.
# Importance-weight degeneracy diagnostics.
#
# Normalising on the log scale with maximum subtraction stops the weights
# overflowing, but successful normalisation says nothing about whether the
# importance sample is informative. A run in which one or two candidates carry
# almost all the weight normalises perfectly well: the retained sample is then
# a handful of repeated points dressed up as `nResample` draws, and every
# interval derived from it is far narrower than the evidence supports.
#
# These are the numbers that make that visible.

# Weights below this fraction of the uniform weight 1/n contribute effectively
# nothing. Used only for the reported count, never to alter the weights.
.sirNegligibleWeightFactor <- 0.01

# Documented warning thresholds. Deliberately loose: they are meant to catch a
# sample that is degenerate rather than merely uneven.
#
# The two ESS thresholds measure DIFFERENT things and have DIFFERENT remedies,
# which is why they are separate warnings rather than one.
#
# .sirEssWarn is an absolute count, and it is the one that bounds what the
# result can support. A retained distribution resting on ~K effectively
# independent points cannot locate its own 2.5th and 97.5th percentiles more
# finely than K draws allow, however many rows `nResample` produced. Absolute
# ESS grows roughly in PROPORTION to the number of samples, so more samples is
# the right remedy here. (Simulated, q = N(0,1) targeting N(1,1): mean ESS of
# 8.2, 22.5, 80.2, 379 and 1842 at n = 16, 50, 200, 1000 and 5000.)
#
# .sirEssFractionWarn is an efficiency ratio, and it is a property of the
# PROPOSAL rather than of the sample size. Asymptotically ESS/n converges to a
# constant fixed by the proposal-target mismatch -- exp(-mu^2) for the normal
# case above -- so drawing more samples does not improve it, and telling a user
# to raise `nSamples` in response to a low ratio is bad advice.
#
# Worse, the ESTIMATE of that ratio is optimistically biased at small n, and
# the bias is largest exactly when the proposal is worst. Same simulation,
# against a true ratio of 0.3679:
#
#     n =   16   ESS/n = 0.514   (+0.146)
#     n =   50           0.450   (+0.082)
#     n =  200           0.401   (+0.033)
#     n = 1000           0.379   (+0.012)
#     n = 5000           0.369   (+0.001)
#
# So a small run reports a flattering efficiency, and raising `nSamples` makes
# the reported ratio FALL as the estimate becomes honest. A user who reads that
# as a regression and responds by drawing more samples again is chasing an
# artefact. The warning says so rather than leaving them to work it out.
# All three are diagnostic CONVENTIONS, not reliability boundaries. Whether a
# given ESS supports a given quantile depends on the weight distribution and
# especially its tail behaviour -- a sample can clear every threshold here and
# still be unreliable if the weights have heavy tails, and a well-behaved one
# below them may be perfectly usable. They are set where they are to catch
# samples that are degenerate rather than merely uneven, and they are quoted in
# the warning so the reader can judge rather than defer.
.sirEssWarn <- 100
.sirEssFractionWarn <- 0.10
.sirMaxWeightWarn <- 0.50

# Below this many samples the efficiency estimate is optimistically biased
# enough to be worth flagging; see the table above. Also a convention -- the
# bias falls off smoothly rather than at a cliff.
.sirEssFractionReliableN <- 200

#' Degeneracy summaries for a set of normalized importance weights
#'
#' @param prob Normalized resampling probabilities, one per candidate. Zero
#'   entries (failed evaluations, non-finite ratios) carry no support and are
#'   dropped rather than counted.
#' @param nSuccessful Number of candidates with a usable OFV, used as the
#'   denominator for `essFraction`.
#' @return A one-row data frame: `ess`, `essFraction`, `maxWeight`,
#'   `perplexity`, `nNonNegligible`.
#' @noRd
.sirWeightDiagnostics <- function(prob, nSuccessful = NULL) {
  p <- prob[is.finite(prob) & prob > 0]
  n <- length(p)
  if (n == 0L) {
    return(data.frame(
      ess = NA_real_,
      essFraction = NA_real_,
      maxWeight = NA_real_,
      perplexity = NA_real_,
      nNonNegligible = 0L
    ))
  }
  # Renormalise defensively: the caller's vector should already sum to one, but
  # these are diagnostics and must not depend on that.
  p <- p / sum(p)

  # Kish's effective sample size: the number of equally weighted draws that
  # would carry the same information.
  ess <- 1 / sum(p^2)
  # Perplexity, exp of the Shannon entropy: the effective number of candidates
  # actually contributing. It penalises a long tail of tiny weights less
  # harshly than ESS does, so the two together say more than either alone.
  entropy <- -sum(p * log(p))
  denom <- if (is.null(nSuccessful) || !is.finite(nSuccessful) || nSuccessful <= 0) {
    n
  } else {
    nSuccessful
  }

  data.frame(
    ess = ess,
    essFraction = ess / denom,
    maxWeight = max(p),
    perplexity = exp(entropy),
    nNonNegligible = sum(p >= .sirNegligibleWeightFactor / n)
  )
}

# Warn when the retained sample rests on too little of the proposal.
#
# Three separate conditions, each with its own remedy, because conflating them
# produced advice that was wrong for two of the three. Every threshold is
# quoted in the message so the user can judge rather than just be alarmed.
.sirWarnWeightDegeneracy <- function(diag, iterNum, nSamples = NULL) {
  if (!is.finite(diag$ess)) {
    return(invisible(FALSE))
  }
  ess <- diag$ess
  essFrac <- diag$essFraction
  maxW <- diag$maxWeight

  # Bound without the leading dot: cli reads {.name} as a style, not a value.
  essWarnAt <- .sirEssWarn                    # nolint: object_usage_linter.
  fracWarnAt <- .sirEssFractionWarn           # nolint: object_usage_linter.
  maxWarnAt <- .sirMaxWeightWarn              # nolint: object_usage_linter.

  bad_ess <- is.finite(ess) && ess < .sirEssWarn
  bad_frac <- is.finite(essFrac) && essFrac < .sirEssFractionWarn
  bad_max <- is.finite(maxW) && maxW > .sirMaxWeightWarn
  if (!bad_ess && !bad_frac && !bad_max) {
    return(invisible(FALSE))
  }

  n <- if (is.null(nSamples) || !is.finite(nSamples)) NA_integer_ else as.integer(nSamples)
  small_n <- !is.na(n) && n < .sirEssFractionReliableN

  msg <- c(
    "Iteration {iterNum}: the importance weights are concentrated on few samples."
  )

  # 1. Absolute effective sample size: what the result can actually support.
  if (bad_ess) {
    msg <- c(msg, "x" = paste0(
      "Effective sample size is {round(ess, 1)} ",
      "(warns below {essWarnAt}): the retained distribution rests on about ",
      "that many independent points, whatever {.arg nResample} says."
    ))
    msg <- c(msg, "i" = paste0(
      "Percentile intervals from so few effective draws are dominated by ",
      "resampling noise. Effective sample size grows roughly in proportion to ",
      "{.arg nSamples}, so more samples raise it."
    ))
  }

  # 2. Efficiency: a property of the proposal, NOT of the sample size.
  if (bad_frac) {
    msg <- c(msg, "x" = paste0(
      "Proposal efficiency is {round(100 * essFrac, 1)}% ",
      "(warns below {round(100 * fracWarnAt)}%): most draws carry ",
      "negligible weight."
    ))
    msg <- c(msg, "i" = paste0(
      "This ratio is a property of the proposal, not of the sample size: ",
      "raising {.arg nSamples} raises the effective sample size but leaves the ",
      "ratio where it is."
    ))
    # Deliberately does NOT say "widen the proposal". ESS/n measures the
    # MAGNITUDE of the proposal-target mismatch, not its direction: a proposal
    # that is too wide scores just as badly as one that is too narrow. For a
    # standard normal target the asymptotic efficiency of a centred normal
    # proposal is sqrt(2 - 1/s^2)/s, which falls monotonically as s grows --
    # 14.1% at SD 10, 7.1% at SD 20, 3.5% at SD 40 -- so widening a proposal
    # that is already too wide makes this number worse. The convergence plot is
    # what distinguishes the two cases, so the user is sent there.
    msg <- c(msg, "i" = paste0(
      "Low efficiency means the proposal is a poor match for the likelihood, ",
      "but not in which direction -- too wide scores as badly as too narrow. ",
      "Check {.code plot(type = \"convergence\")}: widen with the inflation ",
      "controls only if the first iteration's proposal curve sits below the ",
      "reference."
    ))
    if (small_n) {
      msg <- c(msg, "!" = paste0(
        "Measured on {n} sample{?s}, where this ratio is optimistically biased ",
        "-- the true efficiency is likely lower than the figure above."
      ))
    }
  }

  # 3. A single dominating candidate.
  if (bad_max) {
    msg <- c(msg, "x" = paste0(
      "One sample carries {round(100 * maxW, 1)}% of the weight ",
      "(warns above {round(100 * maxWarnAt)}%)."
    ))
  }

  cli::cli_warn(msg)
  invisible(TRUE)
}
