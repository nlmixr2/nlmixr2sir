# The dOFV-vs-chi-square convergence diagnostic.
#
# This is the primary SIR diagnostic and the reason the tool is trusted. Port
# of the dOFV section of PsN R-scripts/sir_default.R.
#
# Three curves per iteration, all as empirical quantiles of dOFV:
#
#   reference : qchisq(q, df = number of estimated parameters)
#   proposal  : over every evaluated sample in the iteration
#   SIR       : over the resampled subset only
#
# Convergence reads as the SIR curve settling onto the reference. A proposal
# curve that sits *below* the reference is too narrow, which SIR cannot fix by
# resampling -- hence the automatic warning.

# PsN's quantile grid: one point per resample, stopping short of 1 so the
# reference chi-square does not run to Inf.
.sirDofvQuantiles <- function(nResample) {
  checkmate::assertCount(nResample, positive = TRUE)
  if (nResample < 3L) {
    cli::cli_abort(c(
      "At least 3 resamples are needed for a dOFV quantile curve.",
      "i" = "The last iteration resampled {nResample}."
    ))
  }
  seq(0, (nResample - 1L) / nResample, length.out = nResample - 1L)
}

# Empirical dOFV quantile curves for every iteration, plus the reference.
.sirDofvCurves <- function(x, nParams = NULL, quant = NULL) {
  iterations <- attr(x, "iterations", exact = TRUE)
  if (is.null(iterations) || length(iterations) == 0L) {
    cli::cli_abort("No stored SIR iterations to build a dOFV curve from.")
  }
  if (is.null(nParams)) {
    nParams <- length(unique(x$param))
  }
  if (is.null(quant)) {
    lastRaw <- iterations[[length(iterations)]]$rawResults
    quant <- .sirDofvQuantiles(sum(lastRaw$resamples > 0L))
  }

  ref <- data.frame(
    iteration = NA_integer_,
    label = "reference",
    type = "reference",
    quantile = quant,
    dOFV = stats::qchisq(quant, df = nParams),
    stringsAsFactors = FALSE
  )

  perIter <- lapply(seq_along(iterations), function(i) {
    raw <- iterations[[i]]$rawResults
    # One row per resample slot, so a vector selected twice counts twice.
    proposalDofv <- .sirProposalRows(raw)$deltaofv
    sirDofv <- raw$deltaofv[raw$resamples > 0L]
    rbind(
      .sirQuantileFrame(proposalDofv, quant, i, "proposal"),
      .sirQuantileFrame(sirDofv, quant, i, "SIR")
    )
  })

  out <- do.call(rbind, c(list(ref), perIter))
  out$label <- factor(
    out$label,
    levels = c("reference", unique(out$label[out$label != "reference"]))
  )
  rownames(out) <- NULL
  out
}

.sirQuantileFrame <- function(dofv, quant, iteration, type) {
  dofv <- dofv[!is.na(dofv)]
  if (length(dofv) == 0L) {
    return(NULL)
  }
  data.frame(
    iteration = as.integer(iteration),
    label = paste0(type, " ", iteration),
    type = type,
    quantile = quant,
    dOFV = unname(stats::quantile(dofv, probs = quant, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

# Resampling noise around a SIR dOFV curve: repeat the weighted resampling many
# times and take the 2.5/97.5 percentiles of the resulting quantile curves.
# PsN uses 2000 replicates; the default here is lower because the cost is borne
# interactively, and it is exposed so it can be raised.
.sirDofvNoise <- function(
  x,
  iteration,
  quant,
  nReplicate = 500L,
  capResampling = 1
) {
  iterations <- attr(x, "iterations", exact = TRUE)
  raw <- iterations[[iteration]]$rawResults
  base <- .sirProposalRows(raw)
  nResample <- sum(raw$resamples > 0L)
  # The normalized resampling probability, not the raw importance ratio.
  # sirCalcWeights() normalizes in log space, so prob_resample is finite by
  # construction, while importance_ratio = exp(log_ir) can overflow to Inf on a
  # strongly favoured candidate. Filtering on the raw ratio therefore discarded
  # exactly the candidate carrying all the weight, and the band then described
  # a different weighted population from the one the resampler drew from.
  prob <- base$prob_resample %||% base$probability_resample
  ok <- !is.na(base$deltaofv) & is.finite(prob) & prob > 0
  base <- base[ok, , drop = FALSE]
  prob <- prob[ok]
  if (nrow(base) < 2L || nResample < 2L) {
    return(NULL)
  }
  # Replicates go through sirResample() itself rather than a local
  # sample.int() call. The local version treated any cap above one as unlimited
  # replacement, while sirResample() expands each candidate into a finite number
  # of slots -- so the band described a different resampler from the one that
  # produced the retained sample, which is the one thing a noise band must not
  # do.
  cap <- max(1L, as.integer(floor(capResampling)))
  usable <- length(prob)
  draws <- min(nResample, usable * cap)
  if (draws < 2L) {
    return(NULL)
  }
  weights <- data.frame(prob_resample = prob / sum(prob))
  dofvMat <- matrix(
    base$deltaofv,
    ncol = 1L,
    dimnames = list(NULL, "deltaofv")
  )

  curves <- vapply(
    seq_len(nReplicate),
    function(i) {
      drawn <- sirResample(
        dofvMat,
        weights,
        m = draws,
        capResampling = cap
      )$samples[, 1L]
      unname(stats::quantile(drawn, probs = quant, na.rm = TRUE))
    },
    numeric(length(quant))
  )
  # vapply() drops to a bare vector when length(quant) == 1, and apply() below
  # needs a dim. Restore it rather than relying on every caller passing several
  # quantiles: the values already come back in the right (quantile, replicate)
  # order, so this only reinstates the shape.
  dim(curves) <- c(length(quant), nReplicate)

  data.frame(
    iteration = as.integer(iteration),
    quantile = quant,
    low = apply(curves, 1L, stats::quantile, probs = 0.025, na.rm = TRUE),
    high = apply(curves, 1L, stats::quantile, probs = 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# The resampling cap the run actually used. Stored on the result by runSIR();
# falls back to PsN's default of one for an object built without it.
.sirEffectiveCap <- function(x) {
  ctl <- attr(x, "control", exact = TRUE)
  cap <- if (is.null(ctl)) NULL else ctl$capResampling
  if (is.null(cap) || !is.finite(cap)) 1 else cap
}

# PsN's automatic warning: if the first iteration's proposal curve falls below
# the reference chi-square for more than a quarter of the quantiles, the
# proposal is too narrow and SIR cannot recover from it by resampling.
.sirProposalTooNarrow <- function(curves, threshold = 0.25) {
  ref <- curves[curves$type == "reference", , drop = FALSE]
  prop <- curves[
    curves$type == "proposal" & curves$iteration == 1L,
    ,
    drop = FALSE
  ]
  if (nrow(ref) == 0L || nrow(prop) == 0L) {
    return(list(warn = FALSE, fraction = NA_real_))
  }
  prop <- prop[match(ref$quantile, prop$quantile), , drop = FALSE]
  frac <- mean(prop$dOFV < ref$dOFV, na.rm = TRUE)
  list(warn = isTRUE(frac > threshold), fraction = frac)
}

.sirWarnProposalTooNarrow <- function(check) {
  if (!isTRUE(check$warn)) {
    return(invisible(FALSE))
  }
  cli::cli_warn(c(
    "The iteration-1 proposal is not entirely above the reference chi-square.",
    "i" = "It falls below for {round(100 * check$fraction)}% of quantiles.",
    "i" = "Consider restarting with an inflated proposal, e.g. {.code runSIRControl(thetaInflation = 1.5, omegaInflation = 1.5, sigmaInflation = 1.5)}."
  ))
  invisible(TRUE)
}

# The convergence plot itself.
.sirConvergencePlot <- function(
  x,
  noise = TRUE,
  nReplicate = 500L,
  warn = TRUE
) {
  curves <- .sirDofvCurves(x)
  check <- .sirProposalTooNarrow(curves)
  if (isTRUE(warn)) {
    .sirWarnProposalTooNarrow(check)
  }

  quant <- curves$quantile[curves$type == "reference"]
  iterations <- attr(x, "iterations", exact = TRUE)

  ribbon <- NULL
  if (isTRUE(noise) && length(iterations) > 0L) {
    # PsN puts a resampling-noise band on the last two iterations, which is
    # where settling should be judged.
    want <- utils::tail(seq_along(iterations), 2L)
    ribbon <- do.call(
      rbind,
      lapply(want, function(i) {
        .sirDofvNoise(
          x,
          i,
          quant,
          nReplicate,
          capResampling = .sirEffectiveCap(x)
        )
      })
    )
  }

  emp <- curves[curves$type != "reference", , drop = FALSE]
  levs <- paste("Iteration", sort(unique(emp$iteration)))
  emp$iterationLabel <- factor(paste("Iteration", emp$iteration), levels = levs)

  # The reference is the same curve in every facet, so it is repeated across
  # them rather than drawn once and lost outside the facetting.
  ref <- curves[curves$type == "reference", , drop = FALSE]
  ref <- do.call(
    rbind,
    lapply(levs, function(l) {
      r <- ref
      r$iterationLabel <- factor(l, levels = levs)
      r
    })
  )

  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$quantile))
  if (!is.null(ribbon) && nrow(ribbon) > 0L) {
    ribbon$iterationLabel <- factor(
      paste("Iteration", ribbon$iteration),
      levels = levs
    )
    p <- p +
      ggplot2::geom_ribbon(
        data = ribbon,
        ggplot2::aes(ymin = .data$low, ymax = .data$high),
        fill = "#6BAED6",
        alpha = 0.25
      )
  }
  p +
    ggplot2::geom_line(
      data = ref,
      ggplot2::aes(y = .data$dOFV),
      color = "grey30",
      linewidth = 0.9,
      linetype = "dashed"
    ) +
    ggplot2::geom_line(
      data = emp,
      ggplot2::aes(y = .data$dOFV, color = .data$type),
      linewidth = 0.7
    ) +
    ggplot2::facet_wrap(stats::as.formula("~ iterationLabel")) +
    ggplot2::scale_color_manual(
      values = c(proposal = "#D94801", SIR = "#2171B5"),
      name = NULL
    ) +
    ggplot2::labs(
      x = "Quantile",
      y = "dOFV",
      caption = "Dashed grey: reference chi-square"
    ) +
    ggplot2::theme_bw()
}
