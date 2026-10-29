## Emitting events on the rxode2 event bus (see rxode2::rxEventListen()), so
## loggers such as nlmixr2log store one SIR result with the fit instead of a
## run per internal fit.  runSIR() enters the bus scope (its internal
## objective evaluations are silent; the workers of nlmixr2utils::.plap()
## enter it themselves) and on exit emits one `fitResult` with a small
## summary.  Nothing happens when rxode2 has no event bus.

#' @noRd
.sirEventBus <- function() {
  exists("rxEventEmit", envir = asNamespace("rxode2"), inherits = FALSE)
}

#' @noRd
.sirEventEnter <- function() {
  if (.sirEventBus()) getExportedValue("rxode2", ".rxEventEnter")()
  invisible()
}

#' Leave runSIR()'s scope and emit its result
#' @noRd
.sirEventExit <- function(result, fit, call) {
  if (!.sirEventBus()) {
    return(invisible())
  }
  .exit <- getExportedValue("rxode2", ".rxEventExit")
  .summary <- if (!is.null(result) && inherits(fit, "nlmixr2FitCore")) {
    tryCatch(.sirEventSummary(result), error = function(e) NULL)
  }
  if (is.null(.summary)) {
    return(.exit())
  }
  .exit("fitResult", fit = fit, result = .summary, kind = "sir", call = call, fun = "runSIR")
}

#' The SIR summary without the resampled parameters or raw results
#'
#' Those stay in the SIR run directory (`outputDir`), which is kept.
#' @noRd
.sirEventSummary <- function(result) {
  .keep <- c("iterationSummary", "covMatrix", "corMatrix", "sdCorMatrix", "outputDir",
             "fitName", "seed", "fingerprint", "class", "names", "row.names")
  .a <- attributes(result)
  attributes(result) <- .a[intersect(names(.a), .keep)]
  result
}
