## Emitting events on the rxode2 event bus (see rxode2::rxEventListen()), so
## loggers such as nlmixr2log store one entry per driver call instead of one
## run per internal fit.  Each driver enters the bus scope (its internal
## nlmixr2() and rxSolve() calls are then silent) and, on exit, emits one
## event chosen from what it returned:
##
## * the input fit itself, changed in place (bootstrapFit), or `update = TRUE`
##   (preconditionFit): `fitUpdate`, a new version of the fit's run;
## * a different fit (linearize, regularmodel, a selected model):
##   `fitComplete` with `source = <kind>`, a run linked to the input fit;
## * anything else: one `fitResult` with a summary in which embedded fits are
##   reduced to their objective function and parameter table.
##
## All of this does nothing when rxode2 has no event bus.

#' @noRd
.extraEventBus <- function() {
  exists("rxEventEmit", envir = asNamespace("rxode2"), inherits = FALSE)
}

#' @noRd
.extraEventEnter <- function() {
  if (.extraEventBus()) getExportedValue("rxode2", ".rxEventEnter")()
  invisible()
}

#' Leave the driver's scope and emit its one event
#'
#' @param result the driver's return value (`returnValue()`; NULL on error)
#' @param fit the input fit
#' @param call the driver's call
#' @param kind the label stored with the result (e.g. "bootstrap")
#' @param fun the driver's name, used as the head of the recorded call
#' @param update always emit `fitUpdate` for `fit` (it was changed in place)
#' @noRd
.extraEventExit <- function(result, fit, call, kind, fun, update = FALSE) {
  if (!.extraEventBus()) {
    return(invisible())
  }
  .exit <- getExportedValue("rxode2", ".rxEventExit")
  if (is.null(result) || !inherits(fit, "nlmixr2FitCore")) {
    return(.exit())
  }
  if (update || (inherits(result, "nlmixr2FitCore") && identical(result$env, fit$env))) {
    return(.exit("fitUpdate", fit = fit, original = fit, name = NULL, what = kind,
                 inPlace = TRUE, call = call, fun = fun))
  }
  if (inherits(result, "nlmixr2FitCore")) {
    return(.exit("fitComplete", fit = result, object = fit, call = call, objName = NULL,
                 source = kind, fun = fun))
  }
  ## never let summarizing break the driver: on failure just leave the scope
  .summary <- tryCatch(.extraEventSummary(result), error = function(e) NULL)
  if (is.null(.summary)) {
    return(.exit())
  }
  .exit("fitResult", fit = fit, result = .summary, kind = kind, call = call, fun = fun)
}

#' A small, fit-free copy of a driver result
#'
#' Fits are replaced by their objective function and parameter table, read
#' from the fit environment without computing anything.
#' @noRd
.extraEventSummary <- function(x, depth = 0L) {
  if (inherits(x, "nlmixr2FitCore")) {
    .env <- if (is.environment(x)) x else tryCatch(x$env, error = function(e) NULL)
    if (!is.environment(.env) && is.list(x)) .env <- .subset2(x, "env")
    if (!is.environment(.env)) {
      return(list(OBJF = NA_real_, parFixedDf = NULL))
    }
    .objDf <- get0("objDf", envir = .env, inherits = FALSE)
    return(list(
      OBJF = if (is.data.frame(.objDf) && nrow(.objDf)) .objDf$OBJF[1] else NA_real_,
      parFixedDf = get0("parFixedDf", envir = .env, inherits = FALSE)
    ))
  }
  if (inherits(x, c("gg", "ggplot", "data.frame"))) {
    return(x)
  }
  if (is.environment(x) || is.function(x)) {
    return(NULL)
  }
  if (is.list(x) && !is.data.frame(x) && depth < 4L) {
    .cls <- class(x)
    x <- lapply(unclass(x), .extraEventSummary, depth = depth + 1L)
    if (!identical(.cls, "list")) attr(x, "nlmixr2extraClass") <- .cls
  }
  x
}
