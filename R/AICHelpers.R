#' Get the model with the minimum Akaike Information Criterion (AIC) value
#'
#' @param ... One or more model fits or lists of model fits
#' @param excludeBoundary Exclude a model from selection if it has a parameter
#'   at its boundary
#' @inheritParams stats::AIC
#' @returns The model fit with the minimum AIC; `NULL` if no model fits are
#'   available after exclusions. If more than one model has the same AIC, the
#'   first is returned.
#' @examples
#' \dontrun{
#' d_noec50 <-
#'   data.frame(
#'     conc = c(rep(0, 10), rep(1:20, each = 10)),
#'     DV = c(rnorm(n = 10, mean = 1, sd = 1e-5), rnorm(n = 200, mean = 5, sd = 1e-5)),
#'     TIME = 0
#'   )
#'
#' modEmax <- function() {
#'   ini({
#'     e0 = 1
#'     emax = 5
#'     ec50 = c(0, 1.1)
#'     addSd = 0.5
#'   })
#'   model({
#'     effect <- e0 + emax*conc/(ec50 + conc)
#'     effect ~ add(addSd)
#'   })
#' }
#'
#' modStep <- function() {
#'   ini({
#'     e0 = 1
#'     emax = 5
#'     addSd = 1e-5
#'   })
#'   model({
#'     effect <- e0 + emax*(conc > 0)
#'     effect ~ add(addSd)
#'   })
#' }
#'
#' modLinear <- function() {
#'   ini({
#'     e0 = 1
#'     slope = 5
#'     addSd = 1
#'   })
#'   model({
#'     effect <- e0 + slope*conc
#'     effect ~ add(addSd)
#'   })
#' }
#'
#' fitEmaxBoundaryIssue <- nlmixr2(modEmax, data = d_noec50, est = "focei", control = list(print = 0))
#' fitStep <- nlmixr2(modStep, data = d_noec50, est = "focei", control = list(print = 0))
#' fitLinear <- nlmixr2(modLinear, data = d_noec50, est = "focei", control = list(print = 0))
#' getMinAICFit(list(fitEmaxBoundaryIssue, fitStep, fitLinear))
#' }
#' @export
getMinAICFit <- function(..., excludeBoundary = TRUE, k = 2) {
  args <- list(...)
  fitList <- list()
  for (currentArg in args) {
    if (is.list(currentArg) && inherits(try(AIC(currentArg), silent = TRUE), "try-error")) {
      fitList <- append(fitList, currentArg)
    } else {
      fitList <- append(fitList, list(currentArg))
    }
  }
  hadFits <- length(fitList) > 0
  if (excludeBoundary) {
    # Remove models at the boundary
    for (idx in rev(seq_along(fitList))) {
      if (isBoundaryFit(fitList[[idx]])) {
        message("Removing model with a parameter at the boundary")
        fitList[[idx]] <- NULL
      }
    }
  }
  if (length(fitList) > 0) {
    fitAIC <- vapply(X = fitList, FUN = AIC, FUN.VALUE = 1)
    ret <- fitList[[which(fitAIC %in% min(fitAIC, na.rm = TRUE))[1]]]
  } else if (!hadFits) {
    warning("No model fits in input list")
    ret <- NULL
  } else {
    warning("No model fits in input after excluding models")
    ret <- NULL
  }
  ret
}

#' Detect if a fit has any parameter at the boundary
#'
#' @param object A model fit to check
#' @returns `TRUE` if the fit is an nlmxir2 fit and it has a parameter at the
#'   boundary; `FALSE` otherwise
#' @export
isBoundaryFit <- function(object) {
  if (inherits(object, "nlmixr2FitCore")) {
    any(grepl(x = object$covMethod, pattern = "Boundary issue"))
  } else {
    FALSE
  }
}

#' Make a list of models tested ready for reporting
#'
#' @param fitList A named list of models that were tested
#' @inheritParams getMinAICFit
#' @param caption The caption for the resulting table
#' @returns A data.frame with a "caption" attribute ready for printing in a
#'   report with `pander::pander()`. The data.frame will have names of
#'   "Description", "AIC", and "dAIC"; if exclusions are applied, the data.frame
#'   will also have an "Exclude" column.
#' @export
listModelsTested <- function(fitList, caption, excludeBoundary = TRUE, k = 2) {
  checkmate::expect_names(names(fitList))
  ret <-
    data.frame(
      Description = names(fitList),
      AIC = vapply(X = fitList, FUN = AIC, FUN.VALUE = 1)
    )
  row.names(ret) <- NULL
  # Default when not calculated
  ret$dAIC <- "-"
  calcdAIC <- rep(!is.na(ret$AIC), nrow(ret))
  if (excludeBoundary) {
    excludes <- vapply(X = fitList, FUN = isBoundaryFit, FUN.VALUE = TRUE)
    if (any(excludes)) {
      ret$Exclude <- ""
      ret$Exclude[excludes] <- "parameter at boundary"
      calcdAIC <- !excludes & !is.na(ret$AIC)
    }
  }
  dAICValues <- ret$AIC[calcdAIC] - min(ret$AIC[calcdAIC], na.rm = TRUE)
  # pretty printing for dAIC
  dAICPretty <- rep(NA_character_, length(dAICValues))
  dAICPretty[dAICValues > 0] <- sprintf("%.3f", dAICValues[dAICValues > 0])
  dAICPretty[dAICValues > 10] <- sprintf("%.2f", dAICValues[dAICValues > 10])
  dAICPretty[dAICValues > 100] <- sprintf("%.1f", dAICValues[dAICValues > 100])
  if (any(dAICValues > 1000)) {
    dAICPretty[dAICValues > 1000] <- ">1000"
  }
  if (any(dAICValues < 0.001)) {
    dAICPretty[dAICValues < 0.001] <- "<0.001"
  }
  if (any(dAICValues %in% 0)) {
    dAICPretty[dAICValues %in% 0] <- "0"
  }

  ret$dAIC[calcdAIC] <- dAICPretty

  attr(ret, "caption") <- paste(caption, "Abbreviations: AIC = Akaike's Information Criterion; dAIC = change from minimum AIC")
  ret
}
