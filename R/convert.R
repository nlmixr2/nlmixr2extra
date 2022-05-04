#' Convert nlmixr compatible data to monolix compatible data (if possible)
#'
#' @param model rxode2 model for conversion
#' @param data Input dataset
#' @param table is the table control; this is mostly to figure out if there are additional columns to keep.
#' @return Monolix compatible dataset
#' @author Matthew L. Fidler
#' @export
#' @examples
#'
nlmixrDataToMonolix <- function(model, data, table=nlmixr2est::tableControl()) {
  model <- rxode2::assertRxUi(model, extra=" to convert the data with 'nlmixrDataToMonolix'")
  rxode2::assertRxUiPrediction(model, extra=" to convert the data with 'nlmixrDataToMonolix'")
  .env <- new.env(parent=emptyenv())
  .env$table <- table
  nlmixr2est::.foceiPreProcessData(data, .env, model)
  .mv <- rxode2::rxModelVars(model)
  .flag <- .mv$flags
  .conv0 <- .Call(`_nlmixr2extra_convertDataBack`, .env$dataSav$ID, .env$dataSav$TIME, .env$dataSav$AMT,
                  .env$dataSav$II, .env$dataSav$EVID, .env$dataSav$CMT,
                  model$predDf$cmt, model$predDf$dvid, .flag["ncmt"], .flag["ka"], length(.mv$state),
                  replaceEvid=5L)
  .conv0
}
