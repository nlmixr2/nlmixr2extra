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
  if (.conv0$hasTinf && .conv0$hasRate) {
    stop("monolix does not support a fixed duration (`tinf`) and rate (`rate`) at the same time",
         call.=FALSE)
  }
  if (.conv0$hasTinf) {
    warning("monolix changes infusion times for `tinf` with bioavailability differently than `nlmixr2`, make sure there is no bioavailibilty changes for this infusion in the model",
            call.=FALSE)
  }
  if (.conv0$turnOffCmt) {
    stop("monolix cannot turn off compartments like `nlmixr2` can, this dataset will not work with monolix",
         call.=FALSE)
  }
  if (.conv0$hasPhantom) {
    stop("transit compartment phantom events are not supported in monolix",
         call.=FALSE)
  }
  if (.conv0$hasReplace) {
    stop("replacement events are not supported in monolix",
         call.=FALSE)
  }
  if (.conv0$hasMult) {
    stop("multiply events are not supported in monolix",
         call.=FALSE)
  }
  if (.conv0$hasSsRate) {
    stop("steady state infusions are not supported in monolix",
         call.=FALSE)
  }
  .df <- .conv0$df
  .new <- .env$dataSav
  .new$EVID <-.df$EVID
  .new$AMT <- .df$AMT
  .new$YTYPE <- .df$DVID
  .new$CMT <- .df$CMT
  if (.conv0$hasSs) {
    .new$SS <- .df$SS
  }
  if (.conv0$hasRate) {
    .new$RATE <- .df$RATE
  }
  if (.conv0$hasTinf) {
    .new$TINF <- .df$TINF
  }
  .new[.df$.nlmixrKeep, ]
}
