#' Make a control step that is quieter and faster
#'
#' @param ctl the control object
#' @return A faster and quieter control object
#' @noRd
setQuietFastControl <- function(ctl) {
  # make estimation steps quieter
  ctl$print <- 0L
  # make estimation steps faster
  ctl$covMethod <- 0L
  ctl$calcTables <- FALSE
  ctl$compress <- FALSE
  ctl
}

#' Attach the data to a figure built without it
#'
#' Figures are built by top-level helpers that call
#' `ggplot2::ggplot(mapping = ...)` without data and receive only column names,
#' titles and flags.  The `aes()` quosures, the layers and the plot keep the
#' frame they were built in, and `serialize()` (so `saveRDS()` and 'targets')
#' writes such a frame out in full with every figure; building in a small frame
#' keeps the fit and the full data out of the figure, which then holds its data
#' only in `$data`.  A builder must evaluate every argument it takes: an
#' argument it never uses stays a promise, which keeps the caller's frame.
#'
#' @param p ggplot built without data
#' @param data data frame for the figure
#' @return `p` with `data` as its data, as `ggplot2::ggplot(data)` stores it
#' @noRd
.plotData <- function(p, data) {
  p$data <- ggplot2::fortify(data)
  p
}
