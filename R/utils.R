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
