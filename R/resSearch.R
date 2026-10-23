#' Exhaustively Search for Residual Error Model
#' @param fit Fitted model object
#' @return list of original model and summary for different residual models searched
#' @author Omar I. Elashkar
#' @export
resSearch <- function(fit) {
  UseMethod("resSearch")
}


#' @export
resSearch.nlmixr2Linearize <- function(fit) {
  # assert FOCE linearization

  getobjDf <- function(fit, name) {
    objDf <- fit$objDf[, c("OBJF", "AIC", "BIC")]
    objDf$search <- name
    objDf
  }

  basefit <- getobjDf(fit, "base fit") # extract type and avoid repeat

  OPRED <- NULL
  add.sd <- NULL
  prop.sd <- NULL

  # scale the residual thetas like residuals, as linearize() does (#133)
  refit <- function(mod) {
    nlmixr(
      mod,
      nlme::getData(fit),
      est = "focei",
      control = nlmixr2est::foceiControl(
        scaleC = .linScaleC(
          mod,
          fit$env$originalFit,
          c(prop.sd = "prop", add.sd = "add")
        )
      )
    )
  }

  linFit.prop <- fit |>
    model(rxR2 <- (prop.sd^2 * OPRED^2)) |>
    ini(prop.sd <- 0.1) |>
    refit()
  linFit.prop <- getobjDf(linFit.prop, "prop")

  linFit.combined2 <- fit |>
    model(rxR2 <- (prop.sd^2 * OPRED^2 + add.sd^2)) |>
    ini(prop.sd <- 0.1) |>
    refit()
  linFit.combined2 <- getobjDf(linFit.combined2, "combined2")

  linFit.combined1 <- fit |>
    model(rxR2 <- (prop.sd * OPRED + add.sd)^2) |>
    ini(prop.sd <- 0.1) |>
    refit()
  linFit.combined1 <- getobjDf(linFit.combined1, "combined1")

  # time-varying

  # transformations

  list(summary = rbind(basefit, linFit.prop, linFit.combined2, linFit.combined1), originalFit = fit)
}
