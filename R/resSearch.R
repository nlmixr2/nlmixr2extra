
#' Exhaustively Search for Residual Error Model
#' @param fit Fitted model object
#' @return list of original model and summary for different residual models searched
#' @author Omar I. Elashkar
#' @export
resSearch <- function(fit){
    UseMethod("resSearch")
}


#' @export
resSearch.nlmixr2Linearize <- function(fit){
  # TODO assert FOCE linearization 
  
  getobjDf <- function(fit, name){
    objDf <- fit$objDf[,c("OBJF", "AIC", "BIC")]
    objDf$search <- name
    objDf
  }
  
  basefit <- getobjDf(fit, "base fit") # TODO extract type and avoid repeat
  
  linFit.prop <- fit |>  
                model(rxR2 <- (prop.sd^2*OPRED^2)) |> ini(prop.sd <- 0.1) |> 
                nlmixr(nlme::getData(fit), est = "focei")
  linFit.prop <- getobjDf(linFit.prop, "prop")
    
  linFit.combined2 <- fit |> 
                    model(rxR2 <- (prop.sd^2*OPRED^2 + add.sd^2)) |> ini(prop.sd <- 0.1) |>
                    nlmixr(nlme::getData(fit), est = "focei")
  linFit.combined2 <- getobjDf(linFit.combined2, "combined2")
  
  linFit.combined1 <- fit |> 
                      model(rxR2 <- (prop.sd*OPRED + add.sd)^2) |> ini(prop.sd <- 0.1) |>
                      nlmixr(nlme::getData(fit), est = "focei")
  linFit.combined1 <- getobjDf(linFit.combined1, "combined1")
  
  
  # TODO time-varying
  
  # TODO transformations
  
  
  list(summary = rbind(basefit, linFit.prop, linFit.combined2, linFit.combined1) ,
       originalFit = fit)
  
}
