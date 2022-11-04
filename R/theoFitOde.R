#' Example single dose Theophylline ODE model 
#'
#' This is a nlmixr2 model that is pre-run so that it can be used in
#' package testing and development.  It is regenerated whenever
#' binaries of `nlmixr2extra` are created.  If there is a binary
#' incompatability between the fit objects, a simple rerun of the
#' installation will fix this nlmixr2 fit object.
#'
#' @name theoFitOde
#' 
#' @format A (modified) data frame with  132 rows and 22 columns.
#' 
#' \describe{
#'  \item{ID}{Patient identifier}
#'  \item{TIME}{Time (hr)}
#'  \item{DV}{Dependent variable (concentration)}
#'  \item{PRED}{Predictions without any between subject variability}
#'  \item{RES}{Population Residual}
#'  \item{WRES}{Weighted Residuals under the FO assumption}
#'  \item{IPRED}{Individual Predictions}
#'  \item{IRES}{Individual Residuals}
#'  \item{IWRES}{Individual Weighted Residuals}
#'  \item{CPRED}{Conditional Prediction under the FOCE assumption}
#'  \item{CRES}{Conditional Residuals under the FOCE assumption}
#'  \item{CWRES}{Conditional Weighted Residuals under the FOCE assumption}
#'  \item{eta.ka}{Between subject changes for ka}
#'  \item{eta.cl}{Between subject changes for v}
#'  \item{depot}{amount in the depot compartment}
#'  \item{center}{amount in the central compartment}
#'  \item{ka}{Individual ka values}
#'  \item{cl}{Individual cl values}
#'  \item{v}{Individual volume of distribution}
#'  \item{tad}{Time after dose}
#'  \item{dosenum}{Dose number}
#' }
#' 
NULL

