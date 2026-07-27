#' Rebuild the stored `theoFitOde` fit
#'
#' Sources `build/build.R`, which re-runs the focei fit and rewrites
#' `data/theoFitOde.rda`.  It is wired into the `theoFitOde` docs with
#' `@eval`, so `devtools::document()` regenerates the fit; keep it next
#' to that documentation so the two are not separated again.  The fit is
#' only rebuilt in a source checkout -- `build/` is `.Rbuildignore`d, so
#' this is a no-op for an installed package.
#'
#' @return "", called for the side effect of rewriting the stored fit
#'
#' @noRd
.buildModel <- function() { # nocov start
  .owd <- getwd()
  on.exit(setwd(.owd))
  try(source(file.path(devtools::package_file(), "build", "build.R")))
  ""
} # nocov end

#' Example single dose Theophylline ODE model
#'
#' This is a nlmixr2 model that is pre-run so that it can be used in
#' package testing and development.  It is regenerated when
#' `devtools::document()` is run on a source checkout.  If there is a
#' binary incompatability between the fit objects, re-documenting (or
#' running `build/build.R`) will fix this nlmixr2 fit object.
#'
#' @name theoFitOde
#'
#' @eval .buildModel()
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

