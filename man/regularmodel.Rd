% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lassocov.R
\name{regularmodel}
\alias{regularmodel}
\title{Regular lasso model}
\usage{
regularmodel(
  fit,
  varsVec,
  covarsVec,
  catvarsVec,
  constraint = 1e-08,
  lassotype = c("regular", "adaptive", "adjusted"),
  stratVar = NULL,
  ...
)
}
\arguments{
\item{fit}{nlmixr2 fit.}

\item{varsVec}{character vector of variables that need to be added}

\item{covarsVec}{character vector of covariates that need to be added}

\item{catvarsVec}{character vector of categorical covariates that need to be added}

\item{constraint}{theta cutoff. below cutoff then the theta will be fixed to zero.}

\item{lassotype}{must be  'regular' , 'adaptive', 'adjusted'}

\item{stratVar}{A variable to stratify on for cross-validation.}

\item{...}{Other parameters to be passed to optimalTvaluelasso}
}
\value{
return fit of the selected lasso coefficients
}
\description{
Regular lasso model
}
\examples{
\dontrun{
one.cmt <- function() {
  ini({
    tka <- 0.45; label("Ka")
    tcl <- log(c(0, 2.7, 100)); label("Cl")
    tv <- 3.45; label("V")
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

d <- nlmixr2data::theo_sd
d$SEX <-0
d$SEX[d$ID<=6] <-1

fit <- nlmixr2(one.cmt, d, est = "saem", control = list(print = 0))
varsVec <- c("ka","cl","v")
covarsVec <- c("WT")
catvarsVec <- c("SEX")

# Model fit with regular lasso coefficients:

lassoDf <- regularmodel(fit,varsVec,covarsVec,catvarsVec)
# Model fit with adaptive lasso coefficients:

lassoDf <- regularmodel(fit,varsVec,covarsVec,catvarsVec,lassotype='adaptive')
# Model fit with adaptive-adjusted lasso coefficients:

lassoDf <- regularmodel(fit,varsVec,covarsVec,catvarsVec, lassotype='adjusted')
}
}
\author{
Vishal Sarsani
}
