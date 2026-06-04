# Stepwise Covariate Model-selection (SCM) method

Stepwise Covariate Model-selection (SCM) method

## Usage

``` r
covarSearchAuto(
  fit,
  varsVec,
  covarsVec,
  pVal = list(fwd = 0.05, bck = 0.01),
  catvarsVec = NULL,
  searchType = c("scm", "forward", "backward"),
  restart = FALSE
)
```

## Arguments

- fit:

  an nlmixr2 'fit' object

- varsVec:

  a list of candidate variables to which the covariates could be added

- covarsVec:

  a list of candidate covariates that need to be tested

- pVal:

  a named list with names 'fwd' and 'bck' for specifying the p-values
  for the forward and backward searches, respectively

- catvarsVec:

  character vector of categorical covariates that need to be added

- searchType:

  one of 'scm', 'forward' and 'backward' to specify the covariate search
  method; default is 'scm'

- restart:

  a boolean that controls if the search should be restarted; default is
  FALSE

## Value

A list summarizing the covariate selection steps and output; This list
has the "summaryTable" for the overall summary of the covariate
selection as well as "resFwd" for the forward selection method and
"resBck" for the backward selection method.

## Author

Vipul Mann, Matthew Fidler, Vishal Sarsani

## Examples

``` r
if (FALSE) { # \dontrun{
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

fit <- nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "saem", control = list(print = 0))
rxode2::.rxWithWd(tempdir(), {# with temporary directory

auto1 <- covarSearchAuto(fit, varsVec = c("ka", "cl"),
    covarsVec = c("WT"))

})

## Note that this didn't include sex, add it to dataset and restart model


d <- nlmixr2data::theo_sd
d$SEX <-0
d$SEX[d$ID<=6] <-1

fit <- nlmixr2(one.cmt, d, est = "saem", control = list(print = 0))

# This would restart if for some reason the search crashed:

rxode2::.rxWithWd(tempdir(), {# with temporary directory

auto2 <- covarSearchAuto(fit, varsVec = c("ka", "cl"), covarsVec = c("WT"),
                catvarsVec= c("SEX"), restart = TRUE)

auto3 <- covarSearchAuto(fit, varsVec = c("ka", "cl"), covarsVec = c("WT"),
                catvarsVec=  c("SEX"), restart = TRUE,
                searchType = "forward")
})
} # }
```
