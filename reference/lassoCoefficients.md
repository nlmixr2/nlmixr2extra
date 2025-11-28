# Return Final lasso coefficients after finding optimal t

Return Final lasso coefficients after finding optimal t

## Usage

``` r
lassoCoefficients(
  fit,
  varsVec,
  covarsVec,
  catvarsVec,
  constraint = 1e-08,
  stratVar = NULL,
  ...
)
```

## Arguments

- fit:

  nlmixr2 fit.

- varsVec:

  character vector of variables that need to be added

- covarsVec:

  character vector of covariates that need to be added

- catvarsVec:

  character vector of categorical covariates that need to be added

- constraint:

  theta cutoff. below cutoff then the theta will be fixed to zero

- stratVar:

  A variable to stratify on for cross-validation

- ...:

  Other parameters to be passed to optimalTvaluelasso

## Value

return data frame of final lasso coefficients

## Author

Vishal Sarsani

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

d <- nlmixr2data::theo_sd
d$SEX <-0
d$SEX[d$ID<=6] <-1

fit <- nlmixr2(one.cmt, d, est = "saem", control = list(print = 0))
varsVec <- c("ka","cl","v")
covarsVec <- c("WT")
catvarsVec <- c("SEX")

# Lasso coefficients:

lassoDf <- lassoCoefficients(fit, varsVec, covarsVec, catvarsVec, constraint=1e-08, stratVar = NULL)
} # }
```
