# Create Horseshoe summary posterior estimates

Create Horseshoe summary posterior estimates

## Usage

``` r
horseshoeSummardf(fit, covarsVec, ...)
```

## Arguments

- fit:

  compiled rxode2 nlmir2 model fit

- covarsVec:

  character vector of covariates that need to be added

- ...:

  other parameters passed to brm(): warmup = 1000, iter = 2000, chains =
  4, cores = 4, control = list(adapt_delta = 0.99, max_treedepth = 15)

## Value

Horse shoe Summary data frame of all covariates

## Author

Vishal Sarsani, Christian Bartels

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
fit <- nlmixr2(one.cmt, d, est = "saem", control = list(print = 0))
covarsVec <- c("WT")

# Horseshoe summary posterior estimates:

#hsDf <- horseshoeSummardf(fit,covarsVec,cores=2)
#brms sometimes may throw a Error in sink(type = “output”)
#Issue Should be fixed by uninstalling and re-installing rstan
} # }
```
