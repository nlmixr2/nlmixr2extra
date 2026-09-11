# Helpers for writing modeling reports

## Introduction

The `nlmixr2` family of packages makes creating reports easier. When you
fit several candidate models, a few common tasks come up over and over
again:

- deciding whether a fit can be trusted (did any parameter land on a
  boundary?),
- choosing the “best” model by an information criterion, and
- summarizing every model that was tested in a single table.

`nlmixr2extra` provides three helpers for these tasks:

- [`isBoundaryFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/isBoundaryFit.md)
  reports whether a fit has a parameter at its boundary.
- [`getMinAICFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/getMinAICFit.md)
  returns the fit with the lowest AIC, optionally excluding boundary
  fits and silently ignoring fits that errored.
- [`listModelsTested()`](https://nlmixr2.github.io/nlmixr2extra/reference/listModelsTested.md)
  builds a report-ready table of every model tested with its AIC and
  change from the minimum AIC (dAIC).

The functions work for a single model or for many models at once.

### Setup

Start with your data and create the candidate models you want to
compare. Here the data follow a step change (a low value at `conc == 0`
and a higher value for any positive concentration), and we fit three
competing structural models.

``` r

library(nlmixr2est)
library(nlmixr2extra)

# Start with your data
set.seed(42)
d_noec50 <-
  data.frame(
    conc = c(rep(0, 10), rep(1:20, each = 10)),
    DV = c(rnorm(n = 10, mean = 1, sd = 1e-5), rnorm(n = 200, mean = 5, sd = 1e-5)),
    TIME = 0
  )

# An Emax model.  Because the data are a step change, ec50 is pushed to its
# lower boundary (0), which makes this fit unreliable.
modEmax <- function() {
  ini({
    e0 = 1
    emax = 5
    ec50 = c(0, 1.1)
    addSd = 0.5
  })
  model({
    effect <- e0 + emax*conc/(ec50 + conc)
    effect ~ add(addSd)
  })
}

# A step-change model
modStep <- function() {
  ini({
    e0 = 1
    emax = 5
    addSd = 1e-5
  })
  model({
    effect <- e0 + emax*(conc > 0)
    effect ~ add(addSd)
  })
}

# A linear model
modLinear <- function() {
  ini({
    e0 = 1
    slope = 5
    addSd = 1
  })
  model({
    effect <- e0 + slope*conc
    effect ~ add(addSd)
  })
}

# Fit the models
fitEmaxBoundaryIssue := nlmixr2est::nlmixr2(modEmax, data = d_noec50, est = "focei", control = list(print = 0))
#> ℹ loading fit from reporting-helpers-fitEmaxBoundaryIssue.zip
#> ℹ loading fit from fitEmaxBoundaryIssue.R
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> ℹ removing unzipped fit files
fitStep := nlmixr2est::nlmixr2(modStep, data = d_noec50, est = "focei", control = list(print = 0))
#> ℹ loading fit from reporting-helpers-fitStep.zip
#> ℹ loading fit from fitStep.R
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> ℹ removing unzipped fit files
fitLinear := nlmixr2est::nlmixr2(modLinear, data = d_noec50, est = "focei", control = list(print = 0))
#> ℹ loading fit from reporting-helpers-fitLinear.zip
#> ℹ loading fit from fitLinear.R
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> ℹ removing unzipped fit files
```

## Detecting boundary issues

A model whose estimate sits on a parameter boundary usually should not
be trusted or selected.
[`isBoundaryFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/isBoundaryFit.md)
returns `TRUE` when a fit has a parameter at its boundary and `FALSE`
otherwise (including for objects that are not `nlmixr2` fits, such as a
fit that errored).

``` r

# The Emax model pushed ec50 to its lower boundary
isBoundaryFit(fitEmaxBoundaryIssue)
#> [1] FALSE

# The step-change model did not have a boundary issue
isBoundaryFit(fitStep)
#> [1] FALSE
```

## Choosing the best model by AIC

[`getMinAICFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/getMinAICFit.md)
returns the fit with the lowest AIC. By default
(`excludeBoundary = TRUE`) it removes any boundary fits before
comparing, and it silently ignores any argument that cannot produce an
AIC (for example, a model that failed to estimate). You can pass fits as
individual arguments or as a list, and it emits a message when a
boundary fit is removed.

``` r

bestFit <- getMinAICFit(fitEmaxBoundaryIssue, fitStep, fitLinear)
```

If every candidate is excluded or has no AIC,
[`getMinAICFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/getMinAICFit.md)
returns `NULL` with a warning, so it is safe to call inside a larger
report-building pipeline.

## Preparing for the report

### Put the models in a named list

By putting the models in a named list, the functions below can build
more parts of the report. The names become the model descriptions in the
summary table.

``` r

allFits <-
  list(
    "Emax model with additive residual error" = fitEmaxBoundaryIssue,
    "Step-change model with additive residual error" = fitStep,
    "Linear model with additive residual error" = fitLinear
  )
```

### Find your best model by AIC

``` r

bestFit <- getMinAICFit(allFits)
```

### Summarize the best model with its equations and parameters

``` r

knit_print(bestFit, inline = FALSE)
```

``` math
\begin{align*}
{effect} & = {e0}+\frac{{emax} {\times} {conc}}{\left({ec50}+{conc}\right)} \\
{effect} & \sim add({addSd})
\end{align*}
```

``` r

pander::pander(bestFit$parFixed, caption = "Model parameters for the best-fit model")
```

|           |  Est.   |   SE    |  %RSE  |  Back-transformed(95%CI)   |
|:---------:|:-------:|:-------:|:------:|:--------------------------:|
|  **e0**   |  0.948  | 0.00350 | 0.369  |    0.948 (0.942, 0.955)    |
| **emax**  |  4.05   | 0.00375 | 0.0925 |     4.05 (4.05, 4.06)      |
| **ec50**  | 0.00166 | 8.09e-4 |  48.7  | 0.00166 (7.41e-5, 0.00324) |
| **addSd** | 0.00998 | 7.18e-4 |  7.19  | 0.00998 (0.00857, 0.0114)  |

Model parameters for the best-fit model {.table style="width:97%;"}

### Summarize all models tested

[`listModelsTested()`](https://nlmixr2.github.io/nlmixr2extra/reference/listModelsTested.md)
returns a data.frame with the model descriptions, their AIC, and the
change from the minimum AIC (dAIC). Models with a boundary issue are
flagged in an `Exclude` column and are left out of the dAIC calculation,
and the returned data.frame carries a `caption` attribute for pretty
printing with
[`pander::pander()`](https://rdrr.io/pkg/pander/man/pander.html).

``` r

pander::pander(
  listModelsTested(allFits, caption = "Listing of all models tested.")
)
```

|                  Description                   |  AIC  |  dAIC  |
|:----------------------------------------------:|:-----:|:------:|
|    Emax model with additive residual error     | -1262 |   0    |
| Step-change model with additive residual error |  392  | \>1000 |
|   Linear model with additive residual error    | 503.8 | \>1000 |

Listing of all models tested. Abbreviations: AIC = Akaike’s Information
Criterion; dAIC = change from minimum AIC {.table style="width:62%;"}
