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
#> Loading required package: nlmixr2data
library(nlmixr2extra)

# Start with your data
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
fitEmaxBoundaryIssue <- nlmixr2est::nlmixr2(modEmax, data = d_noec50, est = "focei", control = list(print = 0))
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling EBE model...
#> ✔ done
#> rxode2 5.1.3 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
#> 
#> Attaching package: 'rxode2'
#> The following objects are masked from 'package:nlmixr2est':
#> 
#>     boxCox, yeoJohnson
#> done
#> → Calculating residuals/tables
#> ✔ done
fitStep <- nlmixr2est::nlmixr2(modStep, data = d_noec50, est = "focei", control = list(print = 0))
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling EBE model...
#> ✔ done
#> calculating covariance matrix
#> done
#> → Calculating residuals/tables
#> ✔ done
fitLinear <- nlmixr2est::nlmixr2(modLinear, data = d_noec50, est = "focei", control = list(print = 0))
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling EBE model...
#> ✔ done
#> calculating covariance matrix
#> done
#> → Calculating residuals/tables
#> ✔ done
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
#> [1] TRUE

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
#> Removing model with a parameter at the boundary
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
#> Removing model with a parameter at the boundary
```

### Summarize the best model with its equations and parameters

``` r

knit_print(bestFit, inline = FALSE)
```

``` math
\begin{align*}
{effect} & = {e0}+{emax} {\times} \left({conc}>{0}\right) \\
{effect} & \sim add({addSd})
\end{align*}
```

``` r

pander::pander(bestFit$parFixed, caption = "Model parameters for the best-fit model")
```

|           |   Est.   |   SE   |  %RSE   | Back-transformed(95%CI)  |
|:---------:|:--------:|:------:|:-------:|:------------------------:|
|  **e0**   |  1.000   | 0.3162 |  31.62  |  1.000 (0.3802, 1.620)   |
| **emax**  |  4.000   | 0.3240 |  8.101  |   4.000 (3.365, 4.635)   |
| **addSd** | 9.866e-6 | 3.203  | 3.247e7 | 9.866e-6 (-6.278, 6.278) |

Model parameters for the best-fit model {.table style="width:96%;"}

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

| Description | AIC | dAIC | Exclude |
|:--:|:--:|:--:|:--:|
| Emax model with additive residual error | -2798 | \- | parameter at boundary |
| Step-change model with additive residual error | 392 | 0 |  |
| Linear model with additive residual error | 503.8 | 111.9 |  |

Listing of all models tested. Abbreviations: AIC = Akaike’s Information
Criterion; dAIC = change from minimum AIC {.table style="width:96%;"}
