# nlmixr2extra

The goal of nlmixr2extra is to provide the tools to help with common
pharmacometric tasks with nlmixr2 models like bootstrapping, covariate
selection etc.

## Installation

You can install the development version of nlmixr2extra from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("nlmixr2/nlmixr2data")
remotes::install_github("nlmixr2/lotri")
remotes::install_github("nlmixr2/rxode2")
remotes::install_github("nlmixr2/nlmixr2est")
remotes::install_github("nlmixr2/nlmixr2extra")
```

## Example of a `bootstrapFit()`

This is a basic example of bootstrapping provided by this package

``` r
library(nlmixr2est)
#> Loading required package: nlmixr2data
library(nlmixr2extra)
# basic example code
# The basic model consists of an ini block that has initial estimates
one.compartment <- function() {
  ini({
    tka <- 0.45; label("Absorption rate, Ka")
    tcl <- 1; label("Clearance, Cl")
    tv <- 3.45; label("Central volumne, V")
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7; label("Additive residual error")
  })
  # and a model block with the error specification and model specification
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(depot) = -ka * depot
    d/dt(center) = ka * depot - cl / v * center
    cp = center / v
    cp ~ add(add.sd)
  })
}

# The fit is performed by the function nlmixr/nlmixr2 specifying the model, data
# and estimate (in a real estimate, nBurn and nEm would be much higher.)
fit <- nlmixr2(one.compartment, theo_sd,  est="saem", saemControl(print=0, nBurn = 10, nEm = 20))
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of saem model...
#> ✔ done
#> → finding duplicate expressions in saem model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in saem model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> ✔ done
#> using C compiler: 'gcc.exe (GCC) 13.2.0'
#> ℹ calculate uninformed etas
#> ℹ done
#> rxode2 3.0.2 using 8 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
#> 
#> Attaching package: 'rxode2'
#> The following objects are masked from 'package:nlmixr2est':
#> 
#>     boxCox, yeoJohnson
#> Calculating covariance matrix
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of saem model...
#> ✔ done
#> → finding duplicate expressions in saem predOnly model 0...
#> → finding duplicate expressions in saem predOnly model 1...
#> → optimizing duplicate expressions in saem predOnly model 1...
#> → finding duplicate expressions in saem predOnly model 2...
#> ✔ done
#> using C compiler: 'gcc.exe (GCC) 13.2.0'
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 5952
#> → compress phiM in nlmixr2 object, save 3712
#> → compress parHistData in nlmixr2 object, save 2456
#> → compress saem0 in nlmixr2 object, save 27920

# In a real bootstrap, nboot would be much higher.
fit2 <- suppressMessages(bootstrapFit(fit, nboot = 5))
fit2
```
