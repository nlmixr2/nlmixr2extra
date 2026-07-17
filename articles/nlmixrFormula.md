# nlmixr2 Algebraic Solutions with Formula

## Introduction

Some models are able to be described by simple algebraic solutions. To
simplify algebraic models, you can use a formula interface similar to
the formula used with the `lme4` library. (The formula interface is
described in more detail below; knowledge of `lme4` is not required to
use the formula interface.)

## Quick start

The simplest, non-trivial model is a linear model, `y = m*x + b`. We
will generate the data for this model and then fit it to show the
simplest use case of the nlmixr2 formula interface.

``` r

library(nlmixr2extra)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
```

``` r

# Simulate the equation y = 3*x + 5 with normally-distributed residual error
# having a standard deviation of 5
withr::with_seed(
  5, # standardize the random seed so that the results are reproducible
  dSim <-
    data.frame(
      x = 1:100,
      y =
        (1:100)*3 +
        5 +
        rnorm(n = 100, mean = 0, sd = 5)
    )
)

ggplot(dSim, aes(x=x, y=y)) + geom_point()
```

![](nlmixrFormula_files/figure-html/quickstart-data-gen-1.png)

To fit this model requires only one line of R code:

``` r

mod <-
  nlmixrFormula(
    y~m*x + b,
    data = dSim,
    start = c(m = 2.5, b = 4, addSd = 2),
    est = "focei"
  )
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
#> Key: U: Unscaled Parameters; X: Back-transformed parameters; G: Gill difference gradient approximation
#> F: Forward difference gradient approximation
#> C: Central difference gradient approximation
#> M: Mixed forward and central difference gradient approximation
#> A: Analytic (forward sensitivity) gradient (fast=TRUE)
#> Unscaled parameters for Omegas=chol(solve(omega));
#> Diagonals are transformed, as specified by foceiControl(diagXform=)
#> 
#> |    #| Function Val. |         m |         b |     addSd |
#> |-----+---------------+-----------+-----------+-----------|
#> |    1|     23328.257 |   -0.5000 |     1.000 |    -1.000 |
#> |    U|               |     2.500 |     4.000 |     2.000 |
#> |    X|               |     2.500 |     4.000 |     2.000 |
#> |    G|      Gill     |-3.500e+04 |    -5281. |-2.309e+04 |
#> |    2|     2192.0220 |    0.3282 |     1.125 |   -0.4537 |
#> |    U|               |     2.831 |     4.500 |     2.546 |
#> |    X|               |     2.831 |     4.500 |     2.546 |
#> |    F|    Forward    |    -7449. |    -1132. |    -1496. |
#> |    3|     2495.1916 |     1.105 |     1.245 |    -1.072 |
#> |    U|               |     3.142 |     4.979 |     1.928 |
#> |    X|               |     3.142 |     4.979 |     1.928 |
#> |    4|     565.67790 |    0.7625 |     1.192 |   -0.6244 |
#> |    U|               |     3.005 |     4.766 |     2.376 |
#> |    X|               |     3.005 |     4.766 |     2.376 |
#> |    F|    Forward    |    -34.29 |    -19.26 |    -246.3 |
#> |    5|     492.12340 |    0.7188 |     1.249 |   -0.1585 |
#> |    U|               |     2.988 |     4.994 |     2.842 |
#> |    X|               |     2.988 |     4.994 |     2.842 |
#> |    F|    Forward    |    -496.5 |    -78.41 |    -129.0 |
#> |    6|     645.12768 |    0.9315 |     1.388 |    0.2385 |
#> |    U|               |     3.073 |     5.551 |     3.239 |
#> |    X|               |     3.073 |     5.551 |     3.239 |
#> |    7|     524.94635 |    0.8367 |     1.267 |   -0.1278 |
#> |    U|               |     3.035 |     5.069 |     2.872 |
#> |    X|               |     3.035 |     5.069 |     2.872 |
#> |    8|     481.66227 |    0.7577 |     1.255 |   -0.1484 |
#> |    U|               |     3.003 |     5.019 |     2.852 |
#> |    X|               |     3.003 |     5.019 |     2.852 |
#> |    F|    Forward    |     38.11 |     2.000 |    -120.7 |
#> |    9|     477.66097 |    0.7455 |     1.254 |   -0.1095 |
#> |    U|               |     2.998 |     5.016 |     2.890 |
#> |    X|               |     2.998 |     5.016 |     2.890 |
#> |   10|     475.66763 |    0.7303 |     1.253 |  -0.06159 |
#> |    U|               |     2.992 |     5.013 |     2.938 |
#> |    X|               |     2.992 |     5.013 |     2.938 |
#> |    F|    Forward    |    -310.4 |    -49.93 |    -109.0 |
#> 
#> |    #| Function Val. |         m |         b |     addSd |
#> |-----+---------------+-----------+-----------+-----------|
#> |   11|     464.93880 |    0.7658 |     1.209 |  0.009492 |
#> |    U|               |     3.006 |     4.836 |     3.009 |
#> |    X|               |     3.006 |     4.836 |     3.009 |
#> |    F|    Forward    |     49.29 |   0.05405 |    -96.08 |
#> |   12|     460.55958 |    0.7361 |     1.269 |   0.07073 |
#> |    U|               |     2.994 |     5.077 |     3.071 |
#> |    X|               |     2.994 |     5.077 |     3.071 |
#> |    F|    Forward    |    -191.0 |    -30.47 |    -88.69 |
#> |   13|     453.92531 |    0.7470 |     1.335 |    0.1323 |
#> |    U|               |     2.999 |     5.342 |     3.132 |
#> |    X|               |     2.999 |     5.342 |     3.132 |
#> |    F|    Forward    |     45.98 |     10.30 |    -80.17 |
#> |   14|     450.88505 |    0.7427 |     1.258 |    0.1798 |
#> |    U|               |     2.997 |     5.032 |     3.180 |
#> |    X|               |     2.997 |     5.032 |     3.180 |
#> |    F|    Forward    |    -125.0 |    -21.38 |    -75.17 |
#> |   15|     448.79046 |    0.7758 |     1.179 |    0.2109 |
#> |    U|               |     3.010 |     4.716 |     3.211 |
#> |    X|               |     3.010 |     4.716 |     3.211 |
#> |    F|    Forward    |     101.2 |     6.398 |    -71.92 |
#> |   16|     443.82733 |    0.7559 |     1.245 |    0.2707 |
#> |    U|               |     3.002 |     4.979 |     3.271 |
#> |    X|               |     3.002 |     4.979 |     3.271 |
#> |   17|     437.85783 |    0.7449 |     1.405 |    0.3880 |
#> |    U|               |     2.998 |     5.618 |     3.388 |
#> |    X|               |     2.998 |     5.618 |     3.388 |
#> |   18|     432.55563 |    0.7119 |     1.646 |    0.5772 |
#> |    U|               |     2.985 |     6.583 |     3.577 |
#> |    X|               |     2.985 |     6.583 |     3.577 |
#> |    F|    Forward    |     130.2 |     41.21 |    -43.40 |
#> |   19|     415.46918 |    0.7883 |    0.8637 |     1.475 |
#> |    U|               |     3.015 |     3.455 |     4.475 |
#> |    X|               |     3.015 |     3.455 |     4.475 |
#> |    F|    Forward    |    -135.1 |    -37.06 |    -7.035 |
#> |   20|     410.69540 |    0.7633 |     1.131 |     1.529 |
#> |    U|               |     3.005 |     4.524 |     4.529 |
#> |    X|               |     3.005 |     4.524 |     4.529 |
#> |    F|    Forward    |    -52.92 |    -14.12 |    -3.790 |
#> 
#> |    #| Function Val. |         m |         b |     addSd |
#> |-----+---------------+-----------+-----------+-----------|
#> |   21|     409.73210 |    0.7492 |     1.288 |     1.599 |
#> |    U|               |     3.000 |     5.151 |     4.599 |
#> |    X|               |     3.000 |     5.151 |     4.599 |
#> |    F|    Forward    |    -4.085 |   -0.7983 |    -1.982 |
#> |   22|     409.65712 |    0.7483 |     1.300 |     1.648 |
#> |    U|               |     2.999 |     5.198 |     4.648 |
#> |    X|               |     2.999 |     5.198 |     4.648 |
#> |    F|    Forward    |    0.7411 |    0.3598 |    -1.024 |
#> |   23|     409.62889 |    0.7489 |     1.296 |     1.697 |
#> |    U|               |     3.000 |     5.186 |     4.697 |
#> |    X|               |     3.000 |     5.186 |     4.697 |
#> |    F|    Forward    |     1.075 |    0.2922 |   -0.1041 |
#> |   24|     409.62838 |    0.7492 |     1.293 |     1.702 |
#> |    U|               |     3.000 |     5.173 |     4.702 |
#> |    X|               |     3.000 |     5.173 |     4.702 |
#> |    F|    Forward    |    0.2724 |   0.06256 | -0.002335 |
#> |   25|     409.62838 |    0.7492 |     1.293 |     1.702 |
#> |    U|               |     3.000 |     5.173 |     4.702 |
#> |    X|               |     3.000 |     5.173 |     4.702 |
#> calculating covariance matrix
#> covType="analytic": a model with no random effects is out of analytic-covariance scope; using the finite-difference covariance instead
#> covType="analytic" not available for this model (out of scope, or the augmented model would not build/solve); using the finite-difference sandwich ("r,s") covariance.
#> done
#> → Calculating residuals/tables
#> ✔ done
```

You can also use a call to `nlmixr` or `nlmixr2`, the arguments are the
same as `nlmixrFormula`:

``` r

mod <-
  nlmixr(
    y ~ m*x + b,
    data = dSim,
    start = c(m = 2.5, b = 4, addSd = 2),
    est = "focei"
  )
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> Key: U: Unscaled Parameters; X: Back-transformed parameters; G: Gill difference gradient approximation
#> F: Forward difference gradient approximation
#> C: Central difference gradient approximation
#> M: Mixed forward and central difference gradient approximation
#> A: Analytic (forward sensitivity) gradient (fast=TRUE)
#> Unscaled parameters for Omegas=chol(solve(omega));
#> Diagonals are transformed, as specified by foceiControl(diagXform=)
#> 
#> |    #| Function Val. |         m |         b |     addSd |
#> |-----+---------------+-----------+-----------+-----------|
#> |    1|     23328.257 |   -0.5000 |     1.000 |    -1.000 |
#> |    U|               |     2.500 |     4.000 |     2.000 |
#> |    X|               |     2.500 |     4.000 |     2.000 |
#> |    G|      Gill     |-3.500e+04 |    -5281. |-2.309e+04 |
#> |    2|     2192.0220 |    0.3282 |     1.125 |   -0.4537 |
#> |    U|               |     2.831 |     4.500 |     2.546 |
#> |    X|               |     2.831 |     4.500 |     2.546 |
#> |    F|    Forward    |    -7449. |    -1132. |    -1496. |
#> |    3|     2495.1916 |     1.105 |     1.245 |    -1.072 |
#> |    U|               |     3.142 |     4.979 |     1.928 |
#> |    X|               |     3.142 |     4.979 |     1.928 |
#> |    4|     565.67790 |    0.7625 |     1.192 |   -0.6244 |
#> |    U|               |     3.005 |     4.766 |     2.376 |
#> |    X|               |     3.005 |     4.766 |     2.376 |
#> |    F|    Forward    |    -34.29 |    -19.26 |    -246.3 |
#> |    5|     492.12340 |    0.7188 |     1.249 |   -0.1585 |
#> |    U|               |     2.988 |     4.994 |     2.842 |
#> |    X|               |     2.988 |     4.994 |     2.842 |
#> |    F|    Forward    |    -496.5 |    -78.41 |    -129.0 |
#> |    6|     645.12768 |    0.9315 |     1.388 |    0.2385 |
#> |    U|               |     3.073 |     5.551 |     3.239 |
#> |    X|               |     3.073 |     5.551 |     3.239 |
#> |    7|     524.94635 |    0.8367 |     1.267 |   -0.1278 |
#> |    U|               |     3.035 |     5.069 |     2.872 |
#> |    X|               |     3.035 |     5.069 |     2.872 |
#> |    8|     481.66227 |    0.7577 |     1.255 |   -0.1484 |
#> |    U|               |     3.003 |     5.019 |     2.852 |
#> |    X|               |     3.003 |     5.019 |     2.852 |
#> |    F|    Forward    |     38.11 |     2.000 |    -120.7 |
#> |    9|     477.66097 |    0.7455 |     1.254 |   -0.1095 |
#> |    U|               |     2.998 |     5.016 |     2.890 |
#> |    X|               |     2.998 |     5.016 |     2.890 |
#> |   10|     475.66763 |    0.7303 |     1.253 |  -0.06159 |
#> |    U|               |     2.992 |     5.013 |     2.938 |
#> |    X|               |     2.992 |     5.013 |     2.938 |
#> |    F|    Forward    |    -310.4 |    -49.93 |    -109.0 |
#> 
#> |    #| Function Val. |         m |         b |     addSd |
#> |-----+---------------+-----------+-----------+-----------|
#> |   11|     464.93880 |    0.7658 |     1.209 |  0.009492 |
#> |    U|               |     3.006 |     4.836 |     3.009 |
#> |    X|               |     3.006 |     4.836 |     3.009 |
#> |    F|    Forward    |     49.29 |   0.05405 |    -96.08 |
#> |   12|     460.55958 |    0.7361 |     1.269 |   0.07073 |
#> |    U|               |     2.994 |     5.077 |     3.071 |
#> |    X|               |     2.994 |     5.077 |     3.071 |
#> |    F|    Forward    |    -191.0 |    -30.47 |    -88.69 |
#> |   13|     453.92531 |    0.7470 |     1.335 |    0.1323 |
#> |    U|               |     2.999 |     5.342 |     3.132 |
#> |    X|               |     2.999 |     5.342 |     3.132 |
#> |    F|    Forward    |     45.98 |     10.30 |    -80.17 |
#> |   14|     450.88505 |    0.7427 |     1.258 |    0.1798 |
#> |    U|               |     2.997 |     5.032 |     3.180 |
#> |    X|               |     2.997 |     5.032 |     3.180 |
#> |    F|    Forward    |    -125.0 |    -21.38 |    -75.17 |
#> |   15|     448.79046 |    0.7758 |     1.179 |    0.2109 |
#> |    U|               |     3.010 |     4.716 |     3.211 |
#> |    X|               |     3.010 |     4.716 |     3.211 |
#> |    F|    Forward    |     101.2 |     6.398 |    -71.92 |
#> |   16|     443.82733 |    0.7559 |     1.245 |    0.2707 |
#> |    U|               |     3.002 |     4.979 |     3.271 |
#> |    X|               |     3.002 |     4.979 |     3.271 |
#> |   17|     437.85783 |    0.7449 |     1.405 |    0.3880 |
#> |    U|               |     2.998 |     5.618 |     3.388 |
#> |    X|               |     2.998 |     5.618 |     3.388 |
#> |   18|     432.55563 |    0.7119 |     1.646 |    0.5772 |
#> |    U|               |     2.985 |     6.583 |     3.577 |
#> |    X|               |     2.985 |     6.583 |     3.577 |
#> |    F|    Forward    |     130.2 |     41.21 |    -43.40 |
#> |   19|     415.46918 |    0.7883 |    0.8637 |     1.475 |
#> |    U|               |     3.015 |     3.455 |     4.475 |
#> |    X|               |     3.015 |     3.455 |     4.475 |
#> |    F|    Forward    |    -135.1 |    -37.06 |    -7.035 |
#> |   20|     410.69540 |    0.7633 |     1.131 |     1.529 |
#> |    U|               |     3.005 |     4.524 |     4.529 |
#> |    X|               |     3.005 |     4.524 |     4.529 |
#> |    F|    Forward    |    -52.92 |    -14.12 |    -3.790 |
#> 
#> |    #| Function Val. |         m |         b |     addSd |
#> |-----+---------------+-----------+-----------+-----------|
#> |   21|     409.73210 |    0.7492 |     1.288 |     1.599 |
#> |    U|               |     3.000 |     5.151 |     4.599 |
#> |    X|               |     3.000 |     5.151 |     4.599 |
#> |    F|    Forward    |    -4.085 |   -0.7983 |    -1.982 |
#> |   22|     409.65712 |    0.7483 |     1.300 |     1.648 |
#> |    U|               |     2.999 |     5.198 |     4.648 |
#> |    X|               |     2.999 |     5.198 |     4.648 |
#> |    F|    Forward    |    0.7411 |    0.3598 |    -1.024 |
#> |   23|     409.62889 |    0.7489 |     1.296 |     1.697 |
#> |    U|               |     3.000 |     5.186 |     4.697 |
#> |    X|               |     3.000 |     5.186 |     4.697 |
#> |    F|    Forward    |     1.075 |    0.2922 |   -0.1041 |
#> |   24|     409.62838 |    0.7492 |     1.293 |     1.702 |
#> |    U|               |     3.000 |     5.173 |     4.702 |
#> |    X|               |     3.000 |     5.173 |     4.702 |
#> |    F|    Forward    |    0.2724 |   0.06256 | -0.002335 |
#> |   25|     409.62838 |    0.7492 |     1.293 |     1.702 |
#> |    U|               |     3.000 |     5.173 |     4.702 |
#> |    X|               |     3.000 |     5.173 |     4.702 |
#> calculating covariance matrix
#> [====
#> covType="analytic": a model with no random effects is out of analytic-covariance scope; using the finite-difference covariance instead
#> covType="analytic" not available for this model (out of scope, or the augmented model would not build/solve); using the finite-difference sandwich ("r,s") covariance.
#> |====|====|====|====|====|====|====|====|====] 0:00:00 
#> done
#> → Calculating residuals/tables
#> ✔ done
```

### Quick start, mixed-effects model

A more typical use case for nlmixr2 is a mixed-effects model. It adds
(subject-level) random effects to the model.

``` r

# Setup the dataset for nonlinear mixed-effects fitting

# Simulate the equation y = 3*x + 5 with normally-distributed residual error
# having a standard deviation of 5
withr::with_seed(
  5, # standardize the random seed so that the results are reproducible
  {
    dSimSetup <-
      data.frame(
        id = rep(1:10, each=10),
        x = rep(1:10, 10),
        y = 1
      )
    dSimNlmePrep <-
      nlmixrFormula(
        y~m*x + b + bRe ~ bRe|id,
        start = c(m=3, b=5, addSd=5),
        data = dSimSetup,
        est = "rxSolve"
      )
  }
)
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments

# Modify the simulated results to be ready for fitting
dSimNlme <-
  dSimNlmePrep |>
  as.data.frame() |>
  select(id, time, x, y=sim)

ggplot(dSimNlme, aes(x=x, y=y, colour=factor(id))) + geom_point()
```

![](nlmixrFormula_files/figure-html/quickstart-data-gen-nlme-1.png)

``` r

mod <-
  nlmixr(
    y~m*x + b + bRe ~ bRe|id,
    start = c(m=3, b=5, addSd=5),
    data = dSimNlme,
    est = "focei"
  )
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate ∂(f)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate ∂(R²)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in inner model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> 
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling inner model...
#> ✔ done
#> → finding duplicate expressions in FD model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling EBE model...
#> ✔ done
#> → compiling events FD model...
#> ✔ done
#> Key: U: Unscaled Parameters; X: Back-transformed parameters; G: Gill difference gradient approximation
#> F: Forward difference gradient approximation
#> C: Central difference gradient approximation
#> M: Mixed forward and central difference gradient approximation
#> A: Analytic (forward sensitivity) gradient (fast=TRUE)
#> Unscaled parameters for Omegas=chol(solve(omega));
#> Diagonals are transformed, as specified by foceiControl(diagXform=)
#> 
#> |    #| Function Val. |         m |         b |     addSd |        o1 |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |    1|     437.57349 |     0.000 |     1.000 |     1.000 |    -1.000 |
#> |    U|               |     3.000 |     5.000 |     5.000 |     1.000 |
#> |    X|               |     3.000 |     5.000 |     5.000 |     1.000 |
#> |    G|      Gill     |     9.554 |     6.246 |    -10.84 |     5.831 |
#> |    2|     452.72409 |   -0.5710 |    0.6267 |     1.648 |    -1.331 |
#> |    U|               |     2.810 |     3.134 |     6.619 |    0.6687 |
#> |    X|               |     2.810 |     3.134 |     6.619 |    0.6687 |
#> |    3|     436.52224 |   -0.1490 |    0.9026 |     1.169 |    -1.091 |
#> |    U|               |     2.950 |     4.513 |     5.422 |    0.9091 |
#> |    X|               |     2.950 |     4.513 |     5.422 |    0.9091 |
#> |    F|    Forward    |    0.7417 |    -12.29 |     5.993 |     3.770 |
#> |    4|     435.21948 |   -0.2334 |     1.114 |     1.122 |    -1.212 |
#> |    U|               |     2.922 |     5.569 |     5.306 |    0.7884 |
#> |    X|               |     2.922 |     5.569 |     5.306 |    0.7884 |
#> |    F|    Forward    |     7.375 |     6.659 |     4.755 |  -0.05658 |
#> |    5|     434.56526 |   -0.4295 |     1.027 |    0.9755 |    -1.241 |
#> |    U|               |     2.857 |     5.133 |     4.939 |    0.7586 |
#> |    X|               |     2.857 |     5.133 |     4.939 |    0.7586 |
#> |    F|    Forward    |     1.294 |    -7.986 |    -7.182 | -0.001548 |
#> |    6|     433.37875 |   -0.6191 |     1.168 |     1.078 |    -1.197 |
#> |    U|               |     2.794 |     5.842 |     5.194 |    0.8031 |
#> |    X|               |     2.794 |     5.842 |     5.194 |    0.8031 |
#> |    F|    Forward    |     1.999 |    -1.449 |     2.547 |    0.7527 |
#> |    7|     433.99510 |   -0.7714 |     1.282 |    0.9054 |    -1.249 |
#> |    U|               |     2.743 |     6.410 |     4.763 |    0.7508 |
#> |    X|               |     2.743 |     6.410 |     4.763 |    0.7508 |
#> |    8|     433.23599 |   -0.6628 |     1.200 |     1.022 |    -1.213 |
#> |    U|               |     2.779 |     6.001 |     5.055 |    0.7866 |
#> |    X|               |     2.779 |     6.001 |     5.055 |    0.7866 |
#> |    F|    Forward    |     2.389 |   0.08358 |    -2.012 |    0.4763 |
#> |    9|     433.15874 |   -0.7227 |     1.198 |     1.072 |    -1.225 |
#> |    U|               |     2.759 |     5.990 |     5.181 |    0.7747 |
#> |    X|               |     2.759 |     5.990 |     5.181 |    0.7747 |
#> |    F|    Forward    |     1.057 |    -2.125 |     2.671 |   -0.8747 |
#> |   10|     433.01948 |   -0.7738 |     1.235 |     1.038 |    -1.193 |
#> |    U|               |     2.742 |     6.177 |     5.094 |    0.8072 |
#> |    X|               |     2.742 |     6.177 |     5.094 |    0.8072 |
#> |    F|    Forward    |     1.373 |   -0.4540 |   -0.6416 |     1.473 |
#> 
#> |    #| Function Val. |         m |         b |     addSd |        o1 |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   11|     432.98141 |   -0.8251 |     1.255 |     1.052 |    -1.248 |
#> |    U|               |     2.725 |     6.274 |     5.129 |    0.7518 |
#> |    X|               |     2.725 |     6.274 |     5.129 |    0.7518 |
#> |    F|    Forward    |     1.034 |   -0.3445 |     1.340 |    -2.465 |
#> |   12|     432.89627 |   -0.8926 |     1.263 |     1.051 |    -1.207 |
#> |    U|               |     2.702 |     6.315 |     5.127 |    0.7925 |
#> |    X|               |     2.702 |     6.315 |     5.127 |    0.7925 |
#> |    F|    Forward    |  -0.03682 |    -1.986 |    0.8318 |    0.5240 |
#> |   13|     432.98973 |   -0.9012 |     1.331 |     1.012 |    -1.220 |
#> |    U|               |     2.700 |     6.654 |     5.031 |    0.7800 |
#> |    X|               |     2.700 |     6.654 |     5.031 |    0.7800 |
#> |   14|     432.87032 |   -0.8922 |     1.286 |     1.041 |    -1.214 |
#> |    U|               |     2.703 |     6.430 |     5.103 |    0.7865 |
#> |    X|               |     2.703 |     6.430 |     5.103 |    0.7865 |
#> |    F|    Forward    |    0.7953 |    0.2577 |   0.07213 |    0.2201 |
#> |   15|     432.86532 |   -0.9158 |     1.278 |     1.039 |    -1.220 |
#> |    U|               |     2.695 |     6.392 |     5.097 |    0.7799 |
#> |    X|               |     2.695 |     6.392 |     5.097 |    0.7799 |
#> |    F|    Forward    |   0.06612 |    -1.282 |  -0.02689 |   -0.1587 |
#> |   16|     432.84939 |   -0.9330 |     1.296 |     1.045 |    -1.218 |
#> |    U|               |     2.689 |     6.482 |     5.113 |    0.7820 |
#> |    X|               |     2.689 |     6.482 |     5.113 |    0.7820 |
#> |    F|    Forward    |    0.3433 |   -0.1879 |    0.5112 |   -0.1240 |
#> |   17|     432.85323 |   -0.9465 |     1.304 |     1.025 |    -1.213 |
#> |    U|               |     2.685 |     6.519 |     5.063 |    0.7868 |
#> |    X|               |     2.685 |     6.519 |     5.063 |    0.7868 |
#> |   18|     432.84592 |   -0.9385 |     1.299 |     1.037 |    -1.216 |
#> |    U|               |     2.687 |     6.497 |     5.092 |    0.7839 |
#> |    X|               |     2.687 |     6.497 |     5.092 |    0.7839 |
#> |    F|    Forward    |    0.3419 |  -0.09513 |   -0.2342 |    0.1077 |
#> |   19|     432.84366 |   -0.9466 |     1.302 |     1.043 |    -1.219 |
#> |    U|               |     2.684 |     6.508 |     5.106 |    0.7814 |
#> |    X|               |     2.684 |     6.508 |     5.106 |    0.7814 |
#> |    F|    Forward    |    0.2579 |   -0.1625 |    0.2968 |   -0.1315 |
#> |   20|     432.84067 |   -0.9549 |     1.305 |     1.038 |    -1.216 |
#> |    U|               |     2.682 |     6.527 |     5.095 |    0.7840 |
#> |    X|               |     2.682 |     6.527 |     5.095 |    0.7840 |
#> |    F|    Forward    |    0.2244 |  -0.09928 |   -0.1442 |   0.09647 |
#> 
#> |    #| Function Val. |         m |         b |     addSd |        o1 |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   21|     432.83885 |   -0.9640 |     1.309 |     1.041 |    -1.218 |
#> |    U|               |     2.679 |     6.544 |     5.103 |    0.7819 |
#> |    X|               |     2.679 |     6.544 |     5.103 |    0.7819 |
#> |    F|    Forward    |    0.1669 |  -0.07973 |    0.1726 |  -0.07722 |
#> |   22|     432.83715 |   -0.9737 |     1.312 |     1.040 |    -1.217 |
#> |    U|               |     2.675 |     6.561 |     5.100 |    0.7828 |
#> |    X|               |     2.675 |     6.561 |     5.100 |    0.7828 |
#> |    F|    Forward    |   0.09230 |  -0.09577 |   0.06487 | -0.005179 |
#> |   23|     432.83715 |   -0.9737 |     1.312 |     1.040 |    -1.217 |
#> |    U|               |     2.675 |     6.561 |     5.100 |    0.7828 |
#> |    X|               |     2.675 |     6.561 |     5.100 |    0.7828 |
#> calculating covariance matrix
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> 
#> covType="analytic" not available for this model (out of scope, or the augmented model would not build/solve); using the finite-difference sandwich ("r,s") covariance.
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> done
#> → Calculating residuals/tables
#> ✔ done

mod
```

``` math
\begin{align*}
{value} & = {m} {\times} {x}+{b}+{bRe} \\
{value} & \sim add({addSd})
\end{align*}
```

In this model, the fixed effects of `m` and `b` were estimated along
with the random effect of `bRe`.

### Quick start, automatic parameter fixed effects

``` r

# Setup the dataset for nonlinear mixed-effects fitting with different types of
# parameters.

# Simulate the equation y = 3*x + 4(when z = 'a') or 6(when z = 'b') with normally-distributed residual error
# having a standard deviation of 5
withr::with_seed(
  5, # standardize the random seed so that the results are reproducible
  {
    dSimSetup <-
      data.frame(
        id = rep(1:10, each=10),
        x = rep(1:10, 10),
        y = 1,
        z = sample(factor(c("a", "b")), size = 100, replace = TRUE)
      )
    # Note that we need to give `start` as a list so that `b` can carry one
    # starting value per factor level of `z` (named-vector `c()` cannot hold
    # multi-element entries).
    dSimNlmePrep <-
      nlmixrFormula(
        y~m*x + b + bRe ~ bRe|id,
        start = list(m=3, b=c(4, 2), addSd=5),
        param = list(b ~ z),
        data = dSimSetup,
        est = "rxSolve"
      )
  }
)
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments

# Modify the simulated results to be ready for fitting
dSimNlme <-
  dSimNlmePrep |>
  as.data.frame() |>
  select(id, time, x, y=sim, z)

ggplot(dSimNlme, aes(x=x, y=y, colour=z)) +
  geom_point() +
  geom_line(aes(group = id), colour = "gray")
```

![](nlmixrFormula_files/figure-html/quickstart-data-gen-param-1.png)

``` r

mod <-
  nlmixr(
    y~m*x + b + bRe ~ bRe|id,
    start = list(m=3, b=c(4, 2), addSd=5),
    param = list(b ~ z),
    data = dSimNlme,
    est = "focei"
  )
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate ∂(f)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate ∂(R²)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in inner model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> 
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling inner model...
#> ✔ done
#> → finding duplicate expressions in FD model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling EBE model...
#> ✔ done
#> → compiling events FD model...
#> ✔ done
#> Key: U: Unscaled Parameters; X: Back-transformed parameters; G: Gill difference gradient approximation
#> F: Forward difference gradient approximation
#> C: Central difference gradient approximation
#> M: Mixed forward and central difference gradient approximation
#> A: Analytic (forward sensitivity) gradient (fast=TRUE)
#> Unscaled parameters for Omegas=chol(solve(omega));
#> Diagonals are transformed, as specified by foceiControl(diagXform=)
#> 
#> |    #| Function Val. |         m |     b.z.b |     b.z.a |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |    1|     416.24007 |     0.000 |    0.5000 |   -0.5000 |     1.000 |
#> |.....................|    -1.000 |...........|...........|...........|
#> |    U|               |     3.000 |     4.000 |     2.000 |     5.000 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |     3.000 |     4.000 |     2.000 |     5.000 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    G|      Gill     |    -7.468 |    -14.93 |     5.862 |     9.291 |
#> |.....................|     1.034 |...........|...........|...........|
#> |    2|     447.77697 |    0.3732 |     1.246 |   -0.7930 |    0.5357 |
#> |.....................|    -1.052 |...........|...........|...........|
#> |    U|               |     3.124 |     6.984 |     1.414 |     3.839 |
#> |.....................|    0.9483 |...........|...........|...........|
#> |    X|               |     3.124 |     6.984 |     1.414 |     3.839 |
#> |.....................|    0.9483 |...........|...........|...........|
#> |    3|     413.85631 |   0.07244 |    0.6448 |   -0.5569 |    0.9099 |
#> |.....................|    -1.010 |...........|...........|...........|
#> |    U|               |     3.024 |     4.579 |     1.886 |     4.775 |
#> |.....................|    0.9900 |...........|...........|...........|
#> |    X|               |     3.024 |     4.579 |     1.886 |     4.775 |
#> |.....................|    0.9900 |...........|...........|...........|
#> |    F|    Forward    |  -0.01479 |    0.1062 |     9.751 |     2.734 |
#> |.....................|   0.06949 |...........|...........|...........|
#> |    4|     412.20435 |   0.07272 |    0.6428 |   -0.7437 |    0.8575 |
#> |.....................|    -1.011 |...........|...........|...........|
#> |    U|               |     3.024 |     4.571 |     1.513 |     4.644 |
#> |.....................|    0.9886 |...........|...........|...........|
#> |    X|               |     3.024 |     4.571 |     1.513 |     4.644 |
#> |.....................|    0.9886 |...........|...........|...........|
#> |    5|     411.09891 |   0.07332 |    0.6384 |    -1.142 |    0.7458 |
#> |.....................|    -1.014 |...........|...........|...........|
#> |    U|               |     3.024 |     4.554 |    0.7157 |     4.364 |
#> |.....................|    0.9858 |...........|...........|...........|
#> |    X|               |     3.024 |     4.554 |    0.7157 |     4.364 |
#> |.....................|    0.9858 |...........|...........|...........|
#> |    F|    Forward    |    -6.768 |    -15.69 |     1.906 |    -13.34 |
#> |.....................|     1.265 |...........|...........|...........|
#> |    6|     409.91678 |    0.1880 |    0.9143 |    -1.600 |     1.011 |
#> |.....................|    -1.041 |...........|...........|...........|
#> |    U|               |     3.063 |     5.657 |   -0.1994 |     5.027 |
#> |.....................|    0.9588 |...........|...........|...........|
#> |    X|               |     3.063 |     5.657 |   -0.1994 |     5.027 |
#> |.....................|    0.9588 |...........|...........|...........|
#> |    F|    Forward    |     4.340 |     6.325 |     2.403 |     15.94 |
#> |.....................|    -3.644 |...........|...........|...........|
#> |    7|     408.25267 | -0.004334 |    0.9956 |    -2.055 |    0.9125 |
#> |.....................|   -0.7115 |...........|...........|...........|
#> |    U|               |     2.999 |     5.982 |    -1.111 |     4.781 |
#> |.....................|     1.289 |...........|...........|...........|
#> |    X|               |     2.999 |     5.982 |    -1.111 |     4.781 |
#> |.....................|     1.289 |...........|...........|...........|
#> |    F|    Forward    |    -2.008 |    -5.176 |    -3.594 |     6.352 |
#> |.....................|   -0.3549 |...........|...........|...........|
#> |    8|     420.44363 |  -0.07216 |     1.285 |    -1.703 |    0.6003 |
#> |.....................|   -0.4673 |...........|...........|...........|
#> |    U|               |     2.976 |     7.140 |   -0.4063 |     4.001 |
#> |.....................|     1.533 |...........|...........|...........|
#> |    X|               |     2.976 |     7.140 |   -0.4063 |     4.001 |
#> |.....................|     1.533 |...........|...........|...........|
#> |    9|     407.89508 |   0.01416 |     1.043 |    -2.022 |    0.8540 |
#> |.....................|   -0.7082 |...........|...........|...........|
#> |    U|               |     3.005 |     6.173 |    -1.044 |     4.635 |
#> |.....................|     1.292 |...........|...........|...........|
#> |    X|               |     3.005 |     6.173 |    -1.044 |     4.635 |
#> |.....................|     1.292 |...........|...........|...........|
#> |    F|    Forward    |     1.766 |     2.716 |    -1.576 |    0.3819 |
#> |.....................|   -0.1594 |...........|...........|...........|
#> |   10|     407.77219 |  -0.03860 |     1.012 |    -1.967 |    0.8471 |
#> |.....................|   -0.6916 |...........|...........|...........|
#> |    U|               |     2.987 |     6.046 |   -0.9343 |     4.618 |
#> |.....................|     1.308 |...........|...........|...........|
#> |    X|               |     2.987 |     6.046 |   -0.9343 |     4.618 |
#> |.....................|     1.308 |...........|...........|...........|
#> |    F|    Forward    |    -1.240 |    -2.774 |    -2.357 |   -0.3235 |
#> |.....................|  -0.05781 |...........|...........|...........|
#> 
#> |    #| Function Val. |         m |     b.z.b |     b.z.a |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |   11|     407.66169 |  -0.05447 |     1.035 |    -1.907 |    0.8537 |
#> |.....................|   -0.6402 |...........|...........|...........|
#> |    U|               |     2.982 |     6.139 |   -0.8131 |     4.634 |
#> |.....................|     1.360 |...........|...........|...........|
#> |    X|               |     2.982 |     6.139 |   -0.8131 |     4.634 |
#> |.....................|     1.360 |...........|...........|...........|
#> |    F|    Forward    |    0.3929 |     1.185 |   -0.7976 |    0.3676 |
#> |.....................|   0.02647 |...........|...........|...........|
#> |   12|     407.61544 | -0.004091 |    0.9929 |    -1.856 |    0.8498 |
#> |.....................|   -0.6235 |...........|...........|...........|
#> |    U|               |     2.999 |     5.971 |   -0.7118 |     4.625 |
#> |.....................|     1.377 |...........|...........|...........|
#> |    X|               |     2.999 |     5.971 |   -0.7118 |     4.625 |
#> |.....................|     1.377 |...........|...........|...........|
#> |    F|    Forward    |    0.3857 |    0.3121 |   -0.5766 |  -0.06291 |
#> |.....................|   0.08571 |...........|...........|...........|
#> |   13|     407.59698 |  -0.06284 |    0.9984 |    -1.799 |    0.8439 |
#> |.....................|   -0.6430 |...........|...........|...........|
#> |    U|               |     2.979 |     5.994 |   -0.5978 |     4.610 |
#> |.....................|     1.357 |...........|...........|...........|
#> |    X|               |     2.979 |     5.994 |   -0.5978 |     4.610 |
#> |.....................|     1.357 |...........|...........|...........|
#> |    F|    Forward    |   -0.6281 |   -0.7060 |   -0.2703 |   -0.6393 |
#> |.....................|    0.1268 |...........|...........|...........|
#> |   14|     407.58271 |  -0.02382 |    0.9823 |    -1.779 |    0.8427 |
#> |.....................|   -0.7077 |...........|...........|...........|
#> |    U|               |     2.992 |     5.929 |   -0.5570 |     4.607 |
#> |.....................|     1.292 |...........|...........|...........|
#> |    X|               |     2.992 |     5.929 |   -0.5570 |     4.607 |
#> |.....................|     1.292 |...........|...........|...........|
#> |    F|    Forward    |   0.09674 |    0.1753 |   0.04940 |   -0.5150 |
#> |.....................|   0.06132 |...........|...........|...........|
#> |   15|     407.58209 |  -0.02066 |    0.9784 |    -1.776 |    0.8510 |
#> |.....................|   -0.7159 |...........|...........|...........|
#> |    U|               |     2.993 |     5.914 |   -0.5514 |     4.628 |
#> |.....................|     1.284 |...........|...........|...........|
#> |    X|               |     2.993 |     5.914 |   -0.5514 |     4.628 |
#> |.....................|     1.284 |...........|...........|...........|
#> |    F|    Forward    |   0.01801 |  -0.05070 |   0.01391 |    0.4726 |
#> |.....................|   0.01250 |...........|...........|...........|
#> |   16|     407.58209 |  -0.02066 |    0.9784 |    -1.776 |    0.8510 |
#> |.....................|   -0.7159 |...........|...........|...........|
#> |    U|               |     2.993 |     5.914 |   -0.5514 |     4.628 |
#> |.....................|     1.284 |...........|...........|...........|
#> |    X|               |     2.993 |     5.914 |   -0.5514 |     4.628 |
#> |.....................|     1.284 |...........|...........|...........|
#> calculating covariance matrix
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> 
#> covType="analytic" not available for this model (out of scope, or the augmented model would not build/solve); using the finite-difference sandwich ("r,s") covariance.
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> done
#> → Calculating residuals/tables
#> ✔ done

mod
```

``` math
\begin{align*}
{b} & = {b.z.b}+{b.z.a} {\times} \left({z}{\equiv}\text{"a"}\right) \\
{value} & = {m} {\times} {x}+{b}+{bRe} \\
{value} & \sim add({addSd})
\end{align*}
```

### Quick start, continuous covariate on a parameter

`param` also accepts numeric (continuous) covariate columns. When `b` is
modelled as a linear function of a continuous covariate `w`, two
parameters are introduced: `pop.b` (the intercept) and `cov_w_b` (the
slope on `w`). The generated model line is `b <- pop.b + cov_w_b * w`.

`start[["b"]]` may be a single value (treated as the intercept; the
slope starts at 0) or `c(intercept, slope)`.

``` r

withr::with_seed(
  5,
  {
    dSimSetup <-
      data.frame(
        id = rep(1:10, each=10),
        x = rep(1:10, 10),
        y = 1,
        w = runif(100, min = 0, max = 10)
      )
    dSimContPrep <-
      nlmixrFormula(
        y~m*x + b + bRe ~ bRe|id,
        start = list(m=3, b=c(4, 0.2), addSd=5),
        param = list(b ~ w),
        data = dSimSetup,
        est = "rxSolve"
      )
  }
)
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments

dSimCont <-
  dSimContPrep |>
  as.data.frame() |>
  select(id, time, x, w, y=sim)
```

``` r

mod <-
  nlmixr(
    y~m*x + b + bRe ~ bRe|id,
    start = list(m=3, b=c(4, 0.2), addSd=5),
    param = list(b ~ w),
    data = dSimCont,
    est = "focei"
  )
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate ∂(f)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate ∂(R²)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in inner model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> 
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling inner model...
#> ✔ done
#> → finding duplicate expressions in FD model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling EBE model...
#> ✔ done
#> → compiling events FD model...
#> ✔ done
#> Key: U: Unscaled Parameters; X: Back-transformed parameters; G: Gill difference gradient approximation
#> F: Forward difference gradient approximation
#> C: Central difference gradient approximation
#> M: Mixed forward and central difference gradient approximation
#> A: Analytic (forward sensitivity) gradient (fast=TRUE)
#> Unscaled parameters for Omegas=chol(solve(omega));
#> Diagonals are transformed, as specified by foceiControl(diagXform=)
#> 
#> |    #| Function Val. |         m |     pop.b |   cov_w_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |    1|     412.28874 |    0.1667 |    0.5833 |    -1.000 |     1.000 |
#> |.....................|   -0.6667 |...........|...........|...........|
#> |    U|               |     3.000 |     4.000 |    0.2000 |     5.000 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |     3.000 |     4.000 |    0.2000 |     5.000 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    G|      Gill     |    -5.350 |    -13.10 |     47.96 |     12.24 |
#> |.....................|    -2.955 |...........|...........|...........|
#> |    2|     3241.9470 |    0.2704 |    0.8374 |    -1.930 |    0.7626 |
#> |.....................|   -0.6094 |...........|...........|...........|
#> |    U|               |     3.035 |     5.016 |    -4.450 |     4.406 |
#> |.....................|     1.057 |...........|...........|...........|
#> |    X|               |     3.035 |     5.016 |    -4.450 |     4.406 |
#> |.....................|     1.057 |...........|...........|...........|
#> |    3|     430.19626 |    0.1770 |    0.6087 |    -1.093 |    0.9763 |
#> |.....................|   -0.6609 |...........|...........|...........|
#> |    U|               |     3.003 |     4.102 |   -0.2650 |     4.941 |
#> |.....................|     1.006 |...........|...........|...........|
#> |    X|               |     3.003 |     4.102 |   -0.2650 |     4.941 |
#> |.....................|     1.006 |...........|...........|...........|
#> |    4|     411.99615 |    0.1678 |    0.5862 |    -1.010 |    0.9973 |
#> |.....................|   -0.6660 |...........|...........|...........|
#> |    U|               |     3.000 |     4.011 |    0.1480 |     4.993 |
#> |.....................|     1.001 |...........|...........|...........|
#> |    X|               |     3.000 |     4.011 |    0.1480 |     4.993 |
#> |.....................|     1.001 |...........|...........|...........|
#> |    F|    Forward    |    -8.058 |    -19.00 |    -8.285 |     12.54 |
#> |.....................|    -1.940 |...........|...........|...........|
#> |    5|     411.77232 |    0.1713 |    0.5945 |    -1.007 |    0.9919 |
#> |.....................|   -0.6652 |...........|...........|...........|
#> |    U|               |     3.002 |     4.045 |    0.1661 |     4.980 |
#> |.....................|     1.001 |...........|...........|...........|
#> |    X|               |     3.002 |     4.045 |    0.1661 |     4.980 |
#> |.....................|     1.001 |...........|...........|...........|
#> |    6|     411.67259 |    0.1759 |    0.6052 |    -1.002 |    0.9848 |
#> |.....................|   -0.6641 |...........|...........|...........|
#> |    U|               |     3.003 |     4.087 |    0.1895 |     4.962 |
#> |.....................|     1.003 |...........|...........|...........|
#> |    X|               |     3.003 |     4.087 |    0.1895 |     4.962 |
#> |.....................|     1.003 |...........|...........|...........|
#> |    F|    Forward    |    -4.829 |    -12.13 |     52.52 |     11.41 |
#> |.....................|    -3.060 |...........|...........|...........|
#> |    7|     410.95731 |    0.1828 |    0.6219 |    -1.016 |    0.9728 |
#> |.....................|   -0.6617 |...........|...........|...........|
#> |    U|               |     3.005 |     4.154 |    0.1225 |     4.932 |
#> |.....................|     1.005 |...........|...........|...........|
#> |    X|               |     3.005 |     4.154 |    0.1225 |     4.932 |
#> |.....................|     1.005 |...........|...........|...........|
#> |    F|    Forward    |    -7.726 |    -18.54 |    -11.69 |     11.24 |
#> |.....................|    -2.017 |...........|...........|...........|
#> |    8|     410.54485 |    0.1910 |    0.6418 |    -1.009 |    0.9607 |
#> |.....................|   -0.6593 |...........|...........|...........|
#> |    U|               |     3.008 |     4.234 |    0.1542 |     4.902 |
#> |.....................|     1.007 |...........|...........|...........|
#> |    X|               |     3.008 |     4.234 |    0.1542 |     4.902 |
#> |.....................|     1.007 |...........|...........|...........|
#> |    F|    Forward    |    -4.927 |    -12.61 |     39.70 |     10.24 |
#> |.....................|    -2.960 |...........|...........|...........|
#> |    9|     409.95434 |    0.1982 |    0.6602 |    -1.021 |    0.9495 |
#> |.....................|   -0.6561 |...........|...........|...........|
#> |    U|               |     3.011 |     4.307 |   0.09644 |     4.874 |
#> |.....................|     1.011 |...........|...........|...........|
#> |    X|               |     3.011 |     4.307 |   0.09644 |     4.874 |
#> |.....................|     1.011 |...........|...........|...........|
#> |    F|    Forward    |    -7.285 |    -17.87 |    -14.33 |     9.950 |
#> |.....................|    -2.116 |...........|...........|...........|
#> |   10|     409.51642 |    0.2063 |    0.6812 |    -1.016 |    0.9388 |
#> |.....................|   -0.6527 |...........|...........|...........|
#> |    U|               |     3.013 |     4.391 |    0.1217 |     4.847 |
#> |.....................|     1.014 |...........|...........|...........|
#> |    X|               |     3.013 |     4.391 |    0.1217 |     4.847 |
#> |.....................|     1.014 |...........|...........|...........|
#> |    F|    Forward    |    -4.746 |    -12.50 |     31.04 |     9.103 |
#> |.....................|    -2.909 |...........|...........|...........|
#> 
#> |    #| Function Val. |         m |     pop.b |   cov_w_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |   11|     409.01288 |    0.2134 |    0.7011 |    -1.026 |    0.9296 |
#> |.....................|   -0.6483 |...........|...........|...........|
#> |    U|               |     3.016 |     4.471 |   0.07008 |     4.824 |
#> |.....................|     1.018 |...........|...........|...........|
#> |    X|               |     3.016 |     4.471 |   0.07008 |     4.824 |
#> |.....................|     1.018 |...........|...........|...........|
#> |    F|    Forward    |    -6.735 |    -16.95 |    -16.21 |     8.836 |
#> |.....................|    -2.219 |...........|...........|...........|
#> |   12|     408.58607 |    0.2208 |    0.7231 |    -1.022 |    0.9211 |
#> |.....................|   -0.6431 |...........|...........|...........|
#> |    U|               |     3.018 |     4.559 |   0.09125 |     4.803 |
#> |.....................|     1.024 |...........|...........|...........|
#> |    X|               |     3.018 |     4.559 |   0.09125 |     4.803 |
#> |.....................|     1.024 |...........|...........|...........|
#> |    F|    Forward    |    -4.345 |    -11.90 |     25.50 |     8.195 |
#> |.....................|    -2.868 |...........|...........|...........|
#> |   13|     408.14986 |    0.2267 |    0.7441 |    -1.031 |    0.9140 |
#> |.....................|   -0.6365 |...........|...........|...........|
#> |    U|               |     3.020 |     4.643 |   0.04445 |     4.785 |
#> |.....................|     1.030 |...........|...........|...........|
#> |    X|               |     3.020 |     4.643 |   0.04445 |     4.785 |
#> |.....................|     1.030 |...........|...........|...........|
#> |    F|    Forward    |    -6.079 |    -15.77 |    -16.77 |     7.996 |
#> |.....................|    -2.292 |...........|...........|...........|
#> |   14|     407.75715 |    0.2322 |    0.7665 |    -1.027 |    0.9072 |
#> |.....................|   -0.6286 |...........|...........|...........|
#> |    U|               |     3.022 |     4.733 |   0.06289 |     4.768 |
#> |.....................|     1.038 |...........|...........|...........|
#> |    X|               |     3.022 |     4.733 |   0.06289 |     4.768 |
#> |.....................|     1.038 |...........|...........|...........|
#> |   15|     407.36920 |    0.2369 |    0.7901 |    -1.030 |    0.9016 |
#> |.....................|   -0.6188 |...........|...........|...........|
#> |    U|               |     3.023 |     4.827 |   0.05119 |     4.754 |
#> |.....................|     1.048 |...........|...........|...........|
#> |    X|               |     3.023 |     4.827 |   0.05119 |     4.754 |
#> |.....................|     1.048 |...........|...........|...........|
#> |   16|     406.18235 |    0.2537 |    0.8728 |    -1.038 |    0.8819 |
#> |.....................|   -0.5844 |...........|...........|...........|
#> |    U|               |     3.029 |     5.158 |   0.01012 |     4.705 |
#> |.....................|     1.082 |...........|...........|...........|
#> |    X|               |     3.029 |     5.158 |   0.01012 |     4.705 |
#> |.....................|     1.082 |...........|...........|...........|
#> |   17|     404.60778 |    0.2930 |     1.067 |    -1.057 |    0.8356 |
#> |.....................|   -0.5036 |...........|...........|...........|
#> |    U|               |     3.042 |     5.935 |  -0.08627 |     4.589 |
#> |.....................|     1.163 |...........|...........|...........|
#> |    X|               |     3.042 |     5.935 |  -0.08627 |     4.589 |
#> |.....................|     1.163 |...........|...........|...........|
#> |    F|    Forward    |     3.300 |     3.302 |     60.68 |     2.150 |
#> |.....................|    -1.658 |...........|...........|...........|
#> |   18|     405.42168 | -0.009269 |     1.155 |    -1.058 |    0.6885 |
#> |.....................|   -0.3922 |...........|...........|...........|
#> |    U|               |     2.941 |     6.287 |  -0.08826 |     4.221 |
#> |.....................|     1.274 |...........|...........|...........|
#> |    X|               |     2.941 |     6.287 |  -0.08826 |     4.221 |
#> |.....................|     1.274 |...........|...........|...........|
#> |   19|     404.24848 |    0.2022 |     1.093 |    -1.065 |    0.7913 |
#> |.....................|   -0.4701 |...........|...........|...........|
#> |    U|               |     3.012 |     6.039 |   -0.1261 |     4.478 |
#> |.....................|     1.197 |...........|...........|...........|
#> |    X|               |     3.012 |     6.039 |   -0.1261 |     4.478 |
#> |.....................|     1.197 |...........|...........|...........|
#> |    F|    Forward    |    -1.318 |    -5.124 |    -9.618 |    -2.862 |
#> |.....................|    -1.229 |...........|...........|...........|
#> |   20|     404.11066 |    0.1481 |     1.154 |    -1.071 |    0.8522 |
#> |.....................|   -0.4297 |...........|...........|...........|
#> |    U|               |     2.994 |     6.284 |   -0.1547 |     4.631 |
#> |.....................|     1.237 |...........|...........|...........|
#> |    X|               |     2.994 |     6.284 |   -0.1547 |     4.631 |
#> |.....................|     1.237 |...........|...........|...........|
#> |    F|    Forward    |    -1.818 |    -5.039 |    -22.26 |     4.217 |
#> |.....................|    -1.168 |...........|...........|...........|
#> 
#> |    #| Function Val. |         m |     pop.b |   cov_w_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |   21|     404.29983 |   0.08454 |     1.188 |    -1.059 |    0.8485 |
#> |.....................|   -0.3475 |...........|...........|...........|
#> |    U|               |     2.973 |     6.418 |  -0.09375 |     4.621 |
#> |.....................|     1.319 |...........|...........|...........|
#> |    X|               |     2.973 |     6.418 |  -0.09375 |     4.621 |
#> |.....................|     1.319 |...........|...........|...........|
#> |   22|     404.35719 |    0.1280 |     1.167 |    -1.059 |    0.8494 |
#> |.....................|   -0.4024 |...........|...........|...........|
#> |    U|               |     2.987 |     6.335 |  -0.09329 |     4.624 |
#> |.....................|     1.264 |...........|...........|...........|
#> |    X|               |     2.987 |     6.335 |  -0.09329 |     4.624 |
#> |.....................|     1.264 |...........|...........|...........|
#> |   23|     404.28896 |    0.1489 |     1.157 |    -1.060 |    0.8502 |
#> |.....................|   -0.4291 |...........|...........|...........|
#> |    U|               |     2.994 |     6.294 |   -0.1014 |     4.626 |
#> |.....................|     1.238 |...........|...........|...........|
#> |    X|               |     2.994 |     6.294 |   -0.1014 |     4.626 |
#> |.....................|     1.238 |...........|...........|...........|
#> |   24|     404.07136 |    0.1483 |     1.155 |    -1.068 |    0.8516 |
#> |.....................|   -0.4295 |...........|...........|...........|
#> |    U|               |     2.994 |     6.287 |   -0.1389 |     4.629 |
#> |.....................|     1.237 |...........|...........|...........|
#> |    X|               |     2.994 |     6.287 |   -0.1389 |     4.629 |
#> |.....................|     1.237 |...........|...........|...........|
#> |    F|    Forward    |   -0.5926 |    -2.381 |     1.565 |     4.182 |
#> |.....................|    -1.193 |...........|...........|...........|
#> |   25|     404.05705 |    0.1487 |     1.157 |    -1.069 |    0.8490 |
#> |.....................|   -0.4287 |...........|...........|...........|
#> |    U|               |     2.994 |     6.293 |   -0.1439 |     4.622 |
#> |.....................|     1.238 |...........|...........|...........|
#> |    X|               |     2.994 |     6.293 |   -0.1439 |     4.622 |
#> |.....................|     1.238 |...........|...........|...........|
#> |   26|     404.04550 |    0.1494 |     1.159 |    -1.071 |    0.8441 |
#> |.....................|   -0.4273 |...........|...........|...........|
#> |    U|               |     2.994 |     6.304 |   -0.1530 |     4.610 |
#> |.....................|     1.239 |...........|...........|...........|
#> |    X|               |     2.994 |     6.304 |   -0.1530 |     4.610 |
#> |.....................|     1.239 |...........|...........|...........|
#> |    F|    Forward    |    -1.369 |    -4.095 |    -15.32 |     3.360 |
#> |.....................|    -1.149 |...........|...........|...........|
#> |   27|     404.01450 |    0.1449 |     1.162 |    -1.068 |    0.8435 |
#> |.....................|   -0.4199 |...........|...........|...........|
#> |    U|               |     2.993 |     6.315 |   -0.1421 |     4.609 |
#> |.....................|     1.247 |...........|...........|...........|
#> |    X|               |     2.993 |     6.315 |   -0.1421 |     4.609 |
#> |.....................|     1.247 |...........|...........|...........|
#> |    F|    Forward    |   -0.5573 |    -2.247 |     1.254 |     3.286 |
#> |.....................|    -1.116 |...........|...........|...........|
#> |   28|     403.99403 |    0.1461 |     1.167 |    -1.071 |    0.8363 |
#> |.....................|   -0.4174 |...........|...........|...........|
#> |    U|               |     2.993 |     6.335 |   -0.1557 |     4.591 |
#> |.....................|     1.249 |...........|...........|...........|
#> |    X|               |     2.993 |     6.335 |   -0.1557 |     4.591 |
#> |.....................|     1.249 |...........|...........|...........|
#> |    F|    Forward    |    -1.264 |    -3.812 |    -14.51 |     2.492 |
#> |.....................|    -1.075 |...........|...........|...........|
#> |   29|     403.96694 |    0.1418 |     1.170 |    -1.069 |    0.8360 |
#> |.....................|   -0.4097 |...........|...........|...........|
#> |    U|               |     2.992 |     6.345 |   -0.1457 |     4.590 |
#> |.....................|     1.257 |...........|...........|...........|
#> |    X|               |     2.992 |     6.345 |   -0.1457 |     4.590 |
#> |.....................|     1.257 |...........|...........|...........|
#> |    F|    Forward    |   -0.5149 |    -2.101 |    0.8364 |     2.445 |
#> |.....................|    -1.041 |...........|...........|...........|
#> |   30|     403.94623 |    0.1432 |     1.175 |    -1.071 |    0.8294 |
#> |.....................|   -0.4069 |...........|...........|...........|
#> |    U|               |     2.992 |     6.367 |   -0.1569 |     4.574 |
#> |.....................|     1.260 |...........|...........|...........|
#> |    X|               |     2.992 |     6.367 |   -0.1569 |     4.574 |
#> |.....................|     1.260 |...........|...........|...........|
#> |    F|    Forward    |   -0.9960 |    -3.177 |    -10.80 |     1.705 |
#> |.....................|    -1.007 |...........|...........|...........|
#> 
#> |    #| Function Val. |         m |     pop.b |   cov_w_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |   31|     403.92729 |    0.1391 |     1.177 |    -1.070 |    0.8293 |
#> |.....................|   -0.3989 |...........|...........|...........|
#> |    U|               |     2.991 |     6.376 |   -0.1495 |     4.573 |
#> |.....................|     1.268 |...........|...........|...........|
#> |    X|               |     2.991 |     6.376 |   -0.1495 |     4.573 |
#> |.....................|     1.268 |...........|...........|...........|
#> |    F|    Forward    |   -0.4424 |    -1.895 |    0.7135 |     1.675 |
#> |.....................|   -0.9689 |...........|...........|...........|
#> |   32|     403.91341 |    0.1406 |     1.184 |    -1.072 |    0.8237 |
#> |.....................|   -0.3957 |...........|...........|...........|
#> |    U|               |     2.991 |     6.402 |   -0.1613 |     4.559 |
#> |.....................|     1.271 |...........|...........|...........|
#> |    X|               |     2.991 |     6.402 |   -0.1613 |     4.559 |
#> |.....................|     1.271 |...........|...........|...........|
#> |    F|    Forward    |   -0.9368 |    -3.000 |    -11.46 |     1.030 |
#> |.....................|   -0.9366 |...........|...........|...........|
#> |   33|     403.89482 |    0.1367 |     1.186 |    -1.071 |    0.8238 |
#> |.....................|   -0.3874 |...........|...........|...........|
#> |    U|               |     2.990 |     6.410 |   -0.1537 |     4.559 |
#> |.....................|     1.279 |...........|...........|...........|
#> |    X|               |     2.990 |     6.410 |   -0.1537 |     4.559 |
#> |.....................|     1.279 |...........|...........|...........|
#> |    F|    Forward    |   -0.3704 |    -1.694 |    0.3090 |     1.025 |
#> |.....................|   -0.9006 |...........|...........|...........|
#> |   34|     403.87689 |    0.1382 |     1.193 |    -1.072 |    0.8194 |
#> |.....................|   -0.3836 |...........|...........|...........|
#> |    U|               |     2.991 |     6.438 |   -0.1602 |     4.549 |
#> |.....................|     1.283 |...........|...........|...........|
#> |    X|               |     2.991 |     6.438 |   -0.1602 |     4.549 |
#> |.....................|     1.283 |...........|...........|...........|
#> |   35|     403.85851 |    0.1420 |     1.210 |    -1.075 |    0.8090 |
#> |.....................|   -0.3744 |...........|...........|...........|
#> |    U|               |     2.992 |     6.507 |   -0.1760 |     4.522 |
#> |.....................|     1.292 |...........|...........|...........|
#> |    X|               |     2.992 |     6.507 |   -0.1760 |     4.522 |
#> |.....................|     1.292 |...........|...........|...........|
#> |    F|    Forward    |   -0.4607 |    -1.995 |    -11.16 |   -0.7404 |
#> |.....................|   -0.8143 |...........|...........|...........|
#> |   36|     403.82344 |    0.1289 |     1.215 |    -1.074 |    0.8099 |
#> |.....................|   -0.3455 |...........|...........|...........|
#> |    U|               |     2.987 |     6.527 |   -0.1692 |     4.525 |
#> |.....................|     1.321 |...........|...........|...........|
#> |    X|               |     2.987 |     6.527 |   -0.1692 |     4.525 |
#> |.....................|     1.321 |...........|...........|...........|
#> |   37|     403.79136 |    0.1064 |     1.223 |    -1.074 |    0.8113 |
#> |.....................|   -0.2962 |...........|...........|...........|
#> |    U|               |     2.980 |     6.558 |   -0.1707 |     4.528 |
#> |.....................|     1.370 |...........|...........|...........|
#> |    X|               |     2.980 |     6.558 |   -0.1707 |     4.528 |
#> |.....................|     1.370 |...........|...........|...........|
#> |   38|     403.76911 |   0.06085 |     1.239 |    -1.075 |    0.8142 |
#> |.....................|   -0.1965 |...........|...........|...........|
#> |    U|               |     2.965 |     6.621 |   -0.1738 |     4.535 |
#> |.....................|     1.470 |...........|...........|...........|
#> |    X|               |     2.965 |     6.621 |   -0.1738 |     4.535 |
#> |.....................|     1.470 |...........|...........|...........|
#> |    F|    Forward    |    -1.618 |    -2.933 |    -15.90 |   -0.5380 |
#> |.....................|   -0.3245 |...........|...........|...........|
#> |   39|     403.77582 |    0.2421 |     1.183 |    -1.076 |    0.8332 |
#> |.....................|    0.1496 |...........|...........|...........|
#> |    U|               |     3.025 |     6.399 |   -0.1780 |     4.583 |
#> |.....................|     1.816 |...........|...........|...........|
#> |    X|               |     3.025 |     6.399 |   -0.1780 |     4.583 |
#> |.....................|     1.816 |...........|...........|...........|
#> |   40|     403.71407 |    0.1490 |     1.212 |    -1.075 |    0.8234 |
#> |.....................|  -0.02816 |...........|...........|...........|
#> |    U|               |     2.994 |     6.513 |   -0.1758 |     4.559 |
#> |.....................|     1.639 |...........|...........|...........|
#> |    X|               |     2.994 |     6.513 |   -0.1758 |     4.559 |
#> |.....................|     1.639 |...........|...........|...........|
#> |    F|    Forward    |   -0.1627 |    -1.462 |    -6.830 |    0.3239 |
#> |.....................|   -0.1651 |...........|...........|...........|
#> 
#> |    #| Function Val. |         m |     pop.b |   cov_w_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |   41|     403.70526 |   0.07172 |     1.293 |    -1.078 |    0.8268 |
#> |.....................|   0.05099 |...........|...........|...........|
#> |    U|               |     2.968 |     6.840 |   -0.1905 |     4.567 |
#> |.....................|     1.718 |...........|...........|...........|
#> |    X|               |     2.968 |     6.840 |   -0.1905 |     4.567 |
#> |.....................|     1.718 |...........|...........|...........|
#> |    F|    Forward    |    0.9706 |     2.512 |     12.52 |    0.6267 |
#> |.....................|   -0.1249 |...........|...........|...........|
#> |   42|     403.67359 |   0.09732 |     1.257 |    -1.077 |    0.8225 |
#> |.....................|   0.09591 |...........|...........|...........|
#> |    U|               |     2.977 |     6.694 |   -0.1827 |     4.556 |
#> |.....................|     1.763 |...........|...........|...........|
#> |    X|               |     2.977 |     6.694 |   -0.1827 |     4.556 |
#> |.....................|     1.763 |...........|...........|...........|
#> |    F|    Forward    |    0.2130 |    0.3607 |     2.061 |    0.1049 |
#> |.....................|   -0.1023 |...........|...........|...........|
#> |   43|     403.66602 |   0.09850 |     1.252 |    -1.076 |    0.8222 |
#> |.....................|    0.1799 |...........|...........|...........|
#> |    U|               |     2.977 |     6.674 |   -0.1816 |     4.556 |
#> |.....................|     1.847 |...........|...........|...........|
#> |    X|               |     2.977 |     6.674 |   -0.1816 |     4.556 |
#> |.....................|     1.847 |...........|...........|...........|
#> |    F|    Forward    |   0.01760 |  -0.08917 |   -0.3348 |  0.004898 |
#> |.....................|  -0.07504 |...........|...........|...........|
#> |   44|     403.65548 |   0.09510 |     1.252 |    -1.076 |    0.8231 |
#> |.....................|    0.3721 |...........|...........|...........|
#> |    U|               |     2.976 |     6.675 |   -0.1816 |     4.558 |
#> |.....................|     2.039 |...........|...........|...........|
#> |    X|               |     2.976 |     6.675 |   -0.1816 |     4.558 |
#> |.....................|     2.039 |...........|...........|...........|
#> |    F|    Forward    |   -0.1049 |   -0.2902 |    -1.512 | -0.009671 |
#> |.....................|  -0.04034 |...........|...........|...........|
#> |   45|     403.64910 |   0.09284 |     1.255 |    -1.076 |    0.8239 |
#> |.....................|    0.5642 |...........|...........|...........|
#> |    U|               |     2.975 |     6.686 |   -0.1822 |     4.560 |
#> |.....................|     2.231 |...........|...........|...........|
#> |    X|               |     2.975 |     6.686 |   -0.1822 |     4.560 |
#> |.....................|     2.231 |...........|...........|...........|
#> |    F|    Forward    |  -0.05584 |   -0.1390 |   -0.7323 |   0.01420 |
#> |.....................|  -0.02362 |...........|...........|...........|
#> |   46|     403.64432 |   0.09169 |     1.257 |    -1.077 |    0.8247 |
#> |.....................|    0.8319 |...........|...........|...........|
#> |    U|               |     2.975 |     6.696 |   -0.1829 |     4.562 |
#> |.....................|     2.499 |...........|...........|...........|
#> |    X|               |     2.975 |     6.696 |   -0.1829 |     4.562 |
#> |.....................|     2.499 |...........|...........|...........|
#> |    F|    Forward    |  0.009493 |   0.02628 |    0.1049 |   0.04159 |
#> |.....................|  -0.01240 |...........|...........|...........|
#> |   47|     403.64432 |   0.09169 |     1.257 |    -1.077 |    0.8247 |
#> |.....................|    0.8319 |...........|...........|...........|
#> |    U|               |     2.975 |     6.696 |   -0.1829 |     4.562 |
#> |.....................|     2.499 |...........|...........|...........|
#> |    X|               |     2.975 |     6.696 |   -0.1829 |     4.562 |
#> |.....................|     2.499 |...........|...........|...........|
#> calculating covariance matrix
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> 
#> covType="analytic" not available for this model (out of scope, or the augmented model would not build/solve); using the finite-difference sandwich ("r,s") covariance.
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> done
#> → Calculating residuals/tables
#> ✔ done

mod
```

``` math
\begin{align*}
{b} & = {pop.b}+{cov\_w\_b} {\times} {w} \\
{value} & = {m} {\times} {x}+{b}+{bRe} \\
{value} & \sim add({addSd})
\end{align*}
```

### Quick start, log-linked parameter

Some parameters are naturally strictly positive (clearances, volumes,
rate constants). Pass `paramLink = c(<parameter> = "log")` to wrap the
linear combination in [`exp()`](https://rdrr.io/r/base/Log.html). The
starting values are then on the log scale.

``` r

nlmixr(
  y~m*x + b + bRe ~ bRe|id,
  start = list(m=3, b=c(log(4), 0.2), addSd=5),
  param = list(b ~ w),
  paramLink = c(b = "log"),
  data = dSimCont,
  est = "focei"
)
```

### Worked example: concentration-QT analysis

A common use of the formula interface is a thorough QT (TQT) or
concentration-QT (C-QT) analysis. The clinical question is whether a 10
msec $`\Delta QT`$ is exceeded at therapeutic drug concentrations. The
standard model has:

- A categorical effect for nominal post-dose time (captures circadian
  and food-related drift)
- A continuous linear effect for drug concentration (the parameter of
  interest)
- A per-subject random intercept (baseline $`\Delta QT`$ offset)

In the formula interface this is one call:

    dQT ~ b + bRe ~ (bRe|id)
    param  = list(b ~ timeF + conc)

`b` is a single intercept parameter; `param` decomposes it into a fixed
effect per `timeF` level *plus* a linear slope on `conc`. The random
per- subject offset `bRe` is added to `b` in the predictor.

#### Simulating the data

We start from a one-compartment oral PK profile (Bateman equation),
scale it so the typical $`C_{max}`$ is about 1000 ng/mL, and apply
per-subject log-normal scaling so individual $`C_{max}`$ values are
log-normally distributed.

``` r

withr::with_seed(42, {
  nSubj    <- 40
  nomTimes <- c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12)  # 10 nominal times

  # Typical 1-compartment oral PK profile (unit Bateman), then scaled
  # so the population Cmax is 1000 ng/mL
  bateman    <- function(t, ka = 2, ke = 0.1) {
    (ka / (ka - ke)) * (exp(-ke * t) - exp(-ka * t))
  }
  cmaxFactor <- 1000 / max(bateman(seq(0, 12, 0.01)))

  # Per-subject log-normal scaling of the concentration profile
  subjScale  <- exp(rnorm(nSubj, mean = 0, sd = 0.4))

  cqtDesign <-
    expand.grid(id = 1:nSubj, timeH = nomTimes)
  cqtDesign$conc  <- bateman(cqtDesign$timeH) * cmaxFactor *
                     subjScale[cqtDesign$id]
  cqtDesign$timeF <- factor(cqtDesign$timeH,
                            levels = as.character(nomTimes))
  cqtDesign$dQT   <- 1  # placeholder; nlmixrFormula(est="rxSolve") fills it in
})
```

Now we simulate `dQT` from the true model. The true slope is 0.01
msec/(ng/mL), which produces about a 10 msec $`\Delta QT`$ at the
population $`C_{max}`$. The categorical time effects encode a typical
circadian drift pattern. The simulation’s `bRe` SD is 1 (the formula
interface’s current default for random-effect starting values).

``` r

cqtSimPrep <-
  nlmixrFormula(
    dQT ~ b + bRe ~ (bRe|id),
    data  = cqtDesign,
    # 10 time-level intercepts (relative to time 0) + slope on conc
    start = list(
      b     = c(0,         # timeF == 0  (pre-dose baseline)
                3, 5, 6, 5, 2, 0, -1, -2, -3,  # post-dose time effects
                0.01),     # slope on conc (msec per ng/mL)
      addSd = 3
    ),
    param = list(b ~ timeF + conc),
    est   = "rxSolve"
  )
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments

cqtData <-
  cqtSimPrep |>
  as.data.frame() |>
  select(id, conc, timeF, dQT = sim)

ggplot(cqtData, aes(x = conc, y = dQT, colour = timeF)) +
  geom_point(alpha = 0.7) +
  labs(x = "Concentration (ng/mL)", y = expression(Delta * "QT (msec)"))
```

![](nlmixrFormula_files/figure-html/cqt-sim-1.png)

#### Fitting the model

Algebraic models with both a factor and a continuous covariate on the
same parameter can be hard for `focei`’s gradient calculation when there
are many factor levels. A robust two-stage strategy is to fit the fixed
effects with `bobyqa` (which is gradient-free and converges fast for
this class of model) and then refine with `focei` using the bobyqa
estimates as starts, which adds the random effect.

``` r

fitFE <-
  nlmixr(
    dQT ~ b,                                       # fixed effects only
    data  = cqtData,
    start = list(
      b     = c(0, 3, 5, 6, 5, 2, 0, -1, -2, -3, 0.01),
      addSd = 3
    ),
    param = list(b ~ timeF + conc),
    est   = "bobyqa"
  )
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → pruning branches (`if`/`else`) of population log-likelihood model...
#> ✔ done
#> → loading llik model into symengine environment...
#> → finding duplicate expressions in population log-likelihood model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in population log-likelihood model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> ✔ done
#> Key: U: Unscaled Parameters; X: Back-transformed parameters; 
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |    1|     1317.9790 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    2|     1317.9790 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    3| 8.8890445e+09 |   -0.1333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               | 2.000e+04 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               | 2.000e+04 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    4|     1336.3469 |   -0.3333 |    0.5333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.600 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.600 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    5|     1316.2619 |   -0.3333 |    0.3333 |    0.9778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     6.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     6.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    6|     1353.0714 |   -0.3333 |    0.3333 |    0.7778 |     1.200 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     7.200 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     7.200 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    7|     1319.5776 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.9778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     6.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     6.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    8|     1317.4572 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.3111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.400 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.400 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    9| 8.8894967e+08 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.1333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 | 2.000e+04 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 | 2.000e+04 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |   10|     1318.4920 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.3556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |   -0.8000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |   -0.8000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   11|     1306.2237 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.5778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -1.600 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -1.600 |    -3.000 |   0.01000 |     3.000 |
#> |   12|     1317.2181 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |   -0.8000 |   -0.3311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -2.400 |   0.01000 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -2.400 |   0.01000 |     3.000 |
#> |   13| 7.5840796e+09 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.1311 |    0.3333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |     20.01 |     3.000 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |     20.01 |     3.000 |
#> |   14|     1267.4243 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.300 |
#> |    X|               |     0.000 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.300 |
#> |   15| 8.8887360e+09 |   -0.5333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.3333 |
#> |    U|               |-2.000e+04 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |    X|               |-2.000e+04 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     0.000 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01000 |     3.000 |
#> |   16| 6.3630137e+09 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3566 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.5298 |    0.5333 |
#> |    U|               |   -0.4040 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -2328. |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |    -19.85 |     3.300 |
#> |    X|               |   -0.4040 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -2328. |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |    -19.85 |     3.300 |
#> |   17| 1.3751194e+08 |   -0.3333 |    0.3635 |    0.8080 |     1.030 |
#> |.....................|    0.8080 |    0.1413 |   -0.3064 |   -0.5254 |
#> |.....................|   -0.7476 |   -0.9698 |   -0.3054 |    0.5696 |
#> |    U|               |  -0.04371 |     3.091 |     5.151 |     6.181 |
#> |.....................|     5.151 |     2.060 |     2698. |   -0.9698 |
#> |.....................|    -1.940 |    -2.909 |     2.577 |     3.354 |
#> |    X|               |  -0.04371 |     3.091 |     5.151 |     6.181 |
#> |.....................|     5.151 |     2.060 |     2698. |   -0.9698 |
#> |.....................|    -1.940 |    -2.909 |     2.577 |     3.354 |
#> |   18| 1.4690412e+08 |   -0.3616 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |    -2828. |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 | -0.002960 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01004 |     3.300 |
#> |    X|               |    -2828. |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 | -0.002960 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01004 |     3.300 |
#> |   19| 1.6038203e+08 |   -0.3317 |    0.3230 |    0.7675 |    0.9897 |
#> |.....................|    0.7675 |    0.1008 |   -0.3431 |   -0.5659 |
#> |.....................|   -0.7881 |    -1.010 |   -0.2992 |    0.5564 |
#> |    U|               |     159.9 |     2.969 |     4.948 |     5.938 |
#> |.....................|     4.948 |     1.979 |    -981.4 |    -1.010 |
#> |.....................|    -2.021 |    -3.031 |     3.198 |     3.335 |
#> |    X|               |     159.9 |     2.969 |     4.948 |     5.938 |
#> |.....................|     4.948 |     1.979 |    -981.4 |    -1.010 |
#> |.....................|    -2.021 |    -3.031 |     3.198 |     3.335 |
#> |   20|     181671.32 |   -0.3333 |    0.3521 |    0.7755 |    0.9977 |
#> |.....................|    0.7755 |    0.1088 |   -0.3358 |   -0.5578 |
#> |.....................|   -0.7801 |    -1.002 |   -0.3301 |    0.5318 |
#> |    U|               |   0.05481 |     3.056 |     4.989 |     5.986 |
#> |.....................|     4.989 |     1.995 |    -251.5 |    -1.002 |
#> |.....................|    -2.005 |    -3.007 |    0.1101 |     3.298 |
#> |    X|               |   0.05481 |     3.056 |     4.989 |     5.986 |
#> |.....................|     4.989 |     1.995 |    -251.5 |    -1.002 |
#> |.....................|    -2.005 |    -3.007 |    0.1101 |     3.298 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   21|     1755535.5 |   -0.3326 |    0.3327 |    0.7785 |     1.001 |
#> |.....................|    0.7785 |    0.1118 |   -0.3262 |   -0.5548 |
#> |.....................|   -0.7771 |   -0.9993 |   -0.3303 |    0.5148 |
#> |    U|               |     78.12 |     2.998 |     5.004 |     6.004 |
#> |.....................|     5.004 |     2.001 |     709.6 |   -0.9993 |
#> |.....................|    -1.999 |    -2.998 |   0.08955 |     3.272 |
#> |    X|               |     78.12 |     2.998 |     5.004 |     6.004 |
#> |.....................|     5.004 |     2.001 |     709.6 |   -0.9993 |
#> |.....................|    -1.999 |    -2.998 |   0.08955 |     3.272 |
#> |   22|     220676.66 |   -0.3333 |    0.3332 |    0.7965 |    0.9976 |
#> |.....................|    0.7753 |    0.1087 |   -0.3369 |   -0.5580 |
#> |.....................|   -0.7802 |    -1.002 |   -0.3304 |    0.5323 |
#> |    U|               |  -0.09543 |     3.000 |     5.094 |     5.985 |
#> |.....................|     4.988 |     1.995 |    -357.8 |    -1.002 |
#> |.....................|    -2.005 |    -3.007 |   0.08022 |     3.298 |
#> |    X|               |  -0.09543 |     3.000 |     5.094 |     5.985 |
#> |.....................|     4.988 |     1.995 |    -357.8 |    -1.002 |
#> |.....................|    -2.005 |    -3.007 |   0.08022 |     3.298 |
#> |   23|     11617273. |   -0.3322 |    0.3346 |    0.7783 |     1.003 |
#> |.....................|    0.7812 |    0.1145 |   -0.3437 |   -0.5521 |
#> |.....................|   -0.7743 |   -0.9966 |   -0.3388 |    0.5207 |
#> |    U|               |     112.4 |     3.004 |     5.003 |     6.021 |
#> |.....................|     5.017 |     2.007 |    -1033. |   -0.9966 |
#> |.....................|    -1.993 |    -2.990 |   -0.7580 |     3.281 |
#> |    X|               |     112.4 |     3.004 |     5.003 |     6.021 |
#> |.....................|     5.017 |     2.007 |    -1033. |   -0.9966 |
#> |.....................|    -1.993 |    -2.990 |   -0.7580 |     3.281 |
#> |   24|     30955.095 |   -0.3333 |    0.3334 |    0.7779 |     1.018 |
#> |.....................|    0.7741 |    0.1075 |   -0.3331 |   -0.5592 |
#> |.....................|   -0.7814 |    -1.004 |   -0.3316 |    0.5337 |
#> |    U|               |   -0.4735 |     3.000 |     5.000 |     6.110 |
#> |.....................|     4.982 |     1.993 |     26.70 |    -1.004 |
#> |.....................|    -2.007 |    -3.011 |  -0.03536 |     3.301 |
#> |    X|               |   -0.4735 |     3.000 |     5.000 |     6.110 |
#> |.....................|     4.982 |     1.993 |     26.70 |    -1.004 |
#> |.....................|    -2.007 |    -3.011 |  -0.03536 |     3.301 |
#> |   25|     8957240.0 |   -0.3325 |    0.3360 |    0.7802 |    0.9969 |
#> |.....................|    0.7759 |    0.1092 |   -0.3161 |   -0.5575 |
#> |.....................|   -0.7797 |    -1.002 |   -0.3390 |    0.5332 |
#> |    U|               |     87.86 |     3.008 |     5.012 |     5.981 |
#> |.....................|     4.990 |     1.996 |     1718. |    -1.002 |
#> |.....................|    -2.004 |    -3.006 |   -0.7800 |     3.300 |
#> |    X|               |     87.86 |     3.008 |     5.012 |     5.981 |
#> |.....................|     4.990 |     1.996 |     1718. |    -1.002 |
#> |.....................|    -2.004 |    -3.006 |   -0.7800 |     3.300 |
#> |   26|     1305.0784 |   -0.3333 |    0.3334 |    0.7778 |     1.000 |
#> |.....................|    0.7957 |    0.1066 |   -0.3333 |   -0.5600 |
#> |.....................|   -0.7822 |    -1.004 |   -0.3311 |    0.5336 |
#> |    U|               |   -0.3104 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.089 |     1.991 |   -0.2462 |    -1.004 |
#> |.....................|    -2.009 |    -3.013 |   0.01186 |     3.300 |
#> |    X|               |   -0.3104 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.089 |     1.991 |   -0.2462 |    -1.004 |
#> |.....................|    -2.009 |    -3.013 |   0.01186 |     3.300 |
#> |   27|     10055567. |   -0.3323 |    0.3308 |    0.7743 |    0.9948 |
#> |.....................|    0.7726 |    0.1068 |   -0.3285 |   -0.5598 |
#> |.....................|   -0.7821 |    -1.004 |   -0.3247 |    0.5470 |
#> |    U|               |     99.91 |     2.992 |     4.983 |     5.969 |
#> |.....................|     4.974 |     1.991 |     484.6 |    -1.004 |
#> |.....................|    -2.009 |    -3.013 |    0.6509 |     3.321 |
#> |    X|               |     99.91 |     2.992 |     4.983 |     5.969 |
#> |.....................|     4.974 |     1.991 |     484.6 |    -1.004 |
#> |.....................|    -2.009 |    -3.013 |    0.6509 |     3.321 |
#> |   28|     72390.052 |   -0.3333 |    0.3340 |    0.7784 |     1.001 |
#> |.....................|    0.7784 |    0.1290 |   -0.3328 |   -0.5606 |
#> |.....................|   -0.7829 |    -1.005 |   -0.3305 |    0.5340 |
#> |    U|               |   0.08308 |     3.002 |     5.003 |     6.004 |
#> |.....................|     5.003 |     2.036 |     55.24 |    -1.005 |
#> |.....................|    -2.010 |    -3.015 |   0.06836 |     3.301 |
#> |    X|               |   0.08308 |     3.002 |     5.003 |     6.004 |
#> |.....................|     5.003 |     2.036 |     55.24 |    -1.005 |
#> |.....................|    -2.010 |    -3.015 |   0.06836 |     3.301 |
#> |   29|     7823083.4 |   -0.3326 |    0.3393 |    0.7836 |     1.006 |
#> |.....................|    0.7837 |    0.1152 |   -0.3269 |   -0.5487 |
#> |.....................|   -0.7709 |   -0.9932 |   -0.3258 |    0.5391 |
#> |    U|               |     76.75 |     3.018 |     5.029 |     6.036 |
#> |.....................|     5.030 |     2.008 |     641.8 |   -0.9932 |
#> |.....................|    -1.986 |    -2.979 |    0.5407 |     3.309 |
#> |    X|               |     76.75 |     3.018 |     5.029 |     6.036 |
#> |.....................|     5.030 |     2.008 |     641.8 |   -0.9932 |
#> |.....................|    -1.986 |    -2.979 |    0.5407 |     3.309 |
#> |   30|     1315.7509 |   -0.3333 |    0.3333 |    0.7778 |    0.9999 |
#> |.....................|    0.7777 |    0.1110 |   -0.3334 |   -0.5392 |
#> |.....................|   -0.7860 |    -1.008 |   -0.3311 |    0.5337 |
#> |    U|               |   -0.4925 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -3.703 |   -0.9837 |
#> |.....................|    -2.016 |    -3.025 |  0.009355 |     3.301 |
#> |    X|               |   -0.4925 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -3.703 |   -0.9837 |
#> |.....................|    -2.016 |    -3.025 |  0.009355 |     3.301 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   31|     10780701. |   -0.3323 |    0.3311 |    0.7746 |    0.9951 |
#> |.....................|    0.7729 |    0.1072 |   -0.3281 |   -0.5602 |
#> |.....................|   -0.7813 |    -1.004 |   -0.3245 |    0.5476 |
#> |    U|               |     101.0 |     2.993 |     4.984 |     5.971 |
#> |.....................|     4.976 |     1.992 |     521.2 |    -1.005 |
#> |.....................|    -2.007 |    -3.011 |    0.6726 |     3.321 |
#> |    X|               |     101.0 |     2.993 |     4.984 |     5.971 |
#> |.....................|     4.976 |     1.992 |     521.2 |    -1.005 |
#> |.....................|    -2.007 |    -3.011 |    0.6726 |     3.321 |
#> |   32|     1268.0456 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7919 |   -0.9859 |   -0.3311 |    0.5333 |
#> |    U|               | 1.312e-07 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |-8.312e-06 |    -1.000 |
#> |.....................|    -2.028 |    -2.958 |   0.01000 |     3.300 |
#> |    X|               | 1.312e-07 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |-8.312e-06 |    -1.000 |
#> |.....................|    -2.028 |    -2.958 |   0.01000 |     3.300 |
#> |   33|     2254474.6 |   -0.3324 |    0.3313 |    0.7763 |    0.9976 |
#> |.....................|    0.7755 |    0.1090 |   -0.3245 |   -0.5573 |
#> |.....................|   -0.7770 |    -1.003 |   -0.3353 |    0.5168 |
#> |    U|               |     90.04 |     2.994 |     4.993 |     5.986 |
#> |.....................|     4.989 |     1.996 |     879.0 |    -1.002 |
#> |.....................|    -1.998 |    -3.008 |   -0.4125 |     3.275 |
#> |    X|               |     90.04 |     2.994 |     4.993 |     5.986 |
#> |.....................|     4.989 |     1.996 |     879.0 |    -1.002 |
#> |.....................|    -1.998 |    -3.008 |   -0.4125 |     3.275 |
#> |   34|     791757.31 |   -0.3325 |    0.3370 |    0.7786 |    0.9976 |
#> |.....................|    0.7773 |    0.1112 |   -0.3382 |   -0.5555 |
#> |.....................|   -0.7733 |    -1.000 |   -0.3297 |    0.5392 |
#> |    U|               |     82.29 |     3.011 |     5.004 |     5.985 |
#> |.....................|     4.997 |     2.000 |    -483.6 |    -1.000 |
#> |.....................|    -1.991 |    -3.001 |    0.1508 |     3.309 |
#> |    X|               |     82.29 |     3.011 |     5.004 |     5.985 |
#> |.....................|     4.997 |     2.000 |    -483.6 |    -1.000 |
#> |.....................|    -1.991 |    -3.001 |    0.1508 |     3.309 |
#> |   35|     1803276.4 |   -0.3361 |    0.3334 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1112 |   -0.3334 |   -0.5555 |
#> |.....................|   -0.7777 |   -0.9999 |   -0.3316 |    0.5337 |
#> |    U|               |    -276.0 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -6.923 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  -0.03664 |     3.301 |
#> |    X|               |    -276.0 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -6.923 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  -0.03664 |     3.301 |
#> |   36|     158743.45 |   -0.3330 |    0.3334 |    0.7776 |    0.9990 |
#> |.....................|    0.7771 |    0.1105 |   -0.3313 |   -0.5560 |
#> |.....................|   -0.7768 |    -1.001 |   -0.3324 |    0.5294 |
#> |    U|               |     36.03 |     3.000 |     4.999 |     5.994 |
#> |.....................|     4.996 |     1.999 |     207.2 |    -1.000 |
#> |.....................|    -1.998 |    -3.002 |   -0.1167 |     3.294 |
#> |    X|               |     36.03 |     3.000 |     4.999 |     5.994 |
#> |.....................|     4.996 |     1.999 |     207.2 |    -1.000 |
#> |.....................|    -1.998 |    -3.002 |   -0.1167 |     3.294 |
#> |   37|     12798.312 |   -0.3334 |    0.3333 |    0.7796 |    0.9998 |
#> |.....................|    0.7776 |    0.1109 |   -0.3336 |   -0.5558 |
#> |.....................|   -0.7781 |    -1.000 |   -0.3308 |    0.5331 |
#> |    U|               |    -8.054 |     3.000 |     5.009 |     5.999 |
#> |.....................|     4.999 |     2.000 |    -26.28 |    -1.000 |
#> |.....................|    -2.001 |    -3.001 |   0.04567 |     3.300 |
#> |    X|               |    -8.054 |     3.000 |     5.009 |     5.999 |
#> |.....................|     4.999 |     2.000 |    -26.28 |    -1.000 |
#> |.....................|    -2.001 |    -3.001 |   0.04567 |     3.300 |
#> |   38|     27669.441 |   -0.3329 |    0.3328 |    0.7785 |     1.001 |
#> |.....................|    0.7780 |    0.1111 |   -0.3332 |   -0.5556 |
#> |.....................|   -0.7793 |   -0.9998 |   -0.3313 |    0.5336 |
#> |    U|               |     47.46 |     2.998 |     5.004 |     6.004 |
#> |.....................|     5.001 |     2.000 |     15.66 |    -1.000 |
#> |.....................|    -2.003 |    -2.999 | -0.006133 |     3.300 |
#> |    X|               |     47.46 |     2.998 |     5.004 |     6.004 |
#> |.....................|     5.001 |     2.000 |     15.66 |    -1.000 |
#> |.....................|    -2.003 |    -2.999 | -0.006133 |     3.300 |
#> |   39|     2359.7142 |   -0.3333 |    0.3351 |    0.7776 |     1.000 |
#> |.....................|    0.7777 |    0.1109 |   -0.3333 |   -0.5558 |
#> |.....................|   -0.7786 |    -1.000 |   -0.3310 |    0.5332 |
#> |    U|               |    -1.378 |     3.005 |     4.999 |     6.001 |
#> |.....................|     4.999 |     2.000 |     2.641 |    -1.000 |
#> |.....................|    -2.002 |    -3.000 |   0.01918 |     3.300 |
#> |    X|               |    -1.378 |     3.005 |     4.999 |     6.001 |
#> |.....................|     4.999 |     2.000 |     2.641 |    -1.000 |
#> |.....................|    -2.002 |    -3.000 |   0.01918 |     3.300 |
#> |   40|     135902.86 |   -0.3331 |    0.3326 |    0.7768 |     1.001 |
#> |.....................|    0.7780 |    0.1112 |   -0.3344 |   -0.5554 |
#> |.....................|   -0.7783 |   -0.9998 |   -0.3304 |    0.5330 |
#> |    U|               |     24.99 |     2.998 |     4.995 |     6.004 |
#> |.....................|     5.001 |     2.000 |    -102.2 |   -0.9999 |
#> |.....................|    -2.001 |    -2.999 |   0.08517 |     3.300 |
#> |    X|               |     24.99 |     2.998 |     4.995 |     6.004 |
#> |.....................|     5.001 |     2.000 |    -102.2 |   -0.9999 |
#> |.....................|    -2.001 |    -2.999 |   0.08517 |     3.300 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   41|     3990.7924 |   -0.3333 |    0.3334 |    0.7778 |     1.001 |
#> |.....................|    0.7775 |    0.1108 |   -0.3333 |   -0.5558 |
#> |.....................|   -0.7772 |    -1.001 |   -0.3313 |    0.5335 |
#> |    U|               |     3.713 |     3.000 |     5.000 |     6.009 |
#> |.....................|     4.999 |     1.999 |    0.4644 |    -1.000 |
#> |.....................|    -1.999 |    -3.003 | -0.006361 |     3.300 |
#> |    X|               |     3.713 |     3.000 |     5.000 |     6.009 |
#> |.....................|     4.999 |     1.999 |    0.4644 |    -1.000 |
#> |.....................|    -1.999 |    -3.003 | -0.006361 |     3.300 |
#> |   42|     107694.11 |   -0.3329 |    0.3336 |    0.7783 |    0.9995 |
#> |.....................|    0.7785 |    0.1117 |   -0.3328 |   -0.5551 |
#> |.....................|   -0.7774 |    -1.001 |   -0.3308 |    0.5339 |
#> |    U|               |     45.48 |     3.001 |     5.002 |     5.997 |
#> |.....................|     5.004 |     2.001 |     53.90 |   -0.9995 |
#> |.....................|    -1.999 |    -3.003 |   0.03796 |     3.301 |
#> |    X|               |     45.48 |     3.001 |     5.002 |     5.997 |
#> |.....................|     5.004 |     2.001 |     53.90 |   -0.9995 |
#> |.....................|    -1.999 |    -3.003 |   0.03796 |     3.301 |
#> |   43|     25002.716 |   -0.3333 |    0.3333 |    0.7778 |    0.9999 |
#> |.....................|    0.7792 |    0.1102 |   -0.3337 |   -0.5564 |
#> |.....................|   -0.7777 |    -1.000 |   -0.3315 |    0.5333 |
#> |    U|               |     5.239 |     3.000 |     5.000 |     5.999 |
#> |.....................|     5.007 |     1.998 |    -32.08 |    -1.001 |
#> |.....................|    -2.000 |    -3.000 |  -0.02847 |     3.300 |
#> |    X|               |     5.239 |     3.000 |     5.000 |     5.999 |
#> |.....................|     5.007 |     1.998 |    -32.08 |    -1.001 |
#> |.....................|    -2.000 |    -3.000 |  -0.02847 |     3.300 |
#> |   44|     31615.575 |   -0.3329 |    0.3335 |    0.7781 |     1.000 |
#> |.....................|    0.7767 |    0.1105 |   -0.3335 |   -0.5561 |
#> |.....................|   -0.7771 |   -0.9989 |   -0.3311 |    0.5335 |
#> |    U|               |     39.82 |     3.001 |     5.002 |     6.001 |
#> |.....................|     4.995 |     1.999 |    -15.83 |    -1.001 |
#> |.....................|    -1.999 |    -2.997 |   0.01260 |     3.300 |
#> |    X|               |     39.82 |     3.001 |     5.002 |     6.001 |
#> |.....................|     4.995 |     1.999 |    -15.83 |    -1.001 |
#> |.....................|    -1.999 |    -2.997 |   0.01260 |     3.300 |
#> |   45|     1436.0098 |   -0.3334 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1125 |   -0.3333 |   -0.5570 |
#> |.....................|   -0.7777 |   -0.9999 |   -0.3311 |    0.5332 |
#> |    U|               |    -3.143 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.003 |     2.176 |    -1.001 |
#> |.....................|    -2.000 |    -3.000 |   0.01532 |     3.300 |
#> |    X|               |    -3.143 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.003 |     2.176 |    -1.001 |
#> |.....................|    -2.000 |    -3.000 |   0.01532 |     3.300 |
#> |   46|     68499.908 |   -0.3330 |    0.3333 |    0.7778 |    0.9993 |
#> |.....................|    0.7770 |    0.1113 |   -0.3346 |   -0.5555 |
#> |.....................|   -0.7781 |    -1.001 |   -0.3317 |    0.5329 |
#> |    U|               |     32.41 |     3.000 |     5.000 |     5.996 |
#> |.....................|     4.996 |     2.000 |    -127.0 |    -1.000 |
#> |.....................|    -2.001 |    -3.002 |  -0.05276 |     3.299 |
#> |    X|               |     32.41 |     3.000 |     5.000 |     5.996 |
#> |.....................|     4.996 |     2.000 |    -127.0 |    -1.000 |
#> |.....................|    -2.001 |    -3.002 |  -0.05276 |     3.299 |
#> |   47|     4502.6532 |   -0.3331 |    0.3335 |    0.7780 |     1.000 |
#> |.....................|    0.7780 |    0.1113 |   -0.3334 |   -0.5552 |
#> |.....................|   -0.7776 |   -0.9997 |   -0.3313 |    0.5328 |
#> |    U|               |     23.33 |     3.001 |     5.001 |     6.002 |
#> |.....................|     5.001 |     2.000 |    -2.949 |   -0.9996 |
#> |.....................|    -2.000 |    -2.999 | -0.009004 |     3.299 |
#> |    X|               |     23.33 |     3.001 |     5.001 |     6.002 |
#> |.....................|     5.001 |     2.000 |    -2.949 |   -0.9996 |
#> |.....................|    -2.000 |    -2.999 | -0.009004 |     3.299 |
#> |   48|     2816.8174 |   -0.3333 |    0.3332 |    0.7777 |    0.9999 |
#> |.....................|    0.7777 |    0.1110 |   -0.3331 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5330 |
#> |    U|               |    0.4947 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     26.88 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01079 |     3.299 |
#> |    X|               |    0.4947 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     26.88 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01079 |     3.299 |
#> |   49|     15282.068 |   -0.3331 |    0.3333 |    0.7779 |     1.000 |
#> |.....................|    0.7776 |    0.1109 |   -0.3333 |   -0.5554 |
#> |.....................|   -0.7779 |   -0.9999 |   -0.3311 |    0.5332 |
#> |    U|               |     27.13 |     3.000 |     5.001 |     6.001 |
#> |.....................|     4.999 |     2.000 |     7.059 |   -0.9999 |
#> |.....................|    -2.000 |    -3.000 |  0.009375 |     3.300 |
#> |    X|               |     27.13 |     3.000 |     5.001 |     6.001 |
#> |.....................|     4.999 |     2.000 |     7.059 |   -0.9999 |
#> |.....................|    -2.000 |    -3.000 |  0.009375 |     3.300 |
#> |   50|     1905.9033 |   -0.3334 |    0.3333 |    0.7779 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3334 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3310 |    0.5333 |
#> |    U|               |    -6.979 |     3.000 |     5.001 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -5.307 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.02165 |     3.300 |
#> |    X|               |    -6.979 |     3.000 |     5.001 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -5.307 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.02165 |     3.300 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   51|     2010.3842 |   -0.3334 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1112 |   -0.3334 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3312 |    0.5332 |
#> |    U|               |    -2.614 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -2.621 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.005898 |     3.300 |
#> |    X|               |    -2.614 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -2.621 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.005898 |     3.300 |
#> |   52|     4274.0783 |   -0.3332 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1112 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5334 |
#> |    U|               |     10.80 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     4.083 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01171 |     3.300 |
#> |    X|               |     10.80 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     4.083 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01171 |     3.300 |
#> |   53|     2036.4211 |   -0.3333 |    0.3334 |    0.7777 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3334 |   -0.5555 |
#> |.....................|   -0.7776 |    -1.000 |   -0.3311 |    0.5332 |
#> |    U|               |     2.867 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -5.071 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01458 |     3.300 |
#> |    X|               |     2.867 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -5.071 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01458 |     3.300 |
#> |   54|     1469.8245 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7777 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7777 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.9834 |     3.000 |     5.000 |     6.001 |
#> |.....................|     5.000 |     2.000 |   -0.5904 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.007330 |     3.300 |
#> |    X|               |   -0.9834 |     3.000 |     5.000 |     6.001 |
#> |.....................|     5.000 |     2.000 |   -0.5904 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.007330 |     3.300 |
#> |   55|     1572.2894 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7777 |    0.1111 |   -0.3334 |   -0.5557 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |    0.5750 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -5.005 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01421 |     3.300 |
#> |    X|               |    0.5750 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -5.005 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01421 |     3.300 |
#> |   56|     1872.6639 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7779 |    0.1110 |   -0.3334 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3312 |    0.5333 |
#> |    U|               |     1.165 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.001 |     2.000 |    -3.891 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.003313 |     3.300 |
#> |    X|               |     1.165 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.001 |     2.000 |    -3.891 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.003313 |     3.300 |
#> |   57|     1281.1571 |   -0.3334 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1110 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |    -1.749 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     1.050 |   -0.9999 |
#> |.....................|    -2.000 |    -3.000 |   0.01185 |     3.300 |
#> |    X|               |    -1.749 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     1.050 |   -0.9999 |
#> |.....................|    -2.000 |    -3.000 |   0.01185 |     3.300 |
#> |   58|     1417.3129 |   -0.3333 |    0.3334 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |    -1.146 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     5.915 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01258 |     3.300 |
#> |    X|               |    -1.146 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |     5.915 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01258 |     3.300 |
#> |   59|     1408.1308 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |    0.3054 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.2681 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01271 |     3.300 |
#> |    X|               |    0.3054 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.2681 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01271 |     3.300 |
#> |   60|     1270.0989 |   -0.3333 |    0.3334 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.5558 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.476 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01016 |     3.300 |
#> |    X|               |   -0.5558 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.476 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01016 |     3.300 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   61|     1308.4439 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |     1.346 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.4625 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01006 |     3.300 |
#> |    X|               |     1.346 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.4625 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01006 |     3.300 |
#> |   62|     1266.6556 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.5466 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.116 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01086 |     3.300 |
#> |    X|               |   -0.5466 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.116 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01086 |     3.300 |
#> |   63|     1318.0955 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |    0.6085 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.8302 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01124 |     3.300 |
#> |    X|               |    0.6085 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.8302 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01124 |     3.300 |
#> |   64|     1266.3067 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.5130 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.8603 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01079 |     3.300 |
#> |    X|               |   -0.5130 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.8603 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01079 |     3.300 |
#> |   65|     1268.5618 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.2610 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.571 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01005 |     3.300 |
#> |    X|               |   -0.2610 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.571 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01005 |     3.300 |
#> |   66|     1266.9184 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.6088 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.023 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01035 |     3.300 |
#> |    X|               |   -0.6088 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.023 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01035 |     3.300 |
#> |   67|     1268.0119 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3576 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.181 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01082 |     3.300 |
#> |    X|               |   -0.3576 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.181 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01082 |     3.300 |
#> |   68|     1265.6817 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3458 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.9467 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01052 |     3.300 |
#> |    X|               |   -0.3458 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.9467 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01052 |     3.300 |
#> |   69|     1267.2488 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.5203 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.212 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01028 |     3.300 |
#> |    X|               |   -0.5203 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.212 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01028 |     3.300 |
#> |   70|     1265.5208 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.2796 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.8161 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01036 |     3.300 |
#> |    X|               |   -0.2796 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.8161 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01036 |     3.300 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   71|     1276.4776 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |  -0.06659 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.252 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01097 |     3.300 |
#> |    X|               |  -0.06659 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |    -1.252 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01097 |     3.300 |
#> |   72|     1272.3755 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.1810 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |  -0.03231 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009530 |     3.300 |
#> |    X|               |   -0.1810 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |  -0.03231 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009530 |     3.300 |
#> |   73|     1265.6821 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.2879 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.3377 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01030 |     3.300 |
#> |    X|               |   -0.2879 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.3377 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01030 |     3.300 |
#> |   74|     1265.5613 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3531 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.9528 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01032 |     3.300 |
#> |    X|               |   -0.3531 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.9528 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01032 |     3.300 |
#> |   75|     1265.8264 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.2473 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.9334 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01020 |     3.300 |
#> |    X|               |   -0.2473 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.9334 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01020 |     3.300 |
#> |   76|     1266.4738 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.1779 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7490 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |    X|               |   -0.1779 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7490 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |   77|     1265.6835 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.2411 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.9019 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01035 |     3.300 |
#> |    X|               |   -0.2411 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.9019 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01035 |     3.300 |
#> |   78|     1265.6737 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.2155 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.8253 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01028 |     3.300 |
#> |    X|               |   -0.2155 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.8253 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01028 |     3.300 |
#> |   79|     1265.3585 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3433 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7299 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |    X|               |   -0.3433 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7299 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |   80|     1265.3591 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3477 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7470 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01039 |     3.300 |
#> |    X|               |   -0.3477 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7470 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01039 |     3.300 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   81|     1265.5519 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3691 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.6327 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01026 |     3.300 |
#> |    X|               |   -0.3691 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.6327 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01026 |     3.300 |
#> |   82|     1265.4163 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3724 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7382 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01032 |     3.300 |
#> |    X|               |   -0.3724 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7382 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01032 |     3.300 |
#> |   83|     1265.3938 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3970 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7935 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01035 |     3.300 |
#> |    X|               |   -0.3970 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7935 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01035 |     3.300 |
#> |   84|     1265.4081 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3595 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7996 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01045 |     3.300 |
#> |    X|               |   -0.3595 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7996 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01045 |     3.300 |
#> |   85|     1265.3724 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3397 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7726 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |    X|               |   -0.3397 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7726 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |   86|     1265.3630 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3600 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7874 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01036 |     3.300 |
#> |    X|               |   -0.3600 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7874 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01036 |     3.300 |
#> |   87|     1265.3485 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3515 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7410 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |    X|               |   -0.3515 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7410 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |   88|     1265.3431 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3618 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7560 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |    X|               |   -0.3618 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7560 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |   89|     1265.3540 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3509 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7616 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |    X|               |   -0.3509 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7616 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |   90|     1265.3525 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3575 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7626 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01036 |     3.300 |
#> |    X|               |   -0.3575 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7626 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01036 |     3.300 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |   91|     1265.3362 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3651 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7532 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01038 |     3.300 |
#> |    X|               |   -0.3651 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7532 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01038 |     3.300 |
#> |   92|     1265.3340 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3648 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7498 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01038 |     3.300 |
#> |    X|               |   -0.3648 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7498 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01038 |     3.300 |
#> |   93|     1265.3302 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3706 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7545 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01038 |     3.300 |
#> |    X|               |   -0.3706 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7545 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01038 |     3.300 |
#> |   94|     1265.3300 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3700 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7553 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01038 |     3.300 |
#> |    X|               |   -0.3700 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7553 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01038 |     3.300 |
#> |   95|     1265.3496 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3550 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7634 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |    X|               |   -0.3550 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7634 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01037 |     3.300 |
#> |   96|     1265.3151 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3834 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7491 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01039 |     3.300 |
#> |    X|               |   -0.3834 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7491 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01039 |     3.300 |
#> |   97|     1265.3008 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3921 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7374 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01040 |     3.300 |
#> |    X|               |   -0.3921 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7374 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01040 |     3.300 |
#> |   98|     1265.2890 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4024 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7367 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01042 |     3.300 |
#> |    X|               |   -0.4024 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7367 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01042 |     3.300 |
#> |   99|     1265.2874 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4117 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7646 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |    X|               |   -0.4117 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7646 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |  100|     1265.2962 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.3999 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7529 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |    X|               |   -0.3999 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7529 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |  101|     1265.2930 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4109 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7819 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |    X|               |   -0.4109 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7819 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |  102|     1265.2905 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4075 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7645 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |    X|               |   -0.4075 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7645 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |  103|     1265.2874 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4104 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7614 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |    X|               |   -0.4104 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7614 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |  104|     1265.2857 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4190 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7553 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01042 |     3.300 |
#> |    X|               |   -0.4190 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7553 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01042 |     3.300 |
#> |  105|     1265.2838 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4195 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7538 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01042 |     3.300 |
#> |    X|               |   -0.4195 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7538 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01042 |     3.300 |
#> |  106|     1265.2867 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4269 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7619 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01042 |     3.300 |
#> |    X|               |   -0.4269 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7619 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01042 |     3.300 |
#> |  107|     1265.2811 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4243 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7594 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |    X|               |   -0.4243 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7594 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |  108|     1265.2843 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4214 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7606 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |    X|               |   -0.4214 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7606 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |  109|     1265.2795 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4254 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7605 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |    X|               |   -0.4254 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7605 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01043 |     3.300 |
#> |  110|     1265.2738 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4276 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7495 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |    X|               |   -0.4276 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7495 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> |  111|     1265.2685 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4287 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7310 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |    X|               |   -0.4287 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7310 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |  112|     1265.2736 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4226 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7331 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01045 |     3.300 |
#> |    X|               |   -0.4226 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7331 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01045 |     3.300 |
#> |  113|     1265.2718 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4247 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7315 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |    X|               |   -0.4247 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7315 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |  114|     1265.2655 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4345 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7308 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01045 |     3.300 |
#> |    X|               |   -0.4345 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7308 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01045 |     3.300 |
#> |  115|     1265.2661 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4305 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7215 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01045 |     3.300 |
#> |    X|               |   -0.4305 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7215 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01045 |     3.300 |
#> |  116|     1265.2705 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4276 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7360 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |    X|               |   -0.4276 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7360 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01044 |     3.300 |
#> |  117|     1265.2621 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01046 |     3.300 |
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01046 |     3.300 |
#> |  118|     1265.2626 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4372 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7217 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01046 |     3.300 |
#> |    X|               |   -0.4372 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7217 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01046 |     3.300 |
#> |  119|     1265.2621 |   -0.3333 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.3333 |   -0.5555 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.5333 |
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01046 |     3.300 |
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01046 |     3.300 |
#> |-----+---------------+-----------+-----------+-----------+-----------|
#> → calculating covariance
#> ✔ done
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling EBE model...
#> ✔ done
#> → loading llik model into symengine environment...
#> → pruning branches (`if`/`else`) of llik full model...
#> ✔ done
#> → finding duplicate expressions in Llik EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in Llik EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling Llik EBE model...
#> ✔ done
#> → Calculating residuals/tables

fitFE
```

``` math
\begin{align*}
{b} & = {b.timeF.0}+{b.timeF.0.5} {\times} \left({timeF}{\equiv}\text{"0.5"}\right)+{b.timeF.1} {\times} \left({timeF}{\equiv}\text{"1"}\right)+{b.timeF.1.5} {\times} \left({timeF}{\equiv}\text{"1.5"}\right)+{b.timeF.2} {\times} \left({timeF}{\equiv}\text{"2"}\right)+{b.timeF.3} {\times} \left({timeF}{\equiv}\text{"3"}\right)+{b.timeF.4} {\times} \left({timeF}{\equiv}\text{"4"}\right)+{b.timeF.6} {\times} \left({timeF}{\equiv}\text{"6"}\right)+{b.timeF.8} {\times} \left({timeF}{\equiv}\text{"8"}\right)+{b.timeF.12} {\times} \left({timeF}{\equiv}\text{"12"}\right)+{cov\_conc\_b} {\times} {conc} \\
{value} & = {b} \\
{value} & \sim add({addSd})
\end{align*}
```

``` r

# Carry the bobyqa estimates forward as starts for focei + random effects
feEst <- setNames(fitFE$iniDf$est, fitFE$iniDf$name)

fitRE <-
  nlmixr(
    dQT ~ b + bRe ~ (bRe|id),
    data  = cqtData,
    start = list(
      b     = c(feEst[1:10], feEst[["cov_conc_b"]]),
      addSd = feEst[["addSd"]]
    ),
    param = list(b ~ timeF + conc),
    est   = "focei"
  )
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate ∂(f)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate ∂(R²)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in inner model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> 
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling inner model...
#> ✔ done
#> → finding duplicate expressions in FD model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling EBE model...
#> ✔ done
#> → compiling events FD model...
#> ✔ done
#> Key: U: Unscaled Parameters; X: Back-transformed parameters; G: Gill difference gradient approximation
#> F: Forward difference gradient approximation
#> C: Central difference gradient approximation
#> M: Mixed forward and central difference gradient approximation
#> A: Analytic (forward sensitivity) gradient (fast=TRUE)
#> Unscaled parameters for Omegas=chol(solve(omega));
#> Diagonals are transformed, as specified by foceiControl(diagXform=)
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |    1|     2231.6817 |   -0.4310 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3310 |    0.4000 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01046 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   0.01046 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    G|      Gill     |     1.132 |    -45.89 |     208.7 |     156.8 |
#> |.....................|     236.8 |     63.34 |     27.26 |    -9.382 |
#> |.....................|    -115.4 |    -174.0 | 2.415e+06 |    -855.0 |
#> |.....................|    -18.95 |...........|...........|...........|
#> |    2| 1.7100212e+11 |   -0.4310 |    0.3334 |    0.7777 |    0.9999 |
#> |.....................|    0.7777 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7777 |   -0.9999 |    -1.331 |    0.4004 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     4.999 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |    -95.63 |     3.301 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     4.999 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |    -95.63 |     3.301 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    3| 1.7101423e+09 |   -0.4310 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.4310 |    0.4000 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |    -9.553 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |    -9.553 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    4|     17082263. |   -0.4310 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3410 |    0.4000 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   -0.9459 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |   -0.9459 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    5|     170861.76 |   -0.4310 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3320 |    0.4000 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  -0.08518 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  -0.08518 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    6|     3700.9266 |   -0.4310 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3311 |    0.4000 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 | 0.0008926 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 | 0.0008926 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    7|     2224.6685 |   -0.4310 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3310 |    0.4000 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009500 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009500 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    M|     Mixed     |    -63.00 |    -61.34 |     186.4 |     126.5 |
#> |.....................|     217.8 |     58.03 |     24.91 |    -9.680 |
#> |.....................|    -124.9 |    -182.9 |-1.009e+06 |    -843.4 |
#> |.....................|   -0.4951 |...........|...........|...........|
#> |    8|     2223.1712 |   -0.4310 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3310 |    0.4000 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009782 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009782 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    F|    Forward    |    -44.04 |    -56.75 |     193.0 |     135.5 |
#> |.....................|     223.5 |     59.61 |     25.62 |    -9.581 |
#> |.....................|    -122.1 |    -180.3 |     2384. |    -844.5 |
#> |.....................|    -10.94 |...........|...........|...........|
#> |    9|     2223.1619 |   -0.4310 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3310 |    0.4000 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009781 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009781 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |   10|     2223.1338 |   -0.4310 |    0.3333 |    0.7778 |     1.000 |
#> |.....................|    0.7778 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3310 |    0.4000 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009781 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009781 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |   11|     2223.0215 |   -0.4310 |    0.3334 |    0.7777 |     1.000 |
#> |.....................|    0.7777 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7778 |    -1.000 |   -0.3310 |    0.4002 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009781 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4395 |     3.000 |     5.000 |     6.000 |
#> |.....................|     5.000 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009781 |     3.300 |
#> |.....................|     1.000 |...........|...........|...........|
#> |   12|     2222.5726 |   -0.4310 |    0.3334 |    0.7776 |    0.9999 |
#> |.....................|    0.7776 |    0.1111 |   -0.4948 |   -0.5556 |
#> |.....................|   -0.7777 |   -0.9999 |   -0.3310 |    0.4006 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4394 |     3.000 |     4.999 |     5.999 |
#> |.....................|     4.999 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009781 |     3.301 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4394 |     3.000 |     4.999 |     5.999 |
#> |.....................|     4.999 |     2.000 |   -0.7266 |    -1.000 |
#> |.....................|    -2.000 |    -3.000 |  0.009781 |     3.301 |
#> |.....................|     1.000 |...........|...........|...........|
#> |   13|     2220.7806 |   -0.4309 |    0.3335 |    0.7772 |    0.9996 |
#> |.....................|    0.7772 |    0.1109 |   -0.4949 |   -0.5555 |
#> |.....................|   -0.7774 |   -0.9995 |   -0.3310 |    0.4023 |
#> |.....................|   -0.1111 |...........|...........|...........|
#> |    U|               |   -0.4392 |     3.000 |     4.997 |     5.998 |
#> |.....................|     4.997 |     2.000 |   -0.7267 |    -1.000 |
#> |.....................|    -1.999 |    -2.998 |  0.009782 |     3.304 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4392 |     3.000 |     4.997 |     5.998 |
#> |.....................|     4.997 |     2.000 |   -0.7267 |    -1.000 |
#> |.....................|    -1.999 |    -2.998 |  0.009782 |     3.304 |
#> |.....................|     1.000 |...........|...........|...........|
#> |   14|     2213.6738 |   -0.4305 |    0.3340 |    0.7757 |    0.9985 |
#> |.....................|    0.7753 |    0.1105 |   -0.4951 |   -0.5555 |
#> |.....................|   -0.7764 |   -0.9980 |   -0.3310 |    0.4093 |
#> |.....................|   -0.1110 |...........|...........|...........|
#> |    U|               |   -0.4384 |     3.002 |     4.989 |     5.991 |
#> |.....................|     4.988 |     1.999 |   -0.7270 |   -0.9999 |
#> |.....................|    -1.997 |    -2.994 |  0.009784 |     3.315 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4384 |     3.002 |     4.989 |     5.991 |
#> |.....................|     4.988 |     1.999 |   -0.7270 |   -0.9999 |
#> |.....................|    -1.997 |    -2.994 |  0.009784 |     3.315 |
#> |.....................|     1.000 |...........|...........|...........|
#> |   15|     2186.1981 |   -0.4291 |    0.3358 |    0.7693 |    0.9941 |
#> |.....................|    0.7680 |    0.1085 |   -0.4959 |   -0.5551 |
#> |.....................|   -0.7724 |   -0.9921 |   -0.3310 |    0.4370 |
#> |.....................|   -0.1106 |...........|...........|...........|
#> |    U|               |   -0.4351 |     3.007 |     4.958 |     5.964 |
#> |.....................|     4.951 |     1.995 |   -0.7281 |   -0.9996 |
#> |.....................|    -1.989 |    -2.976 |  0.009793 |     3.361 |
#> |.....................|     1.000 |...........|...........|...........|
#> |    X|               |   -0.4351 |     3.007 |     4.958 |     5.964 |
#> |.....................|     4.951 |     1.995 |   -0.7281 |   -0.9996 |
#> |.....................|    -1.989 |    -2.976 |  0.009793 |     3.361 |
#> |.....................|     1.000 |...........|...........|...........|
#> |   16|     2090.0745 |   -0.4233 |    0.3433 |    0.7440 |    0.9763 |
#> |.....................|    0.7386 |    0.1007 |   -0.4993 |   -0.5539 |
#> |.....................|   -0.7564 |   -0.9684 |   -0.3310 |    0.5479 |
#> |.....................|   -0.1092 |...........|...........|...........|
#> |    U|               |   -0.4219 |     3.030 |     4.831 |     5.858 |
#> |.....................|     4.804 |     1.979 |   -0.7327 |   -0.9983 |
#> |.....................|    -1.957 |    -2.905 |  0.009831 |     3.544 |
#> |.....................|     1.002 |...........|...........|...........|
#> |    X|               |   -0.4219 |     3.030 |     4.831 |     5.858 |
#> |.....................|     4.804 |     1.979 |   -0.7327 |   -0.9983 |
#> |.....................|    -1.957 |    -2.905 |  0.009831 |     3.544 |
#> |.....................|     1.002 |...........|...........|...........|
#> |   17|     1869.7334 |   -0.4017 |    0.3710 |    0.6497 |    0.9101 |
#> |.....................|    0.6295 |   0.07155 |   -0.5118 |   -0.5492 |
#> |.....................|   -0.6968 |   -0.8804 |   -0.3310 |    0.9605 |
#> |.....................|   -0.1039 |...........|...........|...........|
#> |    U|               |   -0.3729 |     3.113 |     4.359 |     5.460 |
#> |.....................|     4.258 |     1.921 |   -0.7500 |   -0.9936 |
#> |.....................|    -1.838 |    -2.641 |  0.009969 |     4.225 |
#> |.....................|     1.007 |...........|...........|...........|
#> |    X|               |   -0.3729 |     3.113 |     4.359 |     5.460 |
#> |.....................|     4.258 |     1.921 |   -0.7500 |   -0.9936 |
#> |.....................|    -1.838 |    -2.641 |  0.009969 |     4.225 |
#> |.....................|     1.007 |...........|...........|...........|
#> |    M|     Mixed     |    -28.66 |    -30.67 |     107.1 |     72.98 |
#> |.....................|     123.0 |     36.69 |     16.04 |    -5.586 |
#> |.....................|    -71.50 |    -103.5 |-3.145e+04 |    -238.0 |
#> |.....................|    -21.01 |...........|...........|...........|
#> |   18|     1777.8752 |   -0.3178 |    0.4458 |    0.3747 |    0.7292 |
#> |.....................|    0.3172 |  -0.03476 |   -0.5591 |   -0.5340 |
#> |.....................|   -0.5017 |   -0.6044 |   -0.3310 |    0.7929 |
#> |.....................|  0.002176 |...........|...........|...........|
#> |    U|               |   -0.1819 |     3.337 |     2.985 |     4.375 |
#> |.....................|     2.697 |     1.708 |   -0.8150 |   -0.9785 |
#> |.....................|    -1.448 |    -1.813 |   0.01027 |     3.948 |
#> |.....................|     1.113 |...........|...........|...........|
#> |    X|               |   -0.1819 |     3.337 |     2.985 |     4.375 |
#> |.....................|     2.697 |     1.708 |   -0.8150 |   -0.9785 |
#> |.....................|    -1.448 |    -1.813 |   0.01027 |     3.948 |
#> |.....................|     1.113 |...........|...........|...........|
#> |   19|     1601.7054 |  -0.06593 |    0.6703 |   -0.4502 |    0.1866 |
#> |.....................|   -0.6194 |   -0.3537 |   -0.7009 |   -0.4885 |
#> |.....................|   0.08333 |    0.2234 |   -0.3310 |    0.2902 |
#> |.....................|    0.3203 |...........|...........|...........|
#> |    U|               |    0.3912 |     4.011 |    -1.140 |     1.119 |
#> |.....................|    -1.986 |     1.070 |    -1.010 |   -0.9329 |
#> |.....................|   -0.2778 |    0.6702 |   0.01116 |     3.119 |
#> |.....................|     1.431 |...........|...........|...........|
#> |    X|               |    0.3912 |     4.011 |    -1.140 |     1.119 |
#> |.....................|    -1.986 |     1.070 |    -1.010 |   -0.9329 |
#> |.....................|   -0.2778 |    0.6702 |   0.01116 |     3.119 |
#> |.....................|     1.431 |...........|...........|...........|
#> |    M|     Mixed     |     15.82 |     12.60 |     43.41 |     12.33 |
#> |.....................|     37.46 |     77.56 |     41.34 |   -0.9128 |
#> |.....................|    -75.48 |    -70.52 | 2.269e+06 |    -295.7 |
#> |.....................|     9.657 |...........|...........|...........|
#> |   20|     1490.2397 |    0.4255 |    0.5828 |   -0.8239 |    0.4257 |
#> |.....................|   -0.8940 |    -1.616 |    -1.380 |   -0.4692 |
#> |.....................|     1.442 |     1.459 |   -0.3310 |   -0.1175 |
#> |.....................|   0.08807 |...........|...........|...........|
#> |    U|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009211 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009211 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    M|     Mixed     |    -30.17 |    -26.85 |    -102.9 |     63.72 |
#> |.....................|    -67.42 |     52.49 |     50.65 |     4.463 |
#> |.....................|    -71.64 |     21.33 |-2.744e+06 |    -461.2 |
#> |.....................|     29.87 |...........|...........|...........|
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |   21|     1861.3899 |    0.5574 |    0.7519 |   -0.3740 |   -0.5303 |
#> |.....................|   -0.7608 |    -2.879 |    -2.274 |   -0.4898 |
#> |.....................|     2.907 |     1.996 |   -0.3310 |   -0.3758 |
#> |.....................|    0.3996 |...........|...........|...........|
#> |    U|               |     1.810 |     4.256 |   -0.7587 |    -3.182 |
#> |.....................|    -2.693 |    -3.981 |    -3.175 |   -0.9342 |
#> |.....................|     5.369 |     5.989 |   0.01065 |     2.020 |
#> |.....................|     1.511 |...........|...........|...........|
#> |    X|               |     1.810 |     4.256 |   -0.7587 |    -3.182 |
#> |.....................|    -2.693 |    -3.981 |    -3.175 |   -0.9342 |
#> |.....................|     5.369 |     5.989 |   0.01065 |     2.020 |
#> |.....................|     1.511 |...........|...........|...........|
#> |   22|     1533.1526 |    0.4516 |    0.6162 |   -0.7349 |    0.2367 |
#> |.....................|   -0.8677 |    -1.866 |    -1.557 |   -0.4733 |
#> |.....................|     1.732 |     1.566 |   -0.3310 |   -0.1686 |
#> |.....................|    0.1497 |...........|...........|...........|
#> |    U|               |     1.569 |     3.849 |    -2.564 |     1.420 |
#> |.....................|    -3.227 |    -1.954 |    -2.188 |   -0.9177 |
#> |.....................|     3.019 |     4.697 |   0.01095 |     2.362 |
#> |.....................|     1.261 |...........|...........|...........|
#> |    X|               |     1.569 |     3.849 |    -2.564 |     1.420 |
#> |.....................|    -3.227 |    -1.954 |    -2.188 |   -0.9177 |
#> |.....................|     3.019 |     4.697 |   0.01095 |     2.362 |
#> |.....................|     1.261 |...........|...........|...........|
#> |   23|     1544.6382 |    0.4344 |    0.5941 |   -0.7938 |    0.3617 |
#> |.....................|   -0.8851 |    -1.701 |    -1.440 |   -0.4706 |
#> |.....................|     1.540 |     1.495 |   -0.3310 |   -0.1348 |
#> |.....................|    0.1089 |...........|...........|...........|
#> |    U|               |     1.530 |     3.782 |    -2.858 |     2.170 |
#> |.....................|    -3.314 |    -1.623 |    -2.027 |   -0.9150 |
#> |.....................|     2.636 |     4.486 |   0.01100 |     2.418 |
#> |.....................|     1.220 |...........|...........|...........|
#> |    X|               |     1.530 |     3.782 |    -2.858 |     2.170 |
#> |.....................|    -3.314 |    -1.623 |    -2.027 |   -0.9150 |
#> |.....................|     2.636 |     4.486 |   0.01100 |     2.418 |
#> |.....................|     1.220 |...........|...........|...........|
#> |   24|     1552.1000 |    0.4279 |    0.5859 |   -0.8157 |    0.4083 |
#> |.....................|   -0.8916 |    -1.639 |    -1.396 |   -0.4696 |
#> |.....................|     1.469 |     1.469 |   -0.3310 |   -0.1222 |
#> |.....................|   0.09374 |...........|...........|...........|
#> |    U|               |     1.515 |     3.758 |    -2.967 |     2.450 |
#> |.....................|    -3.347 |    -1.500 |    -1.967 |   -0.9140 |
#> |.....................|     2.494 |     4.407 |   0.01102 |     2.438 |
#> |.....................|     1.205 |...........|...........|...........|
#> |    X|               |     1.515 |     3.758 |    -2.967 |     2.450 |
#> |.....................|    -3.347 |    -1.500 |    -1.967 |   -0.9140 |
#> |.....................|     2.494 |     4.407 |   0.01102 |     2.438 |
#> |.....................|     1.205 |...........|...........|...........|
#> |   25|     1554.5222 |    0.4261 |    0.5835 |   -0.8219 |    0.4215 |
#> |.....................|   -0.8935 |    -1.622 |    -1.384 |   -0.4693 |
#> |.....................|     1.449 |     1.462 |   -0.3310 |   -0.1186 |
#> |.....................|   0.08941 |...........|...........|...........|
#> |    U|               |     1.511 |     3.751 |    -2.999 |     2.529 |
#> |.....................|    -3.356 |    -1.465 |    -1.950 |   -0.9138 |
#> |.....................|     2.453 |     4.385 |   0.01102 |     2.444 |
#> |.....................|     1.201 |...........|...........|...........|
#> |    X|               |     1.511 |     3.751 |    -2.999 |     2.529 |
#> |.....................|    -3.356 |    -1.465 |    -1.950 |   -0.9138 |
#> |.....................|     2.453 |     4.385 |   0.01102 |     2.444 |
#> |.....................|     1.201 |...........|...........|...........|
#> |   26|     1555.1223 |    0.4257 |    0.5830 |   -0.8234 |    0.4247 |
#> |.....................|   -0.8939 |    -1.617 |    -1.381 |   -0.4693 |
#> |.....................|     1.444 |     1.460 |   -0.3310 |   -0.1178 |
#> |.....................|   0.08838 |...........|...........|...........|
#> |    U|               |     1.510 |     3.749 |    -3.006 |     2.548 |
#> |.....................|    -3.358 |    -1.457 |    -1.946 |   -0.9137 |
#> |.....................|     2.443 |     4.380 |   0.01103 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.510 |     3.749 |    -3.006 |     2.548 |
#> |.....................|    -3.358 |    -1.457 |    -1.946 |   -0.9137 |
#> |.....................|     2.443 |     4.380 |   0.01103 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   27|     1555.2592 |    0.4256 |    0.5828 |   -0.8238 |    0.4254 |
#> |.....................|   -0.8940 |    -1.616 |    -1.380 |   -0.4692 |
#> |.....................|     1.443 |     1.460 |   -0.3310 |   -0.1176 |
#> |.....................|   0.08814 |...........|...........|...........|
#> |    U|               |     1.510 |     3.748 |    -3.008 |     2.553 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.441 |     4.379 |   0.01103 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.510 |     3.748 |    -3.008 |     2.553 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.441 |     4.379 |   0.01103 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   28|     1555.2899 |    0.4255 |    0.5828 |   -0.8239 |    0.4256 |
#> |.....................|   -0.8940 |    -1.616 |    -1.380 |   -0.4692 |
#> |.....................|     1.442 |     1.459 |   -0.3310 |   -0.1175 |
#> |.....................|   0.08809 |...........|...........|...........|
#> |    U|               |     1.510 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.441 |     4.378 |   0.01103 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.510 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.441 |     4.378 |   0.01103 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   29|     1555.2972 |    0.4255 |    0.5828 |   -0.8239 |    0.4256 |
#> |.....................|   -0.8940 |    -1.616 |    -1.380 |   -0.4692 |
#> |.....................|     1.442 |     1.459 |   -0.3310 |   -0.1175 |
#> |.....................|   0.08808 |...........|...........|...........|
#> |    U|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |   0.01103 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |   0.01103 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   30|     1485.3862 |    0.4255 |    0.5828 |   -0.8239 |    0.4257 |
#> |.....................|   -0.8940 |    -1.616 |    -1.380 |   -0.4692 |
#> |.....................|     1.442 |     1.459 |   -0.3310 |   -0.1175 |
#> |.....................|   0.08807 |...........|...........|...........|
#> |    U|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009777 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009777 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    M|     Mixed     |     42.96 |    -9.712 |    -78.05 |     97.34 |
#> |.....................|    -46.09 |     58.54 |     53.41 |     4.951 |
#> |.....................|    -61.06 |     31.47 | 1.102e+06 |    -455.2 |
#> |.....................|     28.45 |...........|...........|...........|
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |   31|     1484.4516 |    0.4255 |    0.5828 |   -0.8239 |    0.4257 |
#> |.....................|   -0.8940 |    -1.616 |    -1.380 |   -0.4692 |
#> |.....................|     1.442 |     1.459 |   -0.3310 |   -0.1175 |
#> |.....................|   0.08807 |...........|...........|...........|
#> |    U|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    F|    Forward    |     22.16 |    -14.52 |    -85.03 |     87.82 |
#> |.....................|    -52.05 |     56.94 |     52.68 |     4.838 |
#> |.....................|    -64.02 |     28.66 |     6414. |    -454.9 |
#> |.....................|     26.26 |...........|...........|...........|
#> |   32|     1484.4512 |    0.4255 |    0.5828 |   -0.8239 |    0.4257 |
#> |.....................|   -0.8940 |    -1.616 |    -1.380 |   -0.4692 |
#> |.....................|     1.442 |     1.459 |   -0.3310 |   -0.1175 |
#> |.....................|   0.08807 |...........|...........|...........|
#> |    U|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   33|     1484.4502 |    0.4255 |    0.5828 |   -0.8239 |    0.4256 |
#> |.....................|   -0.8940 |    -1.616 |    -1.380 |   -0.4692 |
#> |.....................|     1.442 |     1.459 |   -0.3310 |   -0.1175 |
#> |.....................|   0.08807 |...........|...........|...........|
#> |    U|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   34|     1484.4459 |    0.4255 |    0.5828 |   -0.8239 |    0.4256 |
#> |.....................|   -0.8940 |    -1.616 |    -1.380 |   -0.4692 |
#> |.....................|     1.442 |     1.459 |   -0.3310 |   -0.1175 |
#> |.....................|   0.08807 |...........|...........|...........|
#> |    U|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.509 |     3.748 |    -3.008 |     2.554 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.440 |     4.378 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   35|     1484.4288 |    0.4256 |    0.5828 |   -0.8238 |    0.4255 |
#> |.....................|   -0.8940 |    -1.616 |    -1.380 |   -0.4692 |
#> |.....................|     1.443 |     1.459 |   -0.3310 |   -0.1175 |
#> |.....................|   0.08807 |...........|...........|...........|
#> |    U|               |     1.510 |     3.748 |    -3.008 |     2.553 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.441 |     4.378 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.510 |     3.748 |    -3.008 |     2.553 |
#> |.....................|    -3.359 |    -1.455 |    -1.945 |   -0.9137 |
#> |.....................|     2.441 |     4.378 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   36|     1484.3604 |    0.4256 |    0.5829 |   -0.8237 |    0.4252 |
#> |.....................|   -0.8940 |    -1.617 |    -1.380 |   -0.4693 |
#> |.....................|     1.443 |     1.460 |   -0.3310 |   -0.1176 |
#> |.....................|   0.08808 |...........|...........|...........|
#> |    U|               |     1.510 |     3.749 |    -3.007 |     2.551 |
#> |.....................|    -3.359 |    -1.456 |    -1.945 |   -0.9137 |
#> |.....................|     2.442 |     4.379 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.510 |     3.749 |    -3.007 |     2.551 |
#> |.....................|    -3.359 |    -1.456 |    -1.945 |   -0.9137 |
#> |.....................|     2.442 |     4.379 |  0.009614 |     2.446 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   37|     1484.0875 |    0.4257 |    0.5832 |   -0.8231 |    0.4238 |
#> |.....................|   -0.8939 |    -1.618 |    -1.381 |   -0.4693 |
#> |.....................|     1.445 |     1.460 |   -0.3310 |   -0.1180 |
#> |.....................|   0.08809 |...........|...........|...........|
#> |    U|               |     1.510 |     3.749 |    -3.004 |     2.543 |
#> |.....................|    -3.358 |    -1.459 |    -1.947 |   -0.9137 |
#> |.....................|     2.445 |     4.381 |  0.009614 |     2.445 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.510 |     3.749 |    -3.004 |     2.543 |
#> |.....................|    -3.358 |    -1.459 |    -1.947 |   -0.9137 |
#> |.....................|     2.445 |     4.381 |  0.009614 |     2.445 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   38|     1483.0088 |    0.4262 |    0.5843 |   -0.8207 |    0.4181 |
#> |.....................|   -0.8934 |    -1.625 |    -1.386 |   -0.4694 |
#> |.....................|     1.453 |     1.463 |   -0.3310 |   -0.1195 |
#> |.....................|   0.08812 |...........|...........|...........|
#> |    U|               |     1.511 |     3.753 |    -2.992 |     2.508 |
#> |.....................|    -3.356 |    -1.472 |    -1.953 |   -0.9138 |
#> |.....................|     2.461 |     4.389 |  0.009614 |     2.443 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.511 |     3.753 |    -2.992 |     2.508 |
#> |.....................|    -3.356 |    -1.472 |    -1.953 |   -0.9138 |
#> |.....................|     2.461 |     4.389 |  0.009614 |     2.443 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   39|     1478.9035 |    0.4283 |    0.5887 |   -0.8112 |    0.3954 |
#> |.....................|   -0.8914 |    -1.651 |    -1.405 |   -0.4699 |
#> |.....................|     1.483 |     1.474 |   -0.3310 |   -0.1253 |
#> |.....................|   0.08827 |...........|...........|...........|
#> |    U|               |     1.516 |     3.766 |    -2.945 |     2.372 |
#> |.....................|    -3.346 |    -1.524 |    -1.979 |   -0.9143 |
#> |.....................|     2.522 |     4.421 |  0.009614 |     2.433 |
#> |.....................|     1.199 |...........|...........|...........|
#> |    X|               |     1.516 |     3.766 |    -2.945 |     2.372 |
#> |.....................|    -3.346 |    -1.524 |    -1.979 |   -0.9143 |
#> |.....................|     2.522 |     4.421 |  0.009614 |     2.433 |
#> |.....................|     1.199 |...........|...........|...........|
#> |   40|     1465.9855 |    0.4364 |    0.6063 |   -0.7733 |    0.3046 |
#> |.....................|   -0.8836 |    -1.755 |    -1.480 |   -0.4717 |
#> |.....................|     1.606 |     1.516 |   -0.3310 |   -0.1487 |
#> |.....................|   0.08887 |...........|...........|...........|
#> |    U|               |     1.534 |     3.819 |    -2.755 |     1.828 |
#> |.....................|    -3.307 |    -1.732 |    -2.082 |   -0.9161 |
#> |.....................|     2.768 |     4.549 |  0.009617 |     2.395 |
#> |.....................|     1.200 |...........|...........|...........|
#> |    X|               |     1.534 |     3.819 |    -2.755 |     1.828 |
#> |.....................|    -3.307 |    -1.732 |    -2.082 |   -0.9161 |
#> |.....................|     2.768 |     4.549 |  0.009617 |     2.395 |
#> |.....................|     1.200 |...........|...........|...........|
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |   41|     1458.4244 |    0.4513 |    0.6386 |   -0.7038 |    0.1383 |
#> |.....................|   -0.8692 |    -1.946 |    -1.617 |   -0.4751 |
#> |.....................|     1.832 |     1.595 |   -0.3310 |   -0.1917 |
#> |.....................|   0.08995 |...........|...........|...........|
#> |    U|               |     1.568 |     3.916 |    -2.408 |    0.8298 |
#> |.....................|    -3.235 |    -2.114 |    -2.271 |   -0.9195 |
#> |.....................|     3.219 |     4.784 |  0.009621 |     2.324 |
#> |.....................|     1.201 |...........|...........|...........|
#> |    X|               |     1.568 |     3.916 |    -2.408 |    0.8298 |
#> |.....................|    -3.235 |    -2.114 |    -2.271 |   -0.9195 |
#> |.....................|     3.219 |     4.784 |  0.009621 |     2.324 |
#> |.....................|     1.201 |...........|...........|...........|
#> |    M|     Mixed     |     23.27 |    -5.970 |    -45.30 |    -50.49 |
#> |.....................|    -44.08 |     45.21 |     52.86 |     6.109 |
#> |.....................|    -46.11 |     52.38 |-6.110e+04 |    -494.7 |
#> |.....................|     33.03 |...........|...........|...........|
#> |   42|     1439.6601 |    0.3331 |    0.8146 |   -0.7291 |    0.1861 |
#> |.....................|   -0.9755 |    -2.199 |    -1.951 |   -0.4973 |
#> |.....................|     2.213 |     1.391 |   -0.3310 |   -0.2253 |
#> |.....................|    0.2254 |...........|...........|...........|
#> |    U|               |     1.299 |     4.444 |    -2.534 |     1.117 |
#> |.....................|    -3.766 |    -2.620 |    -2.731 |   -0.9417 |
#> |.....................|     3.981 |     4.174 |  0.009783 |     2.268 |
#> |.....................|     1.336 |...........|...........|...........|
#> |    X|               |     1.299 |     4.444 |    -2.534 |     1.117 |
#> |.....................|    -3.766 |    -2.620 |    -2.731 |   -0.9417 |
#> |.....................|     3.981 |     4.174 |  0.009783 |     2.268 |
#> |.....................|     1.336 |...........|...........|...........|
#> |    M|     Mixed     |    -18.12 |     18.42 |    -58.53 |    -26.36 |
#> |.....................|    -89.90 |     29.95 |     43.97 |     4.349 |
#> |.....................|    -24.94 |     24.25 |-3.028e+05 |    -512.5 |
#> |.....................|     38.92 |...........|...........|...........|
#> |   43|     1430.1295 |    0.5745 |    0.8493 |   -0.8955 |    0.1806 |
#> |.....................|   -0.7822 |    -2.252 |    -2.368 |   -0.5120 |
#> |.....................|     2.447 |     1.259 |   -0.3310 |   -0.2217 |
#> |.....................|    0.4864 |...........|...........|...........|
#> |    U|               |     1.849 |     4.548 |    -3.366 |     1.083 |
#> |.....................|    -2.800 |    -2.727 |    -3.304 |   -0.9564 |
#> |.....................|     4.449 |     3.775 |  0.009341 |     2.274 |
#> |.....................|     1.598 |...........|...........|...........|
#> |    X|               |     1.849 |     4.548 |    -3.366 |     1.083 |
#> |.....................|    -2.800 |    -2.727 |    -3.304 |   -0.9564 |
#> |.....................|     4.449 |     3.775 |  0.009341 |     2.274 |
#> |.....................|     1.598 |...........|...........|...........|
#> |    M|     Mixed     |     23.29 |     23.78 |    -118.2 |    -27.47 |
#> |.....................|    -7.107 |     31.76 |     36.43 |     9.023 |
#> |.....................|    -9.137 |     12.58 |-2.489e+05 |    -507.7 |
#> |.....................|     30.98 |...........|...........|...........|
#> |   44|     1415.7885 |     1.048 |    0.6450 |   -0.6297 |    0.3022 |
#> |.....................|   -0.8770 |    -2.290 |    -3.321 |   -0.5940 |
#> |.....................|     2.685 |     1.237 |   -0.3310 |   -0.1851 |
#> |.....................|     1.143 |...........|...........|...........|
#> |    U|               |     2.926 |     3.935 |    -2.037 |     1.813 |
#> |.....................|    -3.274 |    -2.802 |    -4.616 |    -1.038 |
#> |.....................|     4.925 |     3.709 |  0.008388 |     2.335 |
#> |.....................|     2.254 |...........|...........|...........|
#> |    X|               |     2.926 |     3.935 |    -2.037 |     1.813 |
#> |.....................|    -3.274 |    -2.802 |    -4.616 |    -1.038 |
#> |.....................|     4.925 |     3.709 |  0.008388 |     2.335 |
#> |.....................|     2.254 |...........|...........|...........|
#> |    M|     Mixed     |     131.4 |    -3.491 |    -6.149 |     41.29 |
#> |.....................|    -26.63 |     38.49 |     17.94 |     17.33 |
#> |.....................|     7.504 |     22.91 | 6.338e+05 |    -463.9 |
#> |.....................|     11.41 |...........|...........|...........|
#> |   45|     1390.3976 |    0.8375 |    0.2613 |   -0.6241 |    0.3568 |
#> |.....................|   -0.7204 |    -2.489 |    -4.301 |   -0.8064 |
#> |.....................|     2.640 |     1.471 |   -0.3310 |   -0.1327 |
#> |.....................|     1.808 |...........|...........|...........|
#> |    U|               |     2.447 |     2.784 |    -2.009 |     2.141 |
#> |.....................|    -2.491 |    -3.201 |    -5.965 |    -1.251 |
#> |.....................|     4.835 |     4.412 |  0.008749 |     2.421 |
#> |.....................|     2.919 |...........|...........|...........|
#> |    X|               |     2.447 |     2.784 |    -2.009 |     2.141 |
#> |.....................|    -2.491 |    -3.201 |    -5.965 |    -1.251 |
#> |.....................|     4.835 |     4.412 |  0.008749 |     2.421 |
#> |.....................|     2.919 |...........|...........|...........|
#> |    M|     Mixed     |     25.58 |    -52.67 |    -10.05 |     60.21 |
#> |.....................|     20.16 |     19.88 |    -13.04 |     9.052 |
#> |.....................|     2.440 |     43.22 |-1.248e+06 |    -382.0 |
#> |.....................|     2.319 |...........|...........|...........|
#> |   46|     1351.7125 |     1.035 |    0.4373 |   -0.6433 |    0.2567 |
#> |.....................|   -0.7524 |    -3.319 |    -4.615 |    -1.235 |
#> |.....................|     1.856 |     1.188 |   -0.3310 |  -0.04239 |
#> |.....................|     1.703 |...........|...........|...........|
#> |    U|               |     2.896 |     3.312 |    -2.105 |     1.540 |
#> |.....................|    -2.651 |    -4.861 |    -6.397 |    -1.679 |
#> |.....................|     3.267 |     3.564 |  0.009217 |     2.570 |
#> |.....................|     2.814 |...........|...........|...........|
#> |    X|               |     2.896 |     3.312 |    -2.105 |     1.540 |
#> |.....................|    -2.651 |    -4.861 |    -6.397 |    -1.679 |
#> |.....................|     3.267 |     3.564 |  0.009217 |     2.570 |
#> |.....................|     2.814 |...........|...........|...........|
#> |    M|     Mixed     |     98.81 |     7.232 |     40.46 |     78.88 |
#> |.....................|     60.77 |    -3.401 |    -6.636 |     10.52 |
#> |.....................|    -13.28 |     37.07 | 2.477e+06 |    -249.2 |
#> |.....................|     2.073 |...........|...........|...........|
#> |   47|     1319.9589 |     1.079 |   -0.3061 |   -0.8409 |  -0.01639 |
#> |.....................|   -0.9879 |    -3.512 |    -4.573 |    -1.738 |
#> |.....................|     1.759 |    0.5090 |   -0.3310 |   0.02067 |
#> |.....................|     1.178 |...........|...........|...........|
#> |    U|               |     2.997 |     1.082 |    -3.093 |  -0.09831 |
#> |.....................|    -3.828 |    -5.246 |    -6.339 |    -2.183 |
#> |.....................|     3.073 |     1.527 |  0.009830 |     2.674 |
#> |.....................|     2.289 |...........|...........|...........|
#> |    X|               |     2.997 |     1.082 |    -3.093 |  -0.09831 |
#> |.....................|    -3.828 |    -5.246 |    -6.339 |    -2.183 |
#> |.....................|     3.073 |     1.527 |  0.009830 |     2.674 |
#> |.....................|     2.289 |...........|...........|...........|
#> |    M|     Mixed     |     6.990 |    -41.97 |     22.42 |     14.70 |
#> |.....................|     27.40 |    0.8461 |     2.160 |     8.088 |
#> |.....................|   0.09283 |    -14.25 |-4.051e+05 |    -159.2 |
#> |.....................|     3.511 |...........|...........|...........|
#> |   48|     1311.8358 |    0.9672 |   0.03010 |   -0.8236 |  -0.08934 |
#> |.....................|    -1.024 |    -2.820 |    -4.883 |    -2.459 |
#> |.....................|     1.765 |    0.5314 |   -0.3310 |   0.03354 |
#> |.....................|    0.4439 |...........|...........|...........|
#> |    U|               |     2.742 |     2.090 |    -3.007 |   -0.5360 |
#> |.....................|    -4.009 |    -3.863 |    -6.767 |    -2.904 |
#> |.....................|     3.085 |     1.594 |   0.01004 |     2.695 |
#> |.....................|     1.555 |...........|...........|...........|
#> |    X|               |     2.742 |     2.090 |    -3.007 |   -0.5360 |
#> |.....................|    -4.009 |    -3.863 |    -6.767 |    -2.904 |
#> |.....................|     3.085 |     1.594 |   0.01004 |     2.695 |
#> |.....................|     1.555 |...........|...........|...........|
#> |    M|     Mixed     |     3.307 |    -8.749 |     24.40 |    -16.18 |
#> |.....................|     13.53 |     29.11 |    -6.280 |    -1.848 |
#> |.....................|   -0.4123 |    -14.74 | 5.461e+05 |    -128.2 |
#> |.....................|     13.88 |...........|...........|...........|
#> |   49|     1314.5653 |    0.9891 |    0.4409 |   -0.7020 |   0.06743 |
#> |.....................|   -0.8862 |    -2.928 |    -4.167 |    -3.384 |
#> |.....................|     2.042 |    0.8049 |   -0.3310 |   0.07909 |
#> |.....................|    0.3738 |...........|...........|...........|
#> |    U|               |     2.792 |     3.323 |    -2.399 |    0.4046 |
#> |.....................|    -3.320 |    -4.078 |    -5.780 |    -3.828 |
#> |.....................|     3.639 |     2.415 |  0.008925 |     2.770 |
#> |.....................|     1.485 |...........|...........|...........|
#> |    X|               |     2.792 |     3.323 |    -2.399 |    0.4046 |
#> |.....................|    -3.320 |    -4.078 |    -5.780 |    -3.828 |
#> |.....................|     3.639 |     2.415 |  0.008925 |     2.770 |
#> |.....................|     1.485 |...........|...........|...........|
#> |   50|     1308.5739 |    0.9770 |    0.2149 |   -0.7689 |  -0.01881 |
#> |.....................|   -0.9620 |    -2.869 |    -4.561 |    -2.875 |
#> |.....................|     1.889 |    0.6544 |   -0.3310 |   0.05403 |
#> |.....................|    0.4123 |...........|...........|...........|
#> |    U|               |     2.764 |     2.645 |    -2.733 |   -0.1129 |
#> |.....................|    -3.699 |    -3.960 |    -6.323 |    -3.320 |
#> |.....................|     3.334 |     1.963 |  0.009412 |     2.729 |
#> |.....................|     1.523 |...........|...........|...........|
#> |    X|               |     2.764 |     2.645 |    -2.733 |   -0.1129 |
#> |.....................|    -3.699 |    -3.960 |    -6.323 |    -3.320 |
#> |.....................|     3.334 |     1.963 |  0.009412 |     2.729 |
#> |.....................|     1.523 |...........|...........|...........|
#> |    M|     Mixed     |    -45.94 |    -10.67 |     8.493 |    -27.93 |
#> |.....................|     2.974 |     17.80 |    -4.169 |    -8.256 |
#> |.....................|    -7.613 |    -16.21 |-1.995e+06 |    -109.5 |
#> |.....................|     15.43 |...........|...........|...........|
#> 
#> |    #| Function Val. | b.timeF.0 |b.timeF.0.5 | b.timeF.1 |b.timeF.1.5 |
#> |.....................| b.timeF.2 | b.timeF.3 | b.timeF.4 | b.timeF.6 |
#> |.....................| b.timeF.8 |b.timeF.12 |cov_conc_b |     addSd |
#> |.....................|        o1 |...........|...........|...........|
#> |   51|     1295.9200 |    0.9598 |    0.2363 |   -0.8217 |  -0.05350 |
#> |.....................|    -1.021 |    -3.032 |    -4.357 |    -3.005 |
#> |.....................|     1.834 |    0.6277 |   -0.3310 |   0.07925 |
#> |.....................|  -0.09737 |...........|...........|...........|
#> |    U|               |     2.725 |     2.709 |    -2.997 |   -0.3210 |
#> |.....................|    -3.992 |    -4.286 |    -6.042 |    -3.449 |
#> |.....................|     3.223 |     1.883 |  0.009847 |     2.771 |
#> |.....................|     1.014 |...........|...........|...........|
#> |    X|               |     2.725 |     2.709 |    -2.997 |   -0.3210 |
#> |.....................|    -3.992 |    -4.286 |    -6.042 |    -3.449 |
#> |.....................|     3.223 |     1.883 |  0.009847 |     2.771 |
#> |.....................|     1.014 |...........|...........|...........|
#> |    M|     Mixed     |    -4.933 |     4.716 |     14.01 |    -14.38 |
#> |.....................|     5.050 |     15.96 |     2.914 |    -8.141 |
#> |.....................|    -1.515 |    -9.374 | 2.284e+05 |    -60.47 |
#> |.....................|    -9.150 |...........|...........|...........|
#> |   52|     1290.0330 |     1.051 |    0.1252 |   -0.9004 |  -0.07503 |
#> |.....................|    -1.083 |    -3.365 |    -4.543 |    -2.749 |
#> |.....................|     1.635 |    0.5685 |   -0.3310 |    0.1143 |
#> |.....................|    0.1698 |...........|...........|...........|
#> |    U|               |     2.933 |     2.376 |    -3.391 |   -0.4502 |
#> |.....................|    -4.302 |    -4.952 |    -6.298 |    -3.194 |
#> |.....................|     2.826 |     1.705 |   0.01007 |     2.829 |
#> |.....................|     1.281 |...........|...........|...........|
#> |    X|               |     2.933 |     2.376 |    -3.391 |   -0.4502 |
#> |.....................|    -4.302 |    -4.952 |    -6.298 |    -3.194 |
#> |.....................|     2.826 |     1.705 |   0.01007 |     2.829 |
#> |.....................|     1.281 |...........|...........|...........|
#> |    M|     Mixed     |     14.50 |     6.509 |     12.64 |     2.220 |
#> |.....................|     7.149 |     8.379 |     3.261 |    -2.814 |
#> |.....................|    -1.654 |    -4.416 | 8.872e+05 |    -45.38 |
#> |.....................|     12.56 |...........|...........|...........|
#> |   53|     1288.7485 |     1.064 |   0.09567 |   -0.9735 |   -0.1301 |
#> |.....................|    -1.160 |    -3.478 |    -4.840 |    -2.336 |
#> |.....................|     1.403 |    0.5230 |   -0.3310 |    0.1185 |
#> |.....................|   0.05145 |...........|...........|...........|
#> |    U|               |     2.963 |     2.287 |    -3.757 |   -0.7804 |
#> |.....................|    -4.690 |    -5.178 |    -6.707 |    -2.780 |
#> |.....................|     2.361 |     1.569 |   0.01005 |     2.835 |
#> |.....................|     1.163 |...........|...........|...........|
#> |    X|               |     2.963 |     2.287 |    -3.757 |   -0.7804 |
#> |.....................|    -4.690 |    -5.178 |    -6.707 |    -2.780 |
#> |.....................|     2.361 |     1.569 |   0.01005 |     2.835 |
#> |.....................|     1.163 |...........|...........|...........|
#> |    M|     Mixed     |    -12.96 |     5.927 |    -2.001 |    -13.30 |
#> |.....................|    -8.459 |     5.416 |    -1.223 |     2.186 |
#> |.....................|    -9.471 |    -6.137 |-5.288e+05 |    -37.09 |
#> |.....................|     7.903 |...........|...........|...........|
#> |   54|     1285.7449 |     1.072 |   0.02475 |   -0.9706 |  -0.08340 |
#> |.....................|    -1.128 |    -3.648 |    -4.548 |    -2.634 |
#> |.....................|     1.670 |    0.5860 |   -0.3310 |    0.1657 |
#> |.....................|  0.007691 |...........|...........|...........|
#> |    U|               |     2.981 |     2.074 |    -3.742 |   -0.5004 |
#> |.....................|    -4.527 |    -5.518 |    -6.305 |    -3.078 |
#> |.....................|     2.895 |     1.758 |  0.009912 |     2.913 |
#> |.....................|     1.119 |...........|...........|...........|
#> |    X|               |     2.981 |     2.074 |    -3.742 |   -0.5004 |
#> |.....................|    -4.527 |    -5.518 |    -6.305 |    -3.078 |
#> |.....................|     2.895 |     1.758 |  0.009912 |     2.913 |
#> |.....................|     1.119 |...........|...........|...........|
#> |    M|     Mixed     |    -13.14 |    -3.834 |    -6.292 |    -3.515 |
#> |.....................|    -4.838 |    -2.670 |     3.329 |    -1.013 |
#> |.....................|    -1.049 |    -2.741 |-5.602e+05 |    -8.140 |
#> |.....................|     2.303 |...........|...........|...........|
#> |   55|     1285.0876 |     1.058 |   0.05762 |   -0.9626 |  -0.08390 |
#> |.....................|    -1.122 |    -3.612 |    -4.690 |    -2.449 |
#> |.....................|     1.662 |    0.6066 |   -0.3310 |    0.1759 |
#> |.....................|  -0.03350 |...........|...........|...........|
#> |    U|               |     2.949 |     2.173 |    -3.702 |   -0.5034 |
#> |.....................|    -4.500 |    -5.447 |    -6.501 |    -2.894 |
#> |.....................|     2.880 |     1.820 |   0.01010 |     2.930 |
#> |.....................|     1.078 |...........|...........|...........|
#> |    X|               |     2.949 |     2.173 |    -3.702 |   -0.5034 |
#> |.....................|    -4.500 |    -5.447 |    -6.501 |    -2.894 |
#> |.....................|     2.880 |     1.820 |   0.01010 |     2.930 |
#> |.....................|     1.078 |...........|...........|...........|
#> |    M|     Mixed     |     3.998 |     2.222 |    0.1799 |     2.717 |
#> |.....................|    0.2137 |   -0.4097 |     1.069 |    0.5735 |
#> |.....................|    0.6411 |    0.6238 | 3.452e+05 |   -0.9201 |
#> |.....................|    -2.671 |...........|...........|...........|
#> |   56|     1284.8241 |     1.113 |   0.01677 |   -0.9805 |   -0.1044 |
#> |.....................|    -1.142 |    -3.643 |    -4.813 |    -2.657 |
#> |.....................|     1.609 |    0.5708 |   -0.3310 |    0.1743 |
#> |.....................|  0.001468 |...........|...........|...........|
#> |    U|               |     3.074 |     2.050 |    -3.791 |   -0.6261 |
#> |.....................|    -4.600 |    -5.509 |    -6.670 |    -3.101 |
#> |.....................|     2.773 |     1.712 |   0.01004 |     2.928 |
#> |.....................|     1.113 |...........|...........|...........|
#> |    X|               |     3.074 |     2.050 |    -3.791 |   -0.6261 |
#> |.....................|    -4.600 |    -5.509 |    -6.670 |    -3.101 |
#> |.....................|     2.773 |     1.712 |   0.01004 |     2.928 |
#> |.....................|     1.113 |...........|...........|...........|
#> |    M|     Mixed     |    0.4003 |    0.7889 |   -0.3798 |   -0.1134 |
#> |.....................|   -0.5763 |    0.1752 |    0.2029 |   -0.2980 |
#> |.....................|   0.05393 |    0.1497 | 8.452e+04 |    -2.857 |
#> |.....................|    0.9951 |...........|...........|...........|
#> |   57|     1284.8056 |     1.081 |   0.06167 |   -0.9464 |  -0.07973 |
#> |.....................|    -1.109 |    -3.574 |    -4.727 |    -2.529 |
#> |.....................|     1.683 |    0.6051 |   -0.3310 |    0.1815 |
#> |.....................|  0.001956 |...........|...........|...........|
#> |    U|               |     3.000 |     2.185 |    -3.621 |   -0.4784 |
#> |.....................|    -4.432 |    -5.371 |    -6.551 |    -2.973 |
#> |.....................|     2.922 |     1.815 |  0.009945 |     2.940 |
#> |.....................|     1.113 |...........|...........|...........|
#> |    X|               |     3.000 |     2.185 |    -3.621 |   -0.4784 |
#> |.....................|    -4.432 |    -5.371 |    -6.551 |    -2.973 |
#> |.....................|     2.922 |     1.815 |  0.009945 |     2.940 |
#> |.....................|     1.113 |...........|...........|...........|
#> |    M|     Mixed     |    -2.298 |   -0.2032 |   0.06592 |    -1.283 |
#> |.....................|    0.1484 |    0.2024 |    0.1460 |  -0.06893 |
#> |.....................|   -0.2475 |   -0.8965 |-6.608e+04 |    0.5465 |
#> |.....................|    0.6667 |...........|...........|...........|
#> |   58|     1284.7780 |     1.103 |   0.03780 |   -0.9638 |  -0.09111 |
#> |.....................|    -1.126 |    -3.616 |    -4.780 |    -2.589 |
#> |.....................|     1.648 |    0.5903 |   -0.3310 |    0.1803 |
#> |.....................| -0.007605 |...........|...........|...........|
#> |    U|               |     3.051 |     2.113 |    -3.708 |   -0.5466 |
#> |.....................|    -4.519 |    -5.454 |    -6.625 |    -3.034 |
#> |.....................|     2.852 |     1.771 |  0.009977 |     2.937 |
#> |.....................|     1.104 |...........|...........|...........|
#> |    X|               |     3.051 |     2.113 |    -3.708 |   -0.5466 |
#> |.....................|    -4.519 |    -5.454 |    -6.625 |    -3.034 |
#> |.....................|     2.852 |     1.771 |  0.009977 |     2.937 |
#> |.....................|     1.104 |...........|...........|...........|
#> |    M|     Mixed     |   -0.9487 |   0.06264 |   -0.3461 |   -0.5916 |
#> |.....................|   -0.3676 |  -0.02759 |   0.04549 |  -0.08078 |
#> |.....................|  -0.06953 |   -0.1426 |-2.070e+04 |    0.3853 |
#> |.....................|   -0.2110 |...........|...........|...........|
#> |   59|     1284.7754 |     1.105 |   0.03815 |   -0.9629 |  -0.08998 |
#> |.....................|    -1.125 |    -3.615 |    -4.784 |    -2.587 |
#> |.....................|     1.650 |    0.5917 |   -0.3310 |    0.1804 |
#> |.....................| -0.004075 |...........|...........|...........|
#> |    U|               |     3.057 |     2.114 |    -3.703 |   -0.5399 |
#> |.....................|    -4.515 |    -5.452 |    -6.630 |    -3.031 |
#> |.....................|     2.856 |     1.775 |  0.009975 |     2.938 |
#> |.....................|     1.107 |...........|...........|...........|
#> |    X|               |     3.057 |     2.114 |    -3.703 |   -0.5399 |
#> |.....................|    -4.515 |    -5.452 |    -6.630 |    -3.031 |
#> |.....................|     2.856 |     1.775 |  0.009975 |     2.938 |
#> |.....................|     1.107 |...........|...........|...........|
#> |    M|     Mixed     |   -0.1653 |    0.1191 |  -0.07774 |   -0.1611 |
#> |.....................|   -0.1188 |   0.02259 |  0.003799 |  -0.03973 |
#> |.....................|   0.01731 |   0.01126 |     7391. |    0.3122 |
#> |.....................|    0.1326 |...........|...........|...........|
#> |   60|     1284.7754 |     1.105 |   0.03815 |   -0.9629 |  -0.08998 |
#> |.....................|    -1.125 |    -3.615 |    -4.784 |    -2.587 |
#> |.....................|     1.650 |    0.5917 |   -0.3310 |    0.1804 |
#> |.....................| -0.004075 |...........|...........|...........|
#> |    U|               |     3.057 |     2.114 |    -3.703 |   -0.5399 |
#> |.....................|    -4.515 |    -5.452 |    -6.630 |    -3.031 |
#> |.....................|     2.856 |     1.775 |  0.009975 |     2.938 |
#> |.....................|     1.107 |...........|...........|...........|
#> |    X|               |     3.057 |     2.114 |    -3.703 |   -0.5399 |
#> |.....................|    -4.515 |    -5.452 |    -6.630 |    -3.031 |
#> |.....................|     2.856 |     1.775 |  0.009975 |     2.938 |
#> |.....................|     1.107 |...........|...........|...........|
#> calculating covariance matrix
#> [
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> 
#> covType="analytic" not available for this model (out of scope, or the augmented model would not build/solve); using the finite-difference sandwich ("r,s") covariance.
#> ====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> done
#> → Calculating residuals/tables
#> ✔ done

fitRE
```

``` math
\begin{align*}
{b} & = {b.timeF.0}+{b.timeF.0.5} {\times} \left({timeF}{\equiv}\text{"0.5"}\right)+{b.timeF.1} {\times} \left({timeF}{\equiv}\text{"1"}\right)+{b.timeF.1.5} {\times} \left({timeF}{\equiv}\text{"1.5"}\right)+{b.timeF.2} {\times} \left({timeF}{\equiv}\text{"2"}\right)+{b.timeF.3} {\times} \left({timeF}{\equiv}\text{"3"}\right)+{b.timeF.4} {\times} \left({timeF}{\equiv}\text{"4"}\right)+{b.timeF.6} {\times} \left({timeF}{\equiv}\text{"6"}\right)+{b.timeF.8} {\times} \left({timeF}{\equiv}\text{"8"}\right)+{b.timeF.12} {\times} \left({timeF}{\equiv}\text{"12"}\right)+{cov\_conc\_b} {\times} {conc} \\
{value} & = {b}+{bRe} \\
{value} & \sim add({addSd})
\end{align*}
```

The recovered `cov_conc_b` estimates the slope of concentration on
$`\Delta QT`$. Multiplying by the population $`C_{max}`$ (about 1000
ng/mL) gives the projected $`\Delta QT`$ at $`C_{max}`$ — typically
within ~1 msec of the simulated 10 msec.

## Defining the model equation (the formula)

The model formula has 3 parts: the dependent variable, the equation
predicting the dependent variable, and the optional random effects. A
full equation looks like the below:

    dv~predictors~randomEffects

Any model part may have any parameter name (for example, the dependent
variable does not have to be named `dv`).

### Dependent variable

The dependent variable may be any column in the `data`. It must only be
the column name and no additional modifiers. For example, `log(param)`
is not allowed.

### Predictors

Predictors may be any algebraic equation. In general, there are no
restrictions to the way that the equation may be written.

### Random effects

Random effects are not required for a model; fixed-effect-only models
are possible. When present, fixed and random effects must have different
names.

Random effects are written as a parameter a pipe (`|`) and its grouping
value. From the quick start example above, the parameter is `bRe` and
the grouping value is `id`, so it is written `bRe|id`.

For multiple parameters, they must be setup separately. For example, to
have `mRe` and `bRe` both grouped by `id`, it would be defined as
`~(mRe|id)+(bRe|id)`.

Only one grouping variable is supported per model. Multiple grouping
variables (for inter-occasion variability, for example) are not yet
supported.

### Per-parameter covariate models (`param`)

The `param` argument adds a covariate model to one or more parameters.
It accepts a formula or a list of formulas, each of the form
`<parameter> ~ <covariate columns>`. The right-hand side names one or
more columns from `data`:

- A **factor** column produces one fixed effect per level. The first
  level is the baseline; subsequent levels are estimated as differences
  from the baseline.
- A **numeric** column produces a slope. With only the continuous
  covariate, two parameters are introduced — `pop.<parameter>`
  (intercept) and `cov_<column>_<parameter>` (slope).
- Multiple covariates can be combined on the same parameter with `+`:
  `b ~ z + w` combines factor `z` and continuous `w`. When a factor is
  present, the factor supplies the baseline so `pop.b` is not
  introduced.

When mixing factor and continuous covariates on the same parameter,
`start[[<parameter>]]` must explicitly list one value per factor level
followed by one slope per continuous covariate (in the order they appear
in `param`).

#### `paramLink` (link function per parameter)

By default the assembled linear combination is the parameter itself
(`b <- <linear combination>`). Set `paramLink = c(b = "log")` to wrap
the combination in [`exp()`](https://rdrr.io/r/base/Log.html), so the
parameter is strictly positive on its natural scale.

## Providing the Starting estimates for the model

### Fixed effects

Fixed effects are defined by the `start` argument. The `start` argument
may either be a named vector or a named list. If fixed effects only have
a single starting value, then the two methods are equivalent.
`c(m=3, b=5)` is the same as `list(m=3, b=5)`. If one of the fixed
effects is defined by a factor variable (more on that later), the list
may have multiple starting values such as `list(m=3, b=c(5, 10))`.

For continuous covariates passed via `param`, `start[[<parameter>]]` may
be a single intercept value (the slope defaults to 0) or
`c(intercept, slope)`.

### Random effects

Random effect starting values are currently fixed at 1 without the
ability to modify it.

### Residual error models

The default residual error model is `~ add(addSd)` (additive standard
deviation `addSd`). More complex specifications are supported as long as
every sigma parameter in the residual formula appears in `start`. For
example:

``` r

nlmixr(
  y~m*x + b ~ (bRe|id),
  start = list(m=3, b=5, addSd=1, propSd=0.1),
  residualModel = ~ add(addSd) + prop(propSd),
  data = dSimNlme,
  est = "focei"
)
```
