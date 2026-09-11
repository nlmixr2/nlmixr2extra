# Stabilizing Covariance Estimates with preconditionFit()

## The problem: ill-conditioned covariance matrices

After fitting a nonlinear mixed-effects model with FOCEI, the standard
errors (SEs) of the population parameters are derived from the
variance-covariance matrix of the parameter estimates. This matrix is
assembled from the **R matrix** (second-derivative / Hessian
approximation) and, optionally, the **S matrix** (cross-product of first
derivatives). When the parameters are on very different scales, or when
two parameters are nearly collinear, the R matrix becomes
ill-conditioned: small rounding errors in the gradient calculations are
amplified during the matrix inversion and the resulting covariance
matrix is unreliable.

Symptoms of this problem include:

- `covMethod` reported as `"r"` or `"s"` rather than the preferred
  `"r,s"`
- NaN or negative values on the diagonal of the covariance matrix
- Standard errors that are implausibly large or small relative to the
  estimates
- `%RSE` values that are very large

## The solution: preconditioning

Aoki, Nordgren, and Hooker (2016) showed that a simple **linear
reparameterization** of the fixed-effects parameters can dramatically
improve the condition number of the R matrix without changing the model
or its fit. The idea is to find a matrix **P** such that **P R P^T** is
close to the identity matrix, run FOCEI on the reparameterized model
(which converges to a numerically stable covariance matrix), and then
transform the result back to the original parameterization.

`nlmixr2extra` implements this approach in
[`preconditionFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/preconditionFit.md).

> **Reference:** Aoki Y, Nordgren R, Hooker AC. Preconditioning of
> Nonlinear Mixed Effects Models for Stabilization of
> Variance-Covariance Matrix Computations. *AAPS J.* 2016;18(2):505–518.
> <doi:%5B10.1208/s12248-016-9866-5>\](<https://doi.org/10.1208/s12248-016-9866-5>)

## Basic usage

``` r

library(nlmixr2extra)

oneCompartment <- function() {
  ini({
    tka <- 0.45  # log Ka
    tcl <- 1     # log Cl
    tv  <- 3.45  # log V
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v  ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - cl/v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

fit := nlmixr2(oneCompartment, data = nlmixr2data::theo_sd,
               est = "focei", control = list(print = 0))
fit
```

``` math
\begin{align*}
{ka} & = \exp\left({tka}+{eta.ka}\right) \\
{cl} & = \exp\left({tcl}+{eta.cl}\right) \\
{v} & = \exp\left({tv}+{eta.v}\right) \\
\frac{d \: depot}{dt} & = -{ka} {\times} {depot} \\
\frac{d \: center}{dt} & = {ka} {\times} {depot}-\frac{{cl}}{{v}} {\times} {center} \\
{cp} & = \frac{{center}}{{v}} \\
{cp} & \sim add({add.sd})
\end{align*}
```

If the covariance method reported is not `"r,s"` (or you simply want a
more robust covariance estimate), apply preconditioning:

``` r

invisible(preconditionFit(fit))
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate sensitivities
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate ∂(f)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate ∂(R²)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in inner model...
#> → finding duplicate expressions in EBE model...
#> → compiling inner model...
#> ✔ done
#> → finding duplicate expressions in FD model...
#> → compiling EBE model...
#> ✔ done
#> → compiling events FD model...
#> ✔ done
#> rxode2 5.1.7 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
#> 
#> Attaching package: 'rxode2'
#> The following objects are masked from 'package:nlmixr2est':
#> 
#>     boxCox, yeoJohnson
#> Key: U: Unscaled Parameters; X: Back-transformed parameters; G: Gill difference gradient approximation
#> F: Forward difference gradient approximation
#> C: Central difference gradient approximation
#> M: Mixed forward and central difference gradient approximation
#> A: Analytic (forward sensitivity) gradient (fast=TRUE)
#> Unscaled parameters for Omegas=chol(solve(omega));
#> Diagonals are transformed, as specified by foceiControl(diagXform=)
#> 
#> |    #| Function Val. |nlmixr2Pre_tka |nlmixr2Pre_tcl |nlmixr2Pre_tv |nlmixr2Pre_add.sd |
#> |.....................|        o1 |        o2 |        o3 |...........|
#> |    1|     116.81105 |   -0.9752 |    -1.000 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    2|     116.81098 |   -0.9752 |    -1.000 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    3|     117.16602 |   -0.7752 |    -1.000 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     15.40 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     15.40 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    4|     116.81090 |   -0.9752 |   -0.8000 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.756 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.756 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    5|     355.05484 |   -0.9752 |    -1.000 |     1.200 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1889. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.786 |     1889. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    6|     118.06431 |   -0.9752 |    -1.000 |     1.000 |   -0.5133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1574. |     263.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.786 |     1574. |     263.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    7|     118.02696 |   -0.9752 |    -1.000 |     1.000 |   -0.7133 |
#> |.....................|   -0.7898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.417 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.417 |     1.945 |     2.666 |...........|
#> |    8|     116.96433 |   -0.9752 |    -1.000 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.7890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     2.048 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     2.048 |     2.666 |...........|
#> |    9|     116.82676 |   -0.9752 |    -1.000 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.7880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.741 |...........|
#> |    X|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.741 |...........|
#> |   10|     117.13147 |    -1.175 |    -1.000 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     10.27 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     10.27 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> 
#> |    #| Function Val. |nlmixr2Pre_tka |nlmixr2Pre_tcl |nlmixr2Pre_tv |nlmixr2Pre_add.sd |
#> |.....................|        o1 |        o2 |        o3 |...........|
#> |   11|     116.81517 |   -0.9752 |    -1.200 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.815 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.815 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |   12|     357.44856 |   -0.9752 |    -1.000 |    0.8000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1259. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.786 |     1259. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |   13|     118.45718 |   -0.9752 |    -1.000 |     1.000 |   -0.9133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1574. |     175.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.786 |     1574. |     175.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |   14|     117.97931 |   -0.9752 |    -1.000 |     1.000 |   -0.7133 |
#> |.....................|    -1.190 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.100 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.100 |     1.945 |     2.666 |...........|
#> |   15|     117.03261 |   -0.9752 |    -1.000 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |    -1.189 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.843 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.843 |     2.666 |...........|
#> |   16|     116.87995 |   -0.9752 |    -1.000 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |    -1.188 |...........|
#> |    U|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.591 |...........|
#> |    X|               |     12.83 |    -6.786 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.591 |...........|
#> |   17|     116.81773 |   -0.9803 |   -0.8066 |     1.000 |   -0.6997 |
#> |.....................|   -0.9918 |   -0.9708 |   -0.9248 |...........|
#> |    U|               |     12.77 |    -6.757 |     1575. |     222.9 |
#> |.....................|     1.257 |     1.955 |     2.690 |...........|
#> |    X|               |     12.77 |    -6.757 |     1575. |     222.9 |
#> |.....................|     1.257 |     1.955 |     2.690 |...........|
#> |   18|     116.81210 |   -0.9752 |   -0.7600 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.750 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.750 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |   19|     116.82828 |   -0.9769 |   -0.8042 |     1.000 |   -0.7009 |
#> |.....................|   -0.9921 |   -0.9775 |    -1.035 |...........|
#> |    U|               |     12.81 |    -6.757 |     1575. |     222.6 |
#> |.....................|     1.257 |     1.951 |     2.649 |...........|
#> |    X|               |     12.81 |    -6.757 |     1575. |     222.6 |
#> |.....................|     1.257 |     1.951 |     2.649 |...........|
#> |   20|     116.81430 |   -0.9929 |   -0.7823 |     1.000 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.61 |    -6.754 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.61 |    -6.754 |     1574. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> 
#> |    #| Function Val. |nlmixr2Pre_tka |nlmixr2Pre_tcl |nlmixr2Pre_tv |nlmixr2Pre_add.sd |
#> |.....................|        o1 |        o2 |        o3 |...........|
#> |   21|     116.82036 |   -0.9752 |   -0.8011 |     1.000 |   -0.7090 |
#> |.....................|   -0.9901 |    -1.012 |   -0.9795 |...........|
#> |    U|               |     12.83 |    -6.756 |     1575. |     220.8 |
#> |.....................|     1.258 |     1.934 |     2.669 |...........|
#> |    X|               |     12.83 |    -6.756 |     1575. |     220.8 |
#> |.....................|     1.258 |     1.934 |     2.669 |...........|
#> |   22|     117.46092 |   -0.9752 |   -0.8100 |     1.010 |   -0.7133 |
#> |.....................|   -0.9898 |   -0.9890 |   -0.9880 |...........|
#> |    U|               |     12.83 |    -6.758 |     1590. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |    X|               |     12.83 |    -6.758 |     1590. |     219.9 |
#> |.....................|     1.259 |     1.945 |     2.666 |...........|
#> |   23|     116.81080 |   -0.9679 |   -0.7927 |    0.9997 |   -0.7164 |
#> |.....................|   -0.9880 |   -0.9833 |   -0.9811 |...........|
#> |    U|               |     12.93 |    -6.755 |     1574. |     219.2 |
#> |.....................|     1.260 |     1.948 |     2.669 |...........|
#> |    X|               |     12.93 |    -6.755 |     1574. |     219.2 |
#> |.....................|     1.260 |     1.948 |     2.669 |...........|
#> |   24|     116.81740 |   -0.9675 |   -0.7826 |    0.9997 |   -0.7166 |
#> |.....................|   -0.9781 |   -0.9830 |   -0.9808 |...........|
#> |    U|               |     12.93 |    -6.754 |     1574. |     219.1 |
#> |.....................|     1.268 |     1.949 |     2.669 |...........|
#> |    X|               |     12.93 |    -6.754 |     1574. |     219.1 |
#> |.....................|     1.268 |     1.949 |     2.669 |...........|
#> |   25|     116.81679 |   -0.9684 |   -0.7920 |    0.9996 |   -0.7041 |
#> |.....................|   -0.9947 |   -0.9823 |   -0.9800 |...........|
#> |    U|               |     12.92 |    -6.755 |     1574. |     221.9 |
#> |.....................|     1.255 |     1.949 |     2.669 |...........|
#> |    X|               |     12.92 |    -6.755 |     1574. |     221.9 |
#> |.....................|     1.255 |     1.949 |     2.669 |...........|
#> |   26|     116.96163 |   -0.9704 |   -0.7921 |    0.9947 |   -0.7232 |
#> |.....................|   -0.9979 |   -0.9809 |   -0.9770 |...........|
#> |    U|               |     12.89 |    -6.755 |     1566. |     217.7 |
#> |.....................|     1.252 |     1.950 |     2.670 |...........|
#> |    X|               |     12.89 |    -6.755 |     1566. |     217.7 |
#> |.....................|     1.252 |     1.950 |     2.670 |...........|
#> |   27|     116.81448 |   -0.9654 |   -0.7946 |    0.9996 |   -0.7134 |
#> |.....................|   -0.9882 |   -0.9939 |   -0.9894 |...........|
#> |    U|               |     12.96 |    -6.755 |     1574. |     219.8 |
#> |.....................|     1.260 |     1.943 |     2.666 |...........|
#> |    X|               |     12.96 |    -6.755 |     1574. |     219.8 |
#> |.....................|     1.260 |     1.943 |     2.666 |...........|
#> |   28|     116.81736 |   -0.9682 |   -0.8029 |    0.9997 |   -0.7163 |
#> |.....................|   -0.9782 |   -0.9836 |   -0.9814 |...........|
#> |    U|               |     12.92 |    -6.757 |     1574. |     219.2 |
#> |.....................|     1.268 |     1.948 |     2.669 |...........|
#> |    X|               |     12.92 |    -6.757 |     1574. |     219.2 |
#> |.....................|     1.268 |     1.948 |     2.669 |...........|
#> |   29|     116.81080 |   -0.9732 |   -0.8025 |    0.9996 |   -0.7162 |
#> |.....................|   -0.9913 |   -0.9761 |   -0.9846 |...........|
#> |    U|               |     12.86 |    -6.756 |     1574. |     219.2 |
#> |.....................|     1.257 |     1.952 |     2.667 |...........|
#> |    X|               |     12.86 |    -6.756 |     1574. |     219.2 |
#> |.....................|     1.257 |     1.952 |     2.667 |...........|
#> |   30|     117.08502 |   -0.9632 |   -0.8068 |    0.9930 |   -0.7146 |
#> |.....................|   -0.9960 |   -0.9786 |   -0.9872 |...........|
#> |    U|               |     12.99 |    -6.757 |     1563. |     219.6 |
#> |.....................|     1.254 |     1.951 |     2.666 |...........|
#> |    X|               |     12.99 |    -6.757 |     1563. |     219.6 |
#> |.....................|     1.254 |     1.951 |     2.666 |...........|
#> 
#> |    #| Function Val. |nlmixr2Pre_tka |nlmixr2Pre_tcl |nlmixr2Pre_tv |nlmixr2Pre_add.sd |
#> |.....................|        o1 |        o2 |        o3 |...........|
#> |   31|     116.80912 |   -0.9779 |   -0.7971 |    0.9996 |   -0.7156 |
#> |.....................|   -0.9885 |   -0.9862 |   -0.9784 |...........|
#> |    U|               |     12.80 |    -6.756 |     1574. |     219.4 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |    X|               |     12.80 |    -6.756 |     1574. |     219.4 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |   32|     116.81233 |   -0.9772 |   -0.7855 |    0.9996 |   -0.7161 |
#> |.....................|   -0.9905 |   -0.9940 |   -0.9788 |...........|
#> |    U|               |     12.81 |    -6.754 |     1574. |     219.2 |
#> |.....................|     1.258 |     1.943 |     2.670 |...........|
#> |    X|               |     12.81 |    -6.754 |     1574. |     219.2 |
#> |.....................|     1.258 |     1.943 |     2.670 |...........|
#> |   33|     116.81518 |   -0.9777 |   -0.8072 |    0.9996 |   -0.7255 |
#> |.....................|   -0.9886 |   -0.9863 |   -0.9789 |...........|
#> |    U|               |     12.80 |    -6.757 |     1574. |     217.2 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |    X|               |     12.80 |    -6.757 |     1574. |     217.2 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |   34|     116.81063 |   -0.9799 |   -0.7929 |    0.9996 |   -0.7146 |
#> |.....................|   -0.9873 |   -0.9760 |   -0.9869 |...........|
#> |    U|               |     12.77 |    -6.755 |     1574. |     219.6 |
#> |.....................|     1.261 |     1.952 |     2.667 |...........|
#> |    X|               |     12.77 |    -6.755 |     1574. |     219.6 |
#> |.....................|     1.261 |     1.952 |     2.667 |...........|
#> |   35|     116.80935 |   -0.9793 |   -0.8012 |    0.9996 |   -0.7128 |
#> |.....................|   -0.9892 |   -0.9867 |   -0.9738 |...........|
#> |    U|               |     12.78 |    -6.756 |     1574. |     220.0 |
#> |.....................|     1.259 |     1.947 |     2.671 |...........|
#> |    X|               |     12.78 |    -6.756 |     1574. |     220.0 |
#> |.....................|     1.259 |     1.947 |     2.671 |...........|
#> |   36|     116.84816 |   -0.9774 |   -0.7960 |     1.002 |   -0.7170 |
#> |.....................|   -0.9880 |   -0.9852 |   -0.9770 |...........|
#> |    U|               |     12.80 |    -6.756 |     1577. |     219.1 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |    X|               |     12.80 |    -6.756 |     1577. |     219.1 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |   37|     116.81215 |   -0.9796 |   -0.7984 |    0.9996 |   -0.7135 |
#> |.....................|   -0.9884 |   -0.9871 |   -0.9800 |...........|
#> |    U|               |     12.78 |    -6.756 |     1573. |     219.8 |
#> |.....................|     1.260 |     1.946 |     2.669 |...........|
#> |    X|               |     12.78 |    -6.756 |     1573. |     219.8 |
#> |.....................|     1.260 |     1.946 |     2.669 |...........|
#> |   38|     116.81076 |   -0.9777 |   -0.7973 |    0.9997 |   -0.7152 |
#> |.....................|   -0.9886 |   -0.9879 |   -0.9785 |...........|
#> |    U|               |     12.80 |    -6.756 |     1574. |     219.4 |
#> |.....................|     1.260 |     1.946 |     2.670 |...........|
#> |    X|               |     12.80 |    -6.756 |     1574. |     219.4 |
#> |.....................|     1.260 |     1.946 |     2.670 |...........|
#> |   39|     116.81103 |   -0.9775 |   -0.7966 |    0.9994 |   -0.7159 |
#> |.....................|   -0.9886 |   -0.9848 |   -0.9777 |...........|
#> |    U|               |     12.80 |    -6.756 |     1573. |     219.3 |
#> |.....................|     1.259 |     1.948 |     2.670 |...........|
#> |    X|               |     12.80 |    -6.756 |     1573. |     219.3 |
#> |.....................|     1.259 |     1.948 |     2.670 |...........|
#> |   40|     116.81054 |   -0.9773 |   -0.7968 |    0.9996 |   -0.7150 |
#> |.....................|   -0.9888 |   -0.9860 |   -0.9785 |...........|
#> |    U|               |     12.81 |    -6.756 |     1574. |     219.5 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |    X|               |     12.81 |    -6.756 |     1574. |     219.5 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> 
#> |    #| Function Val. |nlmixr2Pre_tka |nlmixr2Pre_tcl |nlmixr2Pre_tv |nlmixr2Pre_add.sd |
#> |.....................|        o1 |        o2 |        o3 |...........|
#> |   41|     116.81145 |   -0.9782 |   -0.7971 |    0.9996 |   -0.7163 |
#> |.....................|   -0.9889 |   -0.9862 |   -0.9789 |...........|
#> |    U|               |     12.79 |    -6.756 |     1574. |     219.2 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |    X|               |     12.79 |    -6.756 |     1574. |     219.2 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |   42|     116.81172 |   -0.9773 |   -0.7974 |    0.9996 |   -0.7157 |
#> |.....................|   -0.9878 |   -0.9860 |   -0.9786 |...........|
#> |    U|               |     12.81 |    -6.756 |     1574. |     219.3 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |    X|               |     12.81 |    -6.756 |     1574. |     219.3 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |   43|     116.81179 |   -0.9784 |   -0.7965 |    0.9995 |   -0.7155 |
#> |.....................|   -0.9880 |   -0.9863 |   -0.9781 |...........|
#> |    U|               |     12.79 |    -6.756 |     1573. |     219.4 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |    X|               |     12.79 |    -6.756 |     1573. |     219.4 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |   44|     116.80955 |   -0.9780 |   -0.7968 |    0.9996 |   -0.7156 |
#> |.....................|   -0.9884 |   -0.9855 |   -0.9790 |...........|
#> |    U|               |     12.80 |    -6.756 |     1574. |     219.4 |
#> |.....................|     1.260 |     1.947 |     2.669 |...........|
#> |    X|               |     12.80 |    -6.756 |     1574. |     219.4 |
#> |.....................|     1.260 |     1.947 |     2.669 |...........|
#> |   45|     116.81139 |   -0.9783 |   -0.7978 |    0.9999 |   -0.7156 |
#> |.....................|   -0.9887 |   -0.9859 |   -0.9781 |...........|
#> |    U|               |     12.79 |    -6.756 |     1574. |     219.4 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |    X|               |     12.79 |    -6.756 |     1574. |     219.4 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |   46|     116.81065 |   -0.9778 |   -0.7962 |    0.9996 |   -0.7157 |
#> |.....................|   -0.9886 |   -0.9867 |   -0.9785 |...........|
#> |    U|               |     12.80 |    -6.756 |     1574. |     219.3 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |    X|               |     12.80 |    -6.756 |     1574. |     219.3 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |   47|     116.81183 |   -0.9776 |   -0.7973 |    0.9993 |   -0.7160 |
#> |.....................|   -0.9888 |   -0.9866 |   -0.9778 |...........|
#> |    U|               |     12.80 |    -6.756 |     1573. |     219.3 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |    X|               |     12.80 |    -6.756 |     1573. |     219.3 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |   48|     116.81070 |   -0.9775 |   -0.7974 |    0.9996 |   -0.7157 |
#> |.....................|   -0.9887 |   -0.9855 |   -0.9789 |...........|
#> |    U|               |     12.80 |    -6.756 |     1574. |     219.3 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |    X|               |     12.80 |    -6.756 |     1574. |     219.3 |
#> |.....................|     1.259 |     1.947 |     2.670 |...........|
#> |   49|     116.81135 |   -0.9772 |   -0.7965 |    0.9998 |   -0.7159 |
#> |.....................|   -0.9886 |   -0.9859 |   -0.9788 |...........|
#> |    U|               |     12.81 |    -6.756 |     1574. |     219.3 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |    X|               |     12.81 |    -6.756 |     1574. |     219.3 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |   50|     116.81091 |   -0.9779 |   -0.7971 |    0.9996 |   -0.7156 |
#> |.....................|   -0.9885 |   -0.9862 |   -0.9784 |...........|
#> |    U|               |     12.80 |    -6.756 |     1574. |     219.4 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> |    X|               |     12.80 |    -6.756 |     1574. |     219.4 |
#> |.....................|     1.260 |     1.947 |     2.670 |...........|
#> calculating covariance matrix
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> done
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate sensitivities
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate ∂(f)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in inner model...
#> → finding duplicate expressions in EBE model...
#> → compiling inner model...
#> ✔ done
#> → finding duplicate expressions in FD model...
#> → compiling EBE model...
#> ✔ done
#> → compiling events FD model...
#> ✔ done
#> Updated original fit object fit
fit$covMethod
#> [1] "precondition"
```

[`preconditionFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/preconditionFit.md)
modifies the **original fit object in-place**: after the call,
`fit$cov`, `fit$parFixedDf`, and `fit$covMethod` all reflect the
preconditioned covariance. The function returns the preconditioned
covariance matrix invisibly.

## Controlling the amount of re-estimation

The `estType` argument determines how much of the model is re-estimated
on the reparameterized problem:

| `estType` | What is re-estimated | When to use |
|----|----|----|
| `"full"` (default) | Outer *and* inner iterations (theta, ETA) | When you suspect the original estimates may have converged to a local minimum |
| `"posthoc"` | Inner iterations only (ETA, given theta) | When you trust the theta estimates but want a better covariance |
| `"none"` | Nothing – only the covariance matrix is re-computed | Fastest; good for a quick diagnostic check |

``` r

# Shown for reference -- each call re-estimates and modifies `fit` in place, so
# run whichever one you want against a freshly fitted model.
# Only recompute the covariance matrix (fastest)
preconditionFit(fit, estType = "none")

# Fix theta, re-estimate ETAs, then compute covariance
preconditionFit(fit, estType = "posthoc")

# Full re-estimation on the reparameterized model (default)
preconditionFit(fit, estType = "full")
```

## Retrying until convergence: `ntry`

The inner loop of
[`preconditionFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/preconditionFit.md)
iterates the preconditioning until the reparameterized problem achieves
a `"r,s"` covariance (the most reliable type). If that target is not
reached within `ntry` attempts the function stops with an error. The
default is `ntry = 10`.

``` r

# Allow up to 20 attempts before giving up
preconditionFit(fit, ntry = 20L)
```

If
[`preconditionFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/preconditionFit.md)
fails even with more tries, the R matrix itself may be too poorly
determined for preconditioning to help. In that case inspect the model
for identifiability problems or consider likelihood profiling
(`profile(fit)`) instead of Wald-based confidence intervals.

## Switching between covariance estimates

[`preconditionFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/preconditionFit.md)
stores the preconditioned covariance alongside any previously computed
covariances in `fit$covList`. You can inspect what is available and
switch between them with
[`setCov()`](https://nlmixr2.github.io/nlmixr2est/reference/setCov.html):

``` r

# See which covariance estimates are stored
names(fit$covList)
#> [1] "r,s"
```

``` r

# Switch back to the standard r,s covariance
setCov(fit, "r,s")

# Switch to the preconditioned covariance
setCov(fit, "precondition")
```

After
[`setCov()`](https://nlmixr2.github.io/nlmixr2est/reference/setCov.html)
the fit object is updated in-place and the displayed parameter table
(`fit$parFixedDf`) reflects the selected covariance.

## Worked example: comparing SEs before and after preconditioning

``` r

library(nlmixr2extra)

# a second, untouched fit: `fit` above was already preconditioned in place
fitCompare := nlmixr2(oneCompartment, data = nlmixr2data::theo_sd,
                      est = "focei", control = list(print = 0))

# Record the original parameter table
dfBefore <- fitCompare$parFixedDf
cat("Covariance method before:", fitCompare$covMethod, "\n")
#> Covariance method before: r,s

# Apply preconditioning (only recompute covariance, do not re-estimate)
invisible(preconditionFit(fitCompare, estType = "none"))
#> calculating covariance matrix
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
dfAfter <- fitCompare$parFixedDf
cat("Covariance method after:", fitCompare$covMethod, "\n")
#> Covariance method after: precondition

# Compare standard errors
cbind(
  SE_before = dfBefore$SE,
  SE_after  = dfAfter$SE,
  row.names = rownames(dfBefore)
)
#>        SE_before            SE_after              row.names
#> tka    "0.23191882512431"   "0.0216417252282623"  "tka"    
#> tcl    "0.169047370259566"  "0.009228606428201"   "tcl"    
#> tv     "0.0405121558563634" "0.00923909311815786" "tv"     
#> add.sd "0.0730247302042821" "0.0820970850815736"  "add.sd"
```

The point estimates (`Estimate`, `Back-transformed`) and the
random-effects summaries (`BSV(CV%)`, `Shrink(SD)%`) are identical
before and after – only the covariance-derived quantities (SE, %RSE,
confidence intervals) change.

## When to use `preconditionFit()`

Use
[`preconditionFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/preconditionFit.md)
when:

- `fit$covMethod` is not `"r,s"` after a FOCEI fit
- Standard errors look implausibly large or contain NaN
- You want a publication-quality covariance estimate and are willing to
  spend extra computation time

It is not needed when FOCEI already produces a `"r,s"` covariance and
the SEs look reasonable. For models where identifiability is in doubt
(very large %RSE for multiple parameters simultaneously), consider
likelihood profiling rather than trying to stabilize a structurally
unreliable covariance.

## Requirements

[`preconditionFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/preconditionFit.md)
requires:

- A fit produced by FOCEI (or a method that stores an R matrix in
  `fit$R`). SAEM fits do not store an R matrix by default; call
  `fit <- nlmixr2(model, data, est = "focei", ...)` first, or use
  `getVarCov(saemFit)` to trigger covariance computation.
- The `nlmixr2extra` package to be loaded.
