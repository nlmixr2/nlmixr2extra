# Perform linearization of a model fitted using FOCEI

Perform linearization of a model fitted using FOCEI

## Usage

``` r
linearize(
  fit,
  mceta = c(-1, 10, 100, 1000),
  relTol = 0.25,
  focei = NA,
  addEtas = FALSE,
  derivFct = FALSE,
  plot = FALSE,
  est = "focei"
)
```

## Arguments

- fit:

  fit of nonlinear model fitted using any method with at least one eta.
  See details.

- mceta:

  a numeric vector for mceta to try. See details.

- relTol:

  relative deviation tolerance between original and linearized models
  objective functions. Used for switching if focei = NA. See details.

- focei:

  Default is NA for automatic switch from FOCEI to FOCE if failed. See
  details.

- addEtas:

  boolean. If TRUE, add etas on every theta and fix it to small value to
  get derivatives. Default is FALSE.

- derivFct:

  boolean. If TRUE, turn on derivatives for linearization. Default is
  FALSE.

- plot:

  boolean. Print plot of linearized vs original

- est:

  Character. The estimation method used for the linearized model fit.
  Only 'focei' is supported.

## Value

a fit object with subclass nlmixr2linearize

## Details

The function accepts a fit object fitted using any method. However, if
the method is not FOCE+I, the model is going to be evaluated first.

`mceta` vector will be iterated over to find the best linearization if
linearization failed. Escalating to next mceta will depend on the
relative deviation `relTol` of the original and linearized models
objective functions.

If `focei` is set to `NA`, the function will first try to linearize
using FOCEI. If the relative deviation between original and linearized
models objective functions is greater than `relTol`, it will switch to
FOCE where residual linearization is skipped. If `focei` is set to
`TRUE`, the function will use FOCEI linearization with individual and
residual linearization. If `focei` is set to `FALSE`, the function will
use FOCE linearization with residual interaction linearization skipped.

If `derivFct` is set to `TRUE`, the function will use an extra factor
for linearization that might help in stabilization. This might be useful
to try with FOCEI.

`plot` argument can only print ggplot figure with default settings. If a
user wish to capture the plot, one might use
[`linearizePlot()`](https://nlmixr2.github.io/nlmixr2extra/reference/linearizePlot.md)
call.

## Author

Omar I. Elashkar
