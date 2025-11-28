# Perform likelihood profiling on nlmixr2 focei fits

Perform likelihood profiling on nlmixr2 focei fits

## Usage

``` r
# S3 method for class 'nlmixr2FitCore'
profile(
  fitted,
  ...,
  which = NULL,
  method = c("llp", "fixed"),
  control = list()
)
```

## Arguments

- fitted:

  The fit model

- ...:

  ignored

- which:

  The parameter names to perform likelihood profiling on (`NULL`
  indicates all parameters)

- method:

  Method to use for profiling (see the details)

- control:

  Control arguments for the `method`

## Value

A data.frame with one column named `Parameter` indicating the parameter
being fixed on that row, one column for the `OFV` indicating the OFV
estimated for the model at that step, one column named `profileBound`
indicating the estimated value for the profile likelihood and its step
above the minimum profile likelihood value, and columns for each
parameter estimate (or fixed) in the model.

## Log-likelihood profiling

`method = "llp"`

The search will stop when either the OFV is within `ofvtol` of the
desired OFV change or when the parameter is interpolating to more
significant digits than specified in `paramDigits`. The "llp" method
uses the
[`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md)
function. See its help for more details.

## Fixed points

`method = "fixed"`

Estimate the OFV for specific fixed values. The "fixed" method uses the
[`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md)
function. See its help for more details.

## See also

Other Profiling:
[`fixedControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/fixedControl.md),
[`llpControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/llpControl.md),
[`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md),
[`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md),
[`profileNlmixr2FitCoreRet()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileNlmixr2FitCoreRet.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Likelihood profiling takes a long time to run each model multiple times, so
# be aware that running this example may take a few minutes.
oneCmt <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- fixed(log(31.5))
    eta.ka ~ 0.6
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

fit <-
  nlmixr2(
    oneCmt, data = nlmixr2data::theo_sd, est="focei", control = list(print=0)
  )
# profile all parameters
profall <- profile(fit)

# profile a single parameter
proftka <- profile(fit, which = "tka")
} # }
```
