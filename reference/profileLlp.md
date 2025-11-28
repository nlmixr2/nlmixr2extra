# Profile confidence intervals with log-likelihood profiling

Profile confidence intervals with log-likelihood profiling

## Usage

``` r
profileLlp(fitted, which, control)
```

## Arguments

- fitted:

  The fit model

- which:

  Either `NULL` to profile all parameters or a character vector of
  parameters to estimate

- control:

  A list passed to
  [`llpControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/llpControl.md)

## Value

A data.frame with columns named "Parameter" (the parameter name(s) that
were fixed), OFV (the objective function value), and the current
estimate for each of the parameters. In addition, if any boundary is
found, the OFV increase will be indicated by the absolute value of the
"profileBound" column and if that boundary is the upper or lower
boundary will be indicated by the "profileBound" column being positive
or negative, respectively.

## See also

Other Profiling:
[`fixedControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/fixedControl.md),
[`llpControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/llpControl.md),
[`profile.nlmixr2FitCore()`](https://nlmixr2.github.io/nlmixr2extra/reference/profile.nlmixr2FitCore.md),
[`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md),
[`profileNlmixr2FitCoreRet()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileNlmixr2FitCoreRet.md)
