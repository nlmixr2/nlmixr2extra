# Give the output data.frame for a single model for profile.nlmixr2FitCore

Give the output data.frame for a single model for profile.nlmixr2FitCore

## Usage

``` r
profileNlmixr2FitCoreRet(fitted, which, fixedVal)
```

## Arguments

- fitted:

  The fit model

- which:

  The parameter names to perform likelihood profiling on (`NULL`
  indicates all parameters)

- fixedVal:

  The value that `which` is fixed to in case the model does not
  converge.

## Value

A data.frame with columns named "Parameter" (the parameter name(s) that
were fixed), OFV (the objective function value), and the current
estimate for each of the parameters. Omega values are given as their
variances and covariances.

## See also

Other Profiling:
[`fixedControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/fixedControl.md),
[`llpControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/llpControl.md),
[`profile.nlmixr2FitCore()`](https://nlmixr2.github.io/nlmixr2extra/reference/profile.nlmixr2FitCore.md),
[`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md),
[`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md)
