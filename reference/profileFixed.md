# Estimate the objective function values for a model while fixing defined parameter values

Estimate the objective function values for a model while fixing defined
parameter values

## Usage

``` r
profileFixed(fitted, which, control = list())

profileFixedSingle(fitted, which)
```

## Arguments

- fitted:

  The fit model

- which:

  A data.frame with column names of parameters to fix and values of the
  fitted value to fix (one row only).

- control:

  A list passed to
  [`fixedControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/fixedControl.md)
  (currently unused)

## Value

`which` with a column named `OFV` added with the objective function
value of the fitted estimate fixing the parameters in the other columns

## Functions

- `profileFixedSingle()`: Estimate the objective function value for a
  model while fixing a single set of defined parameter values (for use
  in parameter profiling)

## See also

Other Profiling:
[`fixedControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/fixedControl.md),
[`llpControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/llpControl.md),
[`profile.nlmixr2FitCore()`](https://nlmixr2.github.io/nlmixr2extra/reference/profile.nlmixr2FitCore.md),
[`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md),
[`profileNlmixr2FitCoreRet()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileNlmixr2FitCoreRet.md)

Other Profiling:
[`fixedControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/fixedControl.md),
[`llpControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/llpControl.md),
[`profile.nlmixr2FitCore()`](https://nlmixr2.github.io/nlmixr2extra/reference/profile.nlmixr2FitCore.md),
[`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md),
[`profileNlmixr2FitCoreRet()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileNlmixr2FitCoreRet.md)

## Author

Bill Denney (changed by Matt Fidler to take out R 4.1 specific code)
