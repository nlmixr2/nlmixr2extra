# Control options for log-likelihood profiling

Control options for log-likelihood profiling

## Usage

``` r
llpControl(
  ofvIncrease = qchisq(0.95, df = 1),
  rseTheta = 30,
  itermax = 10,
  ofvtol = 0.005,
  paramDigits = 3,
  extrapolateExpand = 1.5
)
```

## Arguments

- ofvIncrease:

  The targetted change in objective function value (3.84 corresponds to
  a Chi-squared test with a 95% confidence interval)

- rseTheta:

  The relative standard error (percent) for the model parameters. It can
  be missing (the default) in which case a default value of 30% will be
  applied. If given as a single number, it will be applied to all
  parameters. If given as a named vector of numbers, it will be applied
  to each named parameter.

- itermax:

  Maximum number of likelihood profiling iterations for each bound
  estimated

- ofvtol:

  The relative tolerance for the objective function being close enough
  to the `ofvIncrease`.

- paramDigits:

  The number of significant digits required for the parameter. When
  interpolation attempts to get smaller than that number of significant
  digits, it will stop.

- extrapolateExpand:

  When extrapolating outside the range previously tested, how far should
  the step occur as a ratio

## Value

A validated list of control options for log-likelihood profiling

## See also

[`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md)

Other Profiling:
[`fixedControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/fixedControl.md),
[`profile.nlmixr2FitCore()`](https://nlmixr2.github.io/nlmixr2extra/reference/profile.nlmixr2FitCore.md),
[`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md),
[`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md),
[`profileNlmixr2FitCoreRet()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileNlmixr2FitCoreRet.md)
