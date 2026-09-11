# Plot a multistart result

Plot a multistart result

## Usage

``` r
# S3 method for class 'nlmixr2Multistart'
plot(x, type = c("waterfall", "parameters"), kBest = 20L, dOfvMax = NULL, ...)
```

## Arguments

- x:

  A `nlmixr2Multistart` object from
  [`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md)

- type:

  `"waterfall"` plots each start's objective function relative to the
  best one, worst to best; `"parameters"` plots how each parameter
  estimate varies across the best starts.

- kBest:

  Number of starts to show in the `"parameters"` plot

- dOfvMax:

  Upper limit for the waterfall's objective-function axis, for zooming
  in when one start is far worse than the rest

- ...:

  ignored

## Value

A ggplot2 object

## See also

Other Multistart:
[`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md),
[`multistartControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistartControl.md)
