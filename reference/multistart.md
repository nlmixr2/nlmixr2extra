# Estimate a model from many starting points

Refits a model from several perturbed sets of initial estimates and
collects the results, so that a fit which settled in a local optimum can
be recognised. Use
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the result
for the objective-function waterfall and the parameter-stability plots.

## Usage

``` r
multistart(object, ...)

# S3 method for class 'nlmixr2FitCore'
multistart(
  object,
  ...,
  data = NULL,
  est = NULL,
  estControl = NULL,
  control = list()
)

# S3 method for class 'rxUi'
multistart(
  object,
  data,
  ...,
  est = "focei",
  estControl = NULL,
  control = list()
)

# S3 method for class '`function`'
multistart(object, data, ...)

# Default S3 method
multistart(object, ...)
```

## Arguments

- object:

  A nlmixr2 fit, a nlmixr2 model function, or a `rxUi` model object

- ...:

  ignored

- data:

  The data to estimate with; taken from `object` when it is a fit

- est:

  The estimation method; taken from `object` when it is a fit

- estControl:

  The control for `est`; taken from `object` when it is a fit, and
  otherwise the method's default

- control:

  A list passed to
  [`multistartControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistartControl.md)

## Value

An object of class `nlmixr2Multistart`, a list with elements `starts`
(one row per candidate starting point), `summary` (one row per estimated
start), `fits`, `best`, and `bestIndex`

## Starting points

The first candidate is always the unperturbed starting point, so the
original fit is always represented in the comparison. Every other
candidate perturbs the unfixed population parameters on the scale the
model is estimated on (which for a mu-referenced parameter is usually
the log scale), and clips the result to the parameter's declared bounds.
Between-subject variances are perturbed multiplicatively so they stay
positive, and the resulting matrix is made positive-definite with
[`lotri::lotriNearPD()`](https://nlmixr2.github.io/lotri/reference/lotriNearPD.html).

## Screening

Fully estimating every candidate is wasteful when many of them start far
from anywhere sensible. With the default `screen = "posthoc"` each
candidate is first evaluated with an empirical-Bayes step only, which
costs a small fraction of a full estimation, and only the best `nFit`
candidates are then fully estimated.

## Parallel estimation

Each estimation already runs across every available thread, so
estimating several starts at once oversubscribes the machine unless the
thread budget is divided. `multistartControl(cores=)` therefore
restricts each worker to a single thread. Whether that is faster than
the serial default depends entirely on the model; a model whose subjects
parallelise well is usually better off left serial. Parallel estimation
uses [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
and is not available on Windows.

## Resuming

Each start is cached to `cacheDir` as it completes, so an interrupted
run resumes where it left off. Increasing `n` on a later call re-uses
the starts already estimated and only estimates the new ones. Pass
`restart = TRUE` to discard the cache, or `cacheDir = NA` to never write
one.

Changing anything that alters what a starting point *is* (`sampling`,
`spread`, `which`, `perturbOmega`, `omegaFold` or `seed`) gives the run
its own cache, so a cached estimation is never re-used for a start it
did not come from. A Latin hypercube is a design over all `n` candidates
at once, so growing an `"lhs"` run moves its earlier starting points;
those starts are detected and re-estimated rather than reported against
the wrong starting point.

## See also

Other Multistart:
[`multistartControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistartControl.md),
[`plot.nlmixr2Multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/plot.nlmixr2Multistart.md)

## Author

Matthew Fidler

## Examples

``` r
if (FALSE) { # \dontrun{
# Every start is a full estimation, so this takes a few minutes.
fit <- nlmixr2extra::theoFitOde

ms <- multistart(fit, control = list(n = 8, spread = 0.3))
ms

# objective function values, best to worst
plot(ms)

# how stable each parameter is across the best starts
plot(ms, "parameters")

# the best fit found, ready to use like any other fit
ms$best
} # }
```
