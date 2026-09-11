# Control options for `multistart()`

Control options for
[`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md)

## Usage

``` r
multistartControl(
  n = 10L,
  nFit = NULL,
  sampling = c("uniform", "lhs", "normal"),
  spread = 0.2,
  which = NULL,
  perturbOmega = TRUE,
  omegaFold = 2,
  around = c("final", "initial"),
  screen = c("posthoc", "none"),
  refitBest = TRUE,
  excludeBoundary = TRUE,
  keepFits = TRUE,
  seed = 1234L,
  cores = 1L,
  cacheDir = NULL,
  restart = FALSE
)
```

## Arguments

- n:

  Number of candidate starting points to generate. The first candidate
  is always the unperturbed starting point, so `n = 1` reproduces the
  original fit.

- nFit:

  Number of candidates to fully estimate after screening. `NULL` (the
  default) estimates every candidate. Ignored when `screen = "none"`.

- sampling:

  How the starting points are drawn around the initial estimates:
  `"uniform"` (the default) draws uniformly within `spread`, `"lhs"`
  uses a Latin hypercube over the same interval so the range is covered
  more evenly, and `"normal"` draws normally with a standard deviation
  of `spread`.

- spread:

  Fractional spread of the perturbation on the estimation scale. A
  parameter with initial estimate `est` is perturbed within
  `est +/- spread*max(abs(est), 1)`.

- which:

  Names of the population parameters to perturb; `NULL` (the default)
  perturbs every unfixed population parameter.

- perturbOmega:

  Should the between-subject variability estimates be perturbed as well?

- omegaFold:

  Fold-range for the `omega` perturbation. A variance `v` is drawn
  within `v/omegaFold` and `v*omegaFold`.

- around:

  When starting from a fit, perturb around the fit's `"final"` estimates
  (the default, which asks "is this a local optimum?") or around the
  `"initial"` estimates the fit started from (which asks "how sensitive
  was this fit to where I started?").

- screen:

  Cheap pre-selection of candidates. `"posthoc"` (the default) evaluates
  the objective function at each candidate with an empirical Bayes step
  only and fully estimates the best `nFit`; `"none"` estimates every
  candidate.

- refitBest:

  Re-run the best start with the full estimation control, so the
  returned fit has the covariance step and tables the exploratory runs
  skip.

- excludeBoundary:

  Should fits with a parameter at a boundary be excluded when picking
  the best start? Matches
  [`getMinAICFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/getMinAICFit.md).

- keepFits:

  Keep each start's fit in the result. Setting this to `FALSE` keeps
  only the summary and the best fit, which is much smaller.

- seed:

  Integer seed. The perturbations are drawn from this seed and each
  start is estimated with its own derived seed.

- cores:

  Number of starts to estimate at once. Each estimation already uses
  every available thread internally, so the default of `1` is usually
  the fastest choice; see the "Parallel estimation" section of
  [`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md).

- cacheDir:

  Directory used to cache the individual starts so that an interrupted
  run can be resumed. `NULL` (the default) derives a name from the
  model; `NA` disables caching.

- restart:

  Discard any cached results and start over.

## Value

A validated list of control options for
[`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md)

## See also

[`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md)

Other Multistart:
[`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md),
[`plot.nlmixr2Multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/plot.nlmixr2Multistart.md)
