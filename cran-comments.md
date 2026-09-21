# nlmixr2extra 5.2.1

This release adds `multistart()` and fixes bugs in `preconditionFit()`,
`linearize()` and `covarSearchAuto()`.

New features:

- `multistart()` re-estimates a model from many perturbed starting points
  to detect fits that settled in a local optimum

Bug fixes:

- `preconditionFit()` works again: it no longer fails on residual names
  containing a `.`, on models with random effects, or when 'nlmixr2est'
  reports a decorated covariance method (`"|r|,|s|"`, `"r,s (full)"`)

- `covarSearchAuto()` selects unit-scale covariates again; since
  'nlmixr2est' 7.0.2 a coefficient added at zero stayed pinned there

- `linearize()` no longer diverges when re-estimating a model with small
  residual error parameters, and works with correlated eta blocks

- The bundled `theoFitOde` fit was regenerated to match the current
  'nlmixr2est'

## Test environments

- Ubuntu 24.04, R 4.6.1 (local), `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

- `checking compilation flags used ... NOTE`: the non-portable flag
  `-mno-omit-leaf-frame-pointer` comes from the distribution r-base
  `Makeconf`, not from anything the package sets.

## Downstream dependencies

We checked all 4 reverse dependencies on CRAN (babelmixr2, nlmixr2,
nlmixr2plot, nlmixr2rpt) with 'revdepcheck', comparing R CMD check results
for the current CRAN version and this release.  No new problems were found.
No exported function was removed and no existing signature changed.
