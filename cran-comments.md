# nlmixr2extra 5.2.0

This release adds reporting helpers for model comparison and fixes
several bugs in the bootstrap and covariate-search code.  It supersedes
the 5.1.1 development version, which was never submitted; its fixes are
included here.

New features:

- Reporting helpers for comparing candidate models: `getMinAICFit()`,
  `listModelsTested()` and `isBoundaryFit()`

Bug fixes:

- `bootstrapFit(stratVar=)` now actually resamples, keeps subject ids
  distinct across strata, and draws whole subjects

- `covarSearchAuto()` no longer errors when a covariate is selected, and
  the forward inclusion test no longer has an inverted sign

- `bootstrapFit()` works for models with a single population parameter,
  a single random effect, or no random effects

- `optimUnisampling()` keeps `N` and `floorT` when it retries internally

- The bundled `theoFitOde` fit was regenerated so it no longer requires
  the `qs2` package to be decoded

## Test environments

- Ubuntu 24.04, R 4.6.1 (local), `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 2 notes

Both notes are properties of the local check environment rather than the
package:

- `checking compilation flags used ... NOTE`: the non-portable flag
  `-mno-omit-leaf-frame-pointer` comes from the distribution r-base
  `Makeconf`, not from anything the package sets.

- `checking HTML version of manual ... NOTE`: HTML validation was
  skipped because no `tidy` command is installed locally.

## Downstream dependencies

The reverse dependencies on CRAN are babelmixr2, nlmixr2, nlmixr2plot
and nlmixr2rpt.  This release only adds exports (`getMinAICFit()`,
`listModelsTested()`, `isBoundaryFit()`); no exported function was
removed and no existing signature changed.
