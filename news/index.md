# Changelog

## nlmixr2extra 3.0.3

- Allow raw fits to be returned (or only the parameters)

## nlmixr2extra 3.0.2

CRAN release: 2025-02-17

- Make sure bootstrapped thetas are named. Fixed issue
  [\#76](https://github.com/nlmixr2/nlmixr2extra/issues/76)

## nlmixr2extra 3.0.1

CRAN release: 2024-10-29

- Remove non-functioning SCM for now
  ([\#71](https://github.com/nlmixr2/nlmixr2extra/issues/71))

## nlmixr2extra 3.0.0

CRAN release: 2024-09-18

- New [`profile()`](https://rdrr.io/r/stats/profile.html) method for
  likelihood profiling (Issue
  [\#1](https://github.com/nlmixr2/nlmixr2extra/issues/1))

## nlmixr2extra 2.0.10

CRAN release: 2024-05-29

- [`bootstrapFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/bootstrapFit.md)
  fixes `se` option (Issue
  [\#66](https://github.com/nlmixr2/nlmixr2extra/issues/66))

## nlmixr2extra 2.0.9

CRAN release: 2024-01-31

- [`bootstrapFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/bootstrapFit.md)
  now will be more careful handling `NA` values so they do not
  completely affect results (Issue
  [\#59](https://github.com/nlmixr2/nlmixr2extra/issues/59))

- [`bootstrapFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/bootstrapFit.md)
  will now only take the correlation of the non-zero diagonals (Issue
  [\#59](https://github.com/nlmixr2/nlmixr2extra/issues/59)).

- New method for
  [`knit_print()`](https://rdrr.io/pkg/knitr/man/knit_print.html) will
  generate model equations for LaTeX reporting automatically.

- Tests are now skipped if they contain linear compartment models that
  need gradients when the gradients are not compiled (as in the case of
  intel c++).

## nlmixr2extra 2.0.8

CRAN release: 2022-10-22

- Use
  [`assignInMyNamespace()`](https://rdrr.io/r/utils/getFromNamespace.html)
  instead of using the global assignment operator for the horseshoe
  prior

- Be specific in version requirements (as requested by CRAN checks)

- Move the `theoFitOde.rda` data build to
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  to reduce CRAN build time (could add more standard models like
  warfarin for package developers which takes way too much time for
  CRAN)

## nlmixr2extra 2.0.7

CRAN release: 2022-10-19

- Fix `cli` issues with the new `cli` 3.4+ release that will allow
  bootstrapping to run again (before `cli` would error, this fixes the
  `donttest` issues on CRAN).

- Fixed step-wise covariate selection to work a bit better with the
  updated UI, thanks to Vishal Sarsani

- Added lasso covariate selection (thanks to Vishal Sarsani)

- Added horseshoe prior covarite selecion (thanks to Vishal Sarsani)

- Added a `NEWS.md` file to track changes to the package.
