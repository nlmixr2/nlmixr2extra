# nlmixr2extra development version

* `bootstrapFit()` now will be more careful handling `NA` values so
  they do not completely affect results (Issue #59)

# nlmixr2extra 2.0.9

* New method for `knit_print()` will generate model equations for LaTeX
  reporting automatically.

# nlmixr2extra 2.0.8

* Use `assignInMyNamespace()` instead of using the global assignment
  operator for the horseshoe prior

* Be specific in version requirements (as requested by CRAN checks)

* Move the `theoFitOde.rda` data build to `devtools::document()` to
  reduce CRAN build time (could add more standard models like warfarin
  for package developers which takes way too much time for CRAN)

# nlmixr2extra 2.0.7

* Fix `cli` issues with the new `cli` 3.4+ release that will allow
  bootstrapping to run again (before `cli` would error, this fixes the
  `donttest` issues on CRAN).

* Fixed step-wise covariate selection to work a bit better with the
  updated UI, thanks to Vishal Sarsani

* Added lasso covariate selection (thanks to Vishal Sarsani)

* Added horseshoe prior covarite selecion (thanks to Vishal Sarsani)

* Added a `NEWS.md` file to track changes to the package.
