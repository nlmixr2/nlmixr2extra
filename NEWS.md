# nlmixr2extra 3.0.0

* New `profile()` method for likelihood profiling (Issue #1)

# nlmixr2extra 2.0.10

* `bootstrapFit()` fixes `se` option (Issue #66)

# nlmixr2extra 2.0.9

* `bootstrapFit()` now will be more careful handling `NA` values so
  they do not completely affect results (Issue #59)

* `bootstrapFit()` will now only take the correlation of the non-zero
  diagonals (Issue #59).

* New method for `knit_print()` will generate model equations for LaTeX
  reporting automatically.

* Tests are now skipped if they contain linear compartment models that
  need gradients when the gradients are not compiled (as in the case
  of intel c++).

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
