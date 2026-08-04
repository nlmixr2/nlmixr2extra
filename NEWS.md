# nlmixr2extra 5.2.0

## New features

- New reporting helpers for comparing candidate models: `getMinAICFit()`
  returns the fit with the lowest AIC, `listModelsTested()` builds a
  `Description`/`AIC`/`dAIC` table ready for `pander::pander()`, and
  `isBoundaryFit()` reports whether a fit has a parameter at its
  boundary.  By default both selection helpers exclude boundary fits.
  See the new "reporting helpers" article.

## Bug fixes

- Fix `bootstrapFit(stratVar=)`, which did not actually resample.  The
  stratified branch called `sample(list(uids), ...)`, and since
  `list(uids)` has length one every draw returned the whole vector of
  subject ids, so the bootstrap datasets did not depend on the seed
  (#99).  Three further problems in the same code are fixed with it:
  the new subject ids restarted at 1 in each stratum, so subjects from
  different strata were merged under a shared id; the sample was split
  across strata by the number of *observations* rather than the number
  of *subjects*, over-weighting strata whose subjects have more
  records; and rounding each stratum up could return more subjects than
  `nSampIndiv` asked for.

- A stratified bootstrap now always draws whole subjects.  When
  `stratVar` changed within a subject, that subject's records were
  split between strata and resampled as separate (partial) subjects;
  each subject is now stratified by its first value, with a warning.

- `nlmixr2extra:::sampling()` now resolves its `uid_colname` default
  before using it.  Called without one it sampled `ncol(data)` subjects
  instead of the number of subjects in the data.  It also accepts a
  tibble, which previously produced a one column tibble where a vector
  of subject ids was expected.

- Fix `covarSearchAuto()` crashing with "wrong arguments for subsetting
  an environment" when a covariate is selected; the best model is now
  re-fit to recover its fit object.  Also corrected the forward
  inclusion test, which had an inverted sign so improving covariates
  were never selected (#103)

- `bootstrapFit()` now works for models with a single estimated
  population parameter, a single random effect, or no random effects at
  all.  Previously the bootstrap summary collapsed 1-row / 1x1 quantile
  arrays to vectors (and could not summarize a `NULL` omega), causing
  `bootstrapFit()` to error with `dim(X) must have a positive length`,
  `incorrect number of dimensions`, or `'data' must be of a vector type,
  was 'NULL'`.  Printing the bootstrap summary of a model with no random
  effects no longer errors either.

- `optimUnisampling()` now keeps `N` and `floorT` when it retries
  internally.  Before, the recursive call reset them to the defaults, so
  asking for a sample size other than 1000, or for un-floored values,
  could silently return 1000 integer samples instead (#97)

- The bundled `theoFitOde` fit was regenerated and can now be read
  without the `qs2` package.  Its `origData` and `parHistData` had been
  serialized with `qs2`, so without that package installed those slots
  could not be decoded and `fit$dataMergeInner()` -- and anything built
  on it, such as the `nlmixr2rpt` figures -- failed.

# nlmixr2extra 5.1.0

- Add focei/foce linearization

- Add formula interface

- Add vignettes on linearization, formula interface, log-likelihood
  profiling and preconditioning.

# nlmixr2extra 5.0.0

* Update internal fit to be nlmixr2est 5.0.0 fit object

# nlmixr2extra 3.0.3

* Allow raw fits to be returned (or only the parameters)

# nlmixr2extra 3.0.2

* Make sure bootstrapped thetas are named.  Fixed issue #76

# nlmixr2extra 3.0.1

* Remove non-functioning SCM for now (#71)

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
