# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**nlmixr2extra** is the "extra support functions" package of the nlmixr2
ecosystem.  It does not fit models itself -- `nlmixr2est` does that -- it
wraps and post-processes fits: bootstrapping, stepwise/lasso/Bayesian
covariate selection, model linearization and IIV/residual-error search,
likelihood profiling, preconditioning, multistart, and reporting helpers.

Almost every entry point takes an `nlmixr2FitCore` fit (or an `rxUi` model)
and repeatedly rebuilds and re-estimates it.  Understanding `rxode2`'s `rxUi`
object and `nlmixr2est::nlmixr2()` is a prerequisite for working here.

## Build and Development Commands

### Install/Build
```r
devtools::install()   # or: R CMD INSTALL .
```

### Document
```r
devtools::document()
```

### Run All Tests
```sh
find src -name "*.so" -o -name "*.o" | xargs rm -f 2>/dev/null; NOT_CRAN=true Rscript -e "devtools::test()"
```

### Run a Single Test File
Filter by the test file name without the `test-` prefix and `.R` suffix:
```sh
NOT_CRAN=true Rscript -e "devtools::test(filter='multistart')"
NOT_CRAN=true Rscript -e "devtools::test(filter='bootstrap')"
NOT_CRAN=true Rscript -e "devtools::test(filter='profile')"
```

### R CMD Check
```r
invisible(lapply(list.files("src", "\\.s?o$", full.names = TRUE), unlink)); devtools::check()
```

Note `tests/testthat.R` deliberately does **not** call `test_check()`, so
`R CMD check` runs no tests.  Nearly every test file starts with
`skip_on_cran()` because the tests fit real models; run them with
`devtools::test()` and `NOT_CRAN=true`.

## Architecture

### The central pattern: rebuild the `rxUi`, refit, collect

Most subsystems follow the same loop -- take a fit, derive a modified `rxUi`,
call `nlmixr2est::nlmixr2()` on it, compare objective functions.  Two shared
helpers make this cheap:

- `R/parsingutil.R` is the **model-mutation engine**.  `buildupatedUI()`
  (note the typo in the exported name -- it is released API) adds or removes
  covariate terms on population parameters; `.expandPopExpr()` /
  `.expandRefMu()` rewrite a population expression given the fit's mu-reference
  and covariate data frames; `buildcovInfo()` turns `varsVec`/`covarsVec` into
  the per-parameter covariate list every search consumes;
  `addCatCovariates()` expands categorical covariates into dummy columns.
  SCM, lasso and the Bayesian selection all go through these.
- `setQuietFastControl()` (`R/utils.R`) is applied to a control object before
  every inner refit: `print = 0`, `covMethod = 0`, `calcTables = FALSE`,
  `compress = FALSE`.  Use it for any new refit loop.

### Subsystems

**Bootstrap** (`R/computingutil.R`, the largest file).  `bootstrapFit()` ->
`modelBootstrap()` -> `sampling()`.  `sampling()` resamples *subject ids*
(`.sampleUid()`), optionally stratified (`.stratSampleSize()`, `.stratProb()`);
`extractVars()` pulls parameters/omega/sigma out of the fit list and
`bootplot()` plots the result.  Fits are cached to a
`nlmixr2BootstrapCache_<fitName>_<md5>/` directory so an interrupted run
resumes; `fit$bootstrapMd5` keys the cache.  The file also holds generic
utilities used elsewhere (`foldgen()` cross-validation folds,
`normalizedData()`, `optimUnisampling()`).

**Stepwise covariate model selection** (`R/SCM.R`).  `covarSearchAuto()` runs
`forwardSearch()` then `backwardSearch()` over `varsVec` x `covarsVec`, caching
to `nlmixr2CovariateSearchCache_<fit>_<digest>/` with `restart=` to resume.

**Lasso covariate selection** (`R/lassocov.R`).  `lassoCoefficients()`,
`adaptivelassoCoefficients()`, `adjustedlassoCoefficients()` and
`regularmodel()`.  Each builds a penalized `rxUi` (`.lassoUicovariate()` /
`.adaptivelassoUicovariate()`), then picks the shrinkage `t` by
cross-validation (`.crossvalidationLasso()`, `.optimalTvaluelasso()`) using
`foldgen()` folds.

**Bayesian covariate selection** (`R/bayesiancovsel.R`).  `horseshoeSummardf()`
and friends; depends on the suggested `brms` package.

**Model linearization** (`R/linearizefocei.R`, `R/rxUiLinearize.R`).
`linearize()` is the pipeline entry: `getDeriv()` extracts FOCEi derivatives
from the fit, `linModGen()` generates the linearized `rxUi`, the result is
refit and gets class `nlmixr2Linearize` (escalating `mceta` until the objective
function agrees within `relTol`; `isLinearizeMatch()` / `linearizePlot()`
check the agreement).  `R/rxUiLinearize.R` supplies the error-model side: the
`linearizeErrorLines` S3 generic and the `.rxGetFctForError*()` variance
factors per residual-error type, registered as the `rxUiGet.linearizeError`
property in `R/zzz.R`.

A linearized fit is cheap to refit, so the searches dispatch on it:
`iivSearch.nlmixr2Linearize()` (`R/iivSearch.R`) enumerates eta combinations
(`iivCombn()`, `addAllEtas()`, `filterEtaMat()`) and `resSearch()`
(`R/resSearch.R`) tries residual-error models.  Both return objects whose
`rerunTopN()` method re-estimates the best candidates on the *original*
non-linear model.  `addCovariate()` also has an `nlmixr2Linearize` method, so
covariate testing can run on the linearized fit.

**Likelihood profiling** (`R/profile.R`).  `profile.nlmixr2FitCore()` is the
S3 entry, dispatching on `control` class: `fixedControl()` -> `profileFixed()`
(evaluate the OFV at fixed values) or `llpControl()` -> `profileLlp()`
(adaptive search for the OFV boundary via `optimProfile()`).  All the pieces
share the `@family Profiling` roxygen tag.

**Multistart** (`R/multistart.R`).  `multistart()` is an S3 generic over
fit/`rxUi`/`function`.  `.msPerturbIni()` builds perturbed `iniDf` starting
points (`"uniform"`, `"lhs"`, `"normal"`, respecting FIXED entries and
bounds), an optional cheap empirical-Bayes screen ranks them, and only the
promising ones are fully estimated.  Every start is cached under
`nlmixr2MultistartCache_<md5>/` as `start`/`screen`/`fit` RDS files
(`.msCacheRead()` / `.msCacheWrite()`), so a run resumes where it stopped --
the md5 covers anything that changes what a starting point *is*.  Returns a
`nlmixr2Multistart` object with `print()`, `plot()` (waterfall and parameter
stability) and `as.data.frame()` methods.

**Formula interface** (`R/nlmixrFormula.R`).  `nlmixrFormula()` assembles a
synthetic nlmixr2 model function from a `brms`-inspired formula.
`.nlmixrFormulaBuild()` drives parse -> data prep -> parameter expansion ->
`ini()`/`model()` assembly; covariate expansion (`.nlmixrFormulaExpandStartParam*()`)
emits `pop.<parameter>` intercepts and `cov_<param>_<start>` slopes.

**Preconditioning** (`R/precondition.R`, `src/preCondInv.cpp`).  The package's
only compiled code: an Rcpp/RcppArmadillo preconditioned matrix inverse used
by `preconditionFit()`.  `src/init.c` hand-maintains the `.Call` registration
table with a hardcoded arity, so adding or changing the arity of an
`[[Rcpp::export]]` function requires both `Rcpp::compileAttributes(".")` and a
hand edit of `src/init.c`.

**Equation printing** (`R/knit_printEquation.R`).  `knit_print` methods for
`rxUi`/`nlmixr2FitCore` that render a model as LaTeX in knitr documents, via
the `extractEqHelper` S3 walk over the model AST.

**Reporting helpers** (`R/AICHelpers.R`).  `getMinAICFit()`,
`listModelsTested()`, `isBoundaryFit()` -- both selection helpers exclude
boundary fits by default.

### Cache directories

Bootstrap, covariate search and multistart all persist intermediate fits to a
`nlmixr2*Cache_*` directory in the working directory, keyed by a digest of the
inputs, and resume from it.  These are listed in `.gitignore` and
`.Rbuildignore`; add any new cache prefix to both.

## R Code Style

Follow the same conventions as `rxode2`:

- **Exported functions**: `camelCase` (`bootstrapFit`, `covarSearchAuto`,
  `multistartControl`)
- **Internal/non-exported functions**: `.camelCase` with a leading dot
  (`.msPerturbIni`, `.expandPopExpr`, `.lassoUicovariate`).  Subsystem-local
  helpers use a short prefix after the dot (`.ms*` for multistart).
- **Local variables inside functions**: `.camelCase` with a leading dot.  Older
  files (`R/computingutil.R`, `R/lassocov.R`, `R/SCM.R`) predate this and use
  bare `snake_case`/`camelCase` locals; match the surrounding file when
  editing, use the dot convention in new code.
- **S3 methods**: `generic.class` (`multistart.rxUi`, `iivSearch.nlmixr2Linearize`).
  Register them in roxygen with `@export`, not by hand in `NAMESPACE`.
- Avoid `snake_case` for new names.
- American English spelling.
- **Never write `pkg:::foo` in package code or in a script.** CodeFactor flags
  every `:::` as a Major Maintainability issue and fails the PR check.  Bind
  the internal once with
  `.foo <- utils::getFromNamespace(".foo", "rxode2")` instead.  Inside
  `tests/testthat/` no qualifier is needed at all -- tests run in the package
  namespace, so call internals by their bare name.
- Do not rename already-exported functions, even obviously misspelled ones
  (`buildupatedUI`, `horseshoeSummardf`): they are released API and reverse
  dependencies rely on them.

## Documentation and Comment Style

- Keep comments and roxygen terse.  Condense multi-line explanations to
  one-liners; state the fact, not the story behind it.  Keep every `@param`,
  `@return`, `@export`, `@family` and `@examples` tag.
- Long-running examples go in `\dontrun{}`.
- `NEWS.md` is organized per version (`# nlmixr2extra X.Y.Z`), user-facing
  changes first under `## New features`, then `## Bug fixes`.  Entries are
  past-tense bullets of a sentence or two, referencing the issue number where
  one exists.
- **ASCII only.  No Unicode anywhere in the repo** (CRAN requirement, and a
  Unicode character in `R/` or `man/` currently fails the check): use `--` for
  em-dashes, `-` for en-dashes, `->` for arrows, straight quotes, `...` for
  ellipses, and spell out Greek letters (`Delta`, not the symbol).  This
  applies to plot labels, roxygen comments, tests, vignettes and `README.md`
  alike.  Check with:
  ```sh
  grep -rnP "[^\x00-\x7F]" R/ src/ tests/ man/ vignettes/ NEWS.md README.md
  ```
