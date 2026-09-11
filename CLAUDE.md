# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Overview

**nlmixr2extra** is the “extra support functions” package of the nlmixr2
ecosystem. It does not fit models itself – `nlmixr2est` does that – it
wraps and post-processes fits: bootstrapping, stepwise/lasso/Bayesian
covariate selection, model linearization and IIV/residual-error search,
likelihood profiling, preconditioning, multistart, and reporting
helpers.

Almost every entry point takes an `nlmixr2FitCore` fit (or an `rxUi`
model) and repeatedly rebuilds and re-estimates it. Understanding
`rxode2`’s `rxUi` object and
[`nlmixr2est::nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)
is a prerequisite for working here.

## Build and Development Commands

### Install/Build

``` r

devtools::install()   # or: R CMD INSTALL .
```

### Document

``` r

devtools::document()
```

### Run All Tests

``` sh
find src -name "*.so" -o -name "*.o" | xargs rm -f 2>/dev/null; NOT_CRAN=true Rscript -e "devtools::test()"
```

### Run a Single Test File

Filter by the test file name without the `test-` prefix and `.R` suffix:

``` sh
NOT_CRAN=true Rscript -e "devtools::test(filter='multistart')"
NOT_CRAN=true Rscript -e "devtools::test(filter='bootstrap')"
NOT_CRAN=true Rscript -e "devtools::test(filter='profile')"
```

### R CMD Check

``` r

invisible(lapply(list.files("src", "\\.s?o$", full.names = TRUE), unlink)); devtools::check()
```

Note `tests/testthat.R` deliberately does **not** call `test_check()`,
so `R CMD check` runs no tests. Nearly every test file starts with
`skip_on_cran()` because the tests fit real models; run them with
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html) and
`NOT_CRAN=true`.

## Architecture

### The central pattern: rebuild the `rxUi`, refit, collect

Most subsystems follow the same loop – take a fit, derive a modified
`rxUi`, call
[`nlmixr2est::nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)
on it, compare objective functions. Two shared helpers make this cheap:

- `R/parsingutil.R` is the **model-mutation engine**.
  [`buildupatedUI()`](https://nlmixr2.github.io/nlmixr2extra/reference/buildupatedUI.md)
  (note the typo in the exported name – it is released API) adds or
  removes covariate terms on population parameters; `.expandPopExpr()` /
  `.expandRefMu()` rewrite a population expression given the fit’s
  mu-reference and covariate data frames;
  [`buildcovInfo()`](https://nlmixr2.github.io/nlmixr2extra/reference/buildcovInfo.md)
  turns `varsVec`/`covarsVec` into the per-parameter covariate list
  every search consumes;
  [`addCatCovariates()`](https://nlmixr2.github.io/nlmixr2extra/reference/addCatCovariates.md)
  expands categorical covariates into dummy columns. SCM, lasso and the
  Bayesian selection all go through these.
- `setQuietFastControl()` (`R/utils.R`) is applied to a control object
  before every inner refit: `print = 0`, `covMethod = 0`,
  `calcTables = FALSE`, `compress = FALSE`. Use it for any new refit
  loop.

### Subsystems

**Bootstrap** (`R/computingutil.R`, the largest file).
[`bootstrapFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/bootstrapFit.md)
-\> `modelBootstrap()` -\> `sampling()`. `sampling()` resamples *subject
ids* (`.sampleUid()`), optionally stratified (`.stratSampleSize()`,
`.stratProb()`); `extractVars()` pulls parameters/omega/sigma out of the
fit list and
[`bootplot()`](https://nlmixr2.github.io/nlmixr2extra/reference/bootplot.md)
plots the result. Fits are cached to a
`nlmixr2BootstrapCache_<fitName>_<md5>/` directory so an interrupted run
resumes; `fit$bootstrapMd5` keys the cache. The file also holds generic
utilities used elsewhere
([`foldgen()`](https://nlmixr2.github.io/nlmixr2extra/reference/foldgen.md)
cross-validation folds,
[`normalizedData()`](https://nlmixr2.github.io/nlmixr2extra/reference/normalizedData.md),
[`optimUnisampling()`](https://nlmixr2.github.io/nlmixr2extra/reference/optimUnisampling.md)).

**Stepwise covariate model selection** (`R/SCM.R`).
[`covarSearchAuto()`](https://nlmixr2.github.io/nlmixr2extra/reference/covarSearchAuto.md)
runs `forwardSearch()` then `backwardSearch()` over `varsVec` x
`covarsVec`, caching to `nlmixr2CovariateSearchCache_<fit>_<digest>/`
with `restart=` to resume.

**Lasso covariate selection** (`R/lassocov.R`).
[`lassoCoefficients()`](https://nlmixr2.github.io/nlmixr2extra/reference/lassoCoefficients.md),
[`adaptivelassoCoefficients()`](https://nlmixr2.github.io/nlmixr2extra/reference/adaptivelassoCoefficients.md),
[`adjustedlassoCoefficients()`](https://nlmixr2.github.io/nlmixr2extra/reference/adjustedlassoCoefficients.md)
and
[`regularmodel()`](https://nlmixr2.github.io/nlmixr2extra/reference/regularmodel.md).
Each builds a penalized `rxUi` (`.lassoUicovariate()` /
`.adaptivelassoUicovariate()`), then picks the shrinkage `t` by
cross-validation (`.crossvalidationLasso()`, `.optimalTvaluelasso()`)
using
[`foldgen()`](https://nlmixr2.github.io/nlmixr2extra/reference/foldgen.md)
folds.

**Bayesian covariate selection** (`R/bayesiancovsel.R`).
[`horseshoeSummardf()`](https://nlmixr2.github.io/nlmixr2extra/reference/horseshoeSummardf.md)
and friends; depends on the suggested `brms` package.

**Model linearization** (`R/linearizefocei.R`, `R/rxUiLinearize.R`).
[`linearize()`](https://nlmixr2.github.io/nlmixr2extra/reference/linearize.md)
is the pipeline entry: `getDeriv()` extracts FOCEi derivatives from the
fit,
[`linModGen()`](https://nlmixr2.github.io/nlmixr2extra/reference/linModGen.md)
generates the linearized `rxUi`, the result is refit and gets class
`nlmixr2Linearize` (escalating `mceta` until the objective function
agrees within `relTol`;
[`isLinearizeMatch()`](https://nlmixr2.github.io/nlmixr2extra/reference/isLinearizeMatch.md)
/
[`linearizePlot()`](https://nlmixr2.github.io/nlmixr2extra/reference/linearizePlot.md)
check the agreement). `R/rxUiLinearize.R` supplies the error-model side:
the `linearizeErrorLines` S3 generic and the `.rxGetFctForError*()`
variance factors per residual-error type, registered as the
`rxUiGet.linearizeError` property in `R/zzz.R`.

A linearized fit is cheap to refit, so the searches dispatch on it:
[`iivSearch.nlmixr2Linearize()`](https://nlmixr2.github.io/nlmixr2extra/reference/iivSearch.md)
(`R/iivSearch.R`) enumerates eta combinations (`iivCombn()`,
[`addAllEtas()`](https://nlmixr2.github.io/nlmixr2extra/reference/addAllEtas.md),
`filterEtaMat()`) and
[`resSearch()`](https://nlmixr2.github.io/nlmixr2extra/reference/resSearch.md)
(`R/resSearch.R`) tries residual-error models. Both return objects whose
[`rerunTopN()`](https://nlmixr2.github.io/nlmixr2extra/reference/rerunTopN.md)
method re-estimates the best candidates on the *original* non-linear
model.
[`addCovariate()`](https://nlmixr2.github.io/nlmixr2extra/reference/addCovariate.md)
also has an `nlmixr2Linearize` method, so covariate testing can run on
the linearized fit.

**Likelihood profiling** (`R/profile.R`).
[`profile.nlmixr2FitCore()`](https://nlmixr2.github.io/nlmixr2extra/reference/profile.nlmixr2FitCore.md)
is the S3 entry, dispatching on `control` class:
[`fixedControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/fixedControl.md)
-\>
[`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md)
(evaluate the OFV at fixed values) or
[`llpControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/llpControl.md)
-\>
[`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md)
(adaptive search for the OFV boundary via `optimProfile()`). All the
pieces share the `@family Profiling` roxygen tag.

**Multistart** (`R/multistart.R`).
[`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md)
is an S3 generic over fit/`rxUi`/`function`. `.msPerturbIni()` builds
perturbed `iniDf` starting points (`"uniform"`, `"lhs"`, `"normal"`,
respecting FIXED entries and bounds), an optional cheap empirical-Bayes
screen ranks them, and only the promising ones are fully estimated.
Every start is cached under `nlmixr2MultistartCache_<md5>/` as
`start`/`screen`/`fit` RDS files (`.msCacheRead()` / `.msCacheWrite()`),
so a run resumes where it stopped – the md5 covers anything that changes
what a starting point *is*. Returns a `nlmixr2Multistart` object with
[`print()`](https://rdrr.io/r/base/print.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) (waterfall and
parameter stability) and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods.

**Formula interface** (`R/nlmixrFormula.R`).
[`nlmixrFormula()`](https://nlmixr2.github.io/nlmixr2extra/reference/nlmixrFormula.md)
assembles a synthetic nlmixr2 model function from a `brms`-inspired
formula. `.nlmixrFormulaBuild()` drives parse -\> data prep -\>
parameter expansion -\>
[`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html)/[`model()`](https://nlmixr2.github.io/rxode2/reference/model.html)
assembly; covariate expansion (`.nlmixrFormulaExpandStartParam*()`)
emits `pop.<parameter>` intercepts and `cov_<param>_<start>` slopes.

**Preconditioning** (`R/precondition.R`, `src/preCondInv.cpp`). The
package’s only compiled code: an Rcpp/RcppArmadillo preconditioned
matrix inverse used by
[`preconditionFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/preconditionFit.md).
`src/init.c` hand-maintains the `.Call` registration table with a
hardcoded arity, so adding or changing the arity of an
`[[Rcpp::export]]` function requires both `Rcpp::compileAttributes(".")`
and a hand edit of `src/init.c`.

**Equation printing** (`R/knit_printEquation.R`). `knit_print` methods
for `rxUi`/`nlmixr2FitCore` that render a model as LaTeX in knitr
documents, via the `extractEqHelper` S3 walk over the model AST.

**Reporting helpers** (`R/AICHelpers.R`).
[`getMinAICFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/getMinAICFit.md),
[`listModelsTested()`](https://nlmixr2.github.io/nlmixr2extra/reference/listModelsTested.md),
[`isBoundaryFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/isBoundaryFit.md)
– both selection helpers exclude boundary fits by default.

### Cache directories

Bootstrap, covariate search and multistart all persist intermediate fits
to a `nlmixr2*Cache_*` directory in the working directory, keyed by a
digest of the inputs, and resume from it. These are listed in
`.gitignore` and `.Rbuildignore`; add any new cache prefix to both.

## Vignettes and the website

`vignettes/` is `.Rbuildignore`d, so the articles build **only for the
pkgdown website**, never during `R CMD check` or on CRAN. That is
deliberate: it lets an article run a real estimation.

Because real fits are slow and the pkgdown runner times out, every
expensive fit is **cached** with the `:=` operator from `nlmixr2save` (a
`Suggests` dependency). Each article opens with

``` r

library(nlmixr2save)
options(nlmixr2save.dir    = "cache",      # "../cache" under vignettes/articles/
        nlmixr2save.prefix = "<article>-",
        nlmixr2save.check  = FALSE)
```

An article renders with its own directory as the working directory, so
`vignettes/*.Rmd` uses `"cache"` and `vignettes/articles/*.Rmd` uses
`"../cache"`. Both then write to the single `vignettes/cache/`.

and then caches each result with `name := nlmixr2(...)`. The first
render fits and writes `vignettes/cache/<prefix><name>.zip` (a fit) or
`.rds` (anything else); later renders load it.
`nlmixr2save.check = FALSE` means the committed cache is trusted and
regenerated only when missing.

**The `:=` target name is the cache key**, so every target must be
unique within an article – do not reuse one variable (`mod`, `fit`,
`ms`) for several different fits.

Populate and refresh the cache with `vignettes/precompute.R`, which
renders each registered article in its own fresh R subprocess (building
many rxode2 models in one long-lived session eventually fails):

``` sh
cd vignettes && Rscript precompute.R           # fit only what is missing
cd vignettes && Rscript precompute.R --clean   # clear cache/ and refit everything
```

When adding a cached article, add its file name to the `vignettes`
vector in `precompute.R`, and commit `vignettes/cache/<prefix>*`
alongside the `.Rmd`.

Keep a chunk `eval=FALSE` only when it is a pure API illustration whose
output is not shown, or when running it would clobber a cached object or
merely repeat work already done earlier in the article.

## R Code Style

Follow the same conventions as `rxode2`:

- **Exported functions**: `camelCase` (`bootstrapFit`,
  `covarSearchAuto`, `multistartControl`)

- **Internal/non-exported functions**: `.camelCase` with a leading dot
  (`.msPerturbIni`, `.expandPopExpr`, `.lassoUicovariate`).
  Subsystem-local helpers use a short prefix after the dot (`.ms*` for
  multistart).

- **Local variables inside functions**: `.camelCase` with a leading dot.
  Older files (`R/computingutil.R`, `R/lassocov.R`, `R/SCM.R`) predate
  this and use bare `snake_case`/`camelCase` locals; match the
  surrounding file when editing, use the dot convention in new code.

- **S3 methods**: `generic.class` (`multistart.rxUi`,
  `iivSearch.nlmixr2Linearize`). Register them in roxygen with
  `@export`, not by hand in `NAMESPACE`.

- Avoid `snake_case` for new names.

- American English spelling.

- **Never write `pkg:::foo` in package code or in a test.** CodeFactor
  flags every `:::` as a Major Maintainability issue and fails the PR
  check. Take a namespace handle once near the top of the file and call
  through it, as `tests/testthat/test-parsing.R` and
  `test-ini-prior-column.R` do:

  ``` r

  .cur <- loadNamespace("nlmixr2extra")
  .cur$.expandPopExpr(...)
  ```

  `utils::getFromNamespace(".foo", "rxode2")` works the same way for a
  single internal from another package.

- Do not rename already-exported functions, even obviously misspelled
  ones (`buildupatedUI`, `horseshoeSummardf`): they are released API and
  reverse dependencies rely on them.

## Documentation and Comment Style

- Keep comments and roxygen terse. Condense multi-line explanations to
  one-liners; state the fact, not the story behind it. Keep every
  `@param`, `@return`, `@export`, `@family` and `@examples` tag.

- Long-running examples go in `\dontrun{}`.

- `NEWS.md` is organized per version (`# nlmixr2extra X.Y.Z`),
  user-facing changes first under `## New features`, then
  `## Bug fixes`. Entries are past-tense bullets of a sentence or two,
  referencing the issue number where one exists.

- **ASCII only. No Unicode anywhere in the repo** (CRAN requirement, and
  a Unicode character in `R/` or `man/` currently fails the check): use
  `--` for em-dashes, `-` for en-dashes, `->` for arrows, straight
  quotes, `...` for ellipses, and spell out Greek letters (`Delta`, not
  the symbol). This applies to plot labels, roxygen comments, tests,
  vignettes and `README.md` alike. Check with:

  ``` sh
  grep -rnP "[^\x00-\x7F]" R/ src/ tests/ man/ vignettes/ NEWS.md README.md
  ```
