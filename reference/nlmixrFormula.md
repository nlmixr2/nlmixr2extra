# A simple formula-based interface for nlmixr2

A simple formula-based interface for nlmixr2

## Usage

``` r
nlmixrFormula(
  object,
  data,
  start,
  param = NULL,
  paramLink = NULL,
  ...,
  residualModel = ~add(addSd)
)
```

## Arguments

- object:

  The formula defining the model (see details)

- data:

  The data to fit

- start:

  A named list of starting estimates. The names define the parameters in
  the model. If a single parameter estimate is desired, it can be given
  here. If a parameter estimate per factor level is desired, either a
  single starting estimate can be given across all factor levels or one
  estimate may be given per factor level. (Specify the factors with the
  `param` argument.)

- param:

  A formula or list of two-sided formulas giving the model used for
  parameters. If a parameter is a simple fixed effect, only, then it
  should not be included here. The right-hand side of a `param` formula
  names one or more columns in `data`: factor columns get a separate
  fixed effect per level; numeric columns get a linear slope. Multiple
  covariates on the same parameter can be combined with `+`, e.g.
  `b ~ z + w`.

- paramLink:

  Optional named character vector mapping parameter names to a link
  function. `"identity"` (the default) emits
  `<param> <- <linear combination>`; `"log"` wraps the linear
  combination in [`exp()`](https://rdrr.io/r/base/Log.html) so the
  parameter is strictly-positive on the natural scale. Only parameters
  that have a `param` entry may appear here.

- ...:

  Arguments passed on to
  [`nlmixr2est::nlmixr`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)

  `est`

  :   estimation method (all methods are shown by \`nlmixr2AllEst()\`).
      Methods can be added for other tools

  `control`

  :   The estimation control object. These are expected to be different
      for each type of estimation method

  `table`

  :   The output table control object (like \`tableControl()\`)

  `save`

  :   Boolean to save a nlmixr2 object in a rds file in the working
      directory. If `NULL`, uses option "nlmixr2.save"

  `envir`

  :   Environment where the nlmixr object/function is evaluated before
      running the estimation routine.

- residualModel:

  The residual model formula to use as a one-sided formula. The default
  is `~ add(addSd)`; richer models such as `~ add(addSd) + prop(propSd)`
  are supported as long as the corresponding sigma parameter names
  appear in `start`.

## Value

The model fit from
[`nlmixr2est::nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)

## Details

The formula is given with different notation than typical formulas. The
formula notation is inspired by and similar to `lme4::nlmer()`. It is a
3-part formula: `dependentVariable~predictorEquation~randomEffects`.

The `dependentVariable` is any variable in the dataset. It may not
include any math; for example, `log(DV)` is not allowed.

The `predictorEquation` is any valid math, and it will be used directly
in the nlmixr2 model.

The `randomEffects` are one or more random effect parameters defined by
putting the parameter in parentheses and putting a vertical bar and the
grouping parameter. Only one grouping parameter is allowed for all
random effects. An example would be `(slope|ID)` to estimate a random
effect parameter named "slope" for each "ID" in the data.

## Examples

``` r
nlmixrFormula(
  height ~ (Asym+AsymRe)+(R0-(Asym+AsymRe))*exp(-exp(lrc)*age) ~ (AsymRe|Seed),
  data = Loblolly,
  start = list(Asym = 103, R0 = -8.5, lrc = -3.3, addSd=1),
  est="focei"
)
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate jacobian
#> → calculate ∂(f)/∂(η)
#> → calculate ∂(R²)/∂(η)
#> → finding duplicate expressions in inner model...
#> → optimizing duplicate expressions in inner model...
#> → finding duplicate expressions in EBE model...
#> → optimizing duplicate expressions in EBE model...
#> → compiling inner model...
#>  
#>  
#> ✔ done
#> → finding duplicate expressions in FD model...
#> → optimizing duplicate expressions in FD model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> → compiling events FD model...
#>  
#>  
#> ✔ done
#> rxode2 5.1.2 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
#> → Calculating residuals/tables
#> ✔ done
#> ── nlmixr² FOCEi (outer: nlminb) ──
#> 
#>           OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> FOCEi 83.99509 248.3768 260.5308      -119.1884        73012.56        75.26153
#> 
#> ── Time (sec value$time): ──
#> 
#>              setup   optimize  covariance preprocess postprocess table compress
#> elapsed 0.03265231 0.02543444 0.008685938      0.039       0.021 0.054    0.001
#>            other
#> elapsed 2.486227
#> 
#> ── Population Parameters (value$parFixed or value$parFixedDf): ──
#> 
#>        Est.     SE %RSE Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> Asym    103   2.53 2.45           103 (98, 108)    2.24     -58.3% 
#> R0    -8.48  0.212  2.5     -8.48 (-8.9, -8.07)                    
#> lrc   -3.25 0.0397 1.22    -3.25 (-3.33, -3.18)                    
#> addSd  0.73                                0.73                    
#>  
#>   Covariance Type (value$covMethod): r,s
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance (value$omega) 
#>     or correlation (value$omegaR; diagonals=SDs)
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in value$shrink 
#>   Information about run found (value$runInfo):
#>    • gradient problems with initial estimate and covariance; see $scaleInfo 
#>    • last objective function was not at minimum, possible problems in optimization 
#>    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#>    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#>   Censoring (value$censInformation): No censoring
#>   Minimization message (value$message):  
#>     relative convergence (4) 
#> 
#> ── Fit Data (object value is a modified tibble): ──
#> # A tibble: 84 × 14
#>   ID     TIME    DV  PRED    RES   WRES IPRED   IRES  IWRES CPRED   CRES  CWRES
#>   <fct> <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>  <dbl>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 301       1  4.51  3.70  0.807  1.05   3.90  0.611  0.838  3.70  0.807  1.05 
#> 2 301       2 10.9  11.1  -0.187 -0.226 11.4  -0.501 -0.687 11.1  -0.187 -0.226
#> 3 301       3 28.7  27.2   1.51   1.48  27.8   0.942  1.29  27.2   1.51   1.48 
#> # ℹ 81 more rows
#> # ℹ 2 more variables: AsymRe <dbl>, age <dbl>
```
