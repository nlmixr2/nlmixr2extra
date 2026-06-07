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
if (FALSE) { # \dontrun{
nlmixrFormula(
  height ~ (Asym+AsymRe)+(R0-(Asym+AsymRe))*exp(-exp(lrc)*age) ~ (AsymRe|Seed),
  data = Loblolly,
  start = list(Asym = 103, R0 = -8.5, lrc = -3.3, addSd=1),
  est="focei"
)
} # }
```
