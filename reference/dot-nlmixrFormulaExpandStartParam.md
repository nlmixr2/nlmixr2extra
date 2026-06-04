# Expand parameters to include their covariate representations, if applicable.

Expand parameters to include their covariate representations, if
applicable.

## Usage

``` r
.nlmixrFormulaExpandStartParam(start, param, paramLink = NULL, data)

.nlmixrFormulaExpandStartParamSingle(
  startName,
  startValue,
  param = NULL,
  link = NULL,
  data = NULL
)
```

## Arguments

- start:

  the starting values for the model

- param:

  The parameter in the model

- paramLink:

  Named character vector giving the link function for each parameter
  that has a `param` entry. Recognised values: `"identity"` (default)
  emits `<param> <- <linear combination>`; `"log"` wraps the linear
  combination in [`exp()`](https://rdrr.io/r/base/Log.html).

- data:

  The dataset

- startName:

  The base name for the parameter

- startValue:

  The initial value for the base parameter

- link:

  Optional link function for this parameter (`"identity"` or `"log"`).
  Defaults to `"identity"`.

## Value

A list with the ini and model parts needed to use the parameter

## Functions

- `.nlmixrFormulaExpandStartParamSingle()`: Expand a single parameter in
  a model using dataset covariates, if applicable
