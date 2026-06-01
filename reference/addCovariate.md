# Add Covariate to Model Fit (Generic)

Add Covariate to Model Fit (Generic)

## Usage

``` r
addCovariate(fit, expr, effect, ref, ...)

# Default S3 method
addCovariate(fit, expr, effect = "power", ref = "median", ...)

# S3 method for class 'rxUi'
addCovariate(fit, expr, effect = "power", ref = "median", ...)

# S3 method for class 'nlmixr2FitCore'
addCovariate(fit, expr, effect = "power", ref = "median", ...)

# S3 method for class 'nlmixr2Linearize'
addCovariate(fit, expr, effect = "power", ...)
```

## Arguments

- fit:

  Model fit object

- expr:

  Expression, or vector or list of expressions. eg. CL ~ WT/70 +
  AGE/80 + ... .

- effect:

  character or list of characters of "linear", "piece_lin", "exp",
  "power". see details

- ref:

  Default normalization for continuous covariates. "mean", "median" or
  NA. Default "median"

- ...:

  other parameters passed to each specific signature

  `effect` and `normaDefault` are only used if covariate is continuous.
  `normaDefault` call will be skipped if what covpar expression is
  normalized by other value

## Author

Omar I. Elashkar
