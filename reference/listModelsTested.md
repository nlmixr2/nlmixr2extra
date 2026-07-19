# Make a list of models tested ready for reporting

Make a list of models tested ready for reporting

## Usage

``` r
listModelsTested(fitList, caption, excludeBoundary = TRUE, k = 2)
```

## Arguments

- fitList:

  A named list of models that were tested

- caption:

  The caption for the resulting table

- excludeBoundary:

  Exclude a model from selection if it has a parameter at its boundary

- k:

  numeric, the *penalty* per parameter to be used; the default `k = 2`
  is the classical AIC.

## Value

A data.frame with a "caption" attribute ready for printing in a report
with [`pander::pander()`](https://rdrr.io/pkg/pander/man/pander.html).
The data.frame will have names of "Description", "AIC", and "dAIC"; if
exclusions are applied, the data.frame will also have an "Exclude"
column.
