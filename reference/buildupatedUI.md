# Build updated from the covariate and variable vector list

Build updated from the covariate and variable vector list

## Usage

``` r
buildupatedUI(ui, varsVec, covarsVec, add = TRUE, indep = FALSE)
```

## Arguments

- ui:

  compiled rxode2 nlmir2 model or fit

- varsVec:

  character vector of variables that need to be added

- covarsVec:

  character vector of covariates that need to be added

- add:

  boolean indicating if the covariate needs to be added or removed

- indep:

  a boolean indicating if the covariates should be added independently,
  or sequentially (append to the previous model). only applicable to
  adding covariate

## Value

updated ui with added covariates

## Author

Vishal Sarsani
