# Stratified cross-validation fold generator function, inspired from the caret

Stratified cross-validation fold generator function, inspired from the
caret

## Usage

``` r
foldgen(data, nfold = 5, stratVar = NULL)
```

## Arguments

- data:

  data frame used in the analysis

- nfold:

  number of k-fold cross validations. Default is 5

- stratVar:

  Stratification Variable. Default is NULL and ID is used for CV

## Value

return data.frame with the fold column attached

## Author

Vishal Sarsani, caret

## Examples

``` r
d <- nlmixr2data::theo_sd
d$SEX <-0
d$SEX[d$ID<=6] <-1

covarsVec <- c("WT")

# Stratified cross-validation data with CMT
df1 <- foldgen(d, nfold=5, stratVar="CMT")

# Stratified cross-validation data with ID (individual)
df2 <- foldgen(d, nfold=5, stratVar=NULL)
```
