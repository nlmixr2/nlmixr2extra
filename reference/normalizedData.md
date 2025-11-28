# Function to return data of normalized covariates

Function to return data of normalized covariates

## Usage

``` r
normalizedData(data, covarsVec, replace = TRUE)
```

## Arguments

- data:

  a dataframe with covariates to normalize

- covarsVec:

  a list of covariate names (parameters) that need to be estimates

- replace:

  replace the original covariate data with normalized data for easier
  updated model.

## Value

data frame with all normalized covariates

## Author

Vishal Sarsani

## Examples

``` r
d <- nlmixr2data::theo_sd
d$SEX <-0
d$SEX[d$ID<=6] <-1

covarsVec <- c("WT")

# Normalized covariate (replaced)
df1 <- normalizedData(d,covarsVec,replace=TRUE)

# Normalized covariate (without replacement)
df2 <- normalizedData(d,covarsVec,replace=FALSE)
```
