# Sample from uniform distribution by optim

Sample from uniform distribution by optim

## Usage

``` r
optimUnisampling(xvec, N = 1000, medValue, floorT = TRUE)
```

## Arguments

- xvec:

  A vector of min,max values . Ex:c(10,20)

- N:

  Desired number of values

- medValue:

  Desired Median

- floorT:

  boolean indicating whether to round up

## Value

Samples with approx desired median.

## Author

Vishal Sarsani

## Examples

``` r
# Simulate 1000 creatine clearance values with median of 71.7 within range of c(6.7,140)
creatCl <- optimUnisampling(xvec=c(6.7,140), N=1000, medValue = 71.7, floorT=FALSE)
```
