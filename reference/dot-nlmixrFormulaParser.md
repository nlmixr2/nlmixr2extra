# Parse the formula to extract the dependent variable, predictor, and random effects

Parse the formula to extract the dependent variable, predictor, and
random effects

## Usage

``` r
.nlmixrFormulaParser(object)
```

## Arguments

- object:

  The formula to parse

## Value

A list with names of "DV", "predictor", and "ranef" which are each part
of the formula broken down into calls

## Examples

``` r
if (FALSE) { # \dontrun{
.nlmixrFormulaParser(a~b+c~(c|id))
} # }
```
