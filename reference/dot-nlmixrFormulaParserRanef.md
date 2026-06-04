# Parse the random effects part of a formula

Parse the random effects part of a formula

## Usage

``` r
.nlmixrFormulaParserRanef(object)
```

## Arguments

- object:

  The formula to parse

## Value

An unnamed list with one element per random effect. The elements each
have names of "ranefVar", "ranefGroup", and "start" indicating the
variable with the random effect, the grouping variable for the random
effect, and the starting value for the standard deviation of the random
effect.

## Examples

``` r
if (FALSE) { # \dontrun{
.nlmixrFormulaParserRanef(str2lang("(c|id)+(d|id2)"))
} # }
```
