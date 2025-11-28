# Extract the equations from an nlmixr2/rxode2 model to produce a 'LaTeX' equation.

Extract the equations from an nlmixr2/rxode2 model to produce a 'LaTeX'
equation.

## Usage

``` r
# S3 method for class 'nlmixr2FitCore'
knit_print(x, ..., output = "equations")

# S3 method for class 'rxUi'
knit_print(x, ...)
```

## Arguments

- x:

  The model to extract equations from

- ...:

  Ignored

- output:

  The type of output to request (currently, just "equations")
