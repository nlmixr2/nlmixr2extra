# Generate LaTeX for if blocks

Generate LaTeX for if blocks

## Usage

``` r
# S3 method for class '`if`'
extractEqHelper(x, ..., inModel, alignment, indent = 0L, firstIf = 0L, name)
```

## Arguments

- x:

  language object to convert to LaTeX

- ...:

  additional arguments passed to methods

- inModel:

  logical; whether the expression is inside a model block

- indent:

  tracks the current indention level

- firstIf:

  tracks if this is the first if at the current indent level

- name:

  optional name for the expression
