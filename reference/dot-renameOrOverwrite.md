# Rename a column in a dataset

If `newName` already exists in `data` and would be overwritten with a
different source column, this is treated as a hard error so that user
data is never silently clobbered.

## Usage

``` r
.renameOrOverwrite(data, newName, oldName)
```

## Arguments

- data:

  The dataset to modify

- newName, oldName:

  The new and old column names

## Value

data with `data[[newName]] <- data[[charOld]]`

## Examples

``` r
if (FALSE) { # \dontrun{
.renameOrOverwrite(data.frame(A=1), newName="B", oldName="A")
} # }
```
