# Add Individual Random Effects and Fix them to Small Value

Add Individual Random Effects and Fix them to Small Value

## Usage

``` r
addAllEtas(ui, fix = FALSE)
```

## Arguments

- ui:

  model with no individual random effects added

- fix:

  If TRUE, the added etas will be fixed (not estimatable). Default is
  FALSE.

## Value

nlmixr2 fit with individual random effects added on all fixed effects
and fixed.

## Author

Omar I. Elashkar
