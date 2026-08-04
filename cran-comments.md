# nlmixr2extra 5.2.0

New features:

- Reporting helpers for comparing candidate models: `getMinAICFit()`,
  `listModelsTested()` and `isBoundaryFit()`

Bug fixes:

- `bootstrapFit(stratVar=)` now actually resamples, keeps subject ids
  distinct across strata, and draws whole subjects

- `covarSearchAuto()` no longer errors when a covariate is selected, and
  the forward inclusion test no longer has an inverted sign

- `bootstrapFit()` works for models with a single population parameter,
  a single random effect, or no random effects

- `optimUnisampling()` keeps `N` and `floorT` when it retries internally

- The bundled `theoFitOde` fit was regenerated so it no longer requires
  the `qs2` package to be decoded
