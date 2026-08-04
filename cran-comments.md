# nlmixr2extra 5.1.1

This is a bug-fix release:

- `covarSearchAuto()` no longer errors when a covariate is selected, and
  the forward inclusion test no longer has an inverted sign

- `bootstrapFit()` works for models with a single population parameter,
  a single random effect, or no random effects

- `optimUnisampling()` keeps `N` and `floorT` when it retries internally

- The bundled `theoFitOde` fit was regenerated so it no longer requires
  the (suggested-elsewhere) `qs2` package to be decoded
