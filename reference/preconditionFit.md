# Linearly re-parameterize the model to be less sensitive to rounding errors

Linearly re-parameterize the model to be less sensitive to rounding
errors

## Usage

``` r
preconditionFit(fit, estType = c("full", "posthoc", "none"), ntry = 10L)
```

## Arguments

- fit:

  A nlmixr2 fit to be preconditioned

- estType:

  Once the fit has been linearly reparameterized, should a "full"
  estimation, "posthoc" estimation or simply a estimation of the
  covariance matrix "none" before the fit is updated

- ntry:

  number of tries before giving up on a pre-conditioned covariance
  estimate

## Value

A nlmixr2 fit object that was preconditioned to stabilize the
variance/covariance calculation

## References

Aoki Y, Nordgren R, Hooker AC. Preconditioning of Nonlinear Mixed
Effects Models for Stabilisation of Variance-Covariance Matrix
Computations. AAPS J. 2016;18(2):505-518. doi:10.1208/s12248-016-9866-5
