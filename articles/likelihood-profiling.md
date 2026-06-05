# Likelihood Profiling for NLME Models

## Introduction

Wald-based confidence intervals (CIs), which are the default output from
`nlmixr2`, rely on the assumption that the parameter estimates are
normally distributed. For parameters near boundaries, or for variance
components of random effects, this assumption often breaks down and the
resulting CIs can be misleading.

**Likelihood profiling** provides a more reliable alternative. Rather
than approximating the likelihood surface with a quadratic, it evaluates
the actual objective function (OFV) at a grid of parameter values,
holding each parameter fixed in turn and re-estimating the rest. A
parameter’s confidence interval is the set of values whose OFV increase
relative to the minimum is below a threshold — by default
$`\chi^2_{0.95,1} \approx 3.84`$.

`nlmixr2extra` provides two profiling methods via the
[`profile()`](https://rdrr.io/r/stats/profile.html) generic:

| Method | Function | Description |
|----|----|----|
| `"llp"` | [`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md) | Adaptive log-likelihood profiling — searches for the exact OFV boundary |
| `"fixed"` | [`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md) | Evaluates the OFV at user-supplied fixed values |

## Quick start

### Define and fit a model

We use the single-dose theophylline dataset (`theo_sd`) from
`nlmixr2data`.

``` r

library(nlmixr2extra)
library(ggplot2)

one_cmt <- function() {
  ini({
    tka  <- log(1.57)   # log absorption rate constant
    tcl  <- log(2.72)   # log clearance
    tv   <- log(31.5)   # log volume (fixed)
    eta.ka ~ 0.6
    eta.cl ~ 0.09
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

fit <- nlmixr2(one_cmt, data = nlmixr2data::theo_sd,
               est = "focei", control = list(print = 0))
fit
```

### Profile all parameters

``` r

# Profiles tka and tcl; eta.ka, eta.cl, and add.sd are also included.
# This may take a few minutes.
prof_all <- profile(fit)
prof_all

# not run
#>      Parameter           OFV         tka        tcl         tv       add.sd       eta.ka     eta.cl profileBound
#> 1        tka  1.302956e+02  0.40707834  1.0277804  3.4301847 7.790360e-01 3.313392e-01 0.11903041           NA
#> 2        tka  1.359081e+02 -0.06205406  1.0357450  3.4201803 7.822078e-01 5.373565e-01 0.11997097           NA
#> 3        tka            NA  0.03246907         NA         NA           NA           NA         NA    -3.841459
#> 4        tka  1.341359e+02  0.03253017  1.0349416  3.4203239 7.813290e-01 4.459202e-01 0.11954391           NA
#> 5        tka  1.340839e+02  0.03536930  1.0339534  3.4211996 7.793176e-01 4.506370e-01 0.11970980           NA
#> 6        tka  1.332169e+02  0.08598386  1.0344147  3.4212832 7.808917e-01 4.048129e-01 0.11966891           NA
#> 7        tka  1.333483e+02  0.73591286  1.0203980  3.4396473 7.820747e-01 4.405237e-01 0.12179708           NA
#> 8        tka  1.340970e+02  0.78149630  1.0155783  3.4395393 7.814438e-01 4.782030e-01 0.12240100           NA
#> 9        tka  1.341362e+02  0.78375611  1.0163098  3.4402083 7.812342e-01 4.813841e-01 0.12266684           NA
#> 10       tka            NA  0.78380398         NA         NA           NA           NA         NA     3.841459
#> 11       tka  1.357760e+02  0.87621075  1.0136331  3.4402359 7.818523e-01 5.617630e-01 0.12254599           NA
#> 12       tcl  1.302956e+02  0.40707834  1.0277804  3.4301847 7.790360e-01 3.313392e-01 0.11903041           NA
#> 13       tcl  1.663811e+02  0.42951804 -0.1566724  3.4457119 7.849234e-01 3.329437e-01 0.59514936           NA
#> 14       tcl            NA          NA  0.7988096         NA           NA           NA         NA    -3.841459
#> 15       tcl  1.341352e+02  0.41822580  0.7988652  3.4390860 7.835247e-01 3.305029e-01 0.17116048           NA
#> 16       tcl  1.341188e+02  0.41746955  0.7994066  3.4384609 7.832626e-01 3.312768e-01 0.17139417           NA
#> 17       tcl  1.340151e+02  0.41805477  0.8030235  3.4381714 7.837632e-01 3.325469e-01 0.17017609           NA
#> 18       tcl  1.334242e+02  0.41788913  0.8242418  3.4379266 7.816312e-01 3.310546e-01 0.15915788           NA
```

### Profile a single parameter

``` r

prof_tka <- profile(fit, which = "tka")

# not run
#>      Parameter           OFV         tka        tcl         tv       add.sd       eta.ka     eta.cl profileBound
#> 1        tka  1.302956e+02  0.40707834  1.0277804  3.4301847 7.790360e-01 3.313392e-01 0.11903041           NA
#> 2        tka  1.359081e+02 -0.06205406  1.0357450  3.4201803 7.822078e-01 5.373565e-01 0.11997097           NA
#> 3        tka            NA  0.03246907         NA         NA           NA           NA         NA    -3.841459
#> 4        tka  1.341359e+02  0.03253017  1.0349416  3.4203239 7.813290e-01 4.459202e-01 0.11954391           NA
#> 5        tka  1.340839e+02  0.03536930  1.0339534  3.4211996 7.793176e-01 4.506370e-01 0.11970980           NA
#> 6        tka  1.332169e+02  0.08598386  1.0344147  3.4212832 7.808917e-01 4.048129e-01 0.11966891           NA
#> 7        tka  1.333483e+02  0.73591286  1.0203980  3.4396473 7.820747e-01 4.405237e-01 0.12179708           NA
#> 8        tka  1.340970e+02  0.78149630  1.0155783  3.4395393 7.814438e-01 4.782030e-01 0.12240100           NA
#> 9        tka  1.341362e+02  0.78375611  1.0163098  3.4402083 7.812342e-01 4.813841e-01 0.12266684           NA
#> 10       tka            NA  0.78380398         NA         NA           NA           NA         NA     3.841459
#> 11       tka  1.357760e+02  0.87621075  1.0136331  3.4402359 7.818523e-01 5.617630e-01 0.12254599           NA
```

## The log-likelihood profiling method (`"llp"`)

[`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md)
uses an adaptive algorithm:

1.  Start from the point estimate and step in both directions.
2.  At each step, fix the target parameter and re-estimate all others.
3.  Stop when the OFV change is within `ofvtol` of the desired
    `ofvIncrease`, or when the next parameter step would not change its
    rounded value to `paramDigits` significant figures.

### Control options

All LLP settings are set via
[`llpControl()`](https://nlmixr2.github.io/nlmixr2extra/reference/llpControl.md):

``` r

# Default settings
ctrl <- llpControl(
  ofvIncrease   = qchisq(0.95, df = 1),  # ~3.84 for 95% CI
  rseTheta      = 30,    # starting step size as % of the estimate
  itermax       = 10,    # max iterations per direction
  ofvtol        = 0.005, # tolerance on OFV difference
  paramDigits   = 3      # significant digits for convergence
)

prof_tka <- profile(fit, which = "tka", control = ctrl)
```

To get a 90% CI instead of 95%, lower the threshold:

``` r

prof_tka_90 <- profile(fit, which = "tka",
                       control = list(ofvIncrease = qchisq(0.90, df = 1)))
```

### Interpreting the output

The returned data frame has one row per model evaluation. Key columns:

- `Parameter` — which parameter was fixed on this row
- `OFV` — the objective function value at that step
- `profileBound` — present on the boundary rows; its absolute value is
  the `ofvIncrease` threshold and its sign indicates lower (`-`) or
  upper (`+`)
- One column per model parameter showing its estimate at each step

``` r

# Show rows near the boundary
prof_tka[!is.na(prof_tka$profileBound), ]
```

The rows where `profileBound` is non-`NA` give the profile CI limits
directly from the `tka` column.

## Plotting the profile

The profile data frame contains everything needed for a profile plot. A
standard presentation shows OFV – OFV_min on the y-axis and the
parameter value on the x-axis, with a horizontal reference line at the
CI threshold.

``` r

# Compute ΔOFV relative to the minimum
ofv_min <- min(prof_tka$OFV, na.rm = TRUE)
prof_tka$dOFV <- prof_tka$OFV - ofv_min

# The two boundary rows (one per direction)
bounds <- prof_tka[!is.na(prof_tka$profileBound), ]

ggplot(prof_tka, aes(x = tka, y = dOFV)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = qchisq(0.95, 1), linetype = "dashed", colour = "red") +
  geom_vline(xintercept = bounds$tka, linetype = "dotted", colour = "blue") +
  labs(
    x     = "tka",
    y     = expression(Delta * "OFV"),
    title = "Likelihood profile for tka",
    subtitle = sprintf(
      "95%% profile CI: [%.3f, %.3f]",
      min(bounds$tka), max(bounds$tka)
    )
  ) +
  theme_bw()
```

A symmetric, approximately parabolic profile indicates the Wald CI is
reliable. An asymmetric or flat profile indicates non-normality; in
these cases the profile CI is more trustworthy.

## The fixed-point method (`"fixed"`)

[`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md)
evaluates the OFV at a user-supplied grid rather than searching
adaptively. It is useful when you want to reproduce a specific grid
(e.g. for a table or a plot comparing multiple models) or when the
adaptive algorithm struggles.

### Evaluating a uniform grid

``` r

# Evaluate tka at seven equally spaced values around the point estimate
tka_est <- fit$theta["tka"]

grid <- data.frame(tka = seq(tka_est - 0.6, tka_est + 0.6, length.out = 7))

prof_grid <- profile(fit, which = grid, method = "fixed")
prof_grid[, c("Parameter", "tka", "OFV")]
```

### Two-parameter joint profile

[`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md)
accepts a data frame with multiple columns. Each row fixes all named
parameters simultaneously, allowing you to trace joint profile contours.

``` r

grid2d <- expand.grid(
  tka = seq(tka_est - 0.4, tka_est + 0.4, length.out = 5),
  tcl = seq(nlmixr2est::fixef(fit)[["tcl"]] - 0.2,
            nlmixr2est::fixef(fit)[["tcl"]] + 0.2, length.out = 5)
)

prof_joint <- profile(fit, which = grid2d, method = "fixed")

ggplot(prof_joint, aes(x = tka, y = tcl,
                       fill = OFV - min(OFV, na.rm = TRUE))) +
  geom_tile() +
  scale_fill_viridis_c(name = expression(Delta * "OFV")) +
  geom_contour(aes(z = OFV - min(OFV, na.rm = TRUE)),
               breaks = qchisq(0.95, 2), colour = "white") +
  labs(title = "Joint profile: tka × tcl",
       subtitle = "White contour = 95% joint confidence region") +
  theme_bw()
```

## Comparing profile CIs with Wald CIs

Profile CIs are wider than Wald CIs when the likelihood surface is
asymmetric. This is common for variance parameters and log-scale
parameters near zero.

``` r

# Wald CI from the fit
wald_ci <- confint(fit)  # uses the covariance matrix

# Extract profile CI from the boundary rows
bound_rows <- prof_all[!is.na(prof_all$profileBound), ]

# For each parameter, the two bound_rows give lower and upper limits
profile_ci <- tapply(
  seq_len(nrow(bound_rows)),
  bound_rows$Parameter,
  function(idx) {
    vals <- bound_rows[idx, bound_rows$Parameter[idx[1]]]
    c(lower = min(vals), upper = max(vals))
  }
)

# Show side-by-side
do.call(rbind, lapply(names(profile_ci), function(p) {
  data.frame(
    parameter   = p,
    wald_lower  = wald_ci[p, 1],
    wald_upper  = wald_ci[p, 2],
    prof_lower  = profile_ci[[p]]["lower"],
    prof_upper  = profile_ci[[p]]["upper"]
  )
}))
```

## Tips

- **Start with `which`**: profiling all parameters is time-consuming.
  Profile only the parameters of primary interest first.
- **Increase `itermax`** if you see “aborted due to too many
  iterations”. The default of 10 is conservative; 20–30 is often
  sufficient for variance parameters.
- **Adjust `rseTheta`**: if the initial step is too large (OFV jumps
  past the target immediately) reduce `rseTheta`. If it is too small
  (many steps needed before the boundary is reached) increase it.
- **Fixed-point fallback**: when
  [`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md)
  fails to converge for a particular parameter,
  [`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md)
  on a coarse grid around the estimate can give a rough CI quickly.
