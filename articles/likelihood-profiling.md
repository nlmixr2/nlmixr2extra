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
```

### Profile a single parameter

``` r

prof_tka <- profile(fit, which = "tka")
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

prof_tka <- profile(fit, which = "tka", control = as.list(ctrl))
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
tka_est <- nlmixr2est::fixef(fit)[["tka"]]

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
