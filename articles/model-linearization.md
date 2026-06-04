# Model Linearization for IIV and Residual Error Search

## Introduction

Fitting a nonlinear mixed-effects (NLME) model for every candidate
random-effect or residual-error structure is computationally expensive.
**Model linearization** replaces the nonlinear model with a first-order
Taylor expansion around the current individual predictions. The
resulting linear model fits in a fraction of the time, allowing
exhaustive search over:

- **Inter-individual variability (IIV) structure** — which parameters
  carry random effects, and which pairs are correlated
  ([`iivSearch()`](https://nlmixr2.github.io/nlmixr2extra/reference/iivSearch.md)).
- **Residual error model** — additive, proportional, or combined
  variance structures
  ([`resSearch()`](https://nlmixr2.github.io/nlmixr2extra/reference/resSearch.md)).

The workflow has three stages:

1.  **Fit** the nonlinear model (any method; at least one $`\eta`$).
2.  **Linearize** the fit
    ([`linearize()`](https://nlmixr2.github.io/nlmixr2extra/reference/linearize.md)).
3.  **Search** IIV or residual structure on the linearized model, then
    refit the top candidates with the original nonlinear model
    ([`rerunTopN()`](https://nlmixr2.github.io/nlmixr2extra/reference/rerunTopN.md)).

## Theory

The FOCEI objective function at the MAP estimates $`\hat\eta_i`$ is:

``` math
\text{OFV}_i = -2\log p(y_i|\hat\eta_i) - 2\log p(\hat\eta_i)
```

A first-order Taylor expansion of $`f(t, \hat\eta_i)`$ and
$`r(t, \hat\eta_i)`$ (the residual variance) around the population
prediction $`\hat\eta_i = 0`$ gives:

``` math
y_i \approx f_0 + \underbrace{\sum_k \frac{\partial f}{\partial\eta_k}\bigg|_{\hat\eta}\!\!(\eta_k - \hat\eta_k)}_{\text{BASE\_TERMS}} + \varepsilon
```

where $`\varepsilon \sim \mathcal{N}(0,\,r_0^2)`$ and the linearized
residual variance is similarly expanded. This is the “FOCE
approximation”. Adding the residual-variance gradient gives the full
“FOCEI approximation”.

Because the linearized model is itself linear in the $`\eta`$
parameters, each candidate structure requires only **one fast FOCEI
fit** rather than a full nonlinear optimization.

## Setup

``` r

library(nlmixr2extra)
library(ggplot2)
```

## Simulate data for demonstration

Define a one-compartment model with additive residual error and IIV with
correlation between CL and V:

``` r

one_cmt <- function() {
  ini({
    tcl <- log(2.7) # Cl
    tv <- log(30) # V
    tka <- log(1.56) #  Ka
    eta.cl + eta.v ~ sd(cor(0.3, 0.99, 0.5))
    eta.ka ~ 0
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
  }

set.seed(42)
ev <- rxode2::et(amt = 300, cmt = "depot") |>
  rxode2::et(time = c(0.25, 0.5, 1, 2, 3, 6, 8, 12, 16, 24))
sim <- rxode2::rxSolve(one_cmt, ev, nSub = 50, addDosing = TRUE)
sim$dv <- sim$sim
sim$id <- sim$sim.id
sim <- sim[, c("id", "time", "amt", "dv", "evid")]
```

## Linearizing a fit

### Define a base model

Define a base model with *no IIV* and additive residual error:

``` r

one_cmt_base <- function() {
  ini({
      tcl <- log(2.7) # Cl
      tv <- log(30) # V
      tka <- log(1.56) #  Ka
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
```

### Base Model Fit

Fit the data to the base model `one_cmt_base` with typical NLME method
(e.g. FOCEI):

``` r

fit <- nlmixr2(one_cmt_base, sim, est = "focei") 
```

### Run Linearization

[`linearize()`](https://nlmixr2.github.io/nlmixr2extra/reference/linearize.md)
takes a fitted `nlmixr2` object and returns an `nlmixr2Linearize` fit.

``` r

# addEtas = TRUE adds fixed etas on every theta before linearizing,
# giving derivatives with respect to every parameter.
fitLin <- linearize(fit, addEtas = TRUE)
```

The function:

1.  If `addEtas = TRUE`, adds small fixed etas on all thetas and
    re-evaluates the model with `maxOuterIterations = 0` to extract
    derivatives.
2.  Builds the linearized model using
    [`linModGen()`](https://nlmixr2.github.io/nlmixr2extra/reference/linModGen.md).
3.  Fits the linear model via FOCEI over a vector of `mceta` values,
    stopping when the relative OFV deviation is within `relTol`.

### Checking linearization quality

After
[`linearize()`](https://nlmixr2.github.io/nlmixr2extra/reference/linearize.md),
confirm that the linear approximation reproduces the original model
well:

``` r

match <- isLinearizeMatch(fitLin)

# OFV agreement (relative tolerance 10%)
match$ofv[[1]]   # TRUE/FALSE if OFV differ by less than tol
match$ofv[[2]]   # all.equal() message if FALSE 

# Omega matrix agreement
match$omega[[1]]  # TRUE/FALSE if Omega matrices differ by less than tol

# Individual eta agreement
match$eta[[1]] # TRUE/FALSE if all individual etas differ by less than tol


# Residual variance agreement
match$err[[1]] # TRUE/FALSE if residual variances differ by less than tol
```

A visual check is available via
[`linearizePlot()`](https://nlmixr2.github.io/nlmixr2extra/reference/linearizePlot.md):

``` r

linearizePlot(fitLin)
```

The plot shows original vs. linearized individual objective values and
etas. Points should lie on the identity line for a good approximation.

### FOCE vs. FOCEI linearization

The `focei` argument controls whether the residual-variance gradient is
included:

| `focei` | Description |
|----|----|
| `NA` (default) | Try FOCEI; switch to FOCE automatically if `relTol` is exceeded |
| `TRUE` | Always use FOCEI (individual + residual linearization) |
| `FALSE` | Always use FOCE (residual interaction linearization skipped) |

FOCEI is more accurate for models with heteroscedastic residuals
(proportional or combined error). FOCE is faster and more stable for
additive-only models.

``` r

# Force FOCE linearization
fitLinFoce <- linearize(fit, addEtas = TRUE, focei = FALSE)

# Force FOCEI
fitLinFocei <- linearize(fit, addEtas = TRUE, focei = TRUE)
```

### Tuning `mceta`

The `mceta` argument is a vector of Monte Carlo eta sample sizes tried
in order. The algorithm stops when the relative OFV deviation is within
`relTol`:

``` r

# Default: try -1, 10, 100, 1000 in order
fitLin <- linearize(fit, addEtas = TRUE, mceta = c(-1, 10, 100, 1000))

# For a difficult model, start with a larger mceta
fitLin <- linearize(fit, addEtas = TRUE, mceta = c(100, 500, 1000))
```

`mceta = -1` uses the exact gradient (no Monte Carlo sampling) and is
the fastest option. Larger values are useful when exact gradients cause
numerical issues.

## IIV structure search

### Overview

[`iivSearch()`](https://nlmixr2.github.io/nlmixr2extra/reference/iivSearch.md)
fits every combination of:

- Which thetas carry a random effect (any subset of those present)
- Which pairs of those etas are correlated

For $`n`$ etas, the number of candidate structures grows as
$`O(2^n \cdot
2^{n(n-1)/2})`$. For $`n = 3`$ this is around 64 models; for $`n = 4`$
it is around 1000. Restrict the search to a manageable number of etas
using
[`addAllEtas()`](https://nlmixr2.github.io/nlmixr2extra/reference/addAllEtas.md)
with `fix = TRUE` to identify which parameters have meaningful
random-effect signal before running the full search.

### Running the search

``` r

iivRes <- iivSearch(fitLin)
```

Progress is shown for each candidate structure. Failed fits
(non-convergence) are stored as `NA` in the summary and do not stop the
search.

### Examining results

``` r

# Print ordered by BIC
print(iivRes)

# The summary data frame
head(iivRes$summary[order(iivRes$summary$BIC), ])
```

The summary contains one row per candidate with columns:

| Column      | Description                           |
|-------------|---------------------------------------|
| `OBJF`      | Objective function value              |
| `AIC`       | Akaike information criterion          |
| `BIC`       | Bayesian information criterion        |
| `search`    | Structure string (see below)          |
| `nParams`   | Number of estimated parameters        |
| `covMethod` | Covariance method (`"r,s"` = success) |

#### Structure string format

The `search` column encodes the IIV structure:

- `eta.cl+eta.v` — both etas present, no correlation
- `eta.cl+eta.v+eta.cl~eta.v` — both etas with correlation between them
- `eta.cl` — only CL carries a random effect

### Refitting top candidates with the original model

The linearized results are used to rank candidates, but final inference
should be based on the original nonlinear model.
[`rerunTopN()`](https://nlmixr2.github.io/nlmixr2extra/reference/rerunTopN.md)
refits the top `n` structures:

``` r

# Refit the 5 best structures with the original nonlinear model
top5 <- rerunTopN(iivRes, n = 5)

# Results ordered by BIC from the nonlinear fits
top5$summary[order(top5$summary$O.BIC), ]
```

The summary column names are prefixed with `O.` (for “original”) to
distinguish them from the linearized model results.

### Visualising the IIV search

``` r

summ <- iivRes$summary
summ <- summ[!is.na(summ$BIC), ]
summ$rank <- rank(summ$BIC)

ggplot(summ, aes(x = rank, y = BIC - min(BIC))) +
  geom_point(aes(colour = covMethod == "r,s"), size = 2) +
  scale_colour_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                      name = "Converged") +
  labs(x = "BIC rank", y = expression(Delta * "BIC"),
       title = "IIV search: BIC relative to best model") +
  theme_bw()
```

## Residual error model search

[`resSearch()`](https://nlmixr2.github.io/nlmixr2extra/reference/resSearch.md)
tests three alternative residual structures on a linearized fit and
returns their OFV, AIC, and BIC for comparison.

The structures tested are:

| Name              | `rxR2` formula                         |
|-------------------|----------------------------------------|
| Base fit          | (whatever was in the linearized model) |
| Proportional      | `prop.sd^2 * OPRED^2`                  |
| Combined (type 2) | `prop.sd^2 * OPRED^2 + add.sd^2`       |
| Combined (type 1) | `(prop.sd * OPRED + add.sd)^2`         |

``` r

# Start from a linearized model; additive-only in this case
fitLinAdd <- linearize(fit, addEtas = TRUE)

isLinearizeMatch(fitLinAdd)  # check quality before searching

resRes <- resSearch(fitLinAdd)

# Compare by BIC
resRes$summary[order(resRes$summary$BIC), ]
```

The returned list also includes `resRes$originalFit` (the linearized fit
used as the base) for reference.

## Full workflow example

``` r

## 1. Fit the base nonlinear model (no IIV, additive error)
fit <- nlmixr2(one_cmt_base, sim, est = "focei")

## 2. Linearize, adding etas on all thetas
fitLin <- linearize(fit, addEtas = TRUE, focei = NA)

## 3. Check linearization quality
match <- isLinearizeMatch(fitLin)
stopifnot(match$ofv[[1]])   # abort if linearization is poor

## 4. Search IIV structure
iivRes <- iivSearch(fitLin)
print(iivRes)   # BIC-ordered summary

## 5. Refit top 3 structures with the original model
top3 <- rerunTopN(iivRes, n = 3)
best_structure <- top3$summary$search[which.min(top3$summary$O.BIC)]
cat("Best IIV structure:", best_structure, "\n")

## 6. Rebuild the model with the best IIV structure and search residual error
# (refit with the chosen IIV structure first)
# ...then:
resRes <- resSearch(fitLin)
resRes$summary[order(resRes$summary$BIC), ]
```

## Tips and troubleshooting

**Linearization fails to converge (`relTol` exceeded)**

- Try larger `mceta` values: `mceta = c(100, 500, 2000)`.
- Switch to `focei = FALSE` (FOCE) which is more numerically stable.
- Increase `relTol` slightly (e.g. `0.30`) for a looser quality
  threshold.

**`iivSearch` produces many `NA` entries**

Non-convergence is common for small or degenerate IIV structures. `NA`
rows are skipped in ranking; they do not indicate a bug. If more than
half the search fails, check whether the starting eta estimates (from
the linearized fit) are reasonable.

**`isLinearizeMatch` fails for omegas**

The omega matrix comparison uses the MAP estimates from the linearized
fit. Small discrepancies in variance components are expected; a
tolerance of 20% (`tol = 0.20`) is often appropriate for omega terms.

**`rerunTopN` is slow**

Each of the `n` calls to
[`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)
is a full nonlinear fit. To reduce runtime, use a coarser
[`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.html)
(e.g. fewer inner iterations) inside
[`rerunTopN()`](https://nlmixr2.github.io/nlmixr2extra/reference/rerunTopN.md)
for screening, then do a clean final fit on the selected structure.
