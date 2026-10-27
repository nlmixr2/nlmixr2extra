# nlmixr2extra events

## Introduction

`rxode2` provides a small event bus: packages *emit* an event when a
user-facing operation finishes, and any package can *listen* for it with
[`rxode2::rxEventListen()`](https://nlmixr2.github.io/rxode2/reference/rxEventListen.html).
The main consumer is
[nlmixr2log](https://github.com/nlmixr2/nlmixr2log), which records every
fit, simulation and result as a run, but a listener can be any function.
See
[`?rxode2::rxEventListen`](https://nlmixr2.github.io/rxode2/reference/rxEventListen.html)
for the bus itself.

The key rule of the bus is that events are only delivered from the
**outermost** operation. Every emitting function enters a shared
operation scope, and an event is dropped when it is emitted inside one.
That matters a great deal for `nlmixr2extra`: its drivers
([`bootstrapFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/bootstrapFit.md),
[`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md),
[`covarSearchAuto()`](https://nlmixr2.github.io/nlmixr2extra/reference/covarSearchAuto.md),
…) run tens or hundreds of
[`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)
fits and
[`rxSolve()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
calls internally. Each driver enters the scope, so its internal fits and
solves are silent, and on the way out it emits **one** event that
describes the whole call. A logger then stores one entry per driver call
instead of one run per internal fit.

Everything here is a no-op when the installed `rxode2` has no event bus,
and nothing is emitted when no listener is registered.

## Which event a driver emits

On exit, a driver looks at what it returned and picks one event:

- the input fit itself, changed in place, gives **`fitUpdate`** – a new
  version of that fit’s run;
- a different fit gives **`fitComplete`** with `source` set to the
  driver’s kind – a new run linked to the input fit (`object`);
- anything else (a data frame, a list of results, a plot) gives one
  **`fitResult`** carrying a fit-free summary of the result.

When a result holds a best fit
([`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md)’s
`$best`), that fit is also announced with `fitComplete`, linked to the
input fit, so it can be kept as its own run.

| Driver | Event(s) | `what` / `source` / `kind` |
|:---|:---|:---|
| [`bootstrapFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/bootstrapFit.md) (default `returnType = "model"`) | `fitUpdate` | `bootstrap` |
| `bootstrapFit(returnType = "fitList"/"modelList")` | `fitResult` | `bootstrap` |
| [`preconditionFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/preconditionFit.md) | `fitUpdate` (always) | `covariance` |
| `multistart(fit)` | `fitResult`, then `fitComplete` for `$best` | `multistart` |
| `multistart(model, data)` | `fitComplete` for `$best`, then `fitResult` | `multistart` |
| [`covarSearchAuto()`](https://nlmixr2.github.io/nlmixr2extra/reference/covarSearchAuto.md) | `fitResult` | `covarSearch` |
| [`resSearch()`](https://nlmixr2.github.io/nlmixr2extra/reference/resSearch.md) | `fitResult` | `resSearch` |
| [`iivSearch()`](https://nlmixr2.github.io/nlmixr2extra/reference/iivSearch.md) | `fitResult` | `iivSearch` |
| [`rerunTopN()`](https://nlmixr2.github.io/nlmixr2extra/reference/rerunTopN.md) | nothing (refits are silent) | – |
| [`linearize()`](https://nlmixr2.github.io/nlmixr2extra/reference/linearize.md) | `fitComplete` | `linearize` |
| [`profile()`](https://rdrr.io/r/stats/profile.html) (fit method) | `fitResult` | `profile` |
| [`profileFixed()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileFixed.md) | `fitResult` | `profile` |
| [`profileLlp()`](https://nlmixr2.github.io/nlmixr2extra/reference/profileLlp.md) | `fitResult` | `profile` |
| [`lassoCoefficients()`](https://nlmixr2.github.io/nlmixr2extra/reference/lassoCoefficients.md) | `fitResult` | `lasso` |
| [`adaptivelassoCoefficients()`](https://nlmixr2.github.io/nlmixr2extra/reference/adaptivelassoCoefficients.md) | `fitResult` | `lasso` |
| [`adjustedlassoCoefficients()`](https://nlmixr2.github.io/nlmixr2extra/reference/adjustedlassoCoefficients.md) | `fitResult` | `lasso` |
| [`regularmodel()`](https://nlmixr2.github.io/nlmixr2extra/reference/regularmodel.md) | `fitComplete` | `lasso` |
| [`bootplot()`](https://nlmixr2.github.io/nlmixr2extra/reference/bootplot.md) | `fitResult` (the plot) | `bootplot` |

The kind is stored in `what` for `fitUpdate`, in `source` for
`fitComplete` and in `kind` for `fitResult`.

## Payloads

Listeners are called as `fun(event, ...)`: the event name, then the
payload fields as named arguments.

### `fitUpdate`

| Field      | Value                                                   |
|:-----------|:--------------------------------------------------------|
| `fit`      | the input fit, as changed by the driver                 |
| `original` | the same fit (it was changed in place)                  |
| `name`     | `NULL`                                                  |
| `what`     | the driver’s kind, e.g. `"bootstrap"` or `"covariance"` |
| `inPlace`  | `TRUE`                                                  |
| `call`     | the driver call                                         |

### `fitComplete`

| Field | Value |
|:---|:---|
| `fit` | the new fit ([`linearize()`](https://nlmixr2.github.io/nlmixr2extra/reference/linearize.md), [`regularmodel()`](https://nlmixr2.github.io/nlmixr2extra/reference/regularmodel.md), or a result’s best fit) |
| `object` | the input fit it came from; `NULL` when the driver started from a model |
| `call` | the driver call |
| `objName` | `NULL` |
| `source` | the driver’s kind, e.g. `"linearize"` |

### `fitResult`

| Field | Value |
|:---|:---|
| `fit` | the input fit (or the best fit when the driver started from a model) |
| `result` | a fit-free summary of what the driver returned (see below) |
| `kind` | the driver’s kind, e.g. `"profile"` |
| `call` | the driver call |

The `result` is built so that a logger can store it cheaply:

- every fit inside the result is replaced by
  `list(OBJF = <objective function>, parFixedDf = <parameter table>)`,
  read from the fit without computing anything;
- data frames and `ggplot` objects are kept as they are;
- functions and environments are dropped;
- lists are summarized up to four levels deep, and the original class of
  a classed list is kept in the `"nlmixr2extraClass"` attribute.

The `call` of every event is normalized by `rxode2`: its head is the
driver’s name, and large values inlined in it are replaced by markers,
so a recorded call never embeds a data set or a fit.

## What is not emitted

- **Internal work.** The fits, refits and solves a driver runs are
  inside its scope and emit nothing – not the bootstrap fits, not the
  multistart candidates, not the profile steps, not the covariate search
  candidates. This includes
  [`bootstrapFit()`](https://nlmixr2.github.io/nlmixr2extra/reference/bootstrapFit.md)
  when it is called by
  [`bootplot()`](https://nlmixr2.github.io/nlmixr2extra/reference/bootplot.md).
- **Nested drivers.** A driver called from inside another driver (or
  from inside any other emitting operation, such as
  [`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html))
  is itself silent; only the outermost call reports.
- **Errors.** When a driver fails, it leaves the scope and emits
  nothing.
- **Non-fit inputs.** A driver whose input is not a fit emits only when
  its result holds a best fit (as `multistart(model, data)` does). So
  [`rerunTopN()`](https://nlmixr2.github.io/nlmixr2extra/reference/rerunTopN.md),
  whose input is an
  [`iivSearch()`](https://nlmixr2.github.io/nlmixr2extra/reference/iivSearch.md)
  result, emits nothing; its refits are simply kept silent.
- **A best fit that is the input.** When a result’s best fit is the
  input fit itself, only the `fitResult` is emitted.
- **Unsummarizable results.** If summarizing the result fails, the
  driver leaves the scope without an event rather than breaking the
  call.

## Listening

A listener is a function `function(event, ...)`. This one records the
event name and the driver’s kind:

``` r

library(nlmixr2extra)

events <- list()
rxode2::rxEventListen("demo", function(event, ...) {
  p <- list(...)
  kind <- if (!is.null(p$what)) p$what else if (!is.null(p$source)) p$source else p$kind
  events[[length(events) + 1L]] <<- data.frame(event = event, kind = kind)
})
```

A fit to drive. It is run inside
[`rxode2::rxEventScope()`](https://nlmixr2.github.io/rxode2/reference/rxEventListen.html)
here so that only the drivers below are recorded (otherwise
[`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)
would report this fit with its own `fitComplete`, when `nlmixr2est` has
the event hooks):

``` r

one.compartment <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- log(31.5)
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

fit <- rxode2::rxEventScope(
  nlmixr2(one.compartment, nlmixr2data::theo_sd, est = "focei",
          control = list(print = 0))
)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> calculating covariance matrix
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> done
```

Now call a few drivers. Each runs many fits internally, but reports once
(or twice, for
[`multistart()`](https://nlmixr2.github.io/nlmixr2extra/reference/multistart.md)’s
best start):

``` r

preconditionFit(fit)
lin <- linearize(fit)
prof <- profile(fit, which = data.frame(tka = log(c(1.4, 1.6))), method = "fixed")
ms <- multistart(fit, control = list(n = 2, spread = 0.3, screen = "none",
                                     cacheDir = NA, refitBest = FALSE))
```

``` r

do.call(rbind, events)
#>         event       kind
#> 1   fitUpdate covariance
#> 2 fitComplete  linearize
#> 3   fitResult    profile
#> 4   fitResult multistart
#> 5 fitComplete multistart
```

The heavier drivers report the same way; for example, this records a
single `fitUpdate` (kind `bootstrap`) for the 100 bootstrap fits:

``` r

bootstrapFit(fit, nboot = 100)
```

When you are done, remove the listener:

``` r

rxode2::rxEventUnlisten("demo")
```

## Full logging

For a complete log of fits, simulations and driver results – with each
driver call stored as one entry linked to the fit it started from – use
[nlmixr2log](https://github.com/nlmixr2/nlmixr2log), which listens to
these events for you.
