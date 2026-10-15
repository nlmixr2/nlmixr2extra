#' Fit one binary dose-response model to counts per dose
#'
#' Doses are scaled by the largest dose (multiplicative only, so Emax and
#' exponential fits match DoseFinding).
#'
#' @param ep A `cuiBinary` endpoint.
#' @param x Events per dose.
#' @param n Patients per dose.
#' @param doses Sorted doses.
#' @return List with `p` (fitted probability per dose) and `fit` (model,
#'   coefficients on the scaled dose, scale) for prediction.
#' @noRd
.cuiDrModel <- function(ep, x, n, doses) {
  .s <- max(abs(doses))
  if (.s == 0) .s <- 1
  .d <- doses / .s
  .fit <- switch(ep$model,
    empirical = list(model = "empirical", p = x / n),
    logitLinear = .cuiDrLogit(x, n, .d, quadratic = FALSE, ep$monotone, ep$penalty),
    logitQuadratic = .cuiDrLogit(x, n, .d, quadratic = TRUE, ep$monotone, ep$penalty),
    emax = ,
    exponential = .cuiDrTwoStage(x, n, .d, ep$model, ep$monotone, ep$boundary, ep$backend),
    stop(sprintf("unknown dose-response model '%s'", ep$model), call. = FALSE)
  )
  .fit$scale <- .s
  .fit$doses <- doses
  if (is.null(.fit$p)) .fit$p <- .cuiDrPredict(.fit, doses)
  list(p = .fit$p, fit = .fit)
}

#' Predict a fitted dose-response model at doses
#' @noRd
.cuiDrPredict <- function(fit, dose) {
  if (fit$model == "empirical") {
    return(fit$p[match(dose, fit$doses)])
  }
  .d <- dose / fit$scale
  # logit models are fit on a centered dose
  if (!is.null(fit$center)) .d <- (.d - fit$center) / fit$spread
  .b <- fit$coef
  .eta <- switch(fit$model,
    logitLinear = .b[1] + .b[2] * .d,
    logitQuadratic = .b[1] + .b[2] * .d + .b[3] * .d^2,
    emax = .b[1] + .b[2] * .d / (.b[3] + .d),
    exponential = .b[1] + .b[2] * (exp(.d / .b[3]) - 1)
  )
  stats::plogis(.eta)
}

#' Binomial negative log-likelihood of logit-scale predictions
#' @noRd
.cuiDrNll <- function(eta, x, n) {
  -sum(x * stats::plogis(eta, log.p = TRUE) + (n - x) * stats::plogis(-eta, log.p = TRUE))
}

#' Logit-linear or logit-quadratic fit by maximum likelihood
#'
#' Monotone fits bound the slope (linear) or the derivative at both ends of
#' the dose range (quadratic: the derivative is linear in dose, so this
#' covers the whole range), allowing a flat fit.  `penalty = TRUE` uses the
#' CUI-MET formulation instead: an exponentiated slope, or a penalty on
#' derivative violations at the observed doses.
#' @noRd
.cuiDrLogit <- function(x, n, d, quadratic, monotone, penalty) {
  .sgn <- switch(monotone, none = 0, increasing = 1, decreasing = -1)
  # center and scale the dose to [-1, 1], so the coefficient bounds never
  # bind for doses far from 0; this keeps the sign of the slope
  .center <- (min(d) + max(d)) / 2
  .spread <- (max(d) - min(d)) / 2
  d <- (d - .center) / .spread
  .lo <- min(d)
  .hi <- max(d)
  .big <- 50
  # map optimizer parameters to (b0, b1[, b2])
  .coef <- function(th) {
    if (!quadratic) {
      if (.sgn != 0 && penalty) return(c(th[1], .sgn * exp(th[2])))
      return(th)
    }
    if (.sgn != 0 && !penalty) {
      # th = (b0, derivative at lo, derivative at hi), both of sign .sgn
      .b2 <- (th[3] - th[2]) / (2 * (.hi - .lo))
      return(c(th[1], th[2] - 2 * .b2 * .lo, .b2))
    }
    th
  }
  .eta <- function(b) b[1] + b[2] * d + if (quadratic) b[3] * d^2 else 0
  .obj <- function(th) {
    .b <- .coef(th)
    .v <- .cuiDrNll(.eta(.b), x, n)
    if (quadratic && .sgn != 0 && penalty) {
      .v <- .v + 1e5 * sum(pmax(-.sgn * (.b[2] + 2 * .b[3] * d), 0))
    }
    .v
  }
  .p0 <- min(max(sum(x) / sum(n), 0.01), 0.99)
  .np <- if (quadratic) 3L else 2L
  .start <- c(stats::qlogis(.p0), rep(0, .np - 1L))
  .lower <- rep(-.big, .np)
  .upper <- rep(.big, .np)
  if (.sgn != 0 && !penalty) {
    .idx <- if (quadratic) 2:3 else 2L
    if (.sgn > 0) .lower[.idx] <- 0 else .upper[.idx] <- 0
  }
  if (.sgn != 0 && penalty && !quadratic) {
    .start[2] <- -5
    .upper[2] <- log(.big)
  }
  if (quadratic && .sgn != 0 && penalty) {
    # the penalty has kinks, so search without derivatives from the exactly
    # constrained fit, which is feasible
    .b <- .cuiDrLogit(x, n, d, quadratic, monotone, penalty = FALSE)$coef
    .o <- stats::optim(.b, .obj, method = "Nelder-Mead", control = list(maxit = 5000))
    return(list(model = "logitQuadratic", coef = .o$par, center = .center, spread = .spread))
  }
  .o <- stats::optim(.start, .obj, method = "L-BFGS-B", lower = .lower, upper = .upper)
  if (.o$convergence > 1) {
    stop(sprintf("dose-response fit did not converge: %s", .o$message), call. = FALSE)
  }
  list(model = if (quadratic) "logitQuadratic" else "logitLinear", coef = .coef(.o$par),
       center = .center, spread = .spread)
}

#' Two-stage Emax or exponential fit (generalized MCP-Mod fitting step)
#'
#' Stage 1 estimates the log-odds per dose and their variances with a
#' no-intercept binomial GLM; stage 2 fits the curve by weighted least
#' squares on the logit scale.
#' @noRd
.cuiDrTwoStage <- function(x, n, d, model, monotone, boundary, backend) {
  .edge <- x == 0 | x == n
  if (any(.edge) && boundary == "error") {
    stop("a dose group has all or no events (boundary = \"error\")", call. = FALSE)
  }
  if (boundary == "correction") {
    .z <- log((x + 0.5) / (n - x + 0.5))
    .var <- 1 / (x + 0.5) + 1 / (n - x + 0.5)
  } else {
    .f <- factor(seq_along(d))
    .g <- suppressWarnings(stats::glm(cbind(x, n - x) ~ .f - 1, family = stats::binomial()))
    .z <- unname(stats::coef(.g))
    .var <- unname(diag(stats::vcov(.g)))
  }
  .bnds <- if (model == "emax") c(0.001, 1.5) else c(0.1, 2)
  if (backend == "DoseFinding") {
    if (!requireNamespace("DoseFinding", quietly = TRUE)) {
      stop("backend = \"DoseFinding\" needs the 'DoseFinding' package", call. = FALSE)
    }
    if (monotone != "none") {
      stop("backend = \"DoseFinding\" does not support 'monotone'", call. = FALSE)
    }
    .m <- DoseFinding::fitMod(d, .z, model = model, S = diag(.var, length(.var)),
                              type = "general", bnds = .bnds)
    .cf <- unname(stats::coef(.m))
    return(list(model = model, coef = .cf))
  }
  .w <- 1 / .var
  .sgn <- switch(monotone, none = 0, increasing = 1, decreasing = -1)
  .basis <- function(sh) if (model == "emax") d / (sh + d) else exp(d / sh) - 1
  # weighted least squares for (e0, effect) given the shape, with an
  # optional sign constraint on the effect
  .wls <- function(sh) {
    .b <- .basis(sh)
    .x <- cbind(1, .b)
    .cf <- tryCatch(solve(crossprod(.x, .w * .x), crossprod(.x, .w * .z))[, 1],
                    error = function(e) c(sum(.w * .z) / sum(.w), 0))
    if (.sgn != 0 && .sgn * .cf[2] < 0) .cf <- c(sum(.w * .z) / sum(.w), 0)
    .r <- .z - (.cf[1] + .cf[2] * .b)
    list(coef = unname(.cf), rss = sum(.w * .r^2))
  }
  # grid then refine, as the shape is often weakly identified
  .grid <- seq(.bnds[1], .bnds[2], length.out = 31)
  .rss <- vapply(.grid, function(sh) .wls(sh)$rss, numeric(1))
  .i <- which.min(.rss)
  .int <- c(.grid[max(1, .i - 1)], .grid[min(length(.grid), .i + 1)])
  .o <- stats::optimize(function(sh) .wls(sh)$rss, .int)
  .sh <- if (.o$objective <= .rss[.i]) .o$minimum else .grid[.i]
  list(model = model, coef = c(.wls(.sh)$coef, .sh))
}
