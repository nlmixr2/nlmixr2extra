# General profiling functions ----

#' Perform likelihood profiling on nlmixr2 focei fits
#'
#'
#' @details
#'
#' # Log-likelihood profiling
#'
#' `method = "llp"`
#'
#' The search will stop when either the OFV is within `ofvtol` of the desired
#' OFV change or when the parameter is interpolating to more significant digits
#' than specified in `paramDigits`.  The "llp" method uses the `profileLlp()`
#' function.  See its help for more details.
#'
#' # Fixed points
#'
#' `method = "fixed"`
#'
#' Estimate the OFV for specific fixed values.  The "fixed" method uses the
#' `profileFixed()` function.  See its help for more details.
#'
#' @param fitted The fit model
#' @param ... ignored
#' @param which The parameter names to perform likelihood profiling on
#'   (`NULL` indicates all parameters)
#' @param method Method to use for profiling (see the details)
#' @param control Control arguments for the `method`
#' @return A data.frame with one column named `Parameter` indicating the
#'   parameter being fixed on that row, one column for the `OFV` indicating the
#'   OFV estimated for the model at that step, one column named `profileBound`
#'   indicating the estimated value for the profile likelihood and its step
#'   above the minimum profile likelihood value, and columns for each parameter
#'   estimate (or fixed) in the model.
#' @importFrom stats coef confint lm
#' @family Profiling
#' @examples
#' \dontrun{
#' # Likelihood profiling takes a long time to run each model multiple times, so
#' # be aware that running this example may take a few minutes.
#' oneCmt <- function() {
#'   ini({
#'     tka <- log(1.57)
#'     tcl <- log(2.72)
#'     tv <- fixed(log(31.5))
#'     eta.ka ~ 0.6
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl)
#'     v <- exp(tv)
#'     cp <- linCmt()
#'     cp ~ add(add.sd)
#'   })
#' }
#'
#' fit <-
#'   nlmixr2(
#'     oneCmt, data = nlmixr2data::theo_sd, est="focei", control = list(print=0)
#'   )
#' # profile all parameters
#' profall <- profile(fit)
#'
#' # profile a single parameter
#' proftka <- profile(fit, which = "tka")
#' }
#' @export
profile.nlmixr2FitCore <- function(
  fitted,
  ...,
  which = NULL,
  method = c("llp", "fixed"),
  control = list()
) {
  method <- match.arg(method)

  if (method == "llp") {
    lifecycle::deprecate_warn(
      when = "5.1.0",
      what = "profile(method = 'llp')",
      with = "runLLP()"
    )
    runLLP(fit = fitted, which = which, control = control, ...)
  } else if (method == "fixed") {
    profileFixed(fitted = fitted, which = which, control = control)
  } else {
    stop("Invalid 'method': ", method) # nocov
  }
}

#' Give the output data.frame for a single model for profile.nlmixr2FitCore
#'
#' @inheritParams profile.nlmixr2FitCore
#' @param fixedVal The value that `which` is fixed to in case the model does not
#'   converge.
#' @return A data.frame with columns named "Parameter" (the parameter name(s)
#'   that were fixed), OFV (the objective function value), and the current
#'   estimate for each of the parameters.  Omega values are given as their
#'   variances and covariances.
#' @family Profiling
profileNlmixr2FitCoreRet <- function(fitted, which, fixedVal) {
  if (inherits(fitted, "try-error")) {
    ret <- data.frame(Parameter = which, OFV = NA_real_, X = fixedVal)
    names(ret)[3] <- which
  } else {
    if (any(names(nlmixr2est::fixef(fitted)) %in% c("Parameter", "OFV"))) {
      cli::cli_abort(
        "Cannot profile a model with a parameter named either 'Parameter' or 'OFV'"
      )
    }
    retOmega <- as.data.frame(t(diag(fitted$omega)))
    for (idxRow in seq_len(nrow(fitted$omega)) - 1) {
      for (idxCol in seq_len(idxRow)) {
        # Only capture the columns when the covariance is nonzero
        if (fitted$omega[idxRow + 1, idxCol] != 0) {
          covColName <-
            sprintf(
              "cov(%s,%s)",
              rownames(fitted$omega)[idxRow + 1],
              colnames(fitted$omega)[idxCol]
            )
          retOmega[[covColName]] <- fitted$omega[idxRow + 1, idxCol]
        }
      }
    }
    ret <-
      cbind(
        data.frame(Parameter = which, OFV = fitted$objective),
        data.frame(t(nlmixr2est::fixef(fitted))),
        retOmega
      )
  }
  rownames(ret) <- NULL
  ret
}

# Fixed parameter estimate profiling ----

#' Estimate the objective function values for a model while fixing defined
#' parameter values
#'
#' @inheritParams profile.nlmixr2FitCore
#' @param which A data.frame with column names of parameters to fix and values
#'   of the fitted value to fix (one row per set of parameters to estimate)
#' @param control A list passed to `fixedControl()` (currently unused)
#' @inherit profileNlmixr2FitCoreRet return
#' @author Bill Denney (changed by Matt Fidler to take out R 4.1 specific code)
#' @family Profiling
#' @export
profileFixed <- function(fitted, which, control = list()) {
  control <- do.call(fixedControl, control)
  checkmate::assert_data_frame(
    which,
    types = "numeric",
    any.missing = FALSE,
    min.rows = 1
  )
  dplyr::bind_rows(lapply(
    X = seq_len(nrow(which)),
    FUN = function(idx, fitted) {
      profileFixedSingle(fitted = fitted, which = which[idx, , drop = FALSE])
    },
    fitted = fitted
  ))
}

#' @describeIn profileFixed Estimate the objective function value for a model
#'   while fixing a single set of defined parameter values (for use in parameter
#'   profiling)
#'
#' @param which A data.frame with column names of parameters to fix and values
#'   of the fitted value to fix (one row only).
#' @returns `which` with a column named `OFV` added with the objective function
#'   value of the fitted estimate fixing the parameters in the other columns
#' @family Profiling
#' @author Bill Denney (changed by Matt Fidler to take out R 4.1 specific code)
#' @export
profileFixedSingle <- function(fitted, which) {
  checkmate::assert_data_frame(
    which,
    types = "numeric",
    any.missing = FALSE,
    nrows = 1
  )
  checkmate::assert_names(
    names(which),
    subset.of = names(nlmixr2est::fixef(fitted))
  )

  cli::cli_inform(
    "Profiling {paste(names(which), '=', unlist(which), collapse = ', ')}"
  )
  # Update the model by fixing all of the parameters
  paramToFix <- lapply(X = which, FUN = function(x) {
    str2lang(sprintf("fixed(%s)", x))
  })
  iniArgs <- append(list(x = fitted), paramToFix)
  modelToFit <- suppressMessages(do.call(rxode2::ini, iniArgs))
  controlFit <- fitted$control
  newFit <-
    try(suppressMessages(
      nlmixr2est::nlmixr2(
        modelToFit,
        est = fitted$est,
        control = setQuietFastControl(fitted$control)
      )
    ))
  ret <- profileNlmixr2FitCoreRet(
    fitted = newFit,
    which = paste(names(which), collapse = ",")
  )
  cli::cli_inform("OFV: {ret$OFV}")
  ret
}

#' Control options for fixed-value likelihood profiling
#'
#' @returns A validated list of control options for fixed-value likelihood
#'   profiling
#' @family Profiling
#' @seealso [profileFixed()]
#' @export
fixedControl <- function() {
  ret <- list()
  class(ret) <- "fixedControl"
  ret
}

# Log-likelihood profiling ----

# Extract SE for all v1-profileable parameters (fixef names) from parFixedDf.
# Returns a named numeric vector; NA where SE is unavailable (fixed params,
# residual-error thetas without a covariance step, or params absent from the df).
llpExtractSE <- function(fit) {
  params <- names(nlmixr2est::fixef(fit))
  parDf <- fit$parFixedDf
  se <- setNames(rep(NA_real_, length(params)), params)
  inDf <- params[params %in% rownames(parDf)]
  if (length(inDf) > 0) {
    se[inDf] <- parDf[inDf, "SE"]
  }
  se
}

#' Profile confidence intervals with log-likelihood profiling
#'
#' @inheritParams profile.nlmixr2FitCore
#' @param control A list passed to `llpControl()`
#' @param which Either `NULL` to profile all parameters or a character vector of
#'   parameters to estimate
#' @return A data.frame with columns named "Parameter" (the parameter name(s)
#'   that were fixed), OFV (the objective function value), and the current
#'   estimate for each of the parameters.  In addition, if any boundary is
#'   found, the OFV increase will be indicated by the absolute value of the
#'   "profileBound" column and if that boundary is the upper or lower boundary
#'   will be indicated by the "profileBound" column being positive or negative,
#'   respectively.
#' @family Profiling
#' @export
profileLlp <- function(fitted, which, control) {
  control <- do.call(llpControl, control)

  if (is.null(which)) {
    which <- names(nlmixr2est::fixef(fitted))
  } else {
    checkmate::assert_subset(
      x = which,
      choices = names(nlmixr2est::fixef(fitted)),
      empty.ok = FALSE
    )
  }

  if (!is.null(names(control$rseTheta))) {
    checkmate::assert_names(
      names(control$rseTheta),
      subset.of = names(nlmixr2est::fixef(fitted))
    )
    checkmate::assert_numeric(
      control$rseTheta,
      lower = 0,
      any.missing = FALSE,
      min.len = 1
    )
  }

  results <- lapply(which, function(w) {
    llpRunOneParameter(fitted = fitted, which = w, control = control)
  })

  raw <- dplyr::bind_rows(lapply(results, `[[`, "raw"))
  rownames(raw) <- NULL
  raw
}

# Side-aware next-value selection for LLP using quadratic interpolation with linear fallback.
#
# paramVals/ofvVals: all tried values so far (may contain NA in ofvVals).
# direction: +1 = upper search, -1 = lower search.
# searchLower/searchUpper: the active one-sided bracket (not used for final clamping).
# hardLower/hardUpper: actual model parameter bounds (used for clamping only).
# Returns the next parameter value to try, or NA_real_ to signal abort.
llpNextParam <- function(
  paramVals,
  ofvVals,
  targetOfv,
  direction,
  mle,
  searchLower,
  searchUpper,
  hardLower,
  hardUpper
) {
  # Remove NAs and restrict to active half-space
  keep <- !is.na(ofvVals) &
    (if (direction < 0) paramVals <= mle else paramVals >= mle)
  x <- paramVals[keep]
  y <- ofvVals[keep]
  if (length(x) < 2L) {
    return(NA_real_)
  }

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  n <- length(x)

  result <- NA_real_

  # PsN-style 3-point window: lower uses 3 leftmost, upper uses 3 rightmost
  if (n >= 3L) {
    wi <- if (direction < 0) seq_len(3L) else seq.int(n - 2L, n)
    wx <- x[wi]
    wy <- y[wi]
    fit <- tryCatch(
      suppressWarnings(lm(wy ~ wx + I(wx^2))),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      coefs <- coef(fit)
      c0 <- unname(coefs[1L])
      c1 <- unname(coefs[2L])
      c2 <- unname(coefs[3L])
      if (!is.finite(c2) || abs(c2) < 1e-10) {
        cli::cli_warn(
          "Flat quadratic profile; falling back to linear interpolation/extrapolation."
        )
      } else {
        disc <- c1^2 - 4 * c2 * (c0 - targetOfv)
        if (is.finite(disc) && disc >= 0) {
          sq <- sqrt(disc)
          roots <- c((-c1 + sq) / (2 * c2), (-c1 - sq) / (2 * c2))
          outer_pt <- if (direction < 0) x[1L] else x[n]
          hs_ok <- if (direction < 0) roots <= mle else roots >= mle
          valid_r <- roots[
            is.finite(roots) &
              roots >= searchLower &
              roots <= searchUpper &
              hs_ok
          ]
          if (length(valid_r) > 0L) {
            result <- valid_r[which.min(abs(valid_r - outer_pt))]
          }
        }
        # negative/non-finite discriminant -> fall through to linear
      }
    }
  }

  # Linear fallback
  if (is.na(result)) {
    if (targetOfv >= min(y) && targetOfv <= max(y)) {
      # Interpolate (no rule = 2 -- never clamp to endpoint values)
      result <- stats::approx(x = y, y = x, xout = targetOfv)$y
    } else {
      # Explicit extrapolation from the two outermost evaluated points in direction
      if (direction < 0) {
        xo <- x[1L]
        yo <- y[1L]
        xi <- x[2L]
        yi <- y[2L]
      } else {
        xo <- x[n]
        yo <- y[n]
        xi <- x[n - 1L]
        yi <- y[n - 1L]
      }
      dxdy <- (xo - xi) / (yo - yi)
      if (!is.finite(dxdy) || dxdy == 0) {
        return(NA_real_)
      }
      result <- xo + (targetOfv - yo) * dxdy
    }
  }

  if (!is.finite(result)) {
    return(NA_real_)
  }

  # Clamp only to hard model bounds (not search brackets)
  margin <- sqrt(.Machine$double.eps)
  min(max(result, hardLower + margin), hardUpper - margin)
}

# Single-parameter, single-direction optimization of likelihood profiling.
# Returns list(raw = <data.frame>, status = <1-row data.frame>).
# searchLower/searchUpper: the active one-sided bracket (MLE on one side).
# hardLower/hardUpper: actual model parameter bounds from iniDf (used for
# near-boundary detection and clamping inside llpNextParam).
optimProfile <- function(
  par,
  fitted,
  optimDf,
  which,
  ofvIncrease,
  direction,
  searchLower,
  searchUpper,
  hardLower,
  hardUpper,
  mle,
  fixFn = NULL,
  itermax = 10,
  ofvtol = 0.005,
  paramDigits = 3
) {
  currentOFVDiff <- Inf
  currentIter <- 0L
  currentPar <- par
  origMinOFV <- min(optimDf$OFV, na.rm = TRUE)
  ret <- optimDf[
    optimDf$Parameter %in%
      which &
      optimDf[[which]] >= searchLower &
      optimDf[[which]] <= searchUpper,
    ,
    drop = FALSE
  ]

  converged <- FALSE
  nearBound <- FALSE
  hitMaxIter <- FALSE
  aborted <- FALSE
  boundEstimate <- NA_real_
  nEval <- 0L
  lastMsg <- ""
  margin <- sqrt(.Machine$double.eps)

  while (itermax > currentIter) {
    currentIter <- currentIter + 1L
    if (is.null(fixFn)) {
      dfWhich <- stats::setNames(data.frame(X = currentPar), nm = which)
      fitResult <- profileFixedSingle(fitted = fitted, which = dfWhich)
    } else {
      fitResult <- fixFn(currentPar)
    }
    nEval <- nEval + 1L

    ret <- dplyr::bind_rows(ret, fitResult)
    ret <- ret[order(ret$OFV), , drop = FALSE]
    currentMinOFV <- min(ret$OFV, na.rm = TRUE)
    if (isTRUE(any(currentMinOFV < origMinOFV))) {
      cli::cli_warn(
        "OFV decreased while profiling {which}; profile may be inaccurate"
      )
      origMinOFV <- NA
    }
    targetOfv <- currentMinOFV + ofvIncrease

    nextPar <- llpNextParam(
      paramVals = ret[[which]],
      ofvVals = ret$OFV,
      targetOfv = targetOfv,
      direction = direction,
      mle = mle,
      searchLower = searchLower,
      searchUpper = searchUpper,
      hardLower = hardLower,
      hardUpper = hardUpper
    )

    if (is.na(nextPar)) {
      aborted <- TRUE
      lastMsg <- "aborted due to numerical difficulties"
      messageProfileComplete(which, direction = direction, lastMsg)
      break
    }

    # Near-boundary: only check hard model bounds, not the MLE-side search bracket
    if (nextPar <= hardLower + margin || nextPar >= hardUpper - margin) {
      nearBound <- TRUE
    }

    if (nextPar %in% ret[[which]]) {
      aborted <- TRUE
      lastMsg <- "aborted due to attempted repeated parameter value estimation"
      messageProfileComplete(which, direction = direction, lastMsg)
      break
    }

    signifSame <- signif(ret[[which]], digits = paramDigits) ==
      signif(nextPar, digits = paramDigits)
    convergedSignif <- any(signifSame & ret[[which]] > nextPar) &&
      any(signifSame & ret[[which]] < nextPar)
    currentOFVDiff <- min(abs(targetOfv - ret$OFV), na.rm = TRUE)
    convergedOfv <- currentOFVDiff <= ofvtol
    converged <- convergedSignif | convergedOfv

    if (converged) {
      if (convergedSignif & convergedOfv) {
        lastMsg <- "complete due to significant digits precision and achieving OFV within specified tolerance"
      } else if (convergedSignif) {
        lastMsg <- "complete due to significant digits precision"
      } else {
        lastMsg <- "complete due to achieving OFV within specified tolerance"
      }
      messageProfileComplete(which, direction = direction, lastMsg)
      boundEstimate <- nextPar
      storeBound <- data.frame(
        Parameter = which,
        X = nextPar,
        profileBound = direction * ofvIncrease
      )
      names(storeBound)[2L] <- which
      ret <- dplyr::bind_rows(storeBound, ret)
      break
    }

    cli::cli_inform(
      "Profiling {which} best distance to target OFV: {round(currentOFVDiff, 4)}"
    )
    currentPar <- nextPar
  }

  if (!converged && !aborted) {
    hitMaxIter <- TRUE
    lastMsg <- "aborted due to too many iterations without convergence"
    messageProfileComplete(which, direction = direction, lastMsg)
  }

  ret <- ret[order(ret[[which]]), , drop = FALSE]
  rownames(ret) <- NULL

  status <- data.frame(
    param = which,
    direction = direction,
    converged = converged,
    nearBound = nearBound,
    hitMaxIter = hitMaxIter,
    message = lastMsg,
    boundEstimate = boundEstimate,
    bestOfvDiff = currentOFVDiff,
    nEval = nEval,
    stringsAsFactors = FALSE
  )

  list(raw = ret, status = status)
}

# Run both profile directions for a single parameter, compute the interval
# ratio, and return list(raw, status).  This is the unit of parallel work in
# Step 7 (parallelisation over parameters).
llpRunOneParameter <- function(fitted, which, control) {
  mle <- nlmixr2est::fixef(fitted)[[which]]
  effectRange <- fitted$iniDf[
    fitted$iniDf$name == which,
    c("lower", "upper"),
    drop = FALSE
  ]
  hardLower <- effectRange$lower
  hardUpper <- effectRange$upper

  seTheta <- llpExtractSE(fitted)

  rseForParam <- control$rseTheta
  if (length(rseForParam) == 1L && is.null(names(rseForParam))) {
    rseForParam <- stats::setNames(rseForParam, which)
  }

  initialDf <- profileNlmixr2FitCoreRet(fitted = fitted, which = which)
  estInitial <- profileNlmixr2FitDataEstInitial(
    estimates = initialDf,
    which = which,
    normq = control$normq,
    rseTheta = rseForParam,
    seTheta = seTheta,
    hardLower = hardLower,
    hardUpper = hardUpper
  )

  # Residual-error-as-theta parameters are absent from fit$cov; use wider threshold
  isResidualError <- tryCatch(
    {
      cov <- fitted$cov
      !is.null(cov) && !which %in% rownames(cov)
    },
    error = function(e) FALSE
  )

  allRaw <- initialDf
  statusList <- list()

  for (direction in c(-1L, 1L)) {
    if (direction == -1L) {
      searchLower <- hardLower
      searchUpper <- mle
      initPar <- estInitial[1L]
    } else {
      searchLower <- mle
      searchUpper <- hardUpper
      initPar <- estInitial[2L]
    }

    res <- optimProfile(
      par = initPar,
      fitted = fitted,
      optimDf = initialDf,
      which = which,
      ofvIncrease = control$ofvIncrease,
      direction = direction,
      searchLower = searchLower,
      searchUpper = searchUpper,
      hardLower = hardLower,
      hardUpper = hardUpper,
      mle = mle,
      itermax = control$itermax,
      ofvtol = control$ofvtol,
      paramDigits = control$paramDigits
    )

    allRaw <- unique(dplyr::bind_rows(allRaw, res$raw))
    statusList[[length(statusList) + 1L]] <- res$status
  }

  status <- dplyr::bind_rows(statusList)

  # Interval ratio warning: flag asymmetric profiles
  lBound <- status$boundEstimate[status$direction == -1L]
  uBound <- status$boundEstimate[status$direction == 1L]
  if (
    length(lBound) == 1L &&
      length(uBound) == 1L &&
      !is.na(lBound) &&
      !is.na(uBound)
  ) {
    lDist <- abs(mle - lBound)
    uDist <- abs(uBound - mle)
    if (lDist > 0 && uDist > 0) {
      ratio <- uDist / lDist
      threshold <- if (isTRUE(isResidualError)) 1.6 else 1.3
      if (ratio > threshold || ratio < 1 / threshold) {
        cli::cli_warn(
          "Asymmetric profile for {.val {which}}: interval ratio {round(ratio, 2)} (threshold {threshold})"
        )
      }
    }
  }

  list(raw = allRaw, status = status)
}

# Omega diagonal profiling (Step O) ----

# Returns names of omega diagonal elements that are not already fixed.
llpOmegaNames <- function(fit) {
  idf <- fit$iniDf
  omega_rows <- idf[
    is.na(idf$ntheta) & !is.na(idf$neta1) & idf$neta1 == idf$neta2,
  ]
  omega_rows$name[!isTRUE(omega_rows$fix)]
}

# Wishart chi-squared approximation for omega diagonal SE.
# Returns a named numeric vector with one entry per profileable omega.
llpOmegaSE <- function(fit) {
  omNames <- llpOmegaNames(fit)
  idf <- fit$iniDf
  n_sub <- tryCatch(.sirNSubjects(fit), error = function(e) NA_integer_)
  se <- setNames(rep(NA_real_, length(omNames)), omNames)
  if (!is.na(n_sub) && n_sub > 1L) {
    for (nm in omNames) {
      omega_kk <- idf[idf$name == nm, "est"]
      se[[nm]] <- sqrt(2 * omega_kk^2 / (n_sub - 1L))
    }
  }
  se
}

# Build a model function (as a parsed R function) that fixes one omega diagonal
# at fixOmegaVal, starting all other parameters from the fit's current estimates.
# Uses eval(parse(text=...)) with eta.X ~ fix(VALUE) so that fix=TRUE is set.
buildFixedOmegaModel <- function(fit, omegaName, fixOmegaVal) {
  idf <- fit$iniDf
  ui <- fit$ui
  thetaRows <- idf[!is.na(idf$ntheta), ]
  omegaRows <- idf[
    is.na(idf$ntheta) & !is.na(idf$neta1) & idf$neta1 == idf$neta2,
  ]

  iniLines <- character(0L)
  for (i in seq_len(nrow(thetaRows))) {
    row <- thetaRows[i, ]
    iniLines <- c(
      iniLines,
      if (isTRUE(row$fix)) {
        sprintf("    %s <- fixed(%s)", row$name, deparse(row$est))
      } else {
        sprintf("    %s <- %s", row$name, deparse(row$est))
      }
    )
  }
  for (i in seq_len(nrow(omegaRows))) {
    row <- omegaRows[i, ]
    nm <- row$name
    iniLines <- c(
      iniLines,
      if (nm == omegaName) {
        sprintf("    %s ~ fix(%s)", nm, deparse(fixOmegaVal))
      } else if (isTRUE(row$fix)) {
        sprintf("    %s ~ fix(%s)", nm, deparse(row$est))
      } else {
        sprintf("    %s ~ %s", nm, deparse(row$est))
      }
    )
  }
  modelLines <- vapply(
    ui$lstExpr,
    function(e) paste0("    ", deparse(e)),
    character(1L)
  )
  fnText <- paste0(
    "function() {\n  ini({\n",
    paste(iniLines, collapse = "\n"),
    "\n  })\n  model({\n",
    paste(modelLines, collapse = "\n"),
    "\n  })\n}"
  )
  eval(parse(text = fnText))
}

# Fix one omega diagonal, re-estimate, return a data.frame row (same format as
# profileFixedSingle / profileNlmixr2FitCoreRet).
profileFixedOmegaSingle <- function(fitted, omegaName, fixOmegaVal, fitData) {
  cli::cli_inform("Profiling {omegaName} = {fixOmegaVal}")
  modelFn <- buildFixedOmegaModel(fitted, omegaName, fixOmegaVal)
  newFit <- try(suppressMessages(
    nlmixr2est::nlmixr2(
      modelFn,
      data = fitData,
      est = fitted$est,
      control = setQuietFastControl(fitted$control)
    )
  ))
  ret <- profileNlmixr2FitCoreRet(
    fitted = newFit,
    which = omegaName,
    fixedVal = fixOmegaVal
  )
  cli::cli_inform("OFV: {ret$OFV}")
  ret
}

# Orchestrate both profile directions for one omega diagonal parameter.
llpRunOneOmega <- function(fitted, which, control) {
  idfRow <- fitted$iniDf[fitted$iniDf$name == which, , drop = FALSE]
  mle <- idfRow$est
  hardLower <- max(0, idfRow$lower) # omega variance is non-negative
  hardUpper <- idfRow$upper

  seOmega <- llpOmegaSE(fitted)
  se <- seOmega[[which]]
  rse <- control$rseTheta
  if (!is.null(names(rse))) {
    rse <- rse[[which]]
  }
  if (is.na(se)) {
    width <- control$normq * (rse / 100) * abs(mle)
  } else {
    width <- control$normq * se
  }
  margin <- sqrt(.Machine$double.eps)
  lowerGuess <- mle - width
  upperGuess <- mle + width
  if (lowerGuess <= hardLower) {
    lowerGuess <- (hardLower - mle) * 0.9 + mle
  }
  if (upperGuess >= hardUpper) {
    upperGuess <- (hardUpper - mle) * 0.9 + mle
  }
  lowerGuess <- max(lowerGuess, hardLower + margin)
  upperGuess <- min(upperGuess, hardUpper - margin)
  estInitial <- c(lowerGuess, upperGuess)

  fitData <- nlmixr2est::getData(fitted)
  fixFn <- function(fixOmegaVal) {
    profileFixedOmegaSingle(
      fitted = fitted,
      omegaName = which,
      fixOmegaVal = fixOmegaVal,
      fitData = fitData
    )
  }

  initialDf <- profileNlmixr2FitCoreRet(fitted = fitted, which = which)
  allRaw <- initialDf
  statusList <- list()

  for (direction in c(-1L, 1L)) {
    searchLower <- if (direction == -1L) hardLower else mle
    searchUpper <- if (direction == -1L) mle else hardUpper
    initPar <- estInitial[if (direction == -1L) 1L else 2L]

    res <- optimProfile(
      par = initPar,
      fitted = fitted,
      optimDf = initialDf,
      which = which,
      ofvIncrease = control$ofvIncrease,
      direction = direction,
      searchLower = searchLower,
      searchUpper = searchUpper,
      hardLower = hardLower,
      hardUpper = hardUpper,
      mle = mle,
      fixFn = fixFn,
      itermax = control$itermax,
      ofvtol = control$ofvtol,
      paramDigits = control$paramDigits
    )

    allRaw <- unique(dplyr::bind_rows(allRaw, res$raw))
    statusList[[length(statusList) + 1L]] <- res$status
  }

  status <- dplyr::bind_rows(statusList)

  lBound <- status$boundEstimate[status$direction == -1L]
  uBound <- status$boundEstimate[status$direction == 1L]
  if (
    length(lBound) == 1L &&
      length(uBound) == 1L &&
      !is.na(lBound) &&
      !is.na(uBound)
  ) {
    lDist <- abs(mle - lBound)
    uDist <- abs(uBound - mle)
    if (lDist > 0 && uDist > 0) {
      ratio <- uDist / lDist
      if (ratio > 1.6 || ratio < 1 / 1.6) {
        cli::cli_warn(
          "Asymmetric profile for {.val {which}}: interval ratio {round(ratio, 2)} (threshold 1.6)"
        )
      }
    }
  }

  list(raw = allRaw, status = status)
}

#' Control options for log-likelihood profiling
#'
#' @param ofvIncrease The targetted change in objective function value (3.84
#'   corresponds to a Chi-squared test with a 95% confidence interval)
#' @param normq Normal quantile used to compute the initial-guess width for each
#'   profile direction: `width = normq * SE` when SE is available, or
#'   `width = normq * rseTheta/100 * |MLE|` as a fallback.  Defaults to 1.96
#'   (PsN convention).
#' @param rseTheta The relative standard error (percent) used as a fallback
#'   initial-guess width when SE is unavailable.  Either a single unnamed
#'   non-negative number applied to all parameters (default 30), or a named
#'   non-negative numeric vector with one entry per parameter.  Name validity
#'   against the fit is checked at profiling time, not here.
#' @param itermax Maximum number of likelihood profiling iterations for each
#'   bound estimated
#' @param ofvtol The relative tolerance for the objective function being close
#'   enough to the `ofvIncrease`.
#' @param paramDigits The number of significant digits required for the
#'   parameter.  When interpolation attempts to get smaller than that number of
#'   significant digits, it will stop.
#' @param extrapolateExpand When extrapolating outside the range previously
#'   tested, how far should the step occur as a ratio
#' @param workers Number of parallel workers for profiling multiple parameters
#'   simultaneously.  `NULL` (default) or `1` runs sequentially without
#'   requiring any optional packages.  A positive integer greater than 1, or
#'   `"auto"`, uses `future` and `future.apply`; both packages must be installed
#'   or profiling aborts with a clear error.
#' @returns A validated list of control options for log-likelihood profiling
#' @family Profiling
#' @seealso [profileLlp()]
#' @export
llpControl <- function(
  ofvIncrease = qchisq(0.95, df = 1),
  normq = 1.96,
  rseTheta = 30,
  itermax = 10,
  ofvtol = 0.005,
  paramDigits = 3,
  extrapolateExpand = 1.5,
  workers = NULL
) {
  checkmate::assert_numeric(
    rseTheta,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    min.len = 1
  )
  isUnnamedScalar <- length(rseTheta) == 1 && is.null(names(rseTheta))
  isNamedVector <- length(rseTheta) >= 1 && !is.null(names(rseTheta))
  if (!isUnnamedScalar && !isNamedVector) {
    cli::cli_abort(
      "{.arg rseTheta} must be a single unnamed number or a named numeric vector"
    )
  }

  if (!is.null(workers) && !identical(workers, "auto")) {
    checkmate::assert_integerish(
      workers,
      lower = 1,
      any.missing = FALSE,
      len = 1,
      null.ok = FALSE,
      coerce = TRUE
    )
    workers <- as.integer(workers)
  }

  ret <-
    list(
      ofvIncrease = checkmate::assert_number(
        ofvIncrease,
        lower = 0,
        finite = TRUE,
        null.ok = FALSE,
        na.ok = FALSE
      ),
      normq = checkmate::assert_number(
        normq,
        lower = 0,
        finite = TRUE,
        null.ok = FALSE,
        na.ok = FALSE
      ),
      rseTheta = rseTheta,
      ofvtol = checkmate::assert_number(
        ofvtol,
        na.ok = FALSE,
        lower = 1e-10,
        upper = 1,
        finite = TRUE,
        null.ok = FALSE
      ),
      itermax = checkmate::assert_integerish(
        itermax,
        lower = 1,
        any.missing = FALSE,
        len = 1,
        null.ok = FALSE,
        coerce = TRUE
      ),
      paramDigits = checkmate::assert_integerish(
        paramDigits,
        lower = 1,
        upper = 10,
        any.missing = FALSE,
        len = 1,
        null.ok = FALSE,
        coerce = TRUE
      ),
      extrapolateExpand = checkmate::assert_number(
        extrapolateExpand,
        lower = 1.001,
        na.ok = FALSE,
        finite = TRUE,
        null.ok = FALSE
      ),
      workers = workers
    )
  class(ret) <- "llpControl"
  ret
}

#' @export
rxUiDeparse.llpControl <- function(object, var) {
  .default <- llpControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}


# Provide the initial estimates going up and down from the initial value.
# Uses normq * SE when seTheta[which] is available (PsN convention); falls back
# to normq * rseTheta/100 * |MLE| otherwise. Out-of-bounds guesses are pulled
# 90% of the way back toward the MLE (PsN soft clamp) rather than hard-clamped.
profileNlmixr2FitDataEstInitial <- function(
  estimates,
  which,
  normq,
  rseTheta,
  seTheta,
  hardLower,
  hardUpper
) {
  checkmate::assert_data_frame(estimates, nrows = 1)
  checkmate::assert_character(which, any.missing = FALSE, len = 1)
  checkmate::assert_subset(which, choices = names(estimates))

  mle <- estimates[[which]]

  if (!is.na(seTheta[which])) {
    width <- normq * unname(seTheta[which])
  } else {
    currentRseTheta <- if (which %in% names(rseTheta)) rseTheta[[which]] else 30
    width <- normq * currentRseTheta / 100 * abs(mle)
  }

  margin <- sqrt(.Machine$double.eps)

  lowerGuess <- mle - width
  if (lowerGuess < hardLower + margin) {
    lowerGuess <- (hardLower - mle) * 0.9 + mle
  }

  upperGuess <- mle + width
  if (upperGuess > hardUpper - margin) {
    upperGuess <- (hardUpper - mle) * 0.9 + mle
  }

  c(lowerGuess, upperGuess)
}

# Parallelization helpers ----

# Creates a unique per-parameter temp directory, runs fun() inside it, then
# restores the original working directory and removes the temp directory.
# This prevents file and compiled-artifact collisions between parallel workers.
llpWithIsolatedFitDir <- function(param, fun) {
  dir <- file.path(
    tempdir(),
    sprintf("llp_%s_%d_%s", param, Sys.getpid(), format(Sys.time(), "%H%M%OS6"))
  )
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  on.exit(
    {
      setwd(old_wd)
      unlink(dir, recursive = TRUE)
    },
    add = TRUE
  )
  setwd(dir)
  fun()
}

# Aborts with a clear installation message if future or future.apply are missing.
# Called before any parallel code path so the user gets an explicit error rather
# than a silent fallback.
.llpCheckParallelDeps <- function() {
  missing_pkgs <- character(0L)
  if (!requireNamespace("future", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "future")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, "future.apply")
  }
  if (length(missing_pkgs) > 0L) {
    cli::cli_abort(c(
      "Parallel LLP requires {.pkg {missing_pkgs}}.",
      "i" = "Install with: {.run install.packages(c({paste0('\"', missing_pkgs, '\"', collapse = ', ')}))}",
      "i" = "Or use {.arg workers = NULL} or {.arg workers = 1L} for sequential profiling."
    ))
  }
}

# Support functions ----

messageProfileComplete <- function(which, direction, msg) {
  dirStr <- if (direction > 0L) {
    "above the point estimate"
  } else {
    "below the point estimate"
  }
  if (startsWith(msg, "complete")) {
    cli::cli_alert_success("Profiling {which} {dirStr} is {msg}")
  } else {
    cli::cli_inform("Profiling {which} {dirStr} is {msg}")
  }
}

# nlmixr2Profile S3 class ----

# Internal constructor: wraps raw profile data.frame in nlmixr2Profile S3 class.
.nlmixr2ProfileNew <- function(
  raw,
  status,
  fit,
  which,
  control,
  fitName = NA_character_
) {
  fixefNames <- names(nlmixr2est::fixef(fit))
  omegaNames <- llpOmegaNames(fit)
  mle <- setNames(
    vapply(
      which,
      function(w) {
        if (w %in% fixefNames) {
          nlmixr2est::fixef(fit)[[w]]
        } else {
          fit$iniDf[fit$iniDf$name == w, "est"]
        }
      },
      numeric(1L)
    ),
    which
  )

  nearBound <- vapply(
    which,
    function(w) {
      any(status$nearBound[status$param == w])
    },
    logical(1L)
  )
  names(nearBound) <- which

  hitMaxIter <- vapply(
    which,
    function(w) {
      any(status$hitMaxIter[status$param == w])
    },
    logical(1L)
  )
  names(hitMaxIter) <- which

  intervalRatio <- vapply(
    which,
    function(w) {
      lb <- status$boundEstimate[status$param == w & status$direction == -1L]
      ub <- status$boundEstimate[status$param == w & status$direction == 1L]
      if (length(lb) != 1L || length(ub) != 1L || is.na(lb) || is.na(ub)) {
        return(NA_real_)
      }
      lDist <- abs(mle[[w]] - lb)
      uDist <- abs(ub - mle[[w]])
      if (lDist == 0 || uDist == 0) {
        return(NA_real_)
      }
      uDist / lDist
    },
    numeric(1L)
  )
  names(intervalRatio) <- which

  attr(raw, "fitName") <- fitName
  attr(raw, "method") <- "llp"
  attr(raw, "ofvIncrease") <- control$ofvIncrease
  attr(raw, "mle") <- mle
  attr(raw, "status") <- status
  attr(raw, "nearBound") <- nearBound
  attr(raw, "hitMaxIter") <- hitMaxIter
  attr(raw, "intervalRatio") <- intervalRatio

  class(raw) <- c("nlmixr2Profile", "data.frame")
  raw
}

#' Profile confidence intervals using log-likelihood profiling
#'
#' @param fit An nlmixr2 fit object (nlmixr2FitCore).
#' @param which Character vector of parameter names to profile.  `NULL`
#'   profiles all v1-profileable fixed-effect parameters.
#' @param control An `llpControl()` object or a named list of arguments to pass
#'   to `llpControl()`.
#' @param ... Ignored.
#' @return An `nlmixr2Profile` object (a `data.frame` subclass).  The
#'   data.frame contents are backward-compatible with the previous plain
#'   data.frame return.  Additional information is available via
#'   `confint()`, `print()`, and `plot()` methods, and via `attr()`:
#'   `status`, `mle`, `nearBound`, `hitMaxIter`, `intervalRatio`,
#'   `ofvIncrease`.
#' @family Profiling
#' @export
runLLP <- function(fit, which = NULL, control = llpControl(), ...) {
  fitName <- tryCatch(deparse(substitute(fit)), error = function(e) {
    NA_character_
  })

  if (!inherits(control, "llpControl")) {
    control <- do.call(llpControl, control)
  }
  checkmate::assert_class(fit, "nlmixr2FitCore")

  fixefNames <- names(nlmixr2est::fixef(fit))
  omegaNames <- llpOmegaNames(fit)

  if (is.null(which)) {
    which <- fixefNames
  } else {
    checkmate::assert_character(which, any.missing = FALSE, min.len = 1L)
    validNames <- c(fixefNames, omegaNames)
    checkmate::assert_subset(x = which, choices = validNames, empty.ok = FALSE)
  }

  if (!is.null(names(control$rseTheta))) {
    checkmate::assert_names(names(control$rseTheta), subset.of = fixefNames)
  }

  whichTheta <- intersect(which, fixefNames)
  whichOmega <- intersect(which, omegaNames)

  workers <- control$workers

  runOne <- function(w) {
    if (w %in% whichOmega) {
      llpRunOneOmega(fitted = fit, which = w, control = control)
    } else {
      llpRunOneParameter(fitted = fit, which = w, control = control)
    }
  }

  if (!is.null(workers) && !identical(workers, 1L)) {
    .llpCheckParallelDeps()
    cli::cli_inform(
      "Profiling {length(which)} parameter{?s} in parallel (workers = {workers})"
    )
    results <- .withWorkerPlan(workers, {
      future.apply::future_lapply(
        which,
        function(w) {
          llpWithIsolatedFitDir(w, function() {
            runOne(w)
          })
        },
        future.seed = TRUE,
        future.packages = "nlmixr2extra"
      )
    })
  } else {
    results <- lapply(seq_along(which), function(i) {
      cli::cli_inform(
        "Profiling parameter {i}/{length(which)}: {.field {which[i]}}"
      )
      runOne(which[i])
    })
  }

  raw <- dplyr::bind_rows(lapply(results, `[[`, "raw"))
  status <- dplyr::bind_rows(lapply(results, `[[`, "status"))
  rownames(raw) <- NULL

  # Keep fixef columns plus any profiled omega columns; drop other omega columns.
  keepCols <- c(
    "Parameter",
    "OFV",
    intersect(c(fixefNames, whichOmega), names(raw)),
    if ("profileBound" %in% names(raw)) "profileBound"
  )
  raw <- raw[, keepCols, drop = FALSE]
  if (!"profileBound" %in% names(raw)) {
    raw$profileBound <- NA_real_
  }

  .nlmixr2ProfileNew(
    raw = raw,
    status = status,
    fit = fit,
    which = which,
    control = control,
    fitName = fitName
  )
}

#' @export
print.nlmixr2Profile <- function(x, ...) {
  mle <- attr(x, "mle")
  if (is.null(mle)) {
    NextMethod()
    return(invisible(x))
  }
  ofvIncrease <- attr(x, "ofvIncrease")
  nearBound <- attr(x, "nearBound")
  hitMaxIter <- attr(x, "hitMaxIter")
  intervalRatio <- attr(x, "intervalRatio")
  params <- names(mle)

  ci <- suppressWarnings(confint(x, parm = params))

  cli::cli_h1("Log-likelihood profile confidence intervals")
  cli::cli_text("{.field dOFV} target: {round(ofvIncrease, 4)}")
  cli::cli_text("")

  col_p <- max(nchar(params), 5L)
  fmt <- paste0("%-", col_p, "s  %10s  %10s  %10s  %6s  %s")
  hdr <- sprintf(fmt, "Param", "Lower CI", "MLE", "Upper CI", "Ratio", "Flags")
  cli::cli_text(hdr)
  cli::cli_text(strrep("-", nchar(hdr)))

  for (p in params) {
    lo <- if (!is.na(ci[p, "lower"])) {
      formatC(ci[p, "lower"], digits = 4L, format = "f")
    } else {
      "      NA"
    }
    hi <- if (!is.na(ci[p, "upper"])) {
      formatC(ci[p, "upper"], digits = 4L, format = "f")
    } else {
      "      NA"
    }
    rat <- if (!is.na(intervalRatio[p])) {
      formatC(intervalRatio[p], digits = 2L, format = "f")
    } else {
      "    NA"
    }
    flags <- character(0L)
    if (isTRUE(nearBound[p])) {
      flags <- c(flags, "near bound")
    }
    if (isTRUE(hitMaxIter[p])) {
      flags <- c(flags, "max iter")
    }
    cli::cli_text(sprintf(
      fmt,
      p,
      lo,
      formatC(mle[[p]], digits = 4L, format = "f"),
      hi,
      rat,
      paste(flags, collapse = ", ")
    ))
  }

  invisible(x)
}

#' @export
confint.nlmixr2Profile <- function(object, parm = NULL, level = 0.95, ...) {
  ofvIncrease <- attr(object, "ofvIncrease")
  expected <- qchisq(level, df = 1L)
  if (abs(ofvIncrease - expected) > 1e-4) {
    cli::cli_warn(c(
      "Requested {level * 100}% CI implies dOFV = {round(expected, 4)},",
      "but this profile used dOFV = {round(ofvIncrease, 4)}."
    ))
  }

  mle <- attr(object, "mle")
  params <- if (is.null(parm)) names(mle) else parm

  bounds <- object[!is.na(object$profileBound), , drop = FALSE]

  mat <- matrix(
    NA_real_,
    nrow = length(params),
    ncol = 2L,
    dimnames = list(params, c("lower", "upper"))
  )

  for (p in params) {
    pb <- bounds[bounds$Parameter == p, , drop = FALSE]
    lo <- pb[pb$profileBound < 0, p, drop = TRUE]
    hi <- pb[pb$profileBound > 0, p, drop = TRUE]
    if (length(lo) > 0L && !is.na(lo[[1L]])) {
      mat[p, "lower"] <- lo[[1L]]
    }
    if (length(hi) > 0L && !is.na(hi[[1L]])) mat[p, "upper"] <- hi[[1L]]
  }

  mat
}

#' @export
plot.nlmixr2Profile <- function(x, parm = NULL, ...) {
  mle <- attr(x, "mle")
  ofvIncrease <- attr(x, "ofvIncrease")
  params <- if (is.null(parm)) names(mle) else parm

  ci <- suppressWarnings(confint(x, parm = params))

  plotList <- lapply(params, function(p) {
    rows <- x[x$Parameter == p & is.na(x$profileBound), , drop = FALSE]
    if (nrow(rows) == 0L) {
      return(NULL)
    }
    minOFV <- min(rows$OFV, na.rm = TRUE)
    data.frame(
      param = p,
      paramVal = rows[[p]],
      dOFV = rows$OFV - minOFV,
      stringsAsFactors = FALSE
    )
  })
  plotData <- dplyr::bind_rows(plotList)

  vlineList <- lapply(params, function(p) {
    xvals <- c(
      MLE = mle[[p]],
      `CI lower` = ci[p, "lower"],
      `CI upper` = ci[p, "upper"]
    )
    xvals <- xvals[!is.na(xvals)]
    data.frame(
      param = p,
      xval = unname(xvals),
      linetype = names(xvals),
      stringsAsFactors = FALSE
    )
  })
  vlineData <- dplyr::bind_rows(vlineList)

  ggplot2::ggplot(
    plotData,
    ggplot2::aes(x = .data[["paramVal"]], y = .data[["dOFV"]])
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_hline(
      yintercept = ofvIncrease,
      linetype = "dashed",
      colour = "grey50"
    ) +
    ggplot2::geom_vline(
      data = vlineData,
      ggplot2::aes(
        xintercept = .data[["xval"]],
        linetype = .data[["linetype"]]
      ),
      colour = "grey30"
    ) +
    ggplot2::scale_linetype_manual(
      values = c("MLE" = "solid", "CI lower" = "dashed", "CI upper" = "dashed")
    ) +
    ggplot2::facet_wrap(~param, scales = "free_x") +
    ggplot2::labs(
      x = "Parameter value",
      y = expression(Delta * "OFV"),
      linetype = NULL
    ) +
    ggplot2::theme_bw()
}
