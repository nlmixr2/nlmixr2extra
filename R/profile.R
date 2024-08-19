#' Perform likelihood profiling on nlmixr2 focei fits
#'
#' The search will stop when either the OFV is within `ofvtol` of the desired
#' OFV change or when the parameter is interpolating to more significant digits
#' than specified in `paramDigits`.
#'
#' @param fitted The fit model
#' @param ... ignored
#' @param which The parameter names to perform likelihood profiling on
#'   (`NULL` indicates all parameters)
#' @param maxpts The maximum number of fits to perform in each direction (higher
#'   or lower values) for each parameter.  The total maximum number of estimates
#'   that may be performed is `maxpts*2*length(which)`.
#' @param ofvIncrease The targetted change in objective function value (3.84
#'   corresponds to a Chi-squared test with a 95% confidence interval)
#' @param normQuantile The quantile of a normal distribution to use for the
#'   initial estimates.
#' @param rseTheta The relative standard error (percent) for the model
#'   parameters.  It can be missing (the default) in which case a default value
#'   of 30% will be applied.  If given as a single number, it will be applied to
#'   all parameters.  If given as a named vector of numbers, it will be applied
#'   to each named parameter.
#' @param ofvtol The relative tolerance for the objective function being close
#'   enough to the `ofvIncrease`.
#' @param paramDigits The number of significant digits required for the
#'   parameter.  When interpolation attempts to get smaller than that number of
#'   significant digits, it will stop.
#' @param ignoreBounds Should the boundary conditions in the model be ignored?
#' @param quiet Print less information from model estimation steps
#' @return A data.frame with one column named `Parameter` indicating the
#'   parameter being fixed on that row, one column for the `OFV` indicating the
#'   OFV estimated for the model at that step, one column named `profileBound`
#'   indicating the estimated value for the profile likelihood and its step
#'   above the minimum profile likelihood value, and columns for each parameter
#'   estimate (or fixed) in the model.
#' @family Profiling
#' @export
profile.nlmixr2FitCore <- function(fitted, ...,
                                   which = NULL, maxpts = 10,
                                   ofvIncrease = qchisq(0.95, df = 1),
                                   normQuantile = qnorm(p = 0.975),
                                   rseTheta,
                                   ignoreBounds = FALSE,
                                   quiet = TRUE,
                                   itermax = 10,
                                   ofvtol = 0.005,
                                   paramDigits = 3) {
  # Input checking
  if (is.null(which)) {
    which <- names(nlmixr2est::fixef(fitted))
  } else {
    checkmate::assert_subset(x = which, choices = names(nlmixr2est::fixef(fitted)), empty.ok = FALSE)
  }

  fittedControl <- setQuietFastControl(fitted$control)
  # fitted$objective sets the ofvType to ensure that we have the same type of
  # OFV (and logLik) later
  fitted$objective
  fittedOfvType <- fitted$ofvType
  # setOfv(currentFit, fittedOfvType)
  # logLik(currentFit)

  checkmate::assert_integerish(maxpts, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assert_number(ofvIncrease, lower = 0, finite = TRUE, null.ok = FALSE, na.ok = FALSE)
  checkmate::assert_integerish(paramDigits, lower = 1, upper = 10, any.missing = FALSE, len = 1, null.ok = FALSE, coerce = TRUE)
  if (missing(rseTheta)) {
    fittedVcov <- stats::vcov(fitted)
    if (is.null(fittedVcov)) {
      # Handle when vcov() is an error
      rseTheta <- 30
      cli::cli_alert_info(paste("covariance is unavailable, using default rseTheta of", rseTheta))
    } else {
      sd_theta <- sqrt(diag(fittedVcov))
      rseTheta <- 100*abs(sd_theta/nlmixr2est::fixef(fitted)[names(sd_theta)])
    }
  }
  if (length(rseTheta) == 1 & is.null(names(rseTheta))) {
    rseTheta <- setNames(rep(rseTheta, length(which)), which)
  } else {
    checkmate::assert_names(names(rseTheta), subset.of = names(nlmixr2est::fixef(fitted)))
    checkmate::assert_numeric(rseTheta, lower = 0, any.missing = FALSE, min.len = 1)
  }
  checkmate::assert_number(ofvtol, na.ok = FALSE, lower = 1e-10, upper = 1, finite = TRUE, null.ok = FALSE)
  checkmate::assert_logical(ignoreBounds, any.missing = FALSE, len = 1)
  checkmate::assert_logical(quiet, any.missing = FALSE, len = 1)

  if (quiet) {
    # TODO: set print = 0 for all the control arguments
  }
  if (length(which) > 1) {
    ret <- data.frame()
    # Make the general case of profiling many parameters to be a concatenation
    # of single cases
    for (currentWhich in which) {
      ret <-
        rbind(
          ret,
          profile(
            fitted = fitted,
            ...,
            which = currentWhich,
            maxpts = maxpts,
            ofvIncrease = ofvIncrease,
            normQuantile = normQuantile,
            rseTheta = rseTheta,
            tol = tol,
            quiet = quiet,
            optimControl = optimControl
          )
        )
    }
  } else {
    ret <- profileNlmixr2FitCoreRet(fitted = fitted, which = which)
    estInitial <-
      profileNlmixr2FitDataEstInitial(
        estimates = ret,
        which = which,
        normQuantile = normQuantile,
        rseTheta = rseTheta
      )
    for (direction in c(-1, 1)) {
      effectRange <- fitted$iniDf[fitted$iniDf$name == which, c("lower", "upper")]
      if (direction == -1) {
        currentInitialEst = estInitial[1]
        currentUpper <- nlmixr2est::fixef(fitted)[[which]]
        currentLower <- effectRange$lower
      } else {
        currentInitialEst = estInitial[2]
        currentLower <- nlmixr2est::fixef(fitted)[[which]]
        currentUpper <- effectRange$upper
      }

      retCurrentDirection <-
        optimProfile(
          par = currentInitialEst,
          fitted = fitted,
          which = which,
          optimDf = ret,
          direction = direction,
          ofvIncrease = ofvIncrease,
          lower = currentLower, upper = currentUpper,
          itermax = itermax,
          ofvtol = ofvtol,
          paramDigits = paramDigits
        )
      ret <- unique(rbindAddCols(ret, retCurrentDirection))
    }
  }
  ret
}

#' Give the output data.frame for a single model for profile.nlmixr2FitCore
#' @inheritParams profile.nlmixr2FitCore
#' @return A data.frame with columns for Parameter (the parameter name), OFV
#'   (the objective function value), and the current estimate for each of the
#'   parameters
#' @noRd
profileNlmixr2FitCoreRet <- function(fitted, which, fixedVal, rowname = 0) {
  if (inherits(fitted, "try-error")) {
    ret <- data.frame(Parameter = which, OFV = NA_real_, X = fixedVal)
    names(ret)[3] <- which
  } else {
    if (any(names(nlmixr2est::fixef(fitted)) %in% c("Parameter", "OFV"))) {
      cli::cli_abort("Cannot profile a model with a parameter named either 'Parameter' or 'OFV'")
    }
    ret <-
      cbind(
        data.frame(Parameter = which, OFV = fitted$objective),
        data.frame(t(nlmixr2est::fixef(fitted)))
      )
  }
  rownames(ret) <- NULL
  ret
}

# Provide the initial estimates going up and down from the initial value
profileNlmixr2FitDataEstInitial <- function(estimates, which, normQuantile, rseTheta) {
  checkmate::assert_data_frame(estimates, nrows = 1)
  checkmate::assert_character(which, any.missing = FALSE, len = 1)
  checkmate::assert_subset(which, choices = names(estimates))
  # use default rseTheta if not available
  if (which %in% names(rseTheta)) {
    currentRseTheta <- rseTheta[which]
  } else {
    currentRseTheta <- 30
  }

  # The -1,1 makes the estimate go up and down
  estimates[[which]] + c(-1, 1)*normQuantile*currentRseTheta/100*abs(estimates[[which]])
}

#' Estimate the objective function value for a model while fixing a single set
#' of defined parameter values (for use in parameter profiling)
#'
#' @inheritParams profile.nlmixr2FitCore
#' @param which A data.frame with column names of parameters to fix and values
#'   of the fitted value to fix (one row only).
#' @returns `which` with a column named `OFV` added with the objective function
#'   value of the fitted estimate fixing the parameters in the other columns
#' @family Profiling
#' @export
profileNlmixr2SingleParam <- function(fitted, which) {
  checkmate::assert_data_frame(which, types = "numeric", any.missing = FALSE, nrows = 1)
  checkmate::assert_names(names(which), subset.of = names(nlmixr2est::fixef(fitted)))

  message("Profiling ", paste(names(which), "=", unlist(which), collapse = ", "))
  # Update the model by fixing all of the parameters
  paramToFix <- lapply(X = which, FUN = \(x) str2lang(sprintf("fixed(%s)", x)))
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
  ret <- profileNlmixr2FitCoreRet(fitted = newFit, which = paste(names(which), collapse = ","))
  message("OFV: ", ret$OFV)
  ret
}

optimProfile <- function(par, fitted, optimDf, which, ofvIncrease, direction, lower = -Inf, upper = Inf,
                         itermax = 10, ofvtol = 0.005, paramDigits = 3,
                         extrapolateExpand = 1.5) {
  currentOFVDiff <- Inf
  currentIter <- 0
  currentPar <- par
  origMinOFV <- min(optimDf$OFV, na.rm = TRUE)
  ret <- optimDf[optimDf$Parameter %in% which & optimDf[[which]] >= lower & optimDf[[which]] <= upper, ]
  while (itermax > currentIter) {
    currentIter <- currentIter + 1
    dfWhich <- stats::setNames(data.frame(X = currentPar), nm = which)
    fitResult <- profileNlmixr2SingleParam(fitted = fitted, which = dfWhich)

    ret <- rbindAddCols(ret, fitResult)
    ret <- ret[order(ret$OFV), , drop = FALSE]
    currentMinOFV <- min(ret$OFV, na.rm = TRUE)
    if (isTRUE(any(currentMinOFV < origMinOFV))) {
      warning("OFV decreased while profiling ", which, " profile may be inaccurate")
      origMinOFV <- NA # Only give the warning once
    }
    targetOfv <- currentMinOFV + ofvIncrease
    if (all(targetOfv > ret$OFV)) { # extrapolate
      # Find the closest two rows for extrapolation
      extrapRows <- ret[!is.na(ret$OFV), ]
      if (nrow(extrapRows) < 2) {
        messageProfileComplete(which, direction = direction, "aborted due to lack of convergence prior to attempted extrapolation")
        currentIter <- Inf
      }
      if (direction > 0) {
        extrapRows <- ret[nrow(ret) + c(-1, 0), ]
      } else {
        extrapRows <- ret[1:2, ]
      }
      currentPar <-
        extrapolateExpand*diff(extrapRows[[which]])/diff(extrapRows$OFV)*(targetOfv - extrapRows$OFV[1]) + extrapRows[[which]][1]
      # ensure that we are within the boundary
      currentPar <- min(max(currentPar, lower), upper)
    } else { # interpolate
      currentPar <- approx(x = ret$OFV, y = ret[[which]], xout = targetOfv)$y
    }
    if (currentPar %in% ret[[which]]) {
      # Don't try to test the same parameter value multiple times (usually the
      # case for hitting a limit)
      messageProfileComplete(which, direction = direction, "aborted due to attempted repeated parameter value estimation")
      currentIter <- Inf
    }
    # Check if the parameter significant digits are sufficiently precise.  This
    # is done by finding values that have the same significant digit rounding
    # where one estimate is above and one is below the next parameter to be
    # estimated.
    signifSame <- signif(ret[[which]], digits = paramDigits) == signif(currentPar, digits = paramDigits)
    signifAbove <- ret[[which]] > currentPar
    signifBelow <- ret[[which]] < currentPar
    if (any(signifSame & signifAbove) && any(signifSame & signifBelow)) {
      messageProfileComplete(which, direction = direction, "complete due to significant digits precision")
      storeBound <- data.frame(Parameter = which, X = currentPar, profileBound = direction*ofvIncrease)
      names(storeBound)[2] <- which
      ret <- rbindAddCols(storeBound, ret)
      currentIter <- Inf
    }

    currentOFVDiff <- min(abs(targetOfv - ret$OFV), na.rm = TRUE)
    if (currentOFVDiff <= ofvtol) {
      messageProfileComplete(which, direction = direction, "complete due to achieving OFV within specified tolerance")
      storeBound <- data.frame(Parameter = which, X = currentPar, profileBound = direction*ofvIncrease)
      names(storeBound)[2] <- which
      ret <- rbindAddCols(storeBound, ret)
      currentIter <- Inf
    } else {
      message("Profiling ", which, " best distance to target OFV: ", currentOFVDiff)
    }
  }
  if (is.finite(currentIter)) {
    messageProfileComplete(which, direction = direction, "aborted due to too many iterations without convergence")
  }

  # Sort the output by parameter value
  ret <- ret[order(ret[[which]]), ]
  # Add extra information to assist with subsequent
  # The output rownames are uninformative and distracting; remove them.
  rownames(ret) <- NULL
  ret
}

messageProfileComplete <- function(which, direction, msg) {
  message(
    "Profiling ",
    which,
    " ",
    ifelse(direction > 0, "above the point estimate", "below the point estimate"),
    " is ",
    msg
  )
}

#' @describeIn profileNlmixr2SingleParam Estimate the objective function values
#'   for a model while fixing defined parameter values (for use in parameter
#'   profiling)
#' @param which A data.frame with column names of parameters to fix and values
#'   of the fitted value to fix (one row per set of parameters to estimate)
#' @export
profileNlmixr2MultiParam <- function(fitted, which) {
  checkmate::assert_data_frame(which, types = "numeric", any.missing = FALSE, min.rows = 1)
  dplyr::bind_rows(lapply(
    X = seq_len(nrow(which)),
    FUN = \(idx, fitted) profileNlmixr2SingleParam(fitted = fitted, which = which[idx, , drop = FALSE]),
    fitted = fitted
  ))
}

# Add more columns to x so that `rbind()` works
rbindAddCols <- function(x, y, .default = NA) {
  extraColsX <- setdiff(names(y), names(x))
  if (length(extraColsX) > 0) {
    x[extraColsX] <- .default
  }
  extraColsY <- setdiff(names(x), names(y))
  if (length(extraColsY) > 0) {
    y[extraColsY] <- .default
  }
  rbind(x, y)
}
