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
profile.nlmixr2FitCore <- function(fitted, ...,
                                   which = NULL,
                                   method = c("llp", "fixed"),
                                   control = list()) {
  method <- match.arg(method)

  if (method == "llp") {
    profileLlp(fitted = fitted, which = which, control = control)
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
      cli::cli_abort("Cannot profile a model with a parameter named either 'Parameter' or 'OFV'")
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
#' @family Profiling
#' @export
profileFixed <- function(fitted, which, control = list()) {
  control <- do.call(fixedControl, control)
  checkmate::assert_data_frame(which, types = "numeric", any.missing = FALSE, min.rows = 1)
  dplyr::bind_rows(lapply(
    X = seq_len(nrow(which)),
    FUN = \(idx, fitted) profileFixedSingle(fitted = fitted, which = which[idx, , drop = FALSE]),
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
#' @export
profileFixedSingle <- function(fitted, which) {
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
  # Validate inputs
  control <- do.call(llpControl, control)

  if (is.null(which)) {
    which <- names(nlmixr2est::fixef(fitted))
  } else {
    checkmate::assert_subset(x = which, choices = names(nlmixr2est::fixef(fitted)), empty.ok = FALSE)
  }

  if (length(which) > 1) {
    ret <- data.frame()
    # Make the general case of profiling many parameters to be a concatenation
    # of single cases
    for (currentWhich in which) {
      ret <-
        dplyr::bind_rows(
          ret,
          profileLlp(
            fitted = fitted,
            which = currentWhich,
            control = control
          )
        )
    }
  } else {
    if (length(control$rseTheta) == 1 & is.null(names(control$rseTheta))) {
      control$rseTheta <- setNames(rep(control$rseTheta, length(which)), which)
    } else {
      checkmate::assert_names(names(control$rseTheta), subset.of = names(nlmixr2est::fixef(fitted)))
      checkmate::assert_numeric(control$rseTheta, lower = 0, any.missing = FALSE, min.len = 1)
    }

    # The original fitted model should be in the output
    ret <- profileNlmixr2FitCoreRet(fitted = fitted, which = which)

    effectRange <- fitted$iniDf[fitted$iniDf$name == which, c("lower", "upper")]
    estInitial <-
      profileNlmixr2FitDataEstInitial(
        estimates = ret,
        which = which,
        ofvIncrease = control$ofvIncrease,
        rseTheta = control$rseTheta,
        lower = effectRange$lower,
        upper = effectRange$upper
      )
    # Find the lower and upper limits
    for (direction in c(-1, 1)) {
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
          ofvIncrease = control$ofvIncrease,
          lower = currentLower, upper = currentUpper,
          itermax = control$itermax,
          ofvtol = control$ofvtol,
          paramDigits = control$paramDigits
        )
      ret <- unique(dplyr::bind_rows(ret, retCurrentDirection))
    }
  }
  ret
}

# Single-parameter, single-direction optimization of likelihood profiling
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
    fitResult <- profileFixedSingle(fitted = fitted, which = dfWhich)

    ret <- dplyr::bind_rows(ret, fitResult)
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
      # When generating the extrapolation rows, always have the first row as the
      # most extreme (farthest out)
      if (direction > 0) {
        extrapRows <- ret[nrow(ret) + c(-1, 0), ]
      } else {
        extrapRows <- ret[1:2, ]
      }
      currentPar <-
        extrapolateExpand*diff(extrapRows[[which]])/diff(extrapRows$OFV)*(targetOfv - extrapRows$OFV[1]) + extrapRows[[which]][1]
      # ensure that we are within the boundaries with a slight margin
      margin <- sqrt(.Machine$double.eps)
      currentPar <- min(max(currentPar, lower + margin), upper - margin)
    } else { # interpolate
      currentPar <- stats::approx(x = ret$OFV, y = ret[[which]], xout = targetOfv)$y
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
    convergedSignif <- any(signifSame & signifAbove) && any(signifSame & signifBelow)

    currentOFVDiff <- min(abs(targetOfv - ret$OFV), na.rm = TRUE)
    convergedOfv <- currentOFVDiff <= ofvtol

    converged <- convergedSignif | convergedOfv

    if (converged) {
      if (convergedSignif & convergedOfv) {
        messageProfileComplete(which, direction = direction, "complete due to significant digits precision and achieving OFV within specified tolerance")
      } else if (convergedSignif) {
        messageProfileComplete(which, direction = direction, "complete due to significant digits precision")
      } else if (convergedOfv) {
        messageProfileComplete(which, direction = direction, "complete due to achieving OFV within specified tolerance")
      } else {
        stop("Unclear convergece criteria, please report a bug") # nocov
      }
      storeBound <- data.frame(Parameter = which, X = currentPar, profileBound = direction*ofvIncrease)
      names(storeBound)[2] <- which
      ret <- dplyr::bind_rows(storeBound, ret)
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

#' Control options for log-likelihood profiling
#'
#' @param ofvIncrease The targetted change in objective function value (3.84
#'   corresponds to a Chi-squared test with a 95% confidence interval)
#' @param rseTheta The relative standard error (percent) for the model
#'   parameters.  It can be missing (the default) in which case a default value
#'   of 30% will be applied.  If given as a single number, it will be applied to
#'   all parameters.  If given as a named vector of numbers, it will be applied
#'   to each named parameter.
#' @param itermax Maximum number of likelihood profiling iterations for each
#'   bound estimated
#' @param ofvtol The relative tolerance for the objective function being close
#'   enough to the `ofvIncrease`.
#' @param paramDigits The number of significant digits required for the
#'   parameter.  When interpolation attempts to get smaller than that number of
#'   significant digits, it will stop.
#' @param extrapolateExpand When extrapolating outside the range previously
#'   tested, how far should the step occur as a ratio
#' @returns A validated list of control options for log-likelihood profiling
#' @family Profiling
#' @seealso [profileLlp()]
#' @export
llpControl <- function(ofvIncrease = qchisq(0.95, df = 1),
                       rseTheta = 30,
                       itermax = 10,
                       ofvtol = 0.005,
                       paramDigits = 3,
                       extrapolateExpand = 1.5) {
  ret <-
    list(
      ofvIncrease = checkmate::assert_number(ofvIncrease, lower = 0, finite = TRUE, null.ok = FALSE, na.ok = FALSE),
      rseTheta = checkmate::assert_number(rseTheta, lower = 0, finite = TRUE, null.ok = FALSE, na.ok = FALSE),
      ofvtol = checkmate::assert_number(ofvtol, na.ok = FALSE, lower = 1e-10, upper = 1, finite = TRUE, null.ok = FALSE),
      itermax = checkmate::assert_integerish(itermax, lower = 1, any.missing = FALSE, len = 1, null.ok = FALSE, coerce = TRUE),
      paramDigits = checkmate::assert_integerish(paramDigits, lower = 1, upper = 10, any.missing = FALSE, len = 1, null.ok = FALSE, coerce = TRUE),
      extrapolateExpand = checkmate::assert_number(extrapolateExpand, lower = 1.001, na.ok = FALSE, finite = TRUE, null.ok = FALSE)
    )
  class(ret) <- "llpControl"
  ret
}

# Provide the initial estimates going up and down from the initial value
profileNlmixr2FitDataEstInitial <- function(estimates, which, ofvIncrease, rseTheta, lower, upper) {
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
  ret <- estimates[[which]] + c(-1, 1)*ofvIncrease*currentRseTheta/100*abs(estimates[[which]])
  # Ensure that the estimate is within the bounds with a slight margin
  margin <- sqrt(.Machine$double.eps)
  pmax(pmin(ret, upper - margin), lower + margin)
}

# Support functions ----

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
