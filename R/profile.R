#' Perform likelihood profiling on nlmixr2 focei fits
#'
#' @param fitted The fit model
#' @param ... ignored
#' @param which The parameter names to perform likelihood profiling on
#'   (\code{NULL} indicates all parameters)
#' @param maxpts The maximum number of fits to perform in each direction (higher
#'   or lower values) for each parameter.  The total maximum number of estimates
#'   that may be performed is \code{maxpts*2*length(which)}.
#' @param ofvIncrease The targetted change in objective function value (3.84
#'   corresponds to a Chi-squared test with a 95% confidence interval)
#' @param normQuantile The quantile of a normal distribution to use for the
#'   initial estimates.
#' @param rse_theta The relative standard error (percent) for the model
#'   parameters.  It can be missing (the default) in which case a default value
#'   of 30% will be applied.  If given as a single number, it will be applied to
#'   all parameters.  If given as a named vector of numbers, it will be applied
#'   to each named parameter.
#' @param tol The relative tolerance for the parameter estimate being close
#'   enough to the \code{ofvIncrease}.  The tolerance is applied to the
#'   parameter estimate and not the the \code{ofvIncrease} to allow for ranges
#'   with a high derivative to complete in a reasonable number of iterations.
#' @param ignoreBounds Should the boundary conditions in the model be ignored?
#' @param quiet Print less information from model estimation steps
#' @return A data.frame with one column for the \code{Parameter} name being
#'   fixed on that row, one column for the \code{OFV}, and columns for each
#'   parameter estimate (or fixed value) in the model.
#' @export
profile.nlmixr2FitCore <- function(fitted, ..., which = NULL, maxpts = 10, ofvIncrease = qchisq(0.95, df = 1), normQuantile = qnorm(p = 0.975), rse_theta, tol = 0.001, ignoreBounds = FALSE, quiet = TRUE) {
  # Input checking
  if (is.null(which)) {
    which <- names(fixef(fitted))
  } else {
    checkmate::assert_subset(x = which, choices = names(fixef(fitted)), empty.ok = FALSE)
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
  if (missing(rse_theta)) {
    fittedVcov <- vcov(fitted)
    if (is.null(fittedVcov)) {
      # Handle when vcov() is an error
      rse_theta <- 30
      cli::cli_alert_info(paste("covariance is unavailable, using default rse_theta of", rse_theta))
    } else {
      sd_theta <- sqrt(diag(fittedVcov))
      rse_theta <- 100*abs(sd_theta/fixef(fitted)[names(sd_theta)])
    }
  }
  if (length(rse_theta) == 1 & is.null(names(rse_theta))) {
    rse_theta <- setNames(rep(rse_theta, length(which)), which)
  } else {
    checkmate::assert_names(names(rse_theta), subset.of = names(fixef(fitted)))
    checkmate::assert_numeric(rse_theta, lower = 0, any.missing = FALSE, min.len = 1)
  }
  checkmate::assert_number(tol, na.ok = FALSE, lower = 1e-10, upper = 1, finite = TRUE, null.ok = FALSE)
  checkmate::assert_logical(ignoreBounds, any.missing = FALSE, len = 1)
  checkmate::assert_logical(quiet, any.missing = FALSE, len = 1)

  if (quiet) {
    # TODO: set print = 0 for all the control arguments
  }
  if (length(which) > 1) {
    ret <- data.frame()
    # Make the general case be a concatenation of single cases
    for (currentWhich in which) {
      ret <-
        rbind(
          ret,
          profile.nlmixr2FitData(fitted = fitted, ..., which = currentWhich, maxpts = maxpts, ofvIncrease = ofvIncrease, normQuantile = normQuantile, rse_theta = rse_theta, tol = tol, quiet = quiet)
        )
    }
  } else {
    ret <- profileNlmixr2FitDataRet(fitted = fitted, which = which)
    estInitial <-
      profileNlmixr2FitDataEstInitial(
        estimates = ret,
        which = which,
        normQuantile = normQuantile,
        rse_theta = rse_theta
      )
    for (direction in c(-1, 1)) {
      found <- FALSE
      boundCol <- ifelse(direction < 0, "lower", "upper")
      if (ignoreBounds) {
        bound <- Inf*direction
      } else {
        bound <- setNames(fitted$ui$iniDf[[boundCol]], fitted$ui$iniDf$name)[which]
      }
      if (abs(bound - ret[[which]][1]/ret[[which]][1]) < tol) {
        found <- TRUE
        cli::cli_warn(
          "parameter %s is close to the %s boundary and likelihood will not be profiled in that direction",
          which, boundCol
        )
      }
      # Confirm that the initial estimate is within the bounds
      if (direction < 0) {
        currentEst <- max(bound, min(estInitial))
      } else {
        currentEst <- min(bound, max(estInitial))
      }
      currentMaxpts <- maxpts
      # Prepare the current parameter to be fixed for each estimation
      fittedFix <-
        suppressMessages(do.call(
          ini, append(list(x=fitted), setNames(list(as.name("fix")), which))
        ))
      while (currentMaxpts > 0 & !found) {
        currentMaxpts <- currentMaxpts - 1
        newEst <- setNames(list(currentEst), which)
        # TODO: How can I detect the current estimation method from the fitted
        # object?  Can you even do likelihood profiling with something other
        # than focei?
        newFit <- try(nlmixr(ini(fittedFix, newEst), est = fitted$est, control = list(print = 0)))
        retRow <- profileNlmixr2FitDataRet(fitted = newFit, which = which, fixedVal = currentEst, rowname = nrow(ret))
        if (inherits(newFit, "try-error")) {
          # Stop for model convergence errors
          extraCols <- setdiff(names(ret), names(retRow))
          retRow <- cbind(retRow, data.frame(t(setNames(rep(NA_character_, length(extraCols)), extraCols))))
          cli::cli_alert_danger(sprintf("Stopping search for %s in %s direction because the fit failed at an estimate of %g", which, boundCol, currentEst))
          found <- TRUE
        } else {
          # Get a new estimate
          currentEst <- profileNlmixr2FitDataNewEst(estimates = ret, which = which, direction = direction, bound = bound, ofvIncrease = ofvIncrease, tol = tol)
          found <- profileNlmixr2FitDatFound(estimates = estimates, newEst = currentEst, bound = bound, direction = direction)
          if (is.character(found)) {
            found <- TRUE
            cli::cli_alert_danger(sprintf("Stopping search for %s in %s direction because %s", which, boundCol, found))
          }
        }
        ret <- rbind(ret, retRow)
        ret <- ret[order(ret[[which]]), ]
        browser()
      }
    }
  }
  ret
}

#' Give the output data.frame for a single model for profile.nlmixr2FitData
#' @inheritParams profile.nlmixr2FitData
#' @return A data.frame with columns for Parameter (the parameter name), OFV
#'   (the objective function value), and the current estimate for each of the
#'   parameters
#' @noRd
profileNlmixr2FitDataRet <- function(fitted, which, fixedVal, rowname = 0) {
  if (inherits(fitted, "try-error")) {
    ret <- data.frame(Parameter = which, OFV = NA_real_, X = fixedVal)
    names(ret)[3] <- which
  } else {
    if (any(names(fixef(fitted)) %in% c("Parameter", "OFV"))) {
      cli::cli_abort("Cannot profile a model with a parameter named either 'Parameter' or 'OFV'")
    }
    ret <-
      cbind(
        data.frame(Parameter = which, OFV = fitted$objDf$OBJF),
        data.frame(t(fixef(fitted)))
      )
  }
  rownames(ret) <- paste(which, rowname)
  ret
}

# Constrain the estimate to be within the bounds, shift away from the bound if
# the boundary is already estimated
profileNlmixr2FitDatEstBound <- function(estimates, which, newEst, bound, direction, tol, atol) {
  # Ensure that the bounds are not exceeded
  if (bound %in% estimates[[which]]) {
    # negative direction ensures that the value is within the bound (rather than
    # outside)
    if (bound == 0) {
      boundAdj <- -direction*atol/2
    } else {
      boundAdj <- bound*(1 - direction*tol/2)
    }
  } else {
    boundAdj <- bound
  }
  if (direction < 0) {
    newEst <- max(newEst, boundAdj)
  } else if (direction > 0) {
    newEst <- min(newEst, boundAdj)
  }
  newEst
}

# Provide the initial estimates going up and down from the initial value
profileNlmixr2FitDataEstInitial <- function(estimates, which, normQuantile, rse_theta) {
  checkmate::assert_data_frame(estimates, nrows = 1)
  checkmate::assert_character(which, any.missing = FALSE, len = 1)
  checkmate::assert_subset(which, choices = names(estimates))
  checkmate::assert_subset(which, choices = names(rse_theta))
  # The -1,1 makes the estimate go up and down
  estimates[[which]] + c(-1, 1)*normQuantile*rse_theta[which]/100*abs(estimates[[which]])
}

#' Generate a new estimate for likelihood profiling, and return a list
#' indicating if the new estimate should not be run (because the current data
#' find the value within the tolerance)
#'
#' @param estimates A data.frame with the parameter estimates
profileNlmixr2FitDataNewEst <- function(estimates, which, direction, bound, ofvIncrease, method = "linapprox") {
  checkmate::assert_data_frame(estimates, min.rows = 2)
  checkmate::assert_subset(which, choices = names(estimates))
  stopifnot(direction != 0)
  checkmate::assert_number(bound, na.ok = FALSE, finite = FALSE, null.ok = FALSE)
  checkmate::assert_number(ofvIncrease, lower = 0, na.ok = FALSE, finite = TRUE, null.ok = FALSE)

  cli::cli_abort("More than one estimation row is required")
  minOFV <- min(estimates$OFV, na.rm = TRUE)
  minRow <- which(estimates$OFV %in% minOFV)
  if (length(minRow) > 1) {
    if (direction < 0) {
      # Choose the lowest minimum for decreasing direction
      minRow <- minRow[1]
    } else {
      # Choose the highest minimum for increasing direction
      minRow <- minRow[length(minRow)]
    }
  }
  maskDirection <- !is.na(estimates$OFV) & estimates[[which]] <= estimates[[which]][minRow]
  dOFV <- estimates$OFV[maskDirection] - minOFV
  # Estimates in the correct direction
  estDir <- estimates[[which]][maskDirection]
  if (length(estDir) == 0) {
    cli::cli_abort("No values detected in the estimation direction, please report a bug") # nocov
  } else if (length(estDir) == 1) {
    browser()
    # Move in the opposite direction by the range
    newEst <- estDir - direction*diff(range(estimates[[which]]))
  } else if (all(dOFV < ofvIncrease)) {
    # Expand the search
    if (length(estDir) == 1) {
      distance <- abs(estimates[[which]][which(rownames(estimates) == "0")] - estimates[[which]][minRow])
    } else {
      distance <- abs(diff(range(estDir)))
    }
    newEst <- estimates[[which]][minRow] + direction*2*distance
  } else {
    # Check that the OFV is monotonic
    monotonic <- length(unique(sign(diff(dOFV)/diff(estDir)))) == 1
    if (!monotonic) {
      newEst <- "OFV is not monotonic"
    } else if (method == "linapprox") {
      newEst <- approx(x = dOFV, y = estDir, xout = ofvIncrease)$y
    } else {
      cli::cli_abort(sprintf("method %s is not implemented", method))
    }
  }
  newEst
}

# Determine if likelihood profiling is complete
profileNlmixr2FitDatFound <- function(estimates, newEst, bound, direction) {
  browser()
  if (direction < 0) {
    maskEstAbove <- estimates[[which]] > newEst
    # Equality allow capturing the bound, if applicable
    maskEstBelow <- estimates[[which]] <= newEst
  } else if (direction > 0) {
    # Equality allow capturing the bound, if applicable
    maskEstAbove <- estimates[[which]] >= newEst
    maskEstBelow <- estimates[[which]] < newEst
  }

  if (any(estimates[[which]] %in% newEst)) {
    found <- TRUE
  } else if (any(maskEstAbove) & any(maskEstBelow)) {
    # If we're interpolating, check the tolerance
    estAbove <- min(estimates[[which]][maskEstAbove])
    estBelow <- min(estimates[[which]][maskEstBelow])

    estCloserToZero <- ifelse(abs(estAbove) < abs(estBelow), estAbove, estBelow)
    if (estAbove == 0 | estBelow == 0) {
      stop("TODO: Handle absolute tolerance when values are at or near zero")
    } else {
      found <- abs((estAbove - estBelow)/estCloserToZero) < tol
    }
  } else {
    # If we're extrapolating, no need to check the tolerance
    found <- FALSE
  }
  if (!found) {
    # TODO: Make the minimum step size half the tolerance
    browser()
  }

}
