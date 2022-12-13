profile.nlmixr2FitData <- function(fitted, ..., which = NULL, maxpts = 10, ofvIncrease = 3.84, normQuantile = 1.96, rse_theta, tol = 0.001, ignoreBounds = FALSE, quiet = TRUE) {
  if (is.null(which)) {
    which <- names(fixef(fitted))
  } else {
    checkmate::assert_subset(x = which, choices = names(fixef(fitted)), empty.ok = FALSE)
  }
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
      currentEst <- ret[[which]][1] + direction*normQuantile*rse_theta[which]/100*abs(ret[[which]][1])
      # Confirm that the initial estimate is within the bounds
      if (direction < 0) {
        currentEst <- max(bound, currentEst)
      } else {
        currentEst <- min(bound, currentEst)
      }
      currentMaxpts <- maxpts
      # Prepare the current parameter to be fixed for each estimation
      fittedFix <- do.call(ini, append(list(x=fitted), setNames(list(as.name("fix")), which)))
      while (currentMaxpts > 0 & !found) {
        currentMaxpts <- currentMaxpts - 1
        newEst <- setNames(list(currentEst), which)
        # TODO: How can I detect the current estimation method from the fitted
        # object?  Can you even do likelihood profiling with something other
        # than focei?
        newFit <- try(nlmixr(ini(fittedFix, newEst), est = "focei", control = list(print = 0)))
        if (inherits(newFit, "try-error")) {
          # Stop trying with failed fits
          cli::cli_alert_danger(sprintf("Stopping search for %s in %s direction because the fit failed at an estimate of %g", which, boundCol, currentEst))
          found <- TRUE
        } else {
          ret <- rbind(ret, profileNlmixr2FitDataRet(fitted = newFit, which = which, rowname = nrow(ret)))
          ret <- ret[order(ret[[which]]), ]
          currentEstList <- profileNlmixr2FitDataNewEst(estimates = ret, which = which, direction = direction, bound = bound, ofvIncrease = ofvIncrease, tol = tol)
          if (is.character(currentEstList$found)) {
            found <- TRUE
            cli::cli_alert_danger(sprintf("Stopping search for %s in %s direction because %s", which, boundCol, currentEstList$found))
          } else {
            currentEst <- currentEstList$newEst
            found <- currentEstList$found
          }
        }
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
profileNlmixr2FitDataRet <- function(fitted, which, rowname = 0) {
  if (any(names(fixef(fitted)) %in% c("Parameter", "OFV"))) {
    cli::cli_abort("Cannot profile a model with a parameter named either 'Parameter' or 'OFV'")
  }
  ret <-
    cbind(
      data.frame(Parameter = which, OFV = fitted$objDf$OBJF),
      data.frame(t(fixef(fitted)))
    )
  rownames(ret) <- as.character(rowname)
  ret
}

#' Generate a new estimate, and return a list indicating if the new estimate
#' should not be run (because the current data find the value within the
#' tolerance)
#'
#' @param estimates A data.frame with the parameter estimates
profileNlmixr2FitDataNewEst <- function(estimates, which, direction, bound, ofvIncrease, tol, method = "linapprox") {
  if (nrow(estimates) == 1) {
    newEst <- estimates[[which]][1] + direction*normQuantile*rse_theta[which]/100*abs(estimates[[which]][1])
  } else {
    minOFV <- min(estimates$OFV, na.rm = TRUE)
    # [1] is required on which in case there are equal OFVs in some scenarios
    minRow <- which(estimates$OFV %in% minOFV)[1]
    maskDirection <- !is.na(estimates$OFV) & estimates[[which]] <= estimates[[which]][minRow]
    dOFV <- estimates$OFV[maskDirection] - minOFV
    # Estimates in the correct direction
    estDir <- estimates[[which]][maskDirection]
    if (all(dOFV < ofvIncrease)) {
      # Expand the search
      if (length(estDir) == 1) {
        distance <- abs(estimates[[which]][which(rownames(estimates) == "0")] - estimates[[which]][minRow])
      } else {
        distance <- abs(diff(range(estDir)))
      }
      newEst <- estimates[[which]][minRow] + direction*2*distance
      found <- FALSE
    } else {
      # Check that the OFV is monotonic
      monotonic <- length(unique(sign(diff(dOFV)/diff(estDir)))) == 1
      if (!monotonic) {
        newEst <- NA
        found <- "OFV is not monotonic"
      } else if (method == "linapprox") {
        newEst <- approx(x = dOFV, y = estDir, xout = ofvIncrease)$y
      } else {
        cli::cli_abort(sprintf("method %s is not implemented", method))
      }
    }
  }

  # Check if the estimation is complete
  if (!is.na(newEst)) {
    maskEstAbove <- estimates[[which]] > newEst
    maskEstBelow <- estimates[[which]] < newEst
    # Ensure that the bounds are not exceeded
    if (direction < 0 & newEst < bound) {
      newEst <- bound
      maskEstBelow <- estimates[[which]] <= newEst
    } else if (direction > 0 & newEst > bound) {
      newEst <- bound
      maskEstAbove <- estimates[[which]] >= newEst
    }

    if (bound %in% estimates[[which]] & newEst == bound) {
      # Avoid the bound, if we've already tried it
      newEst <- bound - direction*(bound*tol)/2
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
  list(
    newEst = newEst,
    found = found
  )
}
