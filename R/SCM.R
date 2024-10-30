#' Stepwise Covariate Model-selection (SCM) method
#'
#' @param fit an nlmixr2 'fit' object
#' @param varsVec a list of candidate variables to which the
#'   covariates could be added
#' @param covarsVec a list of candidate covariates that need to be
#'   tested
#' @param catvarsVec  character vector of categorical covariates that need to be added
#' @param pVal a named list with names 'fwd' and 'bck' for specifying
#'   the p-values for the forward and backward searches, respectively
#' @param searchType one of 'scm', 'forward' and 'backward' to specify
#'   the covariate search method; default is 'scm'
#' @param restart a boolean that controls if the search should be
#'   restarted; default is FALSE
#'
#' @return A list summarizing the covariate selection steps and
#'   output; This list has the "summaryTable" for the overall summary
#'   of the covariate selection as well as "resFwd" for the forward
#'   selection method and "resBck" for the backward selection method.
#'
#' @noRd
#' @author Vipul Mann, Matthew Fidler, Vishal Sarsani
#'
#' @examples
#' \dontrun{
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45; label("Ka")
#'     tcl <- log(c(0, 2.7, 100)); label("Cl")
#'     tv <- 3.45; label("V")
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' fit <- nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "saem", control = list(print = 0))
#' rxode2::.rxWithWd(tempdir(), {# with temporary directory
#'
#' auto1 <- covarSearchAuto(fit, varsVec = c("ka", "cl"),
#'     covarsVec = c("WT"))
#'
#' })
#'
#' ## Note that this didn't include sex, add it to dataset and restart model
#'
#'
#' d <- nlmixr2data::theo_sd
#' d$SEX <-0
#' d$SEX[d$ID<=6] <-1
#'
#' fit <- nlmixr2(one.cmt, d, est = "saem", control = list(print = 0))
#'
#' # This would restart if for some reason the search crashed:
#'
#' rxode2::.rxWithWd(tempdir(), {# with temporary directory
#'
#' auto2 <- covarSearchAuto(fit, varsVec = c("ka", "cl"), covarsVec = c("WT"),
#'                 catvarsVec= c("SEX"), restart = TRUE)
#'
#' auto3 <- covarSearchAuto(fit, varsVec = c("ka", "cl"), covarsVec = c("WT"),
#'                 catvarsVec=  c("SEX"), restart = TRUE,
#'                 searchType = "forward")
#' })
#' }

# unsuccessful runs info store; check for covInformation before resuming
covarSearchAuto <- function(fit,
                            varsVec,
                            covarsVec,
                            pVal = list(fwd = 0.05, bck = 0.01), # diff default vals for fwd and backward
                            catvarsVec=NULL ,#character vector of categorical covariates that need to be added
                            searchType = c("scm", "forward", "backward"),
                            restart = FALSE) {
  if (!is.numeric(AIC(fit))) {
    cli::cli_alert_danger("the 'fit' object needs to have an objective functions value associated with it")
    cli::cli_alert_info("try computing 'AIC(fitobject)' in console to compute and store the corresponding OBJF value")
    stop("aborting...objf value not associated with the current 'fit' object", call. = FALSE)
  }

  ## Update data and covarsvec if categorical variables are provided
  if(!is.null(catvarsVec)){
    covarsVec <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[2]]
    data <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[1]]
  } else {
    covarsVec <- covarsVec
    data <- nlme::getData(fit)
  }

  searchType <- match.arg(searchType)

  if (missing(searchType)) {
    searchType <- "scm"
  }

  if (!all((names(pVal) %in% c("fwd", "bck")))) {
    stop("pVal should be list of two with names  'fwd' and 'bck' ")
  }

  outputDir <-
    paste0("nlmixr2CovariateSearchCache_", as.character(substitute(fit)), "_", digest::digest(fit)) # a new directory with this name will be created

  if (!dir.exists(outputDir)) {
    dir.create(outputDir)
  }
  else if (dir.exists(outputDir) && restart == TRUE) {
    unlink(outputDir, recursive = TRUE, force = TRUE) # unlink any of the previous directories
    dir.create(outputDir) # create a fresh directory
  }

  if (searchType %in% "scm") {
    resFwd <- forwardSearch(varsVec,covarsVec,catvarsVec, fit, pVal$fwd, outputDir = outputDir, restart = restart)
    resBck <- backwardSearch(varsVec,covarsVec,catvarsVec, fitorig = fit, fitupdated = resFwd[[1]], pVal = pVal$bck, reFitCovars = FALSE, outputDir = outputDir, restart = restart)
    summaryTable <- Reduce(rbind, list(resFwd[[2]], resBck[[2]]))

    return(list(summaryTable = summaryTable, resFwd = resFwd, resBck = resBck))
  } else if (searchType %in% "forward") {
    resFwd <- forwardSearch(varsVec,covarsVec,catvarsVec, fit, pVal = pVal$fwd, outputDir = outputDir, restart = restart)
    summaryTable <- Reduce(rbind, list(resFwd[[2]], NULL))

    return(list(summaryTable = summaryTable, resFwd = resFwd, resBck = NULL))
  } else {
    resBck <- backwardSearch(varsVec,covarsVec,catvarsVec, fitorig = fit, pVal = pVal$bck, reFitCovars = TRUE, outputDir = outputDir, restart = restart)
    summaryTable <- Reduce(rbind, list(NULL, resBck[[2]]))

    return(list(summaryTable = summaryTable, resFwd = NULL, resBck = resBck))
  }
}

#' Forward covariate search
#'
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec  character vector of covariates that need to be added
#' @param catvarsVec  character vector of categorical covariates that need to be added
#' @param fit  an nlmixr2 'fit' object
#' @param pVal p-value that should be used for selecting covariates in the forward search
#' @param outputDir the name of the output directory that stores the covariate search result
#' @param restart a boolean that controls if the search should be restarted; default is FALSE
#'
#' @return returns the updated 'fit' object at the end of the forward search and a table of information for all the covariates tested
#' @author Vipul Mann, Matthew Fidler, Vishal Sarsani
#' @noRd
forwardSearch <- function(varsVec,covarsVec,catvarsVec=NULL,fit, pVal = 0.05, outputDir, restart = FALSE) {
  ## Update data and covarsvec if categorical variables are provided
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  } else {
    ui <- fit$finalUiEnv
  }
  if(!is.null(catvarsVec)){
    covarsVec <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[2]]
    data <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[1]]
  } else {
    covarsVec <- covarsVec
    data <- nlme::getData(fit)
  }

  if (missing(outputDir)) {
    stop("please specify output directory to store the results for forward search. aborting ...")
  }

  # construct covInfo
  covInfo <-  buildcovInfo(varsVec,covarsVec)

  resTableComplete <- data.frame(matrix(ncol = 14, nrow = 0))
  cli::cli_h1("starting forward search...")
  stepIdx <- 1

  fnameTablePatternForward <-
    paste0("forward_step_", "[0-9]+", "_table_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )
  fnameFitPatternForward <-
    paste0("forward_step_", "[0-9]+", "_fit_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )
  fnameCompleteTablePatternForward <-
    paste0("forward_step_", "[0-9]+", "_completetable_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )

  fileExistsTab <-
    list.files(paste0("./", outputDir), pattern = fnameTablePatternForward)

  fileExistsFit <-
    list.files(paste0("./", outputDir), pattern = fnameFitPatternForward)

  fileExistsCompleteTable <-
    list.files(paste0("./", outputDir), pattern = fnameCompleteTablePatternForward)

  if (length(fileExistsTab) == 0) {
    restart <- TRUE
  }

  if (!restart) {
    resumeTable <- lapply(fileExistsTab, function(x) {
      readRDS(paste0(outputDir, "/", x))
    })

    resumeTable <- data.table::rbindlist(resumeTable)
    fit <- readRDS(paste0(outputDir, "/", fileExistsFit[[length(fileExistsFit)]]))
    resTableComplete <- readRDS(paste0(outputDir, "/", fileExistsCompleteTable[[length(fileExistsCompleteTable)]]))

    # update covInfo and step idx
    testedCovarVars <- paste0(unlist(resumeTable$covar), unlist(resumeTable$var))

    for (x in testedCovarVars) {
      covInfo[[x]] <- NULL
    }

    stepIdx <- unlist(resumeTable[nrow(resumeTable), ]$step) + 1

    cli::cli_alert_success("loaded forward search data from disk ...")
    cli::cli_alert_success("resuming forward search ...")
  }

  while (length(covInfo) > 0) {
    # forward covariate search
    covSearchRes <- buildupatedUI(ui,varsVec,covarsVec,indep = TRUE,add=TRUE)

    resTable <- lapply(covSearchRes, function(res) {
       xmod <- res
       x <- tryCatch(
         {
           x <-
             suppressWarnings(nlmixr2(xmod,data,fit$est))
           x # to return 'model fit'
         },
         error = function(error_message) {
           print("error fitting the model for the covariate ")
           print(error_message)
         })

      covNames <- utils::tail(res$muRefCovariateDataFrame$covariateParameter,n=1)
      nam_var <- strsplit(covNames,split='_', fixed=TRUE)[[1]][3]
      nam_covar <- strsplit(covNames,split='_', fixed=TRUE)[[1]][2]

      # fwd: if deltObjf <0: pchisq=1-pchisq(-deltObjf, dof), else pchisq=1
      # bck: if deltObjf >0: pchisq=1-pchisq(deltObjf, dof), else pchisq=1

      dObjf <- fit$objf - x$objf
      dof <- length(x$finalUiEnv$ini$est) - length(fit$finalUiEnv$ini$est)
      if (dObjf < 0) {
        pchisqr <- 1 - pchisq(-dObjf, df = dof)
      }
      else {
        pchisqr <- 1
      }

      l1 <- list(step = stepIdx, covar = nam_covar, var = nam_var, objf = x$objf, deltObjf = dObjf, AIC = x$AIC, BIC = x$BIC, numParams = length(x$finalUiEnv$ini$est), qchisqr = qchisq(1 - pVal, dof), pchisqr = pchisqr, included = "no", searchType = "forward")
      l2 <- list(covNames = covNames, covarEffect = x$parFixedDf[covNames, "Estimate"])

      c(l1, l2)
    })

    resTable <- data.frame(do.call(rbind, resTable))
    colnames(resTable) <- c(names(resTable))
    bestRow <- resTable[which.min(resTable$pchisqr), ]

    colnames(resTableComplete) <- colnames(resTable)

    if (bestRow$pchisqr <= pVal) { # should be based on p-value

      # objf function value improved
      resTable[which.min(resTable$pchisqr), "included"] <- "yes"
      bestRow[, "included"] <- "yes"

      cli::cli_h1("best model at step {stepIdx}: ")
      print(bestRow)

      fit <-
        covSearchRes[[which.min(resTable$pchisqr)]][[1]] # extract fit object corresponding to the best model

      covInfo[[paste0(as.character(bestRow$covar), as.character(bestRow$var))]] <- NULL

      cli::cli_h2("excluding {paste0(as.character(bestRow$covar), as.character(bestRow$var))} from list of covariates ...")

      saveRDS(fit, file = paste0(outputDir, "/", "forward_", "step_", stepIdx, "_", "fit", "_", paste0(as.character(bestRow$covar), as.character(bestRow$var)), ".RData"))
      saveRDS(bestRow, file = paste0(outputDir, "/", "forward_", "step_", stepIdx, "_", "table", "_", paste0(as.character(bestRow$covar), as.character(bestRow$var)), ".RData"))
      resTableComplete <- rbind(resTableComplete, resTable)
      saveRDS(resTableComplete, file = paste0(outputDir, "/", "forward_", "step_", stepIdx, "_", "completetable", "_", as.character(bestRow$covar), as.character(bestRow$var), ".RData"))

      stepIdx <- stepIdx + 1
    } else {
      # objf function value did not improve
      cli::cli_h1("objf value did not improve, exiting the search ...")

      resTableComplete <- rbind(resTableComplete, resTable)
      break
    }
  }

  cli::cli_h2(cli::col_red("forward search complete"))

  list(fit, resTableComplete)
}

#' Backward covariate search
#'
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec  character vector of covariates that need to be added
#' @param catvarsVec  character vector of categorical covariates that need to be added
#' @param fitorig the original 'fit' object before forward search
#' @param fitupdated the updated 'fit' object, if any, after the forward search
#' @param pVal p-value that should be used for selecting covariates in the forward search
#' @param reFitCovars if the covariates should be added before performing backward search - useful for directly performing backward search without forward search; default is FALSE
#' @param outputDir the name of the output directory that stores the covariate search result
#' @param restart a boolean that controls if the search should be restarted; default is FALSE
#'
#' @return returns the updated 'fit' object at the end of the backward search and a table of information for all the covariates tested
#' @noRd
#'
#' @author Vipul Mann, Matthew Fidler, Vishal Sarsani
backwardSearch <- function(varsVec,covarsVec,catvarsVec=NULL, fitorig, fitupdated, pVal = 0.01, reFitCovars = FALSE, outputDir, restart = FALSE) {

  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  } else {
    ui <- fit$finalUiEnv
  }
  if(!is.null(catvarsVec)){
    covarsVec <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[2]]
    data <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[1]]
  } else {
    covarsVec <- covarsVec
    data <- nlme::getData(fit)
  }

  if (missing(outputDir)) {
    stop("please specify output directory to store the results for backward search. aborting ...")
  }
  cli::cli_h1("starting backward search...")
  resTableComplete <- data.frame(matrix(ncol = 14, nrow = 0))

  stepIdx <- 1

  if (!missing(fitupdated)) {
    if (all(names(fitupdated$ini$theta) %in% names(fitorig$ini$theta))) {
      cli::cli_alert_warning("no covariates added in the forward search, skipping backward search")
      return(list(fitorig, NULL))
    } else {
      fit <- fitupdated
    }
  }

  if (reFitCovars) {
    xmod <- buildupatedUI(ui,varsVec,covarsVec,indep = FALSE,add=TRUE)
    fitupdated <- suppressWarnings(nlmixr2(xmod,data,fit$est)) # get the last fit object with all covariates added # DOES NOT ADD $ini
    fit <- fitupdated
  }

  fnameTablePatternBackward <-
    paste0("backward_step_", "[0-9]+", "_table_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )
  fnameFitPatternBackward <-
    paste0("backward_step_", "[0-9]+", "_fit_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )

  fnameCompleteTablePatternBackward <-
    paste0("forward_step_", "[0-9]+", "_completetable_", "[a-zA-Z0-9]+", ".RData",
      sep = ""
    )

  fileExistsTab <-
    list.files(paste0("./", outputDir), pattern = fnameTablePatternBackward)

  fileExistsFit <-
    list.files(paste0("./", outputDir), pattern = fnameFitPatternBackward)

  fileExistsCompleteTable <-
    list.files(paste0("./", outputDir), pattern = fnameCompleteTablePatternBackward)

  if (length(fileExistsTab) == 0) {
    restart <- TRUE
  }

  if (!restart) {
    resumeTable <- lapply(fileExistsTab, function(x) {
      readRDS(paste0(outputDir, "/", x))
    })

    resumeTable <- data.table::rbindlist(resumeTable)
    fit <- readRDS(paste0(outputDir, "/", fileExistsFit[[length(fileExistsFit)]]))
    resTableComplete <- readRDS(paste0(outputDir, "/", fileExistsCompleteTable[[length(fileExistsCompleteTable)]]))

    # update covInfo and step idx
    testedCovarVars <- paste0(unlist(resumeTable$covar), unlist(resumeTable$var))

    # construct covInfo
    covInfo <-  buildcovInfo(varsVec,covarsVec)

    for (x in testedCovarVars) {
      covInfo[[x]] <- NULL
    }

    stepIdx <- unlist(resumeTable[nrow(resumeTable), ]$step) + 1

    cli::cli_alert_success("loaded backward search data from disk ...")
    cli::cli_alert_success("resuming backward search ...")
  }

  cli::cli_h2(cli::col_blue("initial function text to remove from:"))
  cli::cli_text(cli::col_red("{fit$fun.txt}"))

  # check if covInfo vars in fit; abort if none of the covariates added in the forward search step

  covSearchRes <- buildupatedUI(ui,varsVec,covarsVec,indep = TRUE,add=TRUE)

  # Now remove covars step by step until the objf fun value...?
  while (length(covInfo) > 0) {
    # Remove covars on by one: if objf val increases retain covar; otherwise (objf val decreases), remove the covar
    # At any stage, retain the one that results in highest increase in objf value; exit if removal of none results in increase

    covSearchRes <- buildupatedUI(ui,varsVec,covarsVec,add=FALSE)

    resTable <- lapply(covSearchRes, function(res) {
      xmod <- res
      x <- tryCatch(
        {
          x <-
            suppressWarnings(nlmixr2(xmod,data,fit$est))
          x # to return 'model fit'
        },
        error = function(error_message) {
          print("error fitting the model for the covariate ")
          print(error_message)
        })

      covNames <- utils::tail(res$muRefCovariateDataFrame$covariateParameter,n=1)
      nam_var <- strsplit(covNames,split='_', fixed=TRUE)[[1]][3]
      nam_covar <- strsplit(covNames,split='_', fixed=TRUE)[[1]][2]
      # fwd: if deltObjf <0: pchisq=1-pchisq(-deltObjf, dof), else pchisq=1
      # bck: if deltObjf >0: pchisq=1-pchisq(deltObjf, dof), else pchisq=1

      dObjf <- x$objf - fit$objf
      dof <- length(fit$finalUiEnv$ini$est) - length(x$finalUiEnv$ini$est)

      if (dObjf > 0) {
        pchisqr <- 1 - pchisq(dObjf, df = dof)
      }
      else {
        pchisqr <- 1
      }

      l1 <- list(step = stepIdx, covar = nam_covar, var = nam_var, objf = x$objf, deltObjf = dObjf, AIC = x$AIC, BIC = x$BIC, numParams = length(x$finalUiEnv$ini$est), qchisqr = qchisq(1 - pVal, dof), pchisqr = pchisqr, included = "", searchType = "backward")
      l2 <- list(covNames = covNames, covarEffect = fit$parFixedDf[covNames, "Estimate"])

      c(l1, l2)
    })

    resTable <- data.frame(do.call(rbind, resTable))
    colnames(resTable) <- c(names(resTable))

    colnames(resTableComplete) <- colnames(resTable)

    bestRow <- resTable[which.min(resTable$pchisqr), ]


    if (bestRow$pchisqr <= pVal) {
      # objf function value increased after removal of covariate: retain the best covariate at this stage, test for the rest

      resTable[which.min(resTable$pchisqr), "included"] <- "yes"
      bestRow[, "included"] <- "yes"


      cli::cli_h1("best model at step {stepIdx}: ")
      print(bestRow)

      fit <-
        covSearchRes[[which.min(resTable$pchisqr)]][[1]] # extract fit object corresponding to the best model
      covInfo[[paste0(as.character(bestRow$covar), as.character(bestRow$var))]] <- NULL

      saveRDS(fit, file = paste0(outputDir, "/", "backward_", "step_", stepIdx, "_", "fit", "_", paste0(as.character(bestRow$covar), as.character(bestRow$var)), ".RData"))
      saveRDS(bestRow, file = paste0(outputDir, "/", "backward_", "step_", stepIdx, "_", "table", "_", paste0(as.character(bestRow$covar), as.character(bestRow$var)), ".RData"))
      cli::cli_h2("retaining {paste0(as.character(bestRow$covar), as.character(bestRow$var))}")
      resTableComplete <- rbind(resTableComplete, resTable)
      saveRDS(resTableComplete, file = paste0(outputDir, "/", "backward_", "step_", stepIdx, "_", "completetable", "_", as.character(bestRow$covar), as.character(bestRow$var), ".RData"))

      stepIdx <- stepIdx + 1
    } else {
      # objf function value did not improve
      cli::cli_h1("objf value did not improve, exiting the search ...")
      resTableComplete <- rbind(resTableComplete, resTable)
      #saveRDS(resTableComplete, file = paste0(outputDir, "/", "backward_", "step_", stepIdx, "_", "completetable", "_", as.character(bestRow$covar), as.character(bestRow$var), ".RData"))

      break
    }
  }

  cli::cli_h2(cli::col_red("backward search complete"))

  list(fit, resTableComplete)
}
