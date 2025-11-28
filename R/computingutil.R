# Function to calculate population standard deviation from group means
#' @param x vector of values to calculate standard deviation.
#' @return population standard deviation
#' @noRd
.sd.p <- function(x) {
  sd(x)*sqrt((length(x)-1)/length(x))
}

.bootstrapEnv <- new.env(parent=emptyenv())
.bootstrapEnv$nSampIndiv <- 0L

#' Function to return pop mean, pop std of a given covariate
#'
#' @param data given data frame
#' @param covariate the covariate that needs popmean,stdev
#'
#' @return list containing values of population mean, standard deviation
#' @noRd
.popMeanStd <- function(data, covariate) {

  checkmate::assertDataFrame(data,col.names = "named")
  checkmate::assertCharacter(covariate,len = 1,any.missing = FALSE )

  if (inherits(try(str2lang(covariate)),"try-error")) {
    stop("`varName` must be a valid R expression",call. = FALSE)
  }

  .new <- intersect(names(data), covariate)
  if (length(.new) == 0L) stop("covariate specified not in original dataset", call.=FALSE)

  #extract Individual ID from data frame

  uidCol <- .idColumn(data)
  # mean by groups (Individual)
  groupMeans <- with(data, ave(get(covariate),get(uidCol), FUN = function(x) mean(x, na.rm = TRUE)))
  # pop mean
  popMean <- mean(groupMeans, na.rm=TRUE)
  # pop std
  popStd <- .sd.p (groupMeans)

  .meanStd <-c(popMean,popStd)
  names(.meanStd) <- c("popmean","popstd")
  .meanStd
}

#' Function to return normalization column of a covariates
#'
#' @param data given a data frame
#' @param covariate the covariate that needs normalization
#'
#' @return data frame with normalized covariate
#' @noRd
.normalizeDf <- function(data, covariate,sub=TRUE) {
  checkmate::assertDataFrame(data,col.names = "named")
  checkmate::assertCharacter(covariate,len = 1,any.missing = FALSE )

  if (inherits(try(str2lang(covariate)),"try-error")) {
    stop("`varName` must be a valid R expression",call. = FALSE)
  }
  .new <- intersect(names(data), covariate)
  if (length(.new) == 0L) stop("covariate specified not in original dataset")

  if (is.factor(data[[covariate]])) {
    data
  } else {
    # Column name for the standardized covariate
    datColNames <- paste0("normalized_", covariate)
    # popMean
    .popMean = .popMeanStd(data,covariate)[[1]]
    # popStdev
    .popStd = .popMeanStd(data,covariate)[[2]]
    # add standardized covariate values to the data frame
    data[,datColNames] <- (data[,covariate]-.popMean)/(.popStd)
    data
  }
}

#' Function to return data of normalized covariates
#'
#' @param data a dataframe with covariates to normalize
#' @param covarsVec a list of covariate names (parameters) that need to be estimates
#' @param replace replace the original covariate data with normalized data for easier updated model.
#'
#' @return data frame with all normalized covariates
#' @author Vishal Sarsani
#' @export
#'
#' @examples
#'
#' d <- nlmixr2data::theo_sd
#' d$SEX <-0
#' d$SEX[d$ID<=6] <-1
#'
#' covarsVec <- c("WT")
#'
#' # Normalized covariate (replaced)
#' df1 <- normalizedData(d,covarsVec,replace=TRUE)
#'
#' # Normalized covariate (without replacement)
#' df2 <- normalizedData(d,covarsVec,replace=FALSE)
normalizedData <- function(data,covarsVec,replace=TRUE) {
  checkmate::assert_character(covarsVec)
  .normalizedDFs <- lapply(covarsVec,.normalizeDf,data=data)

  # final data frame of normalized covariates
  if(replace) {
    .dat <- Reduce(merge,.normalizedDFs)
    dropnormPrefix <- function(x) {
      colnames(x) <- gsub("normalized_", "", colnames(x))
      x
    }
    catCheck <- intersect(covarsVec,names(Filter(is.factor, data)))
    .dat <- cbind(.dat[ , !names(.dat) %in% covarsVec],subset(.dat,select=catCheck))
    .finalDf <- dropnormPrefix(.dat)
  } else {
    .finalDf <- Reduce(merge,.normalizedDFs)
  }
  .finalDf
}

#' Stratified cross-validation fold generator function, inspired from the caret
#'
#' @param data data frame used in the analysis
#' @param nfold number of k-fold cross validations. Default is 5
#' @param stratVar  Stratification Variable. Default is NULL and ID is used for CV
#' @return return data.frame with the fold column attached
#' @author Vishal Sarsani, caret
#' @export
#'
#' @examples
#' d <- nlmixr2data::theo_sd
#' d$SEX <-0
#' d$SEX[d$ID<=6] <-1
#'
#' covarsVec <- c("WT")
#'
#' # Stratified cross-validation data with CMT
#' df1 <- foldgen(d, nfold=5, stratVar="CMT")
#'
#' # Stratified cross-validation data with ID (individual)
#' df2 <- foldgen(d, nfold=5, stratVar=NULL)
foldgen <-  function(data,nfold=5,stratVar=NULL) {
  # check if data frame
  checkmate::assert_data_frame(data,min.cols = 7)
  # check if user want to stratify on a variable , if not default is on individual
  if(!is.null(stratVar)) {
    checkmate::assertCharacter(stratVar,len = 1,any.missing = FALSE)
    stratCheck <- intersect(names(data), stratVar)
    if(!is.null(stratCheck)) {
      y <- data[,stratCheck]
    } else {
      stop(paste0(stratVar, "not in the data to stratify"),
           call.=FALSE)
    }
  } else {
    # extract ID column from the data frame
    ID <- .idColumn(data)
    # Extract list of individuals
    y <- unique(data[,ID])
  }
  ## Group based on magnitudes and sample within groups
  if(is.numeric(y)) {
    ## Group the numeric data based on their magnitudes
    ## and sample within those groups.

    ## When the number of samples is low, we may have
    ## issues further slicing the numeric data into
    ## groups. The number of groups will depend on the
    ## ratio of the number of folds to the sample size.
    ## At most, we will use quantiles. If the sample
    ## is too small, we just do regular unstratified
    ## CV
    cuts <- floor(length(y)/nfold)
    if(cuts < 2) cuts <- 2
    if(cuts > 5) cuts <- 5
    y <- cut(y,
             unique(quantile(y,
                             probs = seq(0, 1, length = cuts),
                             na.rm=TRUE)),
             include.lowest = TRUE)
  }

  if (nfold < length(y)) {
    ## reset levels so that the possible levels and
    ## the levels in the vector are the same
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))

    ## For each class, balance the fold allocation as far
    ## as possible, then resample the remainder.
    ## The final assignment of folds is also randomized.
    for(i in 1:seq_along(numInClass)) {
      ## create a vector of integers from 1:k as many times as possible without
      ## going over the number of samples in the class. Note that if the number
      ## of samples in a class is less than k, nothing is producd here.
      seqVector <- rep(1:nfold, numInClass[i] %/% nfold)
      ## add enough random integers to get  length(seqVector) == numInClass[i]
      if(numInClass[i] %% nfold > 0) {
        seqVector <- c(seqVector, sample(1:nfold, numInClass[i] %% nfold))
      }
      ## shuffle the integers for fold assignment and assign to this classes's data
      foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
    }
  } else {
    foldVector <- seq(along = y)
  }

  out <- split(seq(along = y), foldVector)
  names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), sep = "")
  out <- foldVector

  if(!is.null(stratVar)) {
    out <- cbind(data,fold=out)
  } else {
    indv <- unique(data[,ID])
    out <- merge(data,cbind(ID=indv,fold=out))
  }
  out
}

#' Sample from uniform distribution by optim
#'
#' @param xvec A vector of min,max values . Ex:c(10,20)
#' @param N Desired number of values
#' @param medValue  Desired Median
#' @param floorT boolean indicating whether to round up
#'
#' @return Samples with approx desired median.
#' @author Vishal Sarsani
#' @export
#'
#' @examples
#' # Simulate 1000 creatine clearance values with median of 71.7 within range of c(6.7,140)
#' creatCl <- optimUnisampling(xvec=c(6.7,140), N=1000, medValue = 71.7, floorT=FALSE)
optimUnisampling <- function(xvec,N=1000,medValue,floorT=TRUE) {
  #Function to calculate distance between sampling median and desired
  fun <- function(xvec, N=1000) {
    xmin <- xvec[1]
    xmax <- xvec[2]
    if (floorT) {
      x <- floor(stats::runif(N, xmin, xmax))
    } else {
      x <- stats::runif(N, xmin, xmax)
    }
    xdist <- (median(x)-medValue)^2
    xdist
  }
  # Optimization
  xr <- stats::optim(xvec, fun)
  xrmin <- xr$par[[1]]
  xrmax <- xr$par[[2]]
  sampled <- stats::runif(N, min = xr$par[[1]], max = xr$par[[2]])
  if (xrmin==xvec[1] && xrmax==xvec[2] && floorT) {
    return(floor(sampled))
  }
  else if (xrmin==xvec[1] && xrmax==xvec[2]) {
    return(sampled)
  }
  optimUnisampling(xvec,N=1000,medValue)
}

#' Format confidence bounds for a variable into bracketed notation using string formatting
#'
#' @param var a list of values for the varaible
#' @param confLower the lower bounds for each of the values
#' @param confUpper the upper bounds for each of the values
#' @param sigdig the number of significant digits
#' @author Vipul Mann
#' @noRd
addConfboundsToVar <- function(var, confLower, confUpper, sigdig = 3) {
  res <- lapply(seq_along(var), function(idx) {
    paste0(
      signif(var[idx], sigdig),
      " (",
      signif(confLower[idx], sigdig),
      ", ",
      signif(confUpper[idx], sigdig),
      ")"
    )
  })
  unlist(res)
}

#' Bootstrap nlmixr2 fit
#'
#' Bootstrap input dataset and rerun the model to get confidence bounds and aggregated parameters
#'
#' @param fit the nlmixr2 fit object
#' @param nboot an integer giving the number of bootstrapped models to
#'   be fit; default value is 200
#' @param nSampIndiv an integer specifying the number of samples in
#'   each bootstrapped sample; default is the number of unique
#'   subjects in the original dataset
#' @param stratVar Variable in the original dataset to stratify on;
#'   This is useful to distinguish between sparse and full sampling
#'   and other features you may wish to keep distinct in your
#'   bootstrap
#' @param pvalues a vector of pvalues indicating the probability of
#'   each subject to get selected; default value is NULL implying that
#'   probability of each subject is the same
#' @param restart a boolean that indicates if a previous session has
#'   to be restarted; default value is FALSE
#' @param fitName Name of fit to be saved (by default the variable
#'   name supplied to fit)
#' @param stdErrType This gives the standard error type for the
#'   updated standard errors; The current possibilities are: `"perc"`
#'   which gives the standard errors by percentiles (default), `"sd"`
#'   which gives the standard errors by the using the normal
#'   approximation of the mean with standard devaition, or `"se"`
#'   which uses the normal approximation with standard errors
#'   calculated with `nSampIndiv`
#' @param ci Confidence interval level to calculate.  Default is 0.95
#'   for a 95 percent confidence interval
#' @param plotHist A boolean indicating if a histogram plot to assess
#'   how well the bootstrap is doing.  By default this is turned off
#'   (`FALSE`)
#' @param pvalues a vector of pvalues indicating the probability of
#'   each subject to get selected; default value is `NULL` implying
#'   that probability of each subject is the same
#' @param restart A boolean to try to restart an interrupted or
#'   incomplete boostrap.  By default this is `FALSE`
#' @param fitName is the fit name that is used for the name of the
#'   boostrap files.  By default it is the fit provided though it
#'   could be something else.
#'
#' @param returnType this describes the return type
#'
#' - `model` this is the default and returns an updated model boostrapped
#'   confidence intervals and standard errors updated in the estimates.
#'
#' - `fitList` this returns the list of the full fits from the bootstrap.
#'
#' - `modelList` this returns the model list (which abbreviates the
#'    parameters and messages) used for the bootstrap summary
#'    statistics.
#'
#' @author Vipul Mann, Matthew Fidler
#' @return Nothing, called for the side effects; The original fit is
#'   updated with the bootstrap confidence bands
#' @export
#' @examples
#' \dontrun{
#'
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45; label("Ka")
#'     tcl <- 1; label("Cl")
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
#'
#' withr::with_tempdir({ # Run example in temp dir
#'
#' bootstrapFit(fit, nboot = 5, restart = TRUE) # overwrites any of the existing data or model files
#' bootstrapFit(fit, nboot = 7) # resumes fitting using the stored data and model files
#'
#' # Note this resumes because the total number of bootstrap samples is not 10
#'
#' bootstrapFit(fit, nboot=10)
#'
#' # Note the boostrap standard error and variance/covariance matrix is retained.
#' # If you wish to switch back you can change the covariance matrix by
#'
#' nlmixr2est::setCov(fit, "linFim")
#'
#' # And change it back again
#'
#' nlmixr2est::setCov(fit, "boot10")
#'
#' # This change will affect any simulations with uncertainty in their parameters
#'
#' # You may also do a chi-square diagnostic plot check for the bootstrap with
#' bootplot(fit)
#' })
#' }
bootstrapFit <- function(fit,
                         nboot = 200,
                         nSampIndiv,
                         stratVar,
                         stdErrType = c("perc", "sd", "se"),
                         ci = 0.95,
                         pvalues = NULL,
                         restart = FALSE,
                         plotHist = FALSE,
                         fitName = as.character(substitute(fit)),
                         returnType=c("model", "fitList", "modelList")) {

  stdErrType <- match.arg(stdErrType)
  returnType <- match.arg(returnType)
  checkmate::assertNumeric(ci, lower=0, upper=1, len=1, any.missing=FALSE, null.ok = FALSE)

  if (missing(stratVar)) {
    performStrat <- FALSE
  }
  else {
    if (!(stratVar %in% colnames(nlme::getData(fit)))) {
      cli::cli_alert_danger("{stratVar} not in data")
      stop("aborting ...stratifying variable not in data", call. = FALSE)
    }
    performStrat <- TRUE
  }

  if (is.null(fit$bootstrapMd5)) {
    bootstrapMd5 <- fit$md5
    assign("bootstrapMd5", bootstrapMd5, envir = fit$env)
  }

  if (performStrat) {
    resBootstrap <-
      modelBootstrap(
        fit,
        nboot = nboot,
        nSampIndiv = nSampIndiv,
        stratVar = stratVar,
        pvalues = pvalues,
        restart = restart,
        fitName = fitName
      ) # multiple models

    modelsList <- resBootstrap[[1]]
    fitList <- resBootstrap[[2]]
  }
  else {
    resBootstrap <-
      modelBootstrap(
        fit,
        nboot = nboot,
        nSampIndiv = nSampIndiv,
        pvalues = pvalues,
        restart = restart,
        fitName = fitName
      ) # multiple models
    modelsList <- resBootstrap[[1]]
    fitList <- resBootstrap[[2]]
  }

  if (returnType == "fitList") {
    return(fitList)
  } else if (returnType == "modelList") {
    return(modelsList)
  }


  bootSummary <-
    getBootstrapSummary(modelsList, ci=ci, stdErrType=stdErrType) # aggregate values/summary

  # modify the fit object
  nrws <- nrow(bootSummary$parFixedDf$mean)
  sigdig <- fit$control$sigdigTable

  newParFixedDf <- fit$parFixedDf
  newParFixed <- fit$parFixed

  # Add Estimate_boot
  est <- unname(bootSummary$parFixedDf$mean[1:nrws, 1])
  cLower <- unname(bootSummary$parFixedDf$confLower[1:nrws, 1])
  cUpper <- unname(bootSummary$parFixedDf$confUpper[1:nrws, 1])
  estEst <- est

  estimateBoot <- addConfboundsToVar(est, cLower, cUpper, sigdig)

  # Add SE_boot
  seBoot <- unname(bootSummary$parFixedDf$stdDev[1:nrws, 1])

  # Add Back-transformed
  est <- unname(bootSummary$parFixedDf$mean[1:nrws, 2])
  cLowerBT <- unname(bootSummary$parFixedDf$confLower[1:nrws, 2])
  cUpperBT <- unname(bootSummary$parFixedDf$confUpper[1:nrws, 2])
  backTransformed <-
    addConfboundsToVar(est, cLowerBT, cUpperBT, sigdig)
  estBT <- est

  newParFixedDf["Bootstrap Estimate"] <- estEst
  newParFixedDf["Bootstrap SE"] <- seBoot
  newParFixedDf["Bootstrap %RSE"] <- seBoot / estEst * 100
  newParFixedDf["Bootstrap CI Lower"] <- cLowerBT
  newParFixedDf["Bootstrap CI Upper"] <- cUpperBT
  newParFixedDf["Bootstrap Back-transformed"] <- estBT

  newParFixed["Bootstrap Estimate"] <- estimateBoot
  newParFixed["Bootstrap SE"] <- signif(seBoot, sigdig)
  newParFixed["Bootstrap %RSE"] <-
    signif(seBoot / estEst * 100, sigdig)
  .w <- which(regexpr("^Bootstrap +Back[-]transformed", names(newParFixed)) != -1)
  if (length(.w) >= 1) {
    newParFixed <- newParFixed[, -.w]
  }
  newParFixed[sprintf("Bootstrap Back-transformed(%s%%CI)", ci * 100)] <-
    backTransformed

  # compute bias
  bootParams <- bootSummary$parFixedDf$mean
  origParams <- data.frame(list("Estimate" = fit$parFixedDf$Estimate, "Back-transformed" = fit$parFixedDf$`Back-transformed`))
  bootstrapBiasParfixed <- abs(origParams - bootParams)
  bootstrapBiasOmega <- abs(fit$omega - bootSummary$omega$mean)

  assign("bootBiasParfixed", bootstrapBiasParfixed, envir = fit$env)
  assign("bootBiasOmega", bootstrapBiasOmega, envir = fit$env)

  assign("bootCovMatrix", bootSummary$omega$covMatrix, envir = fit$env)
  assign("bootCorMatrix", bootSummary$omega$corMatrix, envir = fit$env)
  assign("parFixedDf", newParFixedDf, envir = fit$env)
  assign("parFixed", newParFixed, envir = fit$env)
  assign("bootOmegaSummary", bootSummary$omega, envir = fit$env)
  assign("bootSummary", bootSummary, envir = fit$env)

  # plot histogram
  if (plotHist) {

    # compute delta objf values for each of the models
    origData <- nlme::getData(fit)

    if (is.null(fit$bootstrapMd5)) {
      bootstrapMd5 <- fit$md5
      assign("bootstrapMd5", bootstrapMd5, envir = fit$env)
    }

    # already exists
    output_dir <- paste0("nlmixr2BootstrapCache_", fitName, "_", fit$bootstrapMd5)

    deltOBJFloaded <- NULL
    deltOBJF <- NULL
    rxode2::rxProgress(length(fitList))
    cli::cli_h1("Loading/Calculating \u0394 Objective function")
    nlmixr2est::setOfv(fit, "focei") # Make sure we are using focei objective function
    deltOBJF <- lapply(seq_along(fitList), function(i) {
      x <- readRDS(file.path(output_dir, paste0("fitEnsemble_", i, ".rds")))
      .path <- file.path(output_dir, paste0("posthoc_", i, ".rds"))
      if (file.exists(.path)) {
        xPosthoc <- readRDS(.path)
        rxode2::rxTick()
      } else {
        rxode2::rxProgressStop()
        ## rxode2::rxProgressAbort("Starting to posthoc estimates")
        ## Don't calculate the tables
        .msg <- paste0(gettext("Running bootstrap estimates on original data for model index: "), i)
        cli::cli_h1(.msg)
        xPosthoc <- nlmixr2(x,
                            data = origData, est = "posthoc",
                            control = list(calcTables = FALSE, print = 1, compress=FALSE)
                            )
        saveRDS(xPosthoc, .path)
      }
      xPosthoc$objf - fit$objf
    })
    rxode2::rxProgressStop()

    .deltaO <- sort(abs(unlist(deltOBJF)))

    .deltaN <- length(.deltaO)

    .df <- length(fit$ini$est)

    .chisq <- rbind(
      data.frame(
        deltaofv = qchisq(seq(0, 0.99, 0.01), df = .df),
        quantiles = seq(0, 0.99, 0.01),
        Distribution = 1L,
        stringsAsFactors = FALSE
      ),
      data.frame(
        deltaofv = .deltaO,
        quantiles = seq(.deltaN) / .deltaN,
        Distribution = 2L,
        stringsAsFactors = FALSE
      )
    )

    .fdelta <- approxfun(seq(.deltaN) / .deltaN, .deltaO)

    .df2 <- round(mean(.deltaO, na.rm = TRUE))

    .dfD <- data.frame(
      label = paste(c("df\u2248", "df="), c(.df2, .df)),
      Distribution = c(2L, 1L),
      quantiles = 0.7,
      deltaofv = c(.fdelta(0.7), qchisq(0.7, df = .df))
    )

    .dfD$Distribution <- factor(
      .dfD$Distribution, c(1L, 2L),
      c("Reference distribution", "\u0394 objective function")
    )

    .chisq$Distribution <- factor(
      .chisq$Distribution, c(1L, 2L),
      c("Reference distribution", "\u0394 objective function")
    )
    .dataList <- list(
      dfD = .dfD, chisq = .chisq,
      deltaN = .deltaN, df2 = .df2
    )
    assign(".bootPlotData", .dataList, envir = fit$env)
  }
  ## Update covariance estimate
  .nm <- names(fit$theta)[!fit$foceiSkipCov[seq_along(fit$theta)]]
  .nm <- .nm[.nm %in% dimnames(fit$bootSummary$omega$covMatrixCombined)[[1]], drop=FALSE]
  if (length(.nm) == 0) {
    stop("No parameters to update covariance matrix", call.=FALSE)
  }
  .cov <- fit$bootSummary$omega$covMatrixCombined[.nm, .nm]
  .setCov(fit, covMethod = .cov)
  assign("covMethod", paste0("boot", fit$bootSummary$nboot), fit$env)
  invisible(fit)
}

#' Perform bootstrap-sampling from a given dataframe
#'
#' @param data the original dataframe object to sample from for bootstrapping
#'
#' @param nsamp an integer specifying the number of samples in each
#'   bootstrapped sample; default is the number of unique subjects in
#'   the original dataset
#'
#' @param uid_colname a string representing the unique ID of each
#'   subject in the data; default values is 'ID'
#'
#' @param pvalues a vector of pvalues indicating the probability of
#'   each subject to get selected; default value is NULL implying that
#'   probability of each subject is the same
#'
#' @return returns a bootstrap sampled dataframe object
#' @author Vipul Mann, Matthew Fidler
#'
#' @examples
#' sampling(data)
#' sampling(data, 10)
#' @noRd
sampling <- function(data,
                     nsamp=NULL,
                     uid_colname,
                     pvalues = NULL,
                     performStrat = FALSE,
                     stratVar) {
  checkmate::assert_data_frame(data)
  if (is.null(nsamp)) {
    nsamp <- length(unique(data[, uid_colname]))
  }
  else {
    checkmate::assert_integerish(nsamp,
                                 len = 1,
                                 any.missing = FALSE,
                                 lower = 2
                                 )
  }

  if (performStrat && missing(stratVar)) {
    print("stratVar is required for stratifying")
    stop("aborting... stratVar not specified", call. = FALSE)
  }

  checkmate::assert_integerish(nsamp,
                               lower = 2,
                               len = 1,
                               any.missing = FALSE
                               )

  if (missing(uid_colname)) {
    # search the dataframe for a column name of 'ID'
    colNames <- colnames(data)
    colNamesLower <- tolower(colNames)
    if ("id" %in% colNames) {
      uid_colname <- colNames[which("id" %in% colNamesLower)]
    }
    else {
      uid_colname <- "ID"
    }
  }
  else {
    checkmate::assert_character(uid_colname)
  }

  if (performStrat) {
    stratLevels <-
      as.character(unique(data[, stratVar])) # char to access freq. values

    dataSubsets <- lapply(stratLevels, function(x) {
      data[data[, stratVar] == x, ]
    })

    names(dataSubsets) <- stratLevels

    tab <- table(data[stratVar])
    nTab <- sum(tab)

    sampledDataSubsets <- lapply(names(dataSubsets), function(x) {
      dat <- dataSubsets[[x]]

      uids <- unique(dat[, uid_colname])
      uids_samp <- sample(
        list(uids),
        size = ceiling(nsamp * unname(tab[x]) / nTab),
        replace = TRUE,
        prob = pvalues
      )

      sampled_df <-
        data.frame(dat)[0, ] # initialize an empty dataframe with the same col names

      # populate dataframe based on sampled uids
      # new_id = 1
      .env <- environment()
      .env$new_id <- 1
      do.call(rbind, lapply(uids_samp, function(u) {
        data_slice <- dat[dat[, uid_colname] == u, ]
        start <- NROW(sampled_df) + 1
        end <- start + NROW(data_slice) - 1

        data_slice[uid_colname] <-
          .env$new_id # assign a new ID to the sliced dataframe
        .env$new_id <- .env$new_id + 1
        data_slice
      }))
    })
    do.call("rbind", sampledDataSubsets)
  }

  else {
    uids <- unique(data[, uid_colname])
    uids_samp <-
      sample(
        uids,
        size = nsamp,
        replace = TRUE,
        prob = pvalues
      )

    sampled_df <-
      data.frame(data)[0, ] # initialize an empty dataframe with the same col names

    # populate dataframe based on sampled uids
    # new_id = 1
    .env <- environment()
    .env$new_id <- 1

    do.call(rbind, lapply(uids_samp, function(u) {
      data_slice <- data[data[, uid_colname] == u, ]
      start <- NROW(sampled_df) + 1
      end <- start + NROW(data_slice) - 1

      data_slice[uid_colname] <-
        .env$new_id # assign a new ID to the sliced dataframe
      .env$new_id <- .env$new_id + 1
      data_slice
    }))
  }
}

#' Fitting multiple bootstrapped models without aggregaion; called by the function bootstrapFit()
#'
#' @param fit the nlmixr2 fit object
#' @param nboot an integer giving the number of bootstrapped models to be fit; default value is 100
#' @param nSampIndiv an integer specifying the number of samples in each bootstrapped sample; default is the number of unique subjects in the original dataset
#' @param pvalues a vector of pvalues indicating the probability of each subject to get selected; default value is NULL implying that probability of each subject is the same
#' @param restart a boolean that indicates if a previous session has to be restarted; default value is FALSE
#'
#' @return a list of lists containing the different attributed of the fit object for each of the bootstrapped models
#' @author Vipul Mann, Matthew Fidler
#' @examples
#' modelBootstrap(fit)
#' modelBootstrap(fit, 5)
#' modelBootstrap(fit, 5, 20)
#' @noRd
modelBootstrap <- function(fit,
                           nboot = 100,
                           nSampIndiv=NULL,
                           stratVar,
                           pvalues = NULL,
                           restart = FALSE,
                           fitName = "fit") {
  nlmixr2est::assertNlmixrFit(fit)
  if (missing(stratVar)) {
    performStrat <- FALSE
    stratVar <- NULL
  } else {
    performStrat <- TRUE
  }

  data <- nlme::getData(fit)

  .w <- tolower(names(data)) == "id"
  uidCol <- names(data)[.w]

  checkmate::assert_integerish(nboot,
                               len = 1,
                               any.missing = FALSE,
                               lower = 1
                               )

  if (missing(nSampIndiv)) {
    nSampIndiv <- length(unique(data[, uidCol]))
  } else {
    checkmate::assert_integerish(
      nSampIndiv,
      len = 1,
      any.missing = FALSE,
      lower = 2
    )
  }

  # infer the ID column from data
  colNames <- names(data)
  colNamesLower <- tolower(colNames)
  if ("id" %in% colNamesLower) {
    uid_colname <- colNames[which("id" %in% colNamesLower)]
  } else {
    stop("cannot find the 'ID' column! aborting ...", call. = FALSE)
  }

  ui <- fit$finalUiEnv
  fitMeth <- getFitMethod(fit)

  bootData <- vector(mode = "list", length = nboot)

  if (is.null(fit$bootstrapMd5)) {
    bootstrapMd5 <- fit$md5
    assign("bootstrapMd5", bootstrapMd5, envir = fit$env)
  }

  output_dir <-
    paste0("nlmixr2BootstrapCache_", fitName, "_", fit$bootstrapMd5) # a new directory with this name will be created

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  } else if (dir.exists(output_dir) && restart == TRUE) {
    unlink(output_dir, recursive = TRUE, force = TRUE) # unlink any of the previous directories
    dir.create(output_dir) # create a fresh directory
  }

  fnameBootDataPattern <-
    paste0("boot_data_", "[0-9]+", ".rds",
           sep = ""
           )
  fileExists <-
    list.files(paste0("./", output_dir), pattern = fnameBootDataPattern)

  if (length(fileExists) == 0) {
    restart <- TRUE
  }

  if (!restart) {
    # read saved bootData from boot_data files on disk
    if (length(fileExists) > 0) {
      cli::cli_alert_success("resuming bootstrap data sampling using data at {paste0('./', output_dir)}")

      bootData <- lapply(fileExists, function(x) {
        readRDS(paste0("./", output_dir, "/", x, sep = ""))
      })

      startCtr <- length(bootData) + 1
    } else {
      cli::cli_alert_danger(
        cli::col_red(
          "need the data files at {.file {paste0(getwd(), '/', output_dir)}} to resume"
        )
      )
      stop("aborting...resume file missing", call. = FALSE)
    }
  } else {
    startCtr <- 1
  }

  # Generate additional samples (if nboot>startCtr)
  if (nboot >= startCtr) {
    for (mod_idx in startCtr:nboot) {
      bootData[[mod_idx]] <- sampling(
        data,
        nsamp = nSampIndiv,
        uid_colname = uidCol,
        pvalues = pvalues,
        performStrat = performStrat,
        stratVar = stratVar
      )

      # save bootData in curr directory: read the file using readRDS()
      attr(bootData, "randomSeed") <- .Random.seed
      saveRDS(bootData[[mod_idx]],
              file = paste0(
                "./",
                output_dir,
                "/boot_data_",
                mod_idx,
                ".rds"))
    }
  }

  # check if number of samples in stored file is the same as the required number of samples
  fileExists <-
    list.files(paste0("./", output_dir), pattern = fnameBootDataPattern)
  bootData <- lapply(fileExists, function(x) {
    readRDS(paste0("./", output_dir, "/", x, sep = ""))
  })

  currBootData <- length(bootData)

  # Fitting models to bootData now
  .env <- environment()
  fnameModelsEnsemblePattern <-
    paste0("modelsEnsemble_", "[0-9]+",
           ".rds",
           sep = "")
  modFileExists <-
    list.files(paste0("./", output_dir), pattern = fnameModelsEnsemblePattern)

  fnameFitEnsemblePattern <-
    paste0("fitEnsemble_", "[0-9]+",
           ".rds",
           sep = "")
  fitFileExists <- list.files(paste0("./", output_dir), pattern = fnameFitEnsemblePattern)

  if (!restart) {
    if (length(modFileExists) > 0 &&
          (length(fileExists) > 0)) {

      # read bootData and modelsEnsemble files from disk
      cli::cli_alert_success(
        "resuming bootstrap model fitting using data and models stored at {paste0(getwd(), '/', output_dir)}"
      )

      bootData <- lapply(fileExists, function(x) {
        readRDS(paste0("./", output_dir, "/", x, sep = ""))
      })
      modelsEnsembleLoaded <- lapply(modFileExists, function(x) {
        readRDS(paste0("./", output_dir, "/", x, sep = ""))
      })

      fitEnsembleLoaded <- lapply(fitFileExists, function(x) {
        readRDS(paste0("./", output_dir, "/", x, sep = ""))
      })

      .env$mod_idx <- length(modelsEnsembleLoaded) + 1

      currNumModels <- .env$mod_idx - 1

      if (currNumModels > nboot) {
        mod_idx_m1 <- .env$mod_idx-1
        cli::cli_alert_danger(
          cli::col_red(
            "the model file already has {mod_idx_m1} models when max models is {nboot}; using only the first {nboot} model(s)"
          )
        )
        return(list(modelsEnsembleLoaded[1:nboot], fitEnsembleLoaded[1:nboot]))

        # return(modelsEnsembleLoaded[1:nboot])
      } else if (currNumModels == nboot) {
        mod_idx_m1 <- .env$mod_idx-1
        cli::col_red(
          "the model file already has {mod_idx_m1} models when max models is {nboot}; loading from {nboot} models already saved on disk"
        )
        return(list(modelsEnsembleLoaded, fitEnsembleLoaded))

        # return(modelsEnsembleLoaded)
      } else if (currNumModels < nboot) {
        cli::col_red("estimating the additional models ... ")
      }
    }

    else {
      cli::cli_alert_danger(
        cli::col_red(
          "need both the data and the model files at: {paste0(getwd(), '/', output_dir)} to resume"
        )
      )
      stop(
        "aborting...data and model files missing at: {paste0(getwd(), '/', output_dir)}",
        call. = FALSE
      )
    }
  } else {
    .env$mod_idx <- 1
  }

  # get control settings for the 'fit' object and save computation effort by not computing the tables
  .ctl <- setQuietFastControl(fit$control)

  modelsEnsemble <-
    lapply(bootData[.env$mod_idx:nboot], function(boot_data) {
      modIdx <- .env$mod_idx
      cli::cli_h1(paste0("Running nlmixr2 for model index: ",
                         modIdx))

      fit <- tryCatch(
      {
        fit <- suppressWarnings(nlmixr2(ui,
                                        boot_data,
                                        est = fitMeth,
                                        control = .ctl))

        .env$multipleFits <- list(
          # objf = fit$OBJF,
          # aic = fit$AIC,
          omega = fit$omega,
          parFixedDf = fit$parFixedDf[, c("Estimate", "Back-transformed")],
          message = fit$message,
          warnings = fit$warnings)

        fit # to return 'fit'
      },
      error = function(error_message) {
        message("error fitting the model")
        message(error_message)
        message("storing the models as NA ...")
        NA # return NA otherwise (instead of NULL)
      })

      saveRDS(
        .env$multipleFits,
        file = paste0(
          "./",
          output_dir,
          "/modelsEnsemble_",
          .env$mod_idx,
          ".rds"))

      saveRDS(
        fit,
        file = paste0(
          "./",
          output_dir,
          "/fitEnsemble_",
          .env$mod_idx,
          ".rds"
        )
      )

      assign("mod_idx", .env$mod_idx + 1, .env)
    })

  fitEnsemble <- NULL

  if (!restart) {
    modelsEnsemble <- c(modelsEnsembleLoaded, modelsEnsemble)
    fitEnsemble <- c(fitEnsembleLoaded, fitEnsemble)
  }

  modFileExists <-
    list.files(paste0("./", output_dir), pattern = fnameModelsEnsemblePattern)

  modelsEnsemble <- lapply(modFileExists, function(x) {
    readRDS(paste0("./", output_dir, "/", x, sep = ""))
  })

  fitFileExists <- list.files(paste0("./", output_dir), pattern = fnameFitEnsemblePattern)
  fitEnsemble <- lapply(fitFileExists, function(x) {
    readRDS(paste0("./", output_dir, "/", x, sep = ""))
  })

  list(modelsEnsemble, fitEnsemble)
}

#' Get the nlmixr2 method used for fitting the model
#'
#' @param fit the nlmixr2 fit object
#'
#' @return returns a string representing the method used by nlmixr2 for fitting the given model
#'
#' @author Vipul Mann, Matthew Fidler
#' @examples
#' getFitMethod(fit)
#' @noRd
getFitMethod <- function(fit) {
  if (!(inherits(fit, "nlmixr2FitCore"))) {
    stop("'fit' needs to be a nlmixr2 fit", call. = FALSE)
  }
  fit$est
}

#' Extract all the relevant variables from a set of bootstrapped models
#'
#' @param fitlist a list of lists containing information on the multiple bootstrapped models; similar to the output of modelsBootstrap() function
#' @param id a character representing the variable of interest: OBJF, AIC, omega, parFixedDf, method, message, warnings
#'
#' @return returns a vector or list across of the variable of interest from all the fits/bootstrapped models
#'
#' @author Vipul Mann, Matthew Fidler
#' @examples
#' extractVars(fitlist, 1) # returns a vector of OBJF values
#' extractVars(fitlist, 4) # returns a list of dataframes containing parFixedDf values
#' @noRd
extractVars <- function(fitlist, id = "method") {
  if (id == "method") {
    # no lapply for 'method'
    unlist(unname(fitlist[[1]][id]))
  }
  else {
    # if id not equal to 'method'
    res <- lapply(fitlist, function(x) {
      x[[id]]
    })


    if (!(id == "omega" ||
            id == "parFixedDf")) {
      # check if all message strings are empty
      if (id == "message") {
        prev <- TRUE
        for (i in length(res)) {
          status <- (res[[i]] == "") && prev
          prev <- status
        }
        if (status == TRUE) {
          ""
        }
        else {
          # if non-empty 'message'
          unlist(res)
        }
      }

      else {
        # if id does not equal 'message'
        unlist(res)
      }
    }
    else {
      # if id equals 'omega' or 'parFixedDf
      res
    }
  }
}

#' Summarize the bootstrapped fits/models
#'
#' @param fitList a list of lists containing information on the multiple bootstrapped models; similar to the output of modelsBootstrap() function
#' @return returns aggregated quantities (mean, median, standard deviation, and variance) as a list for all the quantities
#' @author Vipul Mann, Matthew Fidler
#' @inheritParams bootstrapFit
#' @examples
#' getBootstrapSummary(fitlist)
#' @noRd
getBootstrapSummary <- function(fitList,
                                ci = 0.95,
                                stdErrType = "perc") {
  checkmate::assertNumeric(ci, len=1, lower=0, upper=1, any.missing=FALSE, null.ok=FALSE)

  quantLevels <-
    c(0.5, (1 - ci)/2, 1 - (1 - ci)/2) # median, (1-ci)/2, 1-(1-ci)/2

  varIds <-
    names(fitList[[1]]) # number of different variables present in fitlist
  summaryList <- lapply(varIds, function(id) {
    # if (!(id %in% c("omega", "parFixedDf", "method", "message", "warnings"))) {
    #   varVec <- extractVars(fitList, id)
    #   mn <- mean(varVec)
    #   median <- median(varVec)
    #   sd <- sd(varVec)
    #
    #   c(
    #     mean = mn,
    #     median = median,
    #     stdDev = sd
    #   )
    # }
    if (id == "omega") {
      # omega estimates
      omegaMatlist <- extractVars(fitList, id)
      varVec <- simplify2array(omegaMatlist)
      mn <- apply(varVec, 1:2, mean, na.rm=TRUE)
      sd <- apply(varVec, 1:2, sd, na.rm=TRUE)

      quants <- apply(varVec, 1:2, function(x) {
        unname(quantile(x, quantLevels, na.rm=TRUE))
      })
      median <- quants[1, , ]
      confLower <- quants[2, , ]
      confUpper <- quants[3, , ]

      if (stdErrType != "perc") {
        confLower <- mn + qnorm(quantLevels[[2]]) * sd
        confUpper <- mn + qnorm(quantLevels[[3]]) * sd
      }

      # computing the covariance and correlation matrices
      # =======================================================
      parFixedOmegaBootVec <- list()
      parFixedlist <- extractVars(fitList, id = "parFixedDf")
      parFixedlistVec <- lapply(parFixedlist, function(x) {
        ret <- x$Estimate
        if (is.null(names(ret))) {
          ret <- stats::setNames(ret, rownames(x))
        }
        ret
      })
      parFixedlistVec <- do.call("rbind", parFixedlistVec)

      omgVecBoot <- list()
      omegaIdx <- seq_along(omegaMatlist)

      omgVecBoot <- lapply(omegaIdx, function(idx) {
        omgMat <- omegaMatlist[[idx]]
        omgVec <- omgMat[lower.tri(omgMat, TRUE)]
        omgVecBoot[[idx]] <- omgVec
      })
      omgVecBoot <- do.call("rbind", omgVecBoot)

      idxName <- 1
      namesList <- list()
      for (nam1 in colnames(omegaMatlist[[1]])) {
        for (nam2 in colnames(omegaMatlist[[1]])) {
          if (nam1 == nam2) {
            if (!(nam1 %in% namesList)) {
              namesList[idxName] <- nam1
              idxName <- idxName + 1
            }
          } else {
            nam <- paste0("(", nam1, ",", nam2, ")")
            namRev <- paste0("(", nam2, ",", nam1, ")")
            if (!(nam %in% namesList | namRev %in% namesList)) {
              namesList[idxName] <- nam
              idxName <- idxName + 1
            }
          }
        }
      }
      colnames(omgVecBoot) <- namesList

      .w <- which(vapply(namesList, function(x) {
        !all(omgVecBoot[, x] == 0)
      }, logical(1), USE.NAMES=FALSE))
      omgVecBoot <- omgVecBoot[, .w]


      parFixedOmegaCombined <- cbind(parFixedlistVec, omgVecBoot)

      covMatrix <- cov(parFixedOmegaCombined)
      w <- which(diag(covMatrix) == 0)
      if (length(w) > 0) {
        d <- dim(covMatrix)[1]
        corMatrix <- matrix(rep(0,d * d), d, d)
        corMatrix[-w, -w] <- cov2cor(covMatrix[-w, -w])
      } else {
        corMatrix <- cov2cor(covMatrix)
      }
      diag(corMatrix) <- sqrt(diag(covMatrix))
      dimnames(corMatrix) <- dimnames(covMatrix)
      lst <- list(
        mean = mn,
        median = median,
        stdDev = sd,
        confLower = confLower,
        confUpper = confUpper,
        covMatrixCombined = covMatrix,
        corMatrixCombined = corMatrix
      )
    }
    else if (id == "parFixedDf") {
      # parameter estimates (dataframe)
      varVec <- extractVars(fitList, id)
      mn <-
        apply(simplify2array(lapply(varVec, as.matrix)), 1:2, mean, na.rm = TRUE)
      sd <-
        apply(simplify2array(lapply(varVec, as.matrix)), 1:2, sd, na.rm = TRUE)

      quants <-
        apply(simplify2array(lapply(varVec, as.matrix)), 1:2, function(x) {
          unname(quantile(x, quantLevels, na.rm = TRUE))
        })

      median <- quants[1, , ]
      confLower <- quants[2, , ]
      confUpper <- quants[3, , ]

      if (stdErrType != "perc") {
        confLower <- mn + qnorm(quantLevels[[2]]) * sd
        confUpper <- mn + qnorm(quantLevels[[3]]) * sd
      }

      lst <- list(
        mean = mn,
        median = median,
        stdDev = sd,
        confLower = confLower,
        confUpper = confUpper
      )
    }

    else {
      # if id equals method, message, or warning
      extractVars(fitList, id)
    }
  })

  names(summaryList) <- varIds
  summaryList$nboot <- length(fitList)
  summaryList$warnings <- unique(summaryList$warnings)
  summaryList$message <- unique(summaryList$message)
  class(summaryList) <- "nlmixr2BoostrapSummary"
  summaryList
}

#' @export
print.nlmixr2BoostrapSummary <- function(x, ..., sigdig = NULL) {
  if (is.null(sigdig)) {
    if (any(names(x) == "sigdig")) {
      sigdig <- x$sigdig
    } else {
      sigdig <- 3
    }
  }

  # objf <- x$objf
  # aic <- x$aic
  message <- x$message
  warnings <- x$warnings

  omega <- x$omega
  parFixedDf <- x$parFixedDf

  nboot <- x$nboot
  cli::cli_h1(
    cli::col_red(
      "Summary of the bootstrap models (nboot: {nboot})"
    )
  )
  cli::cli_li(cli::col_magenta(
    cli::style_bold(
      "Omega matrices: mean, median, standard deviation, and confidence bousnds"
    ),
    cli::col_yellow(" (summary$omega)")
  ))

  lapply(seq_along(omega), function(x) {
    cli::cli_text(cli::col_green(paste0("$", names(omega)[x])))
    print(signif(omega[[x]], sigdig))
  })

  cli::cli_li(cli::col_magenta(
    cli::style_bold(
      "Estimated parameters: mean, median, standard deviation, and confidence bounds"
    ),
    cli::col_yellow(" (summary$parFixedDf)")
  ))
  lapply(seq_along(parFixedDf), function(x) {
    cli::cli_text(cli::col_yellow(paste0("$", names(parFixedDf)[x])))
    print(signif(parFixedDf[[x]], sigdig))
  })

  cli::cli_li(cli::cli_text(
    cli::bg_yellow(cli::style_bold("Messages")),
    cli::col_yellow(" (summary$message)")
  ))
  print(message)

  cli::cli_li(cli::cli_text(
    cli::bg_red(cli::style_bold(cli::col_white("Warnings"))),
    cli::col_yellow(" (summary$warnings)")
  ))
  print(warnings)

  cli::cli_h1("end")
  invisible(x)
}

#' Assign a set of variables to the nlmixr2 fit environment
#'
#' @param namedVars a named list of variables that need to be assigned to the given environment
#' @param fitobject the nlmixr2 fit object that contains its environment information
#' @noRd
#'
assignToEnv <- function(namedVars, fitobject) {
  if (!inherits(fitobject, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit", call. = FALSE)
  }

  if (is.null(names(namedVars))) {
    stop("'namedVars needs to be a named list", call. = FALSE)
  }

  if (length(namedVars) != length(names(namedVars))) {
    stop("'namedVars does not have all the elements named", call. = FALSE)
  }
  env <- fitobject$env
  lapply(names(namedVars), function(x) {
    assign(x, namedVars[[x]], envir = env)
  })
}

#' @title Produce delta objective function for boostrap
#'
#' @param x fit object
#' @param ... other parameters
#' @return Fit traceplot or nothing.
#' @author Vipul Mann,  Matthew L. Fidler
#' @references
#' R Niebecker,  MO Karlsson. (2013)
#' *Are datasets for NLME models large enough for a bootstrap to provide reliable parameter uncertainty distributions?*
#' PAGE 2013.
#' <https://www.page-meeting.org/Abstracts/are-datasets-for-nlme-models-large-enough-for-a-bootstrap-to-provide-reliable-parameter-uncertainty-distributions/>
#' @export
bootplot <- function(x, ...) {
  UseMethod("bootplot")
}

#' @rdname bootplot
#' @export
#' @importFrom ggplot2 .data
bootplot.nlmixr2FitCore <- function(x, ...) {
  .fitName <- as.character(substitute(x))
  if (inherits(x, "nlmixr2FitCore")) {
    if (exists("bootSummary", x$env) & (!exists(".bootPlotData", x$env))) {
      bootstrapFit(x, x$bootSummary$nboot, plotHist = TRUE, fitName = .fitName)
    }
    if (exists(".bootPlotData", x$env)) {
      if (x$bootSummary$nboot != x$env$.bootPlotData$deltaN) {
        bootstrapFit(x, x$bootSummary$nboot, plotHist = TRUE, fitName = .fitName)
      }
      .chisq <- x$env$.bootPlotData$chisq
      .dfD <- x$env$.bootPlotData$dfD
      .deltaN <- x$env$.bootPlotData$deltaN
      .df2 <- x$env$.bootPlotData$df2
      .plot <- ggplot2::ggplot(.chisq, ggplot2::aes(.data$quantiles, .data$deltaofv, color = .data$Distribution)) +
        ggplot2::geom_line() +
        ggplot2::ylab("\u0394 objective function") +
        ggplot2::geom_text(data = .dfD, ggplot2::aes(label = .data$label), hjust = 0) +
        ggplot2::xlab("Distribution quantiles") +
        ggplot2::scale_color_manual(values = c("red", "blue")) +
        rxode2::rxTheme() +
        ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")

      if (requireNamespace("ggtext", quietly = TRUE)) {
        .plot <- .plot +
          ggplot2::theme(
            plot.title = ggtext::element_markdown(),
            legend.position = "none"
          ) +
          ggplot2::labs(
            title = paste0(
              'Bootstrap <span style="color:blue; opacity: 0.2;">\u0394 objective function (', .deltaN,
              " models, df\u2248", .df2, ')</span> vs <span style="color:red; opacity: 0.2;">reference \u03C7\u00B2(df=',
              length(x$ini$est), ")</style>"
            ),
            caption = "\u0394 objective function curve should be on or below the reference distribution curve"
          )
      } else {
        .plot <- ggplot2::labs(
          title = paste0("Distribution of \u0394 objective function values for ", .deltaN, " df=", .df2, " models"),
          caption = "\u0394 objective function curve should be on or below the reference distribution curve"
        )
      }
      .plot
    } else {
      stop("this nlmixr2 object does not include boostrap distribution statics for comparison",
           call. = FALSE
           )
    }
  } else {
    stop("this is not a nlmixr2 object",
         call. = FALSE
         )
  }
}
