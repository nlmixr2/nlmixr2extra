#' Add covariates and lasso string to ui
#'
#' @param ui compiled rxode2 nlmir2 model or fit
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec  character vector of covariates that need to be added
#' @param tvalue float indicating tvalue to be updated for the lasso string
#' @return updated ui with added covariates
#' @noRd
#' @author Vishal Sarsani
.lassoUicovariate <- function(ui,varsVec,covarsVec,tvalue=0.05){

  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$finalUiEnv
  }
  ui <- rxode2::assertRxUi(ui)

  # Extract updated ui with added covariates
  .ui1 <- rxode2::rxUiDecompress(buildupatedUI(ui,varsVec,covarsVec,add = TRUE,indep=FALSE))

  # Extract all the covariate parameters constructed
  .covsparams <- .ui1$muRefCovariateDataFrame$covariateParameter
  checkmate::assertCharacter(.covsparams,min.len = 1,any.missing = FALSE )

  # construct a tvalue string
  tvalueString <- paste0("tvalue <- ",tvalue,"\n")

  # construct an abssum string
  absString <- paste0("abssum <- sum(",paste0(paste0("abs(",.covsparams,")"),collapse="+"),")","\n")

  # construct a ratio string
  ratioString <- paste0("ratio <- ","abssum","/","tvalue","\n")

  # construct a factor string
  factorString <- paste0("factor <- ","exp","(","1-","ratio",")","\n")

  # add lasso factor strings from the covparams to the model function
  .covsfactorString  <- paste0(.covsparams, " * factor")
  .modelfun <- .ui1$funTxt
  for ( i in seq_along(.covsparams)){
    .modelfun  <- gsub(.covsparams[i], .covsfactorString[i], .modelfun)
  }

  # construct the final updated model string
  .newmodelfun <- paste0(tvalueString,absString,ratioString,factorString,.modelfun,sep="")
  .newModel <- eval(parse(text = paste0("quote(model({",paste0(as.character(.newmodelfun),collapse="\n"), "}))")))

  # build ui
  .ini <- .ui1$iniDf
  .ini <- as.expression(lotri::as.lotri(.ini))
  .ini[[1]] <- quote(`ini`)
  # return the updated model
  .getUiFunFromIniAndModel(.ui1, .ini, .newModel)()
}

#' Perform cross validation and return predictive-objective-function value
#'
#' @param data dataset containing all the required information.
#' @param ui compiled rxode2 nlmir2 model or fit
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec  character vector of covariates that need to be added
#' @param tvalue a desired tvalue for constructing lasso constraint
#' @param nfold number of folds for cross-validation
#' @param optcrit criteria for optimization. must be either 'objf' or 'llk'
#' @param estmethod Estimation method for initial nlmixr fit.  must be either 'focei' or 'saem'
#' @return predictive-objective-function value
#' @noRd
#' @author Vishal Sarsani
.crossvalidationLasso <- function(data,ui,varsVec,covarsVec,tvalue=0.10,nfold=5,optcrit='objf',estmethod="focei",adapcoefs=NULL,stratVar=NULL){

  # check if dataframe
  checkmate::assert_data_frame(data,min.cols = 7)
  # check if t-value and nfold are valid
  checkmate::assert_double(tvalue)
  checkmate::assert_int(nfold)
  # List of objective-function values
  ofvList <- list()
  # Build updated ui with given covariate info and add lasso strings
  if(!is.null(adapcoefs)){
    mod <- .adaptivelassoUicovariate(ui,varsVec,covarsVec,tvalue = tvalue,adapcoefs = adapcoefs)
  } else {
    mod <- .lassoUicovariate(ui,varsVec,covarsVec,tvalue = tvalue)
  }

  # Add fold column depending on the stratification
  if (!is.null(stratVar)) {
    dataF <-    foldgen(data,nfold=5,stratVar = stratVar)
  } else {
    dataF <-    foldgen(data,nfold=5)
  }

  #Perform cross-validation
  for (f in seq_len(nfold)) {
    # Training dataset
    trainData <- dataF[dataF$fold!=f,]
    # Testing dataset
    testData <- dataF[dataF$fold==f,]

    ## Stop if test data is empty
    if (dim(testData)[1] == 0) {
      cli::cli_alert_danger("the test dataset is empty. probably not enough individuals to do desired n-fold cross-validation")
      stop("aborting...please re-run by reducing nfold number for the cross-validation", call. = FALSE)
    }

    cli::cli_alert_success("Training and Testing data sets successfully created for cross-validation for fold number {f}")

    # Training Estimation
    fitTrain <- tryCatch(
      {
        fitTrain <-
          suppressWarnings(nlmixr2(mod,trainData ,estmethod))
        fitTrain # to return 'fitTrain'
      },
      error = function(error_message) {
        print("error fitting the model for the training dataset ")
        print(error_message)
      })

    # Testing data model fit with estimates from the training.

    fitTest <- tryCatch(
      {
        fitTest <-
          suppressWarnings(nlmixr2(mod,testData ,"posthoc"))
        fitTest # to return 'fitTest'
      },
      error = function(error_message) {
        print("error fitting the model for the testing dataset ")
        print(error_message)
      })

    # Extract Predictive Objective function value
    if (optcrit=="objf"){
      if (!is.numeric(fitTest$objDf$OBJF)){
        cli::cli_alert_danger("the 'fit' object needs to have an objective functions value associated with it")
        cli::cli_alert_info("try computing 'fit$objDf$OBJF' in console to compute and store the corresponding OBJF value")
        stop("aborting...objf value not associated with the current 'fit' object", call. = FALSE)
      }
      ofvList <- append(ofvList, fitTest$objDf$OBJF[1])
    } else if(optcrit=="llk") {
      if (!is.numeric(fitTest$objDf$`Log-likelihood`)){
        cli::cli_alert_danger("the 'fit' object needs to have an objective functions value associated with it")
        cli::cli_alert_info("try computing 'fit$objDf$$`Log-likelihood`' in console to compute and store the corresponding OBJF value")
        stop("aborting...objf value not associated with the current 'fit' object", call. = FALSE)
      }
      ofvList <- append(ofvList ,-2*fitTest$objDf$`Log-likelihood`[1])
    }
    cli::cli_alert_success("Estimation complete for the fold number : {f}")
  }

  pOFV <- do.call(sum, ofvList)
  pOFV
}

#' Return optimal t-value from the cross validation
#'
#' @param data dataset containing all the required information.
#' @param ui compiled rxode2 nlmir2 model or fit
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec  character vector of covariates that need to be added
#' @param t_start a desired t start value for constructing lasso constraint
#' @param t_stop a desired t stop value for constructing lasso constraint
#' @param t_step a desired t step value for constructing lasso constraint
#' @param convergence either REACHMAX or FIRSTMIN
#' @return Optimal t-value among tvalue range
#' @noRd
#' @author Vishal Sarsani
.optimalTvaluelasso <- function(data,ui,varsVec,covarsVec,t_start=0.01,t_stop=0.99,t_step=0.01,stratVar = NULL,convergence="REACHMAX",...){
  # check if t-start,stop and step values are valid
  checkmate::assert_double(t_start)
  checkmate::assert_double(t_stop)
  checkmate::assert_double(t_step)

  # generate vector of t values
  tvalues <- seq(t_start,t_stop,by=t_step)

  # create a data frame of pOFV values for t values
  firstmin <- Inf
  pofvList <- data.frame()
  for (t in tvalues){

    ofValue <- suppressWarnings(.crossvalidationLasso(data,ui,varsVec,covarsVec,tvalue=t,nfold=5,optcrit='obj',estmethod="focei",adapcoefs=NULL,stratVar = NULL))

    pofv <- data.frame(tvalue=t,POFV=ofValue)
    cli::cli_alert_success("Cross-validation finished for the t-value : {t}")
    pofvList <- rbind(pofvList,pofv)
  }

  # Find optimal t-value
  optimal_t <- pofvList[which.min(pofvList$POFV),]$tvalue
  optimal_t
}

#' Return Final lasso coefficients after finding optimal t
#'
#' @param fit nlmixr2 fit.
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec  character vector of covariates that need to be added
#' @param catvarsVec  character vector of categorical covariates that need to be added
#' @param constraint theta cutoff. below cutoff then the theta will be fixed to zero
#' @param stratVar   A variable to stratify on for cross-validation
#' @param ...  Other parameters to be passed to optimalTvaluelasso
#' @return return data frame of final lasso coefficients
#' @author Vishal Sarsani
#' @export
#'
#' @examples
#' \donttest{
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
#' d <- nlmixr2data::theo_sd
#' d$SEX <-0
#' d$SEX[d$ID<=6] <-1
#'
#' fit <- nlmixr2(one.cmt, d, "focei")
#' varsVec <- c("ka","cl","v")
#' covarsVec <- c("WT")
#' catvarsVec <- c("SEX")
#'
#'
#' # Lasso coefficients:
#'
#' lassoDf <- lassoCoefficients(fit,varsVec,covarsVec,catvarsVec,constraint=1e-08,stratVar = NULL)
#'
#' }
lassoCoefficients <- function(fit,varsVec,covarsVec,catvarsVec,constraint=1e-08,stratVar = NULL,...) {
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  } else {
    ui <- fit$finalUiEnv
  }

  checkmate::assert_character(covarsVec)
  checkmate::assert_character(varsVec)
  checkmate::assert_double(constraint)

  ## Update data and covarsvec if categorical variables are provided
  if (!is.null(catvarsVec)) {
    covarsVec <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[2]]
    data <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[1]]
  } else {
    covarsVec <- covarsVec
    data <- nlme::getData(fit)
  }

  data <- normalizedData(data, covarsVec)

  # construct covInfo
  covInfo <-  buildcovInfo(varsVec, covarsVec)
  # construct covNames
  covNames <- list()
  for ( x in covInfo) {
    covName <- paste0("cov_", x$covariate, "_", x$varName)
    covNames <- append(covNames,covName)
  }

  #Estimated method
  .est <- getFitMethod(fit)
  # Extract optimal t-value
  optTvalue <- .optimalTvaluelasso(data,ui,varsVec,covarsVec,t_start=0.05,t_stop=0.25,t_step=0.05,estmethod=.est,stratVar = NULL,...)

  # Refit model with the optimal t-value
  updatedmod <- .lassoUicovariate(ui,varsVec,covarsVec,tvalue=optTvalue)
  fitobject <- nlmixr2(updatedmod,data,est=.est)

  # Extract covariate estimates
  covEst <- fitobject$parFixedDf[row.names(fitobject$parFixedDf) %in% covNames,"Estimate"]

  # Extract covariate std error
  covStd <- fitobject$parFixedDf[row.names(fitobject$parFixedDf) %in% covNames,"SE"]

  #absolute sum of  lasso THETAs
  abssum <- sum(abs(covEst))
  # ratio
  ratio <- abssum/optTvalue
  # factor
  factor <- exp(1-ratio)

  # Apply lasso constraint
  finalLasso <- as.data.frame(lapply(covEst, function(x) ifelse(abs(x) < constraint , 0, x)),row.names = NULL)

  #Multiply by factor
  finalLasso <- finalLasso *factor
  finalLasso
}

#' Add covariates and adaptive lasso string to ui
#'
#' @param ui compiled rxode2 nlmir2 model or fit
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec  character vector of covariates that need to be added
#' @param tvalue float indicating tvalue to be updated for the lasso string
#' @param adapcoefs dataframe with covariates estimates for the adaptive lasso
#' @return updated ui with added covariates
#' @noRd
#' @author Vishal Sarsani
.adaptivelassoUicovariate <- function(ui,varsVec,covarsVec,tvalue=0.05,adapcoefs=NULL){
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$finalUiEnv
  }
  ui <- rxode2::assertRxUi(ui)

  # Extract updated ui with added covariates
  .ui1 <- buildupatedUI(ui,varsVec,covarsVec)

  # Extract all the covariate parameters constructed
  .covsparams <- .ui1$muRefCovariateDataFrame$covariateParameter
  checkmate::assertCharacter(.covsparams,min.len = 1,any.missing = FALSE )

  #check if adaptive coefficients have values for all covaraite-parameter
  checkmate::assert_data_frame(adapcoefs,ncols = length(varsVec)*length(covarsVec))

  # construct a tvalue string
  tvalueString <- paste0("tvalue <- ",tvalue,"\n")

  # construct an abssum string
  absString <- paste0("abssum <- sum(",paste0(paste0("abs(",.covsparams,")"),collapse="+"),")","\n")

  # construct a ratio string
  ratioString <- paste0("ratio <- ","abssum","/","tvalue","\n")

  # construct a factor string
  factorString <- paste0("factor <- ","exp","(","1-","ratio",")","\n")

  # construct a Adaptive lasso coefficeints
  adaptString <- paste0("AL_",colnames(adapcoefs)," <- ",adapcoefs,"\n",collapse='')

  # add lasso factor strings from the covparams to the model function
  .covsadapString  <-  paste0(.covsparams, " * factor * ","AL_",.covsparams)
  .modelfun <- .ui1$funTxt
  for ( i in seq_along(.covsparams)){
    .modelfun  <- gsub(.covsparams[i], .covsadapString[i], .modelfun)
  }

  # construct the final updated model string
  .newmodelfun <- paste0(tvalueString,absString,ratioString,factorString,adaptString,.modelfun,sep="")
  .newModel <- eval(parse(text = paste0("quote(model({",paste0(as.character(.newmodelfun),collapse="\n"), "}))")))

  # build ui
  .ini <- .ui1$iniDf
  .ini <- as.expression(lotri::as.lotri(.ini))
  .ini[[1]] <- quote(`ini`)
  # return the updated model
  .getUiFunFromIniAndModel(.ui1, .ini, .newModel)()
}

#' Return Adaptive lasso coefficients after finding optimal t
#'
#' @param fit nlmixr2 fit.
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec  character vector of covariates that need to be added
#' @param catvarsVec  character vector of categorical covariates that need to be added
#' @param constraint theta cutoff. below cutoff then the theta will be fixed to zero.
#' @param stratVar   A variable to stratify on for cross-validation.
#' @param ...  Other parameters to be passed to optimalTvaluelasso
#' @return return data frame of final lasso coefficients
#' @author Vishal Sarsani
#' @export
#'
#' @examples
#' \donttest{
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
#' d <- nlmixr2data::theo_sd
#' d$SEX <-0
#' d$SEX[d$ID<=6] <-1
#'
#' fit <-
#'  nlmixr2(
#'    one.cmt, d,
#'    est = "focei", control = nlmixr2est::foceiControl(print = 0)
#'  )
#' varsVec <- c("ka","cl","v")
#' covarsVec <- c("WT")
#' catvarsVec <- c("SEX")
#'
#' # Adaptive Lasso coefficients:
#'
#' lassoDf <- adaptivelassoCoefficients(fit, varsVec, covarsVec, catvarsVec)
#' }
adaptivelassoCoefficients <- function(fit,varsVec,covarsVec,catvarsVec,constraint=1e-08,stratVar = NULL,...) {

  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  else {
    ui <- fit$finalUiEnv
  }

  checkmate::assert_character(covarsVec)
  checkmate::assert_character(varsVec)
  checkmate::assert_double(constraint)

  ## Get initial adaptive coefs from the regular lasso coefficients
  .adapcoefs <- lassoCoefficients(fit,varsVec,covarsVec,catvarsVec,constraint=1e-08,stratVar = NULL,...)

  ## Update data and covarsvec if categorical variables are provided
  if(!is.null(catvarsVec)) {
    covarsVec <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[2]]
    data <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[1]]
  }

  #data
  data <- normalizedData(data,covarsVec)

  # construct covInfo
  covInfo <-  buildcovInfo(varsVec,covarsVec)
  # construct covNames
  covNames <- list()
  for ( x in covInfo) {
    covName <- paste0("cov_", x$covariate, "_", x$varName)
    covNames <- append(covNames,covName)
  }

  #Estimated method
  .est <- getFitMethod(fit)
  # Extract optimal t-value
  optTvalue <- .optimalTvaluelasso(data,ui,varsVec,covarsVec,t_start=0.05,t_stop=0.25,t_step=0.05,estmethod=.est,adapcoefs=.adapcoefs,stratVar=NULL,...)

  # Refit model with the optimal t-value

  updatedmod <- .adaptivelassoUicovariate(ui,varsVec,covarsVec,tvalue=optTvalue,adapcoefs = .adapcoefs)
  fitobject <- nlmixr2(updatedmod,data,est=.est)

  # Extract covariate estimates
  covEst <- fitobject$parFixedDf[row.names(fitobject$parFixedDf) %in% covNames,"Estimate"]

  # Multiply covariate Estimates with adative coefficients
  covestAdap <- covEst*.adapcoefs

  # Extract covariate std error
  covStd <- fitobject$parFixedDf[row.names(fitobject$parFixedDf) %in% covNames,"SE"]

  #absolute sum of  lasso THETAs
  abssum <- sum(abs(covEst))
  # ratio
  ratio <- abssum/optTvalue
  # factor
  factor <- exp(1-ratio)

  # Apply lasso constraint
  finalLasso <- as.data.frame(lapply(covestAdap, function(x) ifelse(abs(x) <= constraint , 0, x)),row.names = NULL)

  finalLasso <- finalLasso *factor
  finalLasso
}

#' Regular lasso model
#'
#' @param fit nlmixr2 fit.
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec  character vector of covariates that need to be added
#' @param catvarsVec  character vector of categorical covariates that need to be added
#' @param constraint theta cutoff. below cutoff then the theta will be fixed to zero.
#' @param lassotype must be  'regular' , 'adaptive', 'adjusted'
#' @param stratVar   A variable to stratify on for cross-validation.
#' @param ...  Other parameters to be passed to optimalTvaluelasso
#' @return return fit of the selected lasso coefficients
#' @author Vishal Sarsani
#' @export
#'
#' @examples
#' \donttest{
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
#' d <- nlmixr2data::theo_sd
#' d$SEX <-0
#' d$SEX[d$ID<=6] <-1
#'
#' fit <- nlmixr2(one.cmt, d, "focei")
#' varsVec <- c("ka","cl","v")
#' covarsVec <- c("WT")
#' catvarsVec <- c("SEX")
#'
#'
#' # Model fit with regular lasso coefficients:
#'
#' lassoDf <- regularmodel(fit,varsVec,covarsVec,catvarsVec)
#' # Model fit with adaptive lasso coefficients:
#'
#' lassoDf <- regularmodel(fit,varsVec,covarsVec,catvarsVec,lassotype='adaptive')
#' # Model fit with adaptive-adjusted lasso coefficients:
#'
#' lassoDf <- regularmodel(fit,varsVec,covarsVec,catvarsVec, lassotype='adjusted')
#' }
regularmodel <- function(fit,varsVec,covarsVec,catvarsVec,constraint=1e-08,
                         lassotype=c("regular", "adaptive", "adjusted"),
                         stratVar = NULL,...) {
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  else {
    ui <- fit$finalUiEnv
  }

  checkmate::assert_character(covarsVec)
  checkmate::assert_character(varsVec)
  checkmate::assert_double(constraint)
  lassotype <- match.arg(lassotype)

  if (lassotype=="regular") {
    .coefValues <- lassoCoefficients(fit,varsVec,covarsVec,catvarsVec,constraint=1e-08,stratVar = NULL,...)
  } else if (lassotype=="adaptive") {
    .coefValues <- adaptivelassoCoefficients(fit,varsVec,covarsVec,catvarsVec,constraint=1e-08,stratVar = NULL,...)
  } else if (lassotype=="adjusted"){
    .coefValues <- adjustedlassoCoefficients(fit,varsVec,covarsVec,catvarsVec,constraint=1e-08,stratVar = NULL,...)
  }

  if(!is.null(catvarsVec)){
    covarsVec <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[2]]
    data <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[1]]
  }

  # Extract updated ui with added covariates
  .ui1 <- buildupatedUI(ui,varsVec,covarsVec,add=TRUE,indep = FALSE)

  # Extract all the covariate parameters constructed
  .covsparams <- .ui1$muRefCovariateDataFrame$covariateParameter
  checkmate::assertCharacter(.covsparams,min.len = 1,any.missing = FALSE)

  # Fix the covariate parameter values
  .modelfun <- .ui1$funTxt
  for ( i in seq_along(.covsparams)){
    .modelfun  <- gsub(.covsparams[i], .coefValues[[i]], .modelfun)
  }

  tvalue <- 0.01

  # construct a tvalue string
  tvalueString <- paste0("tvalue <- ",tvalue,"\n")

  # construct an abssum string
  absString <- paste0("abssum <- sum(",paste0(paste0("abs(",.covsparams,")"),collapse="+"),")","\n")
  for ( i in seq_along(.covsparams)) {
    absString  <- gsub(.covsparams[i], .coefValues[[i]], absString)
  }

  # construct a ratio string
  ratioString <- paste0("ratio <- ","abssum","/","tvalue","\n")

  # construct a factor string
  factorString <- paste0("factor <- ","exp","(","1-","ratio",")","\n")

  # construct the final updated model string
  .newmodelfun <- paste0(tvalueString,absString,ratioString,factorString,.modelfun,sep="")
  .newModel <- eval(parse(text = paste0("quote(model({",paste0(as.character(.newmodelfun),collapse="\n"), "}))")))

  # build ui
  .ini <- .ui1$iniDf
  .ini <- .ini[ ! .ini$name %in% .covsparams, ]
  .ini <- as.expression(lotri::as.lotri(.ini))
  .ini[[1]] <- quote(`ini`)
  # return the updated model
  .ui2 <- rxode2::rxUiDecompress(.getUiFunFromIniAndModel(.ui1, .ini, .newModel)())
  #estimation method
  estmethod=getFitMethod(fit)
  nlmixr2(.ui2,data,est=estmethod)
}

#' Return Adjusted adaptive lasso coefficients after finding optimal t
#'
#' @param fit nlmixr2 fit.
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec  character vector of covariates that need to be added
#' @param catvarsVec  character vector of categorical covariates that need to be added
#' @param constraint theta cutoff. below cutoff then the theta will be fixed to zero.
#' @param stratVar   A variable to stratify on for cross-validation.
#' @param ...  Other parameters to be passed to optimalTvaluelasso
#' @return return data frame of final lasso coefficients
#' @author Vishal Sarsani
#' @export
#'
#' @examples
#' \donttest{
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
#' d <- nlmixr2data::theo_sd
#' d$SEX <-0
#' d$SEX[d$ID<=6] <-1
#'
#' fit <- nlmixr2(one.cmt, d, "focei")
#' varsVec <- c("ka","cl","v")
#' covarsVec <- c("WT")
#' catvarsVec <- c("SEX")
#'
#' # Adaptive Lasso coefficients:
#'
#' lassoDf <- adjustedlassoCoefficients(fit,varsVec,covarsVec,catvarsVec)
#' }
adjustedlassoCoefficients <- function(fit,varsVec,covarsVec,catvarsVec,constraint=1e-08,stratVar = NULL,...) {
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  } else {
    ui <- fit$finalUiEnv
  }

  checkmate::assert_character(covarsVec)
  checkmate::assert_character(catvarsVec)
  checkmate::assert_character(varsVec)
  checkmate::assert_double(constraint)

  if(!is.null(catvarsVec)) {
    covarsVec <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[2]]
    data <- addCatCovariates(nlme::getData(fit),covarsVec = covarsVec,catcovarsVec = catvarsVec)[[1]]
  }

  ## Get updated ui with covariates and without lasso factor
  .mod1 <- buildupatedUI(ui,varsVec,covarsVec,add = TRUE,indep = FALSE)

  # construct covInfo
  covInfo <-  buildcovInfo(varsVec,covarsVec)
  # construct covNames
  covNames <- list()
  for ( x in covInfo) {
    covName <- paste0("cov_", x$covariate, "_", x$varName)
    covNames <- append(covNames,covName)
  }

  data <- normalizedData(data,covarsVec)

  fit1 <- nlmixr2(.mod1,data,est=getFitMethod(fit))

  # Construct AL coefficents
  .adapcoefs <- data.frame(as.list(fit1$parFixedDf[row.names(fit1$parFixedDf) %in% covNames,"Estimate"]/fit1$parFixedDf[row.names(fit1$parFixedDf) %in% covNames,"SE"]))

  # Extract optimal t-value
  optTvalue <- .optimalTvaluelasso(data,ui,varsVec,covarsVec,t_start=0.05,t_stop=0.25,t_step=0.05,estmethod=getFitMethod(fit),adapcoefs=.adapcoefs,stratVar = NULL,...)

  # Refit model with the optimal t-value
  updatedmod <- .adaptivelassoUicovariate(ui,varsVec,covarsVec,tvalue=optTvalue,adapcoefs = .adapcoefs)
  fitobject <- nlmixr2(updatedmod,data,est=getFitMethod(fit))

  # Extract covariate estimates
  covEst <- fitobject$parFixedDf[row.names(fitobject$parFixedDf) %in% covNames,"Estimate"]

  # Multiply covariate Estimates with adative coefficients
  covestAdap <- covEst*.adapcoefs

  # Extract covariate std error
  covStd <- fitobject$parFixedDf[row.names(fitobject$parFixedDf) %in% covNames,"SE"]

  #absolute sum of  lasso THETAs
  abssum <- sum(abs(covEst))
  ratio <- abssum/optTvalue
  factor <- exp(1-ratio)

  # Apply lasso constraint
  finalLasso <- as.data.frame(lapply(covestAdap, function(x) ifelse(abs(x) <= constraint , 0, x)),row.names = NULL)

  finalLasso <- finalLasso *factor
  finalLasso
}
