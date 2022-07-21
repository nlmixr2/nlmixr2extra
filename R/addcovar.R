#' lasso method for predictive covariate model building
#' lasso is suggested over SCM to remedy selection bias ,poor 
#' predictive performance and speed 
#' 


#' Add covariate 
#' 
#' 
#' @param ui compiled rxode2 nlmir2 model or fit  
#' @param varName  the variable name to which the given covariate is to be added
#' @param covariate the covariate that needs string to be constructed
#' @param norm  boolean indicating if the covariate that needs to be added is normalized; default is FALSE

.addCovariate <- function(ui,varName,covariate,norm=FALSE) {
  
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$ui
  }
  ui <- rxode2::assertRxUi(ui)
  
  
  checkmate::assertCharacter(varName,len = 1,any.missing = FALSE )
  checkmate::assertCharacter(covariate,len = 1,any.missing = FALSE )
  
  if (inherits(try(str2lang(varName)),"try-error")){
    stop("`varName` must be a valid R expression",call. = FALSE)
  }
  if (inherits(try(str2lang(covariate)),"try-error")){
    stop("`varName` must be a valid R expression",call. = FALSE)
  }
  .pop <- .getThetaName(ui,varName=varName)
  
  if(norm){
    .cov <- gsub('_normalized', '', paste0("cov_", covariate, "_", varName))
  }
  else {
    .cov <- paste0("cov_", covariate, "_", varName)
  }
  
  
  .covdf <- rbind(ui$muRefCovariateDataFrame,
                    data.frame(theta=.pop,covariate=covariate,covariateParameter=.cov))  
  .split <- ui$getSplitMuModel
  .pars <- c(names(.split$pureMuRef),names(.split$taintMuRef))
  .model <- nlmixr2est::.saemDropMuRefFromModel(ui)
  lapply(.model,.expandRefMu,murefDf=ui$muRefDataFrame,covDf=.covdf,pars=.pars)
}




#' Given model expression, Expand population expression
#'
#' @param x model expression
#' @param murefDf MU referencing data frame
#' @param covDf covariate referencing data frame
#' @param pars theta parameters
#'
#' @return
#' @noRd
.expandRefMu <- function(x,murefDf,covDf,pars) {
  if (is.name(x)) {
    currparam <- as.character(x)
    if (currparam %in% pars){
      return (str2lang(.expandPopExpr(currparam,murefDf,covDf)))
    } 
  }else if(is.call(x)) {
    return(as.call(c(list(x[[1]]), lapply(x[-1], .expandRefMu, murefDf=murefDf, covDf=covDf, pars=pars))))
  }
  x
}
  


#' Expand population expression given the mu reference and covariate reference Data frames
#'
#' @param popParam  population parameter variable
#' @param murefDf MU referencing data frame
#' @param covDf covariate referencing data frame
#' @param factor a  lasso constraint factor already calculated
#'
#' @return expanded expression
#' @noRd
.expandPopExpr <- function(popParam,murefDf,covDf,factor=NULL) {
  
  .par1 <- murefDf[murefDf$theta==popParam,]
  .res <- paste0(.par1$theta,"+",.par1$eta)
  .w <- which(covDf$theta==popParam)
  if(length(.w)>0){
    if (!is.null(factor)){
      ##checkmate::assert_numeric(factor,len = 1)
      .res <- paste0(c(.res,paste0(covDf[.w,"covariateParameter"],"*",covDf[.w,"covariate"],"*",factor)),
                                  collapse="+")}
    else {
      
      .res <- paste0(c(.res,paste0(covDf[.w,"covariateParameter"],"*",covDf[.w,"covariate"])),
                     collapse="+") }
    
  }
  .res
}




#' get the population parameter from variable name
#'
#' @param ui compiled rxode2 nlmir2 model or fi
#' @param varName the variable name to which the given covariate is to be added
#'
#' @return population parameter variable
#' @noRd
.getThetaName <- function(ui,varName){
  .split <-  ui$getSplitMuModel
  if (varName %in% names(.split$pureMuRef)){
    return(varName)
  }
  .w <- which(.split$pureMuRef==varName)
  if(length(.w)==1){
    return(names(.split$pureMuRef)[.w])  
  }
  if (varName %in% names(.split$taintMuRef)){
    return(varName)
  }
  stop("'",varName,"'","has not been found in the model ui",call. = FALSE)
}







#' Given a data frame extract column corresponding to  Individual
#'
#' @param data given data frame 
#' 
#' @return column name of individual
#' @noRd

.idColumn <- function(data) {
  #check if it is a dataframe
  checkmate::assertDataFrame(data,col.names = "named")
  # Extract individual ID from column names
  colNames <- colnames(data)
  colNamesLower <- tolower(colNames)
  if ("id" %in% colNamesLower) {
    uidCol <- colNames[which("id" %in% colNamesLower)]
  } else {
    uidCol <- "ID"
  }
 uidCol
}




# Function to calculate population standard deviation from group means
#' @param x vector of values to calculate standard deviation.
#' @return population standard deviation
#' @noRd
.sd.p=function(x){sd(x)*sqrt((length(x)-1)/length(x))}


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
  
  if (inherits(try(str2lang(covariate)),"try-error")){
    stop("`varName` must be a valid R expression",call. = FALSE)
  }
  
  .new <- intersect(names(data), covariate)
  if (length(.new) == 0L) stop("covariate specified not in original dataset")
  
  #extract Individual ID from data frame
  
  uidCol <- .idColumn(data)
  # mean by groups (Individual)
  groupMeans <- with(data, ave(get(covariate),get(uidCol), FUN = function(x) mean(x, na.rm = TRUE)))
  # pop mean
  popMean <- mean(groupMeans)
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
  
  if (inherits(try(str2lang(covariate)),"try-error")){
    stop("`varName` must be a valid R expression",call. = FALSE)
  }
  .new <- intersect(names(data), covariate)
  if (length(.new) == 0L) stop("covariate specified not in original dataset")

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


#' Function to return data of normalized covariates
#' 
#' 
#' @param fitobject an nlmixr2 'fit' object
#' @param covarsVec a list of covariate names (parameters) that need to be estimates
#' @param replace replace the original covariate data with normalized data for easier updated model.
#' 
 #' @return data frame with all normalized covariates
#' @noRd
.normalizedData <- function(fit,covarsVec,replace=TRUE) {
  
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  
  
  checkmate::assert_character(covarsVec)

  # get data 
  data <- nlme::getData(fit)
  
  ## 
  .normalizedDFs <- lapply(covarsVec,.normalizeDf,data=data)
  
  # final data frame of normalized covariates
  
  if(replace){
    .dat <- Reduce(merge,.normalizedDFs)
    dropnormPrefix <- function(x){ colnames(x) <- gsub("normalized_", "", colnames(x)); x }
    .dat <- .dat[ , !names(.dat) %in% covarsVec]
    .finalDf <- dropnormPrefix(.dat)
  }else{
  
  .finalDf <- Reduce(merge,.normalizedDFs)
  }
  .finalDf
}



#' Build ui from the covariate 
#' 
#' 
#' @param ui compiled rxode2 nlmir2 model or fit  
#' @param varName  the variable name to which the given covariate is to be added
#' @param covariate the covariate that needs string to be constructed
#' @return ui with added covariate
#' @noRd
.builduiCovariate <- function(ui,varName,covariate){
  
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$ui
  }
  ui <- rxode2::assertRxUi(ui)
  
  
  checkmate::assertCharacter(varName,len = 1,any.missing = FALSE )
  checkmate::assertCharacter(covariate,len = 1,any.missing = FALSE )
  
  if (inherits(try(str2lang(varName)),"try-error")){
    stop("`varName` must be a valid R expression",call. = FALSE)
  }
  if (inherits(try(str2lang(covariate)),"try-error")){
    stop("`varName` must be a valid R expression",call. = FALSE)
  }
  
# Add covariate to model expression 
lst <- .addCovariate(ui,varName,covariate,norm=TRUE)
.newModel <- eval(parse(text = paste0("quote(model({",paste0(as.character(lst),collapse="\n"), "}))")))

# Add covariate to initialization 

nthetaLength <- length(which(!is.na(ui$iniDf$ntheta)))
.ini <- ui$iniDf
covName <- paste0("cov_", covariate, "_", varName)


  .ini <- rbind(.ini,data.frame(ntheta=as.integer(nthetaLength+1),neta1=NA_character_,neta2=NA_character_,
                                name=covName,lower=-Inf,est=0,upper=Inf,fix=FALSE,label=NA_character_,
                                backTransform=NA_character_,condition=NA_character_,err=NA_character_)) 
  

# build ui
  .ini <- as.expression(lotri::as.lotri(.ini))
  .ini[[1]] <- quote(`ini`)
  return(.getUiFunFromIniAndModel(ui, .ini, .newModel)())  

}



' Build covInfo list from varsVec and covarsVec
#' 
#' 
#' @param varsVec character vector of variables that need to be added  
#' @param covarsVec  character vector of covariates that need to be added
#' @return covInfo list of covariate info
#' @noRd
.buildcovInfo <- function(varsVec,covarsVec){
  checkmate::assert_character(varsVec,min.len = 1)
  checkmate::assert_character(covarsVec,min.len = 1)
  possiblePerms <- expand.grid(varsVec, covarsVec)
  possiblePerms <-
    list(
      as.character(possiblePerms[[1]]),
      as.character(possiblePerms[[2]])
    )
  names(possiblePerms) <- c("vars", "covars")
  covInfo <- list() # reversivle listVarName!!
  for (item in Map(list, possiblePerms$vars, possiblePerms$covars)) {
    listVarName <- paste0(item[[2]], item[[1]])
    covInfo[[listVarName]] <- list(varName = item[[1]], covariate = item[[2]])
  }
  
  return(covInfo)
}



#' Build updated from the covariate and variable vector list
#' 
#' 
#' @param ui compiled rxode2 nlmir2 model or fit  
#' @param varsVec character vector of variables that need to be added  
#' @param covarsVec  character vector of covariates that need to be added
#' @return updated ui with added covariates
#' @noRd
.buildupatedUI <- function(ui,varsVec,covarsVec){
  
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$ui
  }
  ui <- rxode2::assertRxUi(ui)
  
  
# construct covInfo
 covInfo <-  .buildcovInfo(varsVec,covarsVec)

# check if the covInfo is a list 
 checkmate::assert_list(covInfo)
 
 
 
 # Add covariates one after other
 
covSearchRes <- list()
covsAdded <- list() # to keep track of covariates added and store in a file
covsAddedIdx <- 1
for (x in covInfo) {
  covName <- paste0("cov_", x$covariate, "_", x$varName)
  
  
  if (length(covsAdded) == 0) {
    # covsAdded[[covsAddedIdx]] <- paste0(x$covariate, x$varName)
    covsAdded[[covsAddedIdx]] <- c(x$covariate, x$varName)
    
    # fit2 <-
    #   suppressWarnings(nlmixr2(updatedMod, data, est = getFitMethod(fitobject)))
    
    ui <- tryCatch(
      {
        res <- do.call(.builduiCovariate,c(ui,x))
        res # to return 'ui'
      },
      error = function(error_message) {
        print("error  while SIMULTANEOUSLY ADDING covariates")
        print(error_message)
        print("skipping this covariate")
        return(res) # return NA otherwise (instead of NULL)
      }
    )
    
    covSearchRes[[covsAddedIdx]] <- list(ui, covsAdded[[covsAddedIdx]], covName)
    covsAddedIdx <- covsAddedIdx + 1
  }
  else {
    # covsAdded[[covsAddedIdx]] <-
    #   paste0(covsAdded[[covsAddedIdx - 1]], "_", x$covariate, x$varName)
    covsAdded[[covsAddedIdx]] <-
      c(covsAdded[[covsAddedIdx - 1]], x$covariate, x$varName)
    
    # fit2 <-
    #   suppressWarnings(nlmixr2(updatedMod, data, est = getFitMethod(fitobject)))
    
    ui <- tryCatch(
      {
        res <- do.call(.builduiCovariate,c(ui,x))
        res # to return 'ui'
      },
      error = function(error_message) {
        print("error  while SIMULTANEOUSLY ADDING covariates")
        print(error_message)
        print("skipping this covariate")
        return(res) # return NA otherwise (instead of NULL)
      }
    )
    
    covSearchRes[[covsAddedIdx]] <- list(ui, covsAdded[[covsAddedIdx]], covName)
    covsAddedIdx <- covsAddedIdx + 1
  }
  
  covSearchRes
}


return(covSearchRes[length(covSearchRes)][[1]][[1]])

}





#' Add covariates and lasso string to ui  
#' 
#' 
#' @param ui compiled rxode2 nlmir2 model or fit  
#' @param varsVec character vector of variables that need to be added  
#' @param covarsVec  character vector of covariates that need to be added
#' @param tvalue float indicating tvalue to be updated for the lasso string
#' @return updated ui with added covariates
#' @noRd
.lassoUicovariate <- function(ui,varsVec,covarsVec,tvalue=0.05){
  
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$ui
  }
  ui <- rxode2::assertRxUi(ui)
  
  # Extract updated ui with added covariates 
  .ui1 <- .buildupatedUI(ui,varsVec,covarsVec)
  
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
  return(.getUiFunFromIniAndModel(.ui1, .ini, .newModel)())  
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
#' @author 
#'

.crossvalidationLasso <- function(data,ui,varsVec,covarsVec,tvalue=0.10,nfold=5,optcrit='llk',estmethod="focei",adapcoefs=NULL,stratVar=NULL){
  
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
    
  }
  else{
  mod <- .lassoUicovariate(ui,varsVec,covarsVec,tvalue = tvalue)
  }

  # Add fold column depending on the stratification
  if (!is.null(stratVar))
  {
    dataF <-    foldgen(data,nfold=5,stratVar = stratVar)
  } else {
    dataF <-    foldgen(data,nfold=5)
  }
  
  
  #Perform cross-validation
  for (f in 1:nfold){
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
      
      ofvList <- append(ofvList ,fitTest$objDf$OBJF[1] )
      
    }
    
    else if(optcrit=="llk") {
      
      if (!is.numeric(fitTest$objDf$`Log-likelihood`)){
        
        cli::cli_alert_danger("the 'fit' object needs to have an objective functions value associated with it")
        cli::cli_alert_info("try computing 'fit$objDf$$`Log-likelihood`' in console to compute and store the corresponding OBJF value")
        stop("aborting...objf value not associated with the current 'fit' object", call. = FALSE)
      }
    
    ofvList <- append(ofvList ,-2*fitTest$objDf$`Log-likelihood`[1])
  }
cli::cli_alert_success("Estimation complete for the fold number : {f}")
}
  
# Return the pOFV 
pOFV <- do.call(sum, ofvList)
return(pOFV)  

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
#' @author 
#'

.optimalTvaluelasso <- function(data,ui,varsVec,covarsVec,t_start=0.05,t_stop=0.25,t_step=0.05,stratVar = NULL,convergence="REACHMAX",...){
  
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
    
    ofValue <- suppressWarnings(.crossvalidationLasso(data,ui,varsVec,covarsVec,tvalue=t,nfold=5,optcrit='llk',estmethod="focei",adapcoefs=NULL,stratVar = NULL))
    
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
#' @param constraint theta cutoff. below cutoff then the theta will be fixed to zero
#' @return return data frame of final lasso coefficients
#' @noRd
#' @author 
#'
.lassoCoefficients <- function(fit,varsVec,covarsVec,constraint=1e-08,stratVar = NULL,...) {

  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  else {
      ui <- fit$ui
  }
  
  checkmate::assert_character(covarsVec)
  checkmate::assert_character(varsVec)
  checkmate::assert_double(constraint)
  
  #data
  data <- .normalizedData(fit,covarsVec)
  
  # construct covInfo
  covInfo <-  .buildcovInfo(varsVec,covarsVec)
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
  assign('debug',list(ui,varsVec,covarsVec,tvalue=optTvalue,...),globalenv())
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
  return(finalLasso)
}





#' Add covariates and adaptive lasso string to ui  
#' 
#' 
#' @param ui compiled rxode2 nlmir2 model or fit  
#' @param varsVec character vector of variables that need to be added  
#' @param covarsVec  character vector of covariates that need to be added
#' @param tvalue float indicating tvalue to be updated for the lasso string
#' @param adapcoefs dataframe with covariates estimates for the adaptive lasso
#' @return updated ui with added covariates
#' @noRd
.adaptivelassoUicovariate <- function(ui,varsVec,covarsVec,tvalue=0.05,adapcoefs=NULL){
  
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$ui
  }
  ui <- rxode2::assertRxUi(ui)
  
  # Extract updated ui with added covariates 
  .ui1 <- .buildupatedUI(ui,varsVec,covarsVec)
  
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
  return(.getUiFunFromIniAndModel(.ui1, .ini, .newModel)())  
}




#' Return Adaptive lasso coefficients after finding optimal t
#'
#' @param fit nlmixr2 fit.
#' @param varsVec character vector of variables that need to be added  
#' @param covarsVec  character vector of covariates that need to be added
#' @param constraint theta cutoff. below cutoff then the theta will be fixed to zero.
#' @return return data frame of final lasso coefficients
#' @noRd
#' @author 
#'
.adaptivelassoCoefficients <- function(fit,varsVec,covarsVec,constraint=1e-08,stratVar = NULL,...) {
  
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  else {
    ui <- fit$ui
  }
  
  checkmate::assert_character(covarsVec)
  checkmate::assert_character(varsVec)
  checkmate::assert_double(constraint)
  
  
  ## Get initial adaptive coefs from the regular lasso coefficients
  .adapcoefs=.lassoCoefficients(fit,varsVec,covarsVec,constraint=1e-08,stratVar = NULL,...)
  
  #data
  data <- .normalizedData(fit,covarsVec)
  
  # construct covInfo
  covInfo <-  .buildcovInfo(varsVec,covarsVec)
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
  assign('debug',list(ui,varsVec,covarsVec,tvalue=optTvalue,...),globalenv())
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
  
  #Multiply by factor 
  finalLasso <- finalLasso *factor 
  return(finalLasso)
}



#' Regular lasso model
#'
#' @param fit nlmixr2 fit.
#' @param varsVec character vector of variables that need to be added  
#' @param covarsVec  character vector of covariates that need to be added
#' @param constraint theta cutoff. below cutoff then the theta will be fixed to zero.
#' @param lassotype must be  'regular' , 'adaptive', 'adjusted'
#' @return return fit of the regular coefficients
#' @noRd
#' @author 
#'
.regularmodel <- function(fit,varsVec,covarsVec,constraint=1e-08,lassotype='regular',stratVar = NULL,...) {
  
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  else {
    ui <- fit$ui
  }
  
  checkmate::assert_character(covarsVec)
  checkmate::assert_character(varsVec)
  checkmate::assert_double(constraint)
  
  
  
  if (lassotype=="regular") {
.coefValues <- .lassoCoefficients(fit,varsVec,covarsVec,constraint=1e-08,stratVar = NULL,...)
  }
  else if (lassotype=="adaptive") {
.coefValues <- .adaptivelassoCoefficients(fit,varsVec,covarsVec,constraint=1e-08,stratVar = NULL,...) 
  }
  
  else if (lassotype=="adaptive"){
    
  .coefValues <- .adjustedlassoCoefficients(fit,varsVec,covarsVec,constraint=1e-08,stratVar = NULL,...)    
  }

# Extract updated ui with added covariates 
.ui1 <- .buildupatedUI(ui,varsVec,covarsVec)

# Extract all the covariate parameters constructed 
.covsparams <- .ui1$muRefCovariateDataFrame$covariateParameter
checkmate::assertCharacter(.covsparams,min.len = 1,any.missing = FALSE )


# Fix the covariate parameter values
.modelfun <- .ui1$funTxt
for ( i in seq_along(.covsparams)){
  .modelfun  <- gsub(.covsparams[i], .coefValues[[i]], .modelfun)
}

# construct the final updated model string 

.newmodelfun <- paste0(tvalueString,absString,ratioString,factorString,.modelfun,sep="")
.newModel <- eval(parse(text = paste0("quote(model({",paste0(as.character(.newmodelfun),collapse="\n"), "}))")))

# build ui
.ini <- .ui1$iniDf
.ini <- .ini[ ! .ini$name %in% .covsparams, ]
.ini <- as.expression(lotri::as.lotri(.ini))
.ini[[1]] <- quote(`ini`)
# return the updated model
.ui2 <- .getUiFunFromIniAndModel(.ui1, .ini, .newModel)()
#estimation method
estmethod=getFitMethod(fit)
return(nlmixr2(ui2,data,est=estmethod)) 
}


#' Return Adjusted adaptive lasso coefficients after finding optimal t
#'
#' @param fit nlmixr2 fit.
#' @param varsVec character vector of variables that need to be added  
#' @param covarsVec  character vector of covariates that need to be added
#' @param constraint theta cutoff. below cutoff then the theta will be fixed to zero.
#' @return return data frame of final lasso coefficients
#' @noRd
#' @author 
#'
.adjustedlassoCoefficients <- function(fit,varsVec,covarsVec,constraint=1e-08,stratVar = NULL,...) {
  
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  else {
    ui <- fit$ui
  }
  
  checkmate::assert_character(covarsVec)
  checkmate::assert_character(varsVec)
  checkmate::assert_double(constraint)
  
  
  ## Get updated ui with covariates and without lasso factor 
  .mod1 <- .buildupatedUI(ui,varsVec,covarsVec)
  
  # construct covInfo
  covInfo <-  .buildcovInfo(varsVec,covarsVec)
  # construct covNames   
  covNames <- list()    
  for ( x in covInfo) {
    covName <- paste0("cov_", x$covariate, "_", x$varName)
    covNames <- append(covNames,covName)
  }
  
  #data
  data <- .normalizedData(fit,covarsVec)
  
  #Fit
  fit1 <- nlmixr2(.mod1,data,est=getFitMethod(fit))
  
  # Construct AL coefficents
  .adapcoefs=data.frame(as.list(fit1$parFixedDf[row.names(fit1$parFixedDf) %in% covNames,"Estimate"]/fit1$parFixedDf[row.names(fit1$parFixedDf) %in% covNames,"SE"]))
  
  # Extract optimal t-value
  optTvalue <- .optimalTvaluelasso(data,ui,varsVec,covarsVec,t_start=0.05,t_stop=0.25,t_step=0.05,estmethod=getFitMethod(fit),adapcoefs=.adapcoefs,stratVar = NULL,...)   
  
  # Refit model with the optimal t-value 
  assign('debug',list(ui,varsVec,covarsVec,tvalue=optTvalue,...),globalenv())
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
  # ratio 
  ratio <- abssum/optTvalue
  # factor 
  factor <- exp(1-ratio)
  
  
  # Apply lasso constraint 
  finalLasso <- as.data.frame(lapply(covestAdap, function(x) ifelse(abs(x) <= constraint , 0, x)),row.names = NULL)
  
  #Multiply by factor 
  finalLasso <- finalLasso *factor 
  return(finalLasso)
}


#' Stratified cross-validation fold generator function inspired from the caret
#' @param data data frame used in the analysis
#' @param nfold number of k-fold cross validations. Default is 5 
#' @param stratVar  Stratification Variable. Default is NULL and ID is used for CV
#' @return return dataframe with the fold column attached
#' @noRd
#' @author 


foldgen <-  function(data,nfold=5,stratVar=NULL){

  
  # check if data frame
  checkmate::assert_data_frame(data,min.cols = 7)


# check if user want to stratify on a variable , if not default is on individual

if(!is.null(stratVar)){
  checkmate::assertCharacter(stratVar,len = 1,any.missing = FALSE )
  stratCheck <- intersect(names(data), stratVar)
  if(!is.null(stratCheck)){
    y <- data[,stratCheck]
  }
  else {
    stop(paste0(stratVar, "not in the data to stratify"))
  }
} else {
  
  # extract ID column from the data frame 
  ID <- .idColumn(data)
  # Extract list of individuals
  y <- unique(data[,ID])
}
## Group based on magnitudes and sample within groups 

if(is.numeric(y))
{
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
  y <- cut(
    y, 
    unique(
      quantile(y,
               probs =
                 seq(0, 1, length = cuts))), 
    include.lowest = TRUE)
}

if(nfold < length(y))
{
  ## reset levels so that the possible levels and 
  ## the levels in the vector are the same
  y <- factor(as.character(y))
  numInClass <- table(y)
  foldVector <- vector(mode = "integer", length(y))
  
  ## For each class, balance the fold allocation as far 
  ## as possible, then resample the remainder.
  ## The final assignment of folds is also randomized. 
  for(i in 1:length(numInClass))
  {
    ## create a vector of integers from 1:k as many times as possible without 
    ## going over the number of samples in the class. Note that if the number 
    ## of samples in a class is less than k, nothing is producd here.
    seqVector <- rep(1:nfold, numInClass[i] %/% nfold)
    ## add enough random integers to get  length(seqVector) == numInClass[i]
    if(numInClass[i] %% nfold > 0) seqVector <- c(seqVector, sample(1:nfold, numInClass[i] %% nfold))
    ## shuffle the integers for fold assignment and assign to this classes's data
    foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
  }
} else foldVector <- seq(along = y)


out <- split(seq(along = y), foldVector)
names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), sep = "")
out <- foldVector

if(!is.null(stratVar)){
out <- cbind(data,fold=out)
} else {
indv <- unique(data[,ID])
out <- merge(data,cbind(ID=indv,fold=out))  
}
out
}







  

