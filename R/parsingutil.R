
#' Add covariate 
#' 
#' 
#' @param ui compiled rxode2 nlmir2 model or fit  
#' @param varName  the variable name to which the given covariate is to be added
#' @param covariate the covariate that needs string to be constructed
#' @author Matthew Fidler, Vishal Sarsani
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
#' @noRd
#' @author Matthew Fidler
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
#' @author Matthew Fidler
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
#' 
#' @author Matthew Fidler
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
#' @author Matthew Fidler, Vishal Sarsani

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

#' Build ui from the covariate 
#' 
#' 
#' @param ui compiled rxode2 nlmir2 model or fit  
#' @param varName  the variable name to which the given covariate is to be added
#' @param covariate the covariate that needs string to be constructed
#' @return ui with added covariate
#' @noRd
#' @author Matthew Fidler, Vishal Sarsani
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
#' @author Matthew Fidler, Vishal Sarsani
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
#' @param indep a boolean indicating if the covariates should be added independently, or
#'  sequentially (append to the previous model)
#' @return updated ui with added covariates
#' @noRd
#' @author Matthew Fidler, Vishal Sarsani
.buildupatedUI <- function(ui,varsVec,covarsVec,indep=FALSE){
  
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
  
  if (indep) {
    for (i in seq_along(covInfo)) {
      x <- covInfo[[i]]
      covName <- paste0("cov_", x$covariate, "_", x$varName)
      
      reassignVars <- rownames(fit$parFixedDf)[fit$parFixedDf$Estimate != fit$parFixedDf[, "Back-transformed"] & fit$parFixedDf[, "Back-transformed"] == 0]
      
      res <- do.call(.builduiCovariate,c(ui,x))
      if (length(reassignVars) > 0) {
        ini2 <- as.data.frame(res$ini)
        
        for (r in reassignVars) {
          ini2[ini2$name == r, "est"] <- 1.0
          cli::cli_alert_warning("reasssigned initial value for {r} to 1.0")
        }
        
        class(ini2) <- c("nlmixr2Bounds", "data.frame")
        res$ini <- ini2
      }
      
      covSearchRes[[i]] <- list(ui,c(x$covariate, x$varName), covName)
    }
    return(covSearchRes)
  }
  
  else {
    ## Add all at once
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
}

#' Make dummy variable cols and updated covarsVec
#' @param data data frame used in the analysis
#' @param nfold number of k-fold cross validations. Default is 5 
#' @param stratVar  Stratification Variable. Default is NULL and ID is used for CV
#' @return return dataframe with the fold column attached
#' @noRd
#' @author Matthew Fidler, Vishal Sarsani

.addCatCovariates <- function(data,covarsVec,catcovarsVec) {
  
  # check for valid inputs
  checkmate::assert_data_frame(data,min.cols = 7)
  checkmate::assert_character(covarsVec)
  checkmate::assert_character(catcovarsVec)
  #create new catcovarsvec
  newcatvars <- character(0)
  for (col in catcovarsVec) {
    if (is.factor(data[[col]])) {
      uniqueVals <- levels(data[[col]])
      if (any(is.na(data[[col]]))) {
        uniqueVals <- c(uniqueVals, NA)
      }
      # Else by order values appear.
    } else {
      uniqueVals <- unique(data[[col]])
    }
    uniqueVals <- as.character(uniqueVals)
    
    # Remove NA values and first dummy
    uniqueVals <- uniqueVals[!is.na(uniqueVals)][-1]
    
    
    for (uniqueValue in uniqueVals) {
      colname= paste0(col,"_",uniqueValue)
      data[,colname] <- as.factor(match(data[[col]],"0",nomatch=0))
      newcatvars <- c(newcatvars ,colname)
    }
  }
  
  # Remove original categorical variables
  updatedData <- data[ , !(names(data) %in% catcovarsVec)]
  #Update entire covarsvec with added categorical variables
  updcovarsVec <- c(covarsVec,newcatvars)
  
  return(list(updatedData,updcovarsVec))
  
}












