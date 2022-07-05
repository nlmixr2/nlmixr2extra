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

.addCovariate <- function(ui,varName,covariate,norm=FALSE,...) {
  
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
  lapply(.model,.expandRefMu,murefDf=ui$muRefDataFrame,covDf=.covdf,pars=.pars,...)
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
.expandRefMu <- function(x,murefDf,covDf,pars,...) {
  if (is.name(x)) {
    currparam <- as.character(x)
    if (currparam %in% pars){
      return (str2lang(.expandPopExpr(currparam,murefDf,covDf,...)))
    } 
  }else if(is.call(x)) {
    return(as.call(c(list(x[[1]]), lapply(x[-1], .expandRefMu, murefDf=murefDf, covDf=covDf, pars=pars,...))))
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

.normalizeDf <- function(data, covariate) {
  
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
#' 
 #' @return data frame with all normalized covariates
#' @noRd
.testCovariate <- function(fit,covarsVec) {
  
  if (!inherits(fitobject, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  
  
  checkmate::assert_character(covarsVec)

  # get data 
  data <- nlme::getData(fitobject)
  
  ## 
  .normalizedDFs <- lapply(covarsVec,.normalizeDf,data=data)
  
  # final data frame of normalized covariates
  
  .finalDf <- Reduce(merge,.normalizedDFs)
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
lst <- .addCovariate(ui,varName,covariate)
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
 
 
 
 # Add covariates to 
 
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








