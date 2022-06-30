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
#' @norm  boolean indicating if the covariate that needs to be added is normalized; default is FALSE

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
#' @param addfactor a boolean indicating if the lasso constraint factor needs to be added; default is FALSE
#' @param factor a int of lasso constraint factor already calculated
#'
#' @return expanded expression
#' @noRd
.expandPopExpr <- function(popParam,murefDf,covDf,addfactor = FALSE,factor) {
  
  .par1 <- murefDf[murefDf$theta==popParam,]
  .res <- paste0(.par1$theta,"+",.par1$eta)
  .w <- which(covDf$theta==popParam)
  if(length(.w)>0){
    if (addfactor){
      checkmate::assert_numeric(factor,len = 1)
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


## recursive factor string 

#1) var and covar list
#2) mapply to add covar 
#3) Intilaize covariates 
#4) Recursion





