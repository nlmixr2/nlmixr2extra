
#' lasso method for predictive covariate model building
#' lasso is suggested over SCM to remedy selection bias ,poor 
#' predictive performance and speed 
#' 
#' Initializing covariates before estimation
#' 
#' 
#' @param ui compiled rxode2 nlmir2 model or fit  
#' @param varName  the variable name to which the given covariate is to be added
#' @param covariate the covariate that needs string to be constructed 
addCovariate <- function(ui,varName,covariate) {

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


  
}



#' get the population parameter from variable name
#'
#' @param ui compiled rxode2 nlmir2 model or fi
#' @param varName the variable name to which the given covariate is to be added
#'
#' @return population parameter variable
#' @noRd
.getThetaName <- function(ui,varName){
  .split <-  ui1$getSplitMuModel
  if (varName %in% names(.split$pureMuRef)){
    return(varName)
  }
  .w <- which(.split$pureMuRef==varName)
  if(length(.w)==1){
    return(names(.split$pureMuRef)[.w])  
  }
}




#'
#' @param fitobject an nlmixr2 'fit' object
#' @param fstring a string giving the entire expression for the model function string
#' @param covNames  a list of covariate names (parameters) that need to be estimates
#' @param initialEst the initial estimate for the covariate parameters to be estimated; default is 0
#' @param initialEstLB a lower bound for the covariate parameters to be estimated; default is -Inf
#' @param initialEstUB an upper bound for the covariate parameters to be estimated; default is Inf
#'
#' @return updated model object with the modified initial values
#' @export
#' @author Vipul Mann, Matthew Fidler
#'
initializeCovars <- function(fitobject,
                             fstring,
                             covNames,
                             initialEst,
                             initialEstLB,
                             initialEstUB) {
  updatedMod <- paste0("model(fitobject,{", fstring, "})")
  updatedMod <- eval(parse(text = updatedMod))

  ini2 <- updatedMod$iniDf
  for (covName in covNames) {
    ini2[ini2$name == covName, "est"] <- initialEst
    ini2[ini2$name == covName, "lower"] <- initialEstLB
    ini2[ini2$name == covName, "upper"] <- initialEstUB
  }

  class(ini2) <- c("nlmixr2Bounds", "data.frame")
  updatedMod$ini <- ini2

  updatedMod
}
#' 
getFitMethod <- function(fit) {
  if (!(inherits(fit, "nlmixr2FitCore"))) {
    stop("'fit' needs to be a nlmixr2 fit", call. = FALSE)
  }
  fit$est

}
#'
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

      if (item[[2]] %in% catCovariates) {
        covInformation[[listVarName]]$categorical <- TRUE
      }
      else {
        covInformation[[listVarName]]$categorical <- FALSE
      }

      if (listVarName %in% names(covInformation)) {
        covInfo[[listVarName]] <- c(list(varName = item[[1]], covariate = item[[2]]), covInformation[[listVarName]])
      }
      else {
        covInfo[[listVarName]] <- list(varName = item[[1]], covariate = item[[2]])
      }
    }
#' Prepare normalized and factor covariate strings for the model
#'
#' @param data a dataframe containing the dataset that needs to be used
#' @param covariate the covariate that needs string to be constructed 
#' @param varName the variable name to which the given covariate is to be added
#' @param addfactor a boolean indicating if the lasso constraint factor needs to be added; default is FALSE
#' @param factor a int of lasso constraint factor already calculated
#'
#' @return a list comprising the covariate factor string, and covariate names
#' @export
#' @author 
#'



normcovariateString <-  function(data,covariate,varName,addfactor=FALSE,factor=NULL){
    
# Extract individual ID from column names
colNames <- colnames(data)
  colNamesLower <- tolower(colNames)
  if ("id" %in% colNamesLower) {
    uidCol <- colNames[which("id" %in% colNamesLower)]
  } else {
    uidCol <- "ID"
  }
# Function to calculate population standard deviation from group means 
sd.p=function(x){sd(x)*sqrt((length(x)-1)/length(x))}
# mean by groups (Individual)
groupMeans <- with(data, ave(get(covariate),get(uidCol), FUN = function(x) mean(x, na.rm = TRUE)))
# pop mean
popMean <- mean(groupMeans)
 # pop std
popStd <- sd.p (groupMeans)

# Prepare variable-covariate names                            
covNames <- paste0("cov_", covariate, "_", varName)
 # Add Normalization term                 
standardization <- paste0("(",covariate,"-",popMean,")","/",popStd) 
# Prepare Covariate Name modification 
if(addfactor){
covNameMod <- paste0(covNames," * ",standardization," * ",factor)
}
                             else{
covNameMod <- paste0(covNames," * ",standardization)
}
 list( covNameMod, covNames)
}   
#' Create a lasso model fit for the cross-validation
#'
#' @param covInfo a list containing information about each variable-covariate pair
#' @param fitobject an nlmixr2 'fit' object
#' @param addfactor a boolean indicating if the lasso constraint factor needs to be added; default is FALSE
#' @param factor a int of lasso constraint factor already calculated
#'
#' @return a fit  comprising the updated model with normlization and lasso constraints. 
#' @export
#' @author 

setupLassomodel <- function(covInfo, fitobject,addfactor=FALSE,factor=NULL) {
    
# covariate factor strings to be added
covfactorStrings <- list()
# covariate name list
covNamesList <- list()
# get data from fit object
data <- nlme::getData(fitobject)
for (f in 1:length(covInfo)) {
    
    varName <- covInfo[[f]]$varName
    covariate <- covInfo[[f]]$covariate
    x <- normcovariateString(data,covariate,varName,addfactor,factor)
    covfactorStrings <- append(covfactorStrings,x[[1]])
    covNamesList <- append(covNamesList,x[[2]])
}
# Extract fitobject string
funstring <- fitobject$ui$fun.txt
funstringSplit <-  unlist(strsplit(funstring, split = "\\\n"))
# Generate a variate vector from the covInfo
varsTemp <- list()
for ( i in 1:length(covInfo)){
  varsTemp <- append(varsTemp,covInfo[[i]]$varName) 
}
varsList <- unique(varsTemp)   
    
# construct covariate factor string
fstringList <- list()
for (var in varsList){
 idx <- grep(paste0("(?<!\\w)", var), funstringSplit, perl = TRUE)    
covNameMod <- paste(grep(var,covfactorStrings, perl = TRUE,value = TRUE),collapse="+")  
isLog <-grepl("\\bexp *\\(",funstringSplit[[idx]], perl = TRUE) 
fstringList <- append(fstringList,addCovariate(
        funstringSplit[[idx]],
        var,
        covNameMod,
        names(fitobject$theta),
        isLog
      )  )      
}
fstring_new <- paste(fstringList,collapse = "\n")
# Initialize Covariates 
updatedMod  <- initializeCovars(fitobject,fstring_new,covNamesList,initialEst=0,initialEstLB=-Inf,initialEstUB=Inf)
reassignVars <- rownames(fitobject$parFixedDf)[fitobject$parFixedDf$Estimate != fitobject$parFixedDf[, "Back-transformed"] & fitobject$parFixedDf[, "Back-transformed"] == 0]
if (length(reassignVars) > 0) {
    ini2 <- as.data.frame(updatedMod$ini)
    for (r in reassignVars) {
          ini2[ini2$name == r, "est"] <- 1.0
          cli::cli_alert_warning("reasssigned initial value for {r} to 1.0")
        }
        class(ini2) <- c("nlmixr2Bounds", "data.frame")
        updatedMod$ini <- ini2
      }
# fit the model
 fit2 <- tryCatch(
      {
        fit2 <-
          suppressWarnings(nlmixr2(updatedMod, data, est = getFitMethod(fitobject)))
        fit2 # to return 'fit2'
      },
      error = function(error_message) {
        print("error fitting the model after REMOVING the covariates")
        print(error_message)
        print("returning the previous fitobject")
        return(fitobject) # return NA otherwise (instead of NULL)
      }
    )
#print (fit2$ui)
return(fit2)
    }
#' Return lasso constraint for a specified t-value
#'
#' @param fitobject an nlmixr2 'fit' object
#' @param covNames  a list of covariate names (parameters) that need to be estimates
#' @param tvalue a desired tvalue for constructing lasso constraint
#' @return updated model object with the modified initial values
#' @export
#' @author 
#'
lassoConstraint <- function(fitobject,
                             covInfo,
                             tvalue=0.5) {
# covariate name list
covNames <- list()
# get data from fit object
data <- nlme::getData(fitobject)
# construct covNames
for (f in 1:length(covInfo)) {
    varName <- covInfo[[f]]$varName
    covariate <- covInfo[[f]]$covariate
    x <- normcovariateString(data,covariate,varName)
    covNames <- append(covNames,x[[2]])
}
#absolute sum of  lasso THETAs
abssum <- sum(abs(fitobject$parFixedDf[row.names(fitobject$parFixedDf) %in% covNames,"Estimate"]))
# ratio 
ratio <- abssum/tvalue
# factor 
factor <- exp(1-ratio)
return(factor)
}
#' Perform cross validation and return predictive-objective-function value
#'
#' @param data dataset containing all the required information.
#' @param mod an model code for nlmixr without covariates for the initial model.
#' @param covInfo a list containing information about each variable-covariate pair
#' @param tvalue a desired tvalue for constructing lasso constraint
#' @param nfold number of folds for cross-validation
#' @return predictive-objective-function value
#' @export
#' @author 
#'

crossvalidationLasso <- function(data,mod,covInfo,tvalue=0.10,nfold=2){

# Extract list of individuals
indv <- unique(data$ID)
# Fold splits   
foldSplits <- cut(1:length(indv), breaks = nfold, labels = FALSE)
# Randomize splits
foldSplitsR <- sample(foldSplits, length(foldSplits)) 
# data frame of mapping of Indvidual to fold splits
foldInddf <- data.frame(ID=indv,fold=foldSplitsR )
# List of objective-function values
ofvList <- list()
for (f in 1:nfold){
# Training dataset 
trainIndv  <- (foldInddf[foldInddf$fold!=f,])$ID
trainData <- data[data$ID %in% trainIndv, ]
# Testing dataset
testIndv  <- (foldInddf[foldInddf$fold==f,])$ID
testData <- data[data$ID %in% testIndv, ]

# Training Estimation   
fit1Train <-  nlmixr2(mod,trainData ,"focei")
# Setup  model without lasso constraint 
fit2Train <- setupLassomodel(covInfo,fit1Train,addfactor=FALSE,factor=NULL)
# calculate lasso constraint from the estimates 
lassoFactor <- lassoConstraint(fit2Train,covInfo,tvalue)
# Testing Initial Model
fit1Test <-  nlmixr2(mod,testData ,"focei")
# Setup model with constraint lasso constraint from the training
fit2Test <- setupLassomodel(covInfo,fit1Test,addfactor=TRUE,factor=lassoFactor)
# Store OFV from the test dataset 
ofvList <- append(ofvList ,fit2Test$objDf$OBJF )   
}

# Return the pOFV 
pOFV <- do.call(sum, ofvList)
return(pOFV)
}
#' Run cross validation and return final lasso coefficients
#'
#' @param data dataset containing all the required information.
#' @param mod an model code for nlmixr without covariates for the initial model.
#' @param covInfo a list containing information about each variable-covariate pair
#' @param t_start a desired t start value for constructing lasso constraint
#' @param t_stop a desired t stop value for constructing lasso constraint
#' @param t_step a desired t step value for constructing lasso constraint
#' @param nfold number of folds for cross-validation
#' @return predictive-objective-function value
#' @export
#' @author 
#'

runLasso <- function(data,mod,covInfo,t_start=0.05,t_stop=0.15,t_step=0.05,nfold=2){

# generate vector of t values
tvalues <- seq(t_start,t_stop,by=t_step)    

 # create a data frame of pOFV values for t values
pofvList <- data.frame()
for (t in tvalues){
 pofv <- data.frame(tvalue=t,POFV=suppressWarnings(suppressMessages(crossvalidationLasso(data,mod,covInfo,tvalue=t,nfold))))
pofvList <- rbind(pofvList,pofv)
}
# optimal t value
 optimal_t <- pofvList[which.min(pofvList$POFV),]$tvalue
    
# use optimal t-value with the regular model
fit1 <- setupLassomodel(covInfo,fit,addfactor=FALSE,factor=NULL)
# lasso factor from tvalue
lassoFactor <- lassoConstraint(fit1,covInfo,optimal_t)
# final model
finalfit <- setupLassomodel(covInfo,fit,addfactor=TRUE,factor=lassoFactor)
## compute final lasso coefficients 
# covariate name list
covNames <- list()
# get data from fit object
data <- nlme::getData(finalfit)
# construct covNames
for (f in 1:length(covInfo)) {
    varName <- covInfo[[f]]$varName
    covariate <- covInfo[[f]]$covariate
    x <- normcovariateString(data,covariate,varName)
    covNames <- append(covNames,x[[2]])
}
# final lasso coefficients
flc <- finalfit$parFixedDf[row.names(finalfit$parFixedDf) %in% covNames,"Estimate"]*lassoFactor
return(list(flc,lassoFactor))   
}