
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
#' @param indep a boolean indicating if the covariates should be added independently, or
#'  sequentially (append to the previous model)
#' @return updated ui with added covariates
#' @noRd
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
#' @author 

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


#' Add covariate expression to a function string
#'
#' @param funstring a string giving the expression that needs to be modified
#' @param varName the variable to which the given string corresponds to in the model expression
#' @param covariate the covariate expression that needs to be added (at the appropriate place)
#' @param theta a list of names of the 'theta' parameters in the 'fit' object
#' @param isLog a boolean signifying the presence of log-transformation in the funstring
#' @return returns the modified string with the covariate added to function string
#' @author Vipul Mann, Matthew Fidler
#' @export

oldaddCovariate <-
  function(funstring,
           varName,
           covariate,
           theta,
           isLog) {
    f <- function(x, isCov = FALSE) {
      if (is.atomic(x)) {
        return(x)
      }
      else if (is.name(x)) {
        if (isCov) {
          if (any(as.character(x) == theta) && regexpr("^cov\\_", as.character(x)) == -1) {
            return(eval(parse(
              text = paste0("quote(", x, "+", covariate, ")")
            )))
          }
        }
        return(x)
      } else if (is.pairlist(x)) {
        return(x)
      } else if (is.call(x)) {
        if (identical(x[[1]], quote(`{`))) {
          return(paste(unlist(lapply(x[-1], function(x) {
            gsub(" +", "", paste(deparse1(f(x, isCov)), collapse = ""))
          })), collapse = "\n"))
        }
        else if (identical(x[[1]], quote(`~`)) ||
                 identical(x[[1]], quote(`=`)) ||
                 identical(x[[1]], quote(`<-`))) {
          if (length(x[[2]]) == 1) {
            if (as.character(x[[2]]) == varName) {
              isCov <- TRUE
            }
          }
          return(as.call(lapply(x, f, isCov = isCov)))
        }
        else {
          return(as.call(lapply(x, f, isCov = isCov)))
        }
      }
    }
    
    f(eval(parse(text = paste0(
      "quote({", funstring, "})"
    ))))
    
    f2 <- function(varName) {
      funstringLhsRhs <- strsplit(funstring, "(<-|=)")[[1]]
      funstringRhs <- funstringLhsRhs[2]
      funstringLhs <- funstringLhsRhs[1]
      
      expr <-
        paste0("(", funstringRhs, ")", "*", "(", covariate, ")")
      expr <-
        gsub(" ", "", expr, perl = TRUE) # remove white spaces from the above string
      funstringLhs <-
        gsub(" ", "", funstringLhs, perl = TRUE) # remove white spaces from the above string
      return(paste0(funstringLhs, "<-", expr))
    }
    
    if (!isLog) {
      f2(varName)
    }
    else {
      f(eval(parse(text = paste0(
        "quote({", funstring, "})"
      ))))
    }
  }

















#' Remove covariate expression from a function string
#'
#' @param funstring a string giving the expression that needs to be modified
#' @param varName the variable to which the given string corresponds to in the model expression
#' @param covariate the covariate expression that needs to be removed (from the appropriate place)
#' @param theta a list of names of the 'theta' parameters in the 'fit' object
#'
#' @return returns the modified string with the covariate removed from the function string
#'
#' @author Vipul Mann, Matthew Fidler
#' @export
#'
removeCovariate <- function(funstring, varName, covariate, theta) {
  covariate <-
    gsub(" ", "", covariate, perl = TRUE) # remove white spaces from the above string
  covariateSplit <- strsplit(covariate, "\\*|\\+")[[1]]
  
  f <- function(x, isCov = FALSE) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
      return(x)
    } else if (is.pairlist(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`{`))) {
        return(paste(unlist(lapply(x[-1], function(x) {
          gsub(" +", "", paste(deparse1(f(x, isCov)), collapse = ""))
        })), collapse = "\n"))
      }
      else if (identical(x[[1]], quote(`~`)) ||
               identical(x[[1]], quote(`=`)) ||
               identical(x[[1]], quote(`<-`))) {
        if (length(x[[2]]) == 1) {
          if (as.character(x[[2]]) == varName) {
            isCov <- TRUE
          }
        }
        return(as.call(lapply(x, f, isCov = isCov)))
      }
      else if (isCov && identical(x[[1]], quote(`+`))) {
        # Case 1: + centered_factor*cov_factor
        if (length(x[[3]]) > 1 &&
            identical(quote(`*`), x[[3]][[1]])) {
          if (as.character(x[[3]][[2]]) %in% covariateSplit &&
              as.character(x[[3]][[3]]) %in% covariateSplit) {
            # return(x[[2]])
            if (length(x[[2]]) == 1) {
              return(x[[2]])
            }
            return(as.call(lapply(x[[2]], f, isCov = isCov)))
          }
        }
        return(as.call(lapply(x, f, isCov = isCov)))
      }
      
      else if (isCov && identical(x[[1]], quote(`*`))) {
        # Case 2: *(centered_factor*cov_factor)
        if (length(x[[3]]) > 1 &&
            identical(x[[3]][[1]], quote(`(`))) {
          covExpr <- x[[3]][[2]]
          
          if (length(covExpr) > 1 &&
              identical(quote(`*`), covExpr[[1]])) {
            if (as.character(covExpr[[2]]) %in% covariateSplit &&
                as.character(covExpr[[3]]) %in% covariateSplit) {
              # return(x[[2]])
              return(as.call(lapply(x[[2]], f, isCov = isCov)))
            }
          }
          
          else if (length(covExpr) > 1 &&
                   identical(quote(`+`), covExpr[[1]])) {
            if (length(covExpr[[3]]) > 1 &&
                identical(quote(`*`), covExpr[[3]][[1]])) {
              if (as.character(covExpr[[3]][[2]]) %in% covariateSplit &&
                  as.character(covExpr[[3]][[3]]) %in% covariateSplit) {
                # return(x[[2]])
                return(as.call(lapply(x[[2]], f, isCov = isCov)))
              }
            }
          }
        }
        
        return(as.call(lapply(x, f, isCov = isCov)))
      }
      
      
      else {
        return(as.call(lapply(x, f, isCov = isCov)))
      }
    }
  }
  
  nch <- 0
  ret <- funstring
  
  while (nch != nchar(ret)) {
    nch <- nchar(ret)
    
    ret <- f(eval(parse(text = paste0(
      "quote({", ret, "})"
    ))))
  }
  
  return(ret)
}


#' Adding covariate to a given variable in an nlmixr2 model expression
#'
#' @param fitobject an nlmixr2 'fit' object
#' @param varName a string giving the variable name to which covariate needs to be added
#' @param covariate a string giving the covariate name; must be present in the data used for 'fit'
#' @param norm the kind of normalization to be used while normalizing covariates; para
#' @param norm_type a string defining operator to be used for transforming covariates using 'norm'; must be one among 'mul', 'div', 'sub', 'add'
#' @param categorical a boolean indicating if the 'covariate' is categorical
#' @param isHS a boolean indicating if 'covariate' is of Hockey-stick kind
#' @param initialEst the initial estimate for the covariate parameters to be estimated; default is 0
#' @param initialEstLB a lower bound for the covariate parameters to be estimated; default is -Inf
#' @param initialEstUB an upper bound for the covariate parameters to be estimated; default is Inf
#'
#' @return a list with the updated model expression and data with columns corresponding to normalized covariate(s) appended
#' @author Vipul Mann, Matthew Fidler
#' @export
#'
addCovVar <- function(fitobject,
                      varName,
                      covariate,
                      norm = c("median", "mean", "autoscale"),
                      norm_type = c("mul", "div", "sub", "add", "autoscale"),
                      categorical = FALSE,
                      isHS = FALSE,
                      initialEst = 0,
                      initialEstLB = -Inf,
                      initialEstUB = Inf) {
  if (!inherits(fitobject, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  
  if (inherits(norm, "numeric")) {
    checkmate::assert_numeric(
      norm,
      len = 1,
      lower = .Machine$double.eps,
      null.ok = FALSE,
      any.missing = FALSE
    )
  }
  else {
    norm <- match.arg(norm)
  }
  
  norm_type <- match.arg(norm_type)
  
  if (missing(norm_type)) {
    norm_type <- "sub"
  }
  norm_ops <-
    list(
      "div" = `/`,
      "mul" = `*`,
      "sub" = `-`,
      "add" = `+`,
      "autoscale" = c(`-`, `/`)
    )
  normOp <- norm_ops[[norm_type]]
  
  funstring <- fitobject$ui$fun.txt
  funstringSplit <-
    unlist(strsplit(funstring, split = "\\\n")) # split the string at \n
  
  # # look for the string that needs to be modified
  # idx <- grep(varName, funstringSplit)
  idx <-
    grep(paste0("(?<!\\w)", varName), funstringSplit, perl = TRUE) # varName not preceded by any other character
  
  if (length(idx) >= 2) {
    print(funstring)
    stop("cannot update more than one variable at a time")
  }
  
  # Get the normalization value from data
  data <- nlme::getData(fitobject)
  
  # search the dataframe for a column name of 'ID'
  colNames <- colnames(data)
  colNamesLower <- tolower(colNames)
  if ("id" %in% colNamesLower) {
    uidCol <- colNames[which("id" %in% colNamesLower)]
  }
  else {
    uidCol <- "ID"
  }
  
  if (covariate %in% colnames(data)) {
    
    # infer categorical covariates if not specified
    if (missing(categorical) && is.factor(data[, covariate])) {
      categorical <- TRUE
      cli::cli_alert_info("treating {covariate} as categorical variable ...")
    }
    
    if (inherits(norm, "numeric")) {
      normValVec <- norm
    }
    else {
      if (norm %in% "autoscale") {
        normValVecMean <- mean(data[, covariate])
        normValVecSd <- sd(data[, covariate])
        
        if (normValVecSd == 0) {
          # if (!categorical){  # print warning only for non-categorical variables
          #   cli::cli_alert_warning('normalization value for subject ID: {x} is zero; using with 1...')
          # }
          normValVecSd <- 1
        }
        
        normValVec <- list(normValVecMean, normValVecSd)
      }
      
      else if (norm %in% "mean") {
        # mean of the mean values
        uids <- unlist(unname(data[, uidCol]))
        normValVec <- lapply(uids, function(x) {
          datSlice <- data[data[, uidCol] == x, ]
          normVal <- mean(unlist(datSlice[covariate]))
          if (normVal == 0) {
            if (!categorical && !all(unlist(datSlice[covariate]) == 0)) { # print warning only for non-categorical variables
              cli::cli_alert_warning("normalization value for subject ID: {x} is zero; using with 1...")
            }
            normVal <- 1
          }
          
          normVal
        })
        
        normValVec <- mean(unlist(normValVec))
      }
      else {
        # mean of the median values
        uids <- unlist(unname(data[, uidCol]))
        normValVec <- lapply(uids, function(x) {
          datSlice <- data[data[, uidCol] == x, ]
          normVal <- median(unlist(datSlice[covariate]))
          if (normVal == 0) {
            if (!categorical && !all(unlist(datSlice[covariate]) == 0)) { # print warning only for non-categorical variables
              cli::cli_alert_warning("normalization value for subject ID: {x} is zero; using with 1...")
            }
            normVal <- 1
          }
          normVal
        })
        normValVec <- mean(unlist(normValVec))
      }
    }
    
    # check if log transform required
    isLog <-
      grepl("\\bexp *\\(", funstringSplit[[idx]], perl = TRUE)
    
    res <- performNorm(
      data = data,
      covariate = covariate,
      varName = varName,
      normOp = normOp,
      normValVec = normValVec,
      isLog = isLog,
      isCat = categorical,
      isHS = isHS
    )
    
    data <- res[[1]]
    covNameMod <- res[[2]]
    covNames <- res[[3]]
    
    # cli::cli_h1('theta: {names(fitobject$theta)}')
    # cli::cli_h1('previous value: {funstringSplit[[idx]]}')
    funstringSplit[[idx]] <-
      oldaddCovariate(
        funstringSplit[[idx]],
        varName,
        covNameMod,
        names(fitobject$theta),
        isLog
      )
  }
  
  cli::cli_alert_success("added {covNameMod} to {varName}'s equation in the model")
  cli::cli_alert_success("updated value: {funstringSplit[[idx]]}")
  
  updatedMod <- initializeCovars(
    fitobject,
    funstringSplit[[idx]],
    covNames,
    initialEst,
    initialEstLB,
    initialEstUB
  )
  
  list(updatedMod, data, covNames) # return updated model and updated data
}


#' Remove covariate from function string
#'
#' Function to remove covariates from a given variable's equation in the function string text
#'
#' @param fitobject an nlmixr2 'fit' object
#' @param varName a string giving the variable name to which covariate needs to be added
#' @param covariate a string giving the covariate name; must be present in the data used for 'fit'
#' @param categorical a boolean to represent if the covariate to be added is categorical
#' @param isHS a boolean to represent if the covariate to be added is hockey-stick normalized
#'
#' @return returns a list containing the updated model and the parameter names for the covariates added
#'
#' @author Vipul Mann, Matthew Fidler
#'
#' @export
removeCovVar <- function(fitobject,
                         varName,
                         covariate,
                         categorical = FALSE,
                         isHS = FALSE) {
  if (!inherits(fitobject, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  
  funstring <- fitobject$ui$fun.txt
  funstringSplit <-
    unlist(strsplit(funstring, split = "\\\n")) # split the string at \n
  
  # # look for the string that needs to be modified
  idx <-
    grep(paste0("(?<!\\w)", varName), funstringSplit, perl = TRUE) # varName not preceded by any other character
  
  if (length(idx) >= 2) {
    print(funstring)
    stop("cannot remove more than one covariate at a time")
  }
  
  if (!isHS && !categorical) {
    # not HS, not CAT
    covNameMod <-
      paste0(
        paste0("centered_", covariate),
        "*",
        paste0("cov_", covariate, "_", varName)
      )
    covNames <- paste0("cov_", covariate, "_", varName)
  }
  
  else if (isHS) {
    # HS
    prefix <- paste0("centered_", covariate, "_")
    prefix2 <- paste0("cov_", covariate, "_")
    s <- c("lower", "upper")
    covModExpr <-
      paste0(paste0(prefix, s), "*", paste0(prefix2, s, "_", varName))
    covNameMod <- paste(covModExpr, collapse = "+")
    
    covNames <- paste0(prefix2, s, "_", varName)
  }
  
  else if (categorical) {
    # CAT
    prefix <- paste0("categorical_", covariate, "_")
    prefix2 <- paste0("cov_", covariate, "_")
    
    v <- unlist(nlme::getData(fitobject)[, covariate])
    s <- head(sort(unique(v)), -1) # remove the last column
    
    covModExpr <-
      paste0(paste0(prefix, s), "*", paste0(prefix2, s, "_", varName))
    covNameMod <- paste(covModExpr, collapse = "+")
    
    covNames <- paste0(prefix2, s, "_", varName)
  }
  
  else {
    stop("aborting...unknown covariate type", call. = FALSE)
  }
  
  funstringSplit[[idx]] <-
    removeCovariate(
      funstringSplit[[idx]],
      varName,
      covNameMod,
      names(fitobject$theta)
    )
  
  
  cli::cli_alert_success("removed {covNameMod} from {varName}'s equation in the model")
  cli::cli_alert_success("updated function text: {funstringSplit[[idx]]}")
  
  updatedMod <- paste0("model(fitobject,{", funstringSplit[[idx]], "})")
  updatedMod <- eval(parse(text = updatedMod))
  
  list(updatedMod, covNames)
}

#' Perform normalization of the covariate
#'
#' @param data a dataframe consisting the covariates added
#' @param covariate a string giving the covariate name; must be present in the data used for 'fit'
#' @param varName the variable name to which the covariate is being added
#' @param normOp an operator indicating the kind transformation to be done on the covariate
#' @param normValVec a numeric value to be used for normalization of the covariate
#' @param isLog a boolean indicating the presence of log-transformation in the funstring; default is FALSE
#' @param isCat a boolean indicating if the covariate is categorical; default is FALSE
#' @param isHS a boolean indicating if the covariate is of Hockey-stick kind; default is FALSE
#'
#' @return a list comprising the update dataframe, the expression for covariate, and a list of covariate names
#' @export
#' @author Vipul Mann, Matthew Fidler
performNorm <- function(data,
                        covariate,
                        varName,
                        normOp,
                        normValVec,
                        isLog = FALSE,
                        isCat = FALSE,
                        isHS = FALSE) {
  if (!(isCat)) {
    # not categorical variable
    if (!(isHS)) {
      # not hockey stick
      datColNames <- paste0("centered_", covariate)
      
      
      if (length(normOp) > 1) {
        if (normValVec[[2]] == 0) {
          normValVec <- 1
          cli::cli_alert_warning("replacing the normalization value from 0 to 1")
        }
        
        data[, datColNames] <-
          normOp[[1]](unname(unlist(data[, covariate])), normValVec[[1]])
        data[, datColNames] <-
          normOp[[2]](unname(unlist(data[, datColNames])), normValVec[[2]])
      }
      else {
        if (normValVec == 0) {
          normValVec <- 1
          cli::cli_alert_warning("replacing the normalization value from 0 to 1")
        }
        
        data[, datColNames] <-
          normOp(unname(unlist(data[, covariate])), normValVec)
      }
      
      covNameMod1 <- datColNames
      covNameParam1 <- paste0("cov_", covariate, "_", varName)
      covNameMod <- paste0(covNameMod1, "*", covNameParam1)
      covNames <- covNameParam1
    }
    
    else {
      res <- makeHockeyStick(data, covariate = covariate, varName = varName)
      
      data <- res[[1]]
      covModExpr <- res[[2]]
      covNames <- res[[3]]
      datColNames <- res[[4]]
      covNameMod <- paste(covModExpr, collapse = "+")
      return(list(data, covNameMod, covNames))
    }
  }
  
  else {
    # categorical variable
    # datColNames <- paste0("categorical_", covariate)
    res <- makeDummies(data, covariate = covariate, varName = varName)
    data <- res[[1]]
    covModExpr <- res[[2]]
    covNames <- res[[3]]
    datColNames <- res[[4]]
    
    covNameMod <- paste(covModExpr, collapse = "+")
    
    return(list(data, covNameMod, covNames))
  }
  
  if (isLog) {
    if (varName %in% c("cl")) {
      # with 0.75 prefactor
      for (datColName in datColNames) {
        if (!all(is.finite(log(data[, datColName])))) {
          stop("non-finite values encountered in log-normalization. aborting...", call. = FALSE)
        }
        
        # for loop to handle both non-categorical and categorical vars
        data[, datColName] <- 0.75 * log(data[, datColName])
      }
    }
    else {
      for (datColName in datColNames) {
        if (!all(is.finite(log(data[, datColName])))) {
          stop("non-finite values encountered in log-normalization. aborting...", call. = FALSE)
        }
        
        data[, datColName] <- log(data[, datColName])
      }
    }
  }
  
  list(data, covNameMod, covNames)
}

#' Initializing covariates before estimation
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

#' Creating Hockey-stick covariates
#'
#' @param data a dataframe containing the dataset that needs to be used
#' @param covariate the covariate that needs to be converted to hockey-stick; must be present in the data
#' @param varName the variable name to which the given covariate is to be added
#'
#' @return a list of updated data with covariates added, an expression that needs to be added to the model expression, the list of covariate names, and the column names corresponding to the hockey-stick covariates
#' @export
#' @author Vipul Mann, Matthew Fidler
makeHockeyStick <- function(data, covariate, varName) {
  v <- unlist(data[, covariate])
  
  prefix <- paste0("centered_", covariate, "_")
  s <- c("lower", "upper")
  
  # create two columns for below and above the median
  med <- median(v)
  d <- list(v * 1L * (v < med), v * 1L * (v >= med))
  
  names(d) <- paste0(prefix, s)
  newdat <- cbind(data, d)
  
  prefix2 <- paste0("cov_", covariate, "_")
  covNames <- paste0(prefix2, s, "_", varName)
  
  covModExpr <- paste0(names(d), "*", covNames)
  list(newdat, covModExpr, covNames, colnames(d))
}


#' Create categorical covariates
#'
#' @param data a dataframe containing the dataset that needs to be used
#' @param covariate the covariate that needs to be converted to categorical; must be present in the data
#' @param varName the variable name to which the given covariate is to be added
#'
#' @return a list of updated data with covariates added, an expression that needs to be added to the model expression, the list of covariate names, and the column names corresponding to the categorical covariates
#' @export
#' @author Vipul Mann, Matthew Fidler
#'
makeDummies <- function(data, covariate, varName) {
  v <- unlist(data[, covariate])
  
  prefix <- paste0("categorical_", covariate, "_")
  s <- head(sort(unique(v)), -1) # remove the last column
  
  d <- outer(v, s, function(v, s) {
    1L * (v == s)
  })
  colnames(d) <- paste0(prefix, s)
  newdat <- cbind(data, d)
  
  prefix2 <- paste0("cov_", covariate, "_")
  covNames <- paste0(prefix2, s, "_", varName)
  
  covModExpr <- paste0(colnames(d), "*", covNames)
  list(newdat, covModExpr, covNames, colnames(d))
}


#' Removing multiple covariates
#'
#' @param covInfo a list containing information about each variable-covariate pair
#' @param fitobject an nlmixr2 'fit' object
#'
#' @return a list with the updated fit object, the variable-covariate pair string, and the parameter names for the corresponding covariates removed
#' @export
#'
#' @author Vipul Mann, Matthew Fidler
#'
removeCovMultiple <- function(covInfo, fitobject) {
  covSearchRes <- list() # list to store fitobjects during the search
  
  # removing multiple covariates (independently)
  .env <- environment()
  .env$covSearchRes <- list()
  lapply(1:length(covInfo), function(idx) {
    x <- covInfo[[idx]]
    
    res <- do.call(removeCovVar, c(list(fitobject), x))
    updatedMod <- res[[1]]
    data <- nlme::getData(fitobject)
    covNames <- res[[2]]
    
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
    
    # fit2 <-
    #   suppressWarnings(nlmixr2(updatedMod, data, est = getFitMethod(fitobject)))
    
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
    
    AIC(fit2)
    
    .env$covSearchRes[[idx]] <- list(fit2, c(x$covariate, x$varName), covNames)
  })
  
  covSearchRes <- .env$covSearchRes
}


#' Add multiple covariates to a given model, sequentially or all at once
#'
#' @param covInfo a list containing information about each variable-covariate pair
#' @param fitobject an nlmixr2 'fit' object
#' @param indep a boolean indicating if the covariates should be added independently, or sequentially (append to the previous model); default is TRUE
#'
#' @return A list of fitobject searched
#' @export
#' @keywords internal
#' @author Vipul Mann, Matthew Fidler
#'
addCovMultiple <- function(covInfo, fitobject, indep = TRUE) {
  covSearchRes <- list() # list to store fitobjects during the search
  
  # adding covariates independent of each other
  .env <- environment()
  .env$covSearchRes <- list()
  if (indep) {
    lapply(1:length(covInfo), function(idx) {
      x <- covInfo[[idx]]
      res <- do.call(addCovVar, c(list(fitobject), x))
      updatedMod <- res[[1]]
      data <- res[[2]]
      covNames <- res[[3]]
      
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
      
      # fit2 <-
      #   suppressWarnings(nlmixr2(updatedMod, data, est = getFitMethod(fitobject)))
      
      fit2 <- tryCatch(
        {
          fit2 <-
            suppressWarnings(nlmixr2(updatedMod, data, est = getFitMethod(fitobject)))
          fit2 # to return 'fit2'
        },
        error = function(error_message) {
          print("error fitting the model after ADDING the covariate")
          print(error_message)
          print("returning the previous fitobject")
          return(fitobject) # return NA otherwise (instead of NULL)
        }
      )
      
      AIC(fit2)
      
      .env$covSearchRes[[idx]] <- list(fit2, c(x$covariate, x$varName), covNames)
    })
    
    covSearchRes <- .env$covSearchRes
  }
  
  # adding covariates one after the other, appending to the previous model
  else {
    covsAdded <-
      list() # to keep track of covariates added and store in a file
    covsAddedIdx <- 1
    for (x in covInfo) {
      res <- do.call(addCovVar, c(list(fitobject), x))
      updatedMod <- res[[1]]
      data <- res[[2]]
      covNames <- res[[3]]
      
      if (length(covsAdded) == 0) {
        # covsAdded[[covsAddedIdx]] <- paste0(x$covariate, x$varName)
        covsAdded[[covsAddedIdx]] <- c(x$covariate, x$varName)
        
        # fit2 <-
        #   suppressWarnings(nlmixr2(updatedMod, data, est = getFitMethod(fitobject)))
        
        fit2 <- tryCatch(
          {
            fit2 <-
              suppressWarnings(nlmixr2(updatedMod, data, est = getFitMethod(fitobject)))
            fit2 # to return 'fit2'
          },
          error = function(error_message) {
            print("error fitting the model while SIMULTANEOUSLY ADDING covariates")
            print(error_message)
            print("skipping this covariate")
            return(fitobject) # return NA otherwise (instead of NULL)
          }
        )
        
        covSearchRes[[covsAddedIdx]] <- list(fit2, covsAdded[[covsAddedIdx]], covNames)
        covsAddedIdx <- covsAddedIdx + 1
      }
      else {
        # covsAdded[[covsAddedIdx]] <-
        #   paste0(covsAdded[[covsAddedIdx - 1]], "_", x$covariate, x$varName)
        covsAdded[[covsAddedIdx]] <-
          c(covsAdded[[covsAddedIdx - 1]], x$covariate, x$varName)
        
        # fit2 <-
        #   suppressWarnings(nlmixr2(updatedMod, data, est = getFitMethod(fitobject)))
        
        fit2 <- tryCatch(
          {
            fit2 <-
              suppressWarnings(nlmixr2(updatedMod, data, est = getFitMethod(fitobject)))
            fit2 # to return 'fit2'
          },
          error = function(error_message) {
            print("error fitting the model while SIMULTANEOUSLY ADDING covariates")
            print(error_message)
            print("skipping this covariate")
            return(fitobject) # return NA otherwise (instead of NULL)
          }
        )
        
        covSearchRes[[covsAddedIdx]] <- list(fit2, covsAdded[[covsAddedIdx]], covNames)
        covsAddedIdx <- covsAddedIdx + 1
      }
      
      print(fit2$fun.txt)
      
      fitobject <- fit2
    }
  }
  
  covSearchRes
}





