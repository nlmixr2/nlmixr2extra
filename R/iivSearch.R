#' Automated Inter-Individual Variability Search
#' @param fit a model fit
#' @param sortBy either "BIC" or "OBJ". Default is "BIC"
#' @param mceta an integer specifying mceta during fitting.
#' @param ... parameters for sub functions
#' 
#' @author Omar I. Elashkar
#' @export 
iivSearch <- function(fit, ...){
    UseMethod("iivSearch")
}


#'@rdname iivSearch
#'@export 
iivSearch.default <- function(fit, ...){
    stop("Currently supports only linearized models")
}


#'@rdname iivSearch
#'@export 
iivSearch.nlmixr2Linearize <- function(fit, sortBy = "BIC", mceta=5, ...){
    checkmate::assertChoice(sortBy, c("OBJ", "BIC"))
    checkmate::assertInteger(mceta, lower = -1, upper = 1000)
  
    if(!hasAnyEta(fit)){stop("Use `addEtas=TRUE` while linearizing this model")}
    # get eta names
    etaAll <- fit$iniDf[!is.na(fit$iniDf$neta1), ]
    etaNoCorr <- etaAll[etaAll$neta1 == etaAll$neta2,]
    etaNames <- etaNoCorr$name

    iivSpace <- iivCombn(etaNames)

    varCovMat <- cov(fit$eta[,-1])
    res <- lapply(cli::cli_progress_along(iivSpace, "IIV Search"), function(x){
        x <- iivSpace[x]
        message(x)
        omegaMat <- filterEtaMat(varCovMat, x)
        newMod <- fit %>% rxode2::ini(omegaMat) 
        noCorrSpace <- unlist(strsplit(x, "\\+"))
        noCorrSpace <- grep("~", noCorrSpace, invert = TRUE, value = TRUE)
        etaToRemove <- setdiff(etaNames, noCorrSpace)
        # for(i in etaToRemove){ # deriviatives are to be kept not removed
        #     newMod <- eval(str2lang(paste0("model(newMod, base_", i, " = 0)")))
        #     newMod <- eval(str2lang(paste0("model(newMod, err_", i, " = 0)")))
        # }
        unfixedIni <- newMod$iniDf 
        unfixedIni$fix[!is.na(unfixedIni$neta1) & !(unfixedIni$name %in% etaToRemove)] <- FALSE
        unfixedIni$upper[!is.na(unfixedIni$ntheta) &  unfixedIni$fix == FALSE] <-  unfixedIni$est[!is.na(unfixedIni$ntheta) &  unfixedIni$fix == FALSE]*2  # set upper to residual to 2x
        rxode2::ini(newMod) <- unfixedIni
        tryCatch({
            fit <- nlmixr(newMod, nlme::getData(fit), est = "focei",
                control = nlmixr2est::foceiControl(mceta=mceta, 
                                                    print = 10,
                                                    maxInnerIterations = 9999,
                                                    etaMat = as.matrix(fit$eta[,c(noCorrSpace)])))
          
          summ <- fit$objDf[, c("OBJF", "AIC", "BIC")]
          summ$search <- x
          summ$nParams <- nrow(fit$iniDf[fit$iniDf$est != 0,]) # including res + etas
          summ$covMethod  <- fit$covMethod
          summ$outerOptTxt <- fit$control$outerOptTxt
          
          list(fit = fit, objDf = summ, search = x)
        }, error = function(y){
            return(list(fit = NULL, objDf = NA, search = x))
          }
        )
    })

    objDfAll <- do.call(rbind, lapply(res, function(x) x$objDf))

    if(any(duplicated(objDfAll$BIC))){
        warning("Multiple models have the same BIC. Consider using a different search space.")
    }
    

    finalRes <- list(res = res, summary = objDfAll, linearFit = fit)
    class(finalRes) <- "linIIVSearch"
    finalRes
}

#' Filter covariance matrix
#' @param oMat covariance matrix 
#' @param filterStr single string with pattern of "+" between entities and "~" for within entity correlation
#' @author Omar I. Elashkar
#' @noRd 
filterEtaMat <- function(oMat, filterStr, minIni=0.1){
        stopifnot(length(filterStr) == 1)
        stopifnot(isSymmetric(oMat))

        resMat <- oMat
        resMat[] <- 0 
        elem <- unlist(strsplit(filterStr, "\\+"))
        for(i in elem){
            elem2 <- unlist(strsplit(i, "\\~"))
            x <- elem2[1]
            y <- ifelse(length(elem2) == 1, elem2[1], elem2[2])
            resMat[x,y] <- oMat[x, y]
            resMat[y, x] <- oMat[y, x]  # Ensure symmetry

        }
        resMat[resMat < minIni & resMat !=0] <- minIni

        stopifnot(isSymmetric(resMat))
        resMat <- lotri::lotriNearPD(resMat)
        
        etaToRemove <-  setdiff(colnames(oMat), elem)
        if(length(etaToRemove > 0)){
          resMat[etaToRemove, ] <- 0
          resMat[,etaToRemove ] <- 0
          resMat[etaToRemove, etaToRemove ] <- 0
        }

        resMat
    }

#' Check if a model has any individual variability (eta) parameters
#' @return TRUE if any eta parameters are present, FALSE otherwise
#' @author Omar I. Elashkar
#' @noRd
hasAnyEta <- function(ui){
    etadf <- ui$iniDf[!is.na(ui$iniDf$neta1),]
    # thetadf <- ui$iniDf[is.na(ui$iniDf$neta1) & is.na(ui$ini$condition),]
    
    nrow(etadf) > 0
}

#' Check if a model has unfixed eta parameters
#' @param ui rxode2 ui or fit
#' @return TRUE if any eta parameters are unfixed, FALSE otherwise
#' This function is used to ensure that the model is suitable for linearization intended for IIV search.
#' @author Omar I. Elashkar
#' @noRd
hasUnFixedEta <- function(ui){
    etadf <- ui$iniDf[!is.na(ui$iniDf$neta1) & !ui$iniDf$fix,]
    nrow(etadf) > 0
}


#' Generate all combinations of input parameters
#' @param params character vector of parameter names
#' @author Omar I. Elashkar
#' @return char vector of possible combinations separated by "+". Correlations are marked by "~".
#' @noRd 
iivCombn <- function(params) {
    # Generate all combinations of input parameters
    all_combinations <- unlist(lapply(1:length(params), function(i) {
    combn(params, i, FUN = function(x) paste(x, collapse = "+"), simplify = TRUE)
    }))

    # Generate results with/without all possible corr per all available etas
    all_results <- unlist(lapply(all_combinations, function(combo) {
        individual <- strsplit(combo, "\\+")[[1]]
        if (length(individual) == 1) {
            return(combo)
        }
        # Generate all pairwise correlations
        all_corr <- combn(individual, 2, FUN = function(x) paste(x, collapse = "~"), simplify = TRUE)
        
        # Generate subsets of correlations
        corr_subsets <- unlist(lapply(0:length(all_corr), function(k) {
            combn(all_corr, k, FUN = function(x) paste(x, collapse = "+"), simplify = TRUE)
        }))
        
        # Combine the base combination with each correlation subset
        paste0(combo, "+", corr_subsets)
    }))

    all_results
}

#' Add Individual Random Effects and Fix them to Small Value
#' @param ui model with no individual random effects added
#' @param fix If TRUE, the added etas will be fixed (not estimatable). Default is FALSE.
#' @return nlmixr2 fit with individual random effects added on all fixed effects and fixed.
#' @author Omar I. Elashkar
#' @export 
addAllEtas <- function(ui, fix = FALSE){ 
    if(inherits(ui, "nlmixr2FitCore")){ui <- ui$ui}
    if(inherits(ui, "function")){ui <- ui()}
    
    if(!is.null(ui$noMuEtas)){stop("Some Etas are not mu referenced")}
    if(!inherits(ui, "rxUi")){stop("`ui` must be a model or fit")}

    # ui$singleTheta 
    # ui$muRefDataFrame

    thetaNames <- ui$iniDf[is.na(ui$iniDf$neta1) & is.na(ui$ini$condition), "name"]
    freeTheta <- setdiff(thetaNames, ui$muRefDataFrame$theta)
    if(length(freeTheta) > 0){
      newMod <- nlmixr2lib::addEta(ui, freeTheta)
    } else{
      newMod <- ui
    }

    # fix all random effects to 10^-3
    iniDf <- newMod$iniDf
    iniDf$fix[!is.na(iniDf$neta1)] <- fix
    iniDf$est[!is.na(iniDf$neta1)] <- 0.001
    rxode2::ini(newMod) <- iniDf

    newMod

}

#' Print Summary Table For Linearized IIV Search
#' @param x IIV search results object
#' @param ... Other arguments passed to print
#' @author Omar I. Elashkar
#' @export
print.linIIVSearch <- function(x, ...){
  df <- x$summary
  df[order(df$BIC), ]
}

#' Rerun Top N Original Models From A Search
#' @param x a Search object
#' @param n number of models to rerun
#' @param ... Other Parameters
#' @author Omar I. Elashkar
#' @export
rerunTopN <- function(x, ...){
  UseMethod("rerunTopN")
}
  
#' @export
rerunTopN.default <- function(x, ...){
  stop("Method not supported")
}


#' @rdname rerunTopN
#' @export
rerunTopN.linIIVSearch <- function(x, n = 5, ...){
  
  oObj <- x$linearFit$env$originalFit$objDf$OBJF
  
  passed <- x$summary 
  passed <- passed[passed$OBJF >= oObj/2,] # more than that, likely combination
  
  topCand <- passed[order(passed$BIC), ]$search[1:n]
  
  cli::cli_alert_info("Starting refitting original models")
  cli::cli_alert_info("Starting refit: {topCand}")
  stopifnot(!any(is.na(topCand)))
  
  origData <- nlme::getData(x$linearFit$env$originalFit)
  nlui <- x$linearFit$env$originalFit$finalUi
  nlui <- addAllEtas(nlui)
  resOrig <- lapply(cli::cli_progress_along(topCand), function(i){
    i <- topCand[i]
    cli::cli_alert_info("Fitting Original Model with {x}")
    omegaMat <- filterEtaMat(cov(x$linearFit$eta[,-1]), i)
    iniDf <- nlui$iniDf
    iniDf[!is.na(iniDf$neta1), "fix"] <- FALSE 
    rxode2::ini(nlui) <- iniDf
    nlui <- nlui %>% rxode2::ini(omegaMat)
    nlmixr(nlui, origData, est="focei", control = nlmixr2est::foceiControl())
  })
  
  resOrigSummary <- do.call(rbind, lapply(resOrig, function(x) x$objDf))
  colnames(resOrigSummary) <- paste0("O.", colnames(resOrigSummary))
  resOrigSummary$search <- topCand
  
  list(summary = resOrigSummary, fits = resOrig)
}
