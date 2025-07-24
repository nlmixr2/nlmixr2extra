#' Automated Inter-Individual Variability Search
#' @param fit a model fit
#' 
#' @author Omar Elashkar
#' @export 
iivSearch <- function(fit){
    UseMethod("iivSearch")
}

#'@export 
iivSearch.default <- function(fit){
    stop("Currently supports only linearized models")
}


#'@export 
iivSearch.nlmixr2Linearize <- function(fit, sortBy = "BIC"){
    if(hasUnFixedEta(fit)){stop("This model has unfixed IIV and not suitable for this procedure")}
    # get eta names
    etaAll <- fit$iniDf[!is.na(fit$iniDf$neta1), ]
    etaNoCorr <- etaAll[etaAll$neta1 == etaAll$neta2,]
    etaNames <- etaNoCorr$name

    iivSpace <- iivCombn(etaNames)

    varCovMat <- cov(fit$eta[,-1])
    res <- lapply(iivSpace, function(x){
        omegaMat <- filterEtaMat(varCovMat, x)
        newMod <- fit %>% ini(omegaMat) 
        noCorrSpace <- unlist(strsplit(x, "-"))
        noCorrSpace <- grep(",", noCorrSpace, invert = TRUE, value = TRUE)
        etaToRemove <- setdiff(etaNames, noCorrSpace)
        for(i in etaToRemove){
            newMod <- eval(str2lang(paste0("model(newMod, base_", i, " = 0)")))
            newMod <- eval(str2lang(paste0("model(newMod, err_", i, " = 0)")))
        }
        unfixedIni <- newMod$iniDf 
        unfixedIni$fix[!is.na(unfixedIni$neta1)] <- FALSE
        ini(newMod) <- unfixedIni

        fit <- nlmixr(newMod, nlme::getData(fit), est = "focei")
        summ <- fit$objDf[, c("OBJF", "AIC", "BIC")]
        summ$search <- x
        summ$nParams <- nrow(fit$iniDf[fit$iniDf$est != 0,])
        list(fit = fit, objDf = summ, search = x)
    })

    objDfAll <- do.call(rbind, lapply(res, function(x) x$objDf))

    if(any(duplicated(objDfAll$BIC))){
        warning("Multiple models have the same BIC. Consider using a different search space.")
    }
    finalVarCov <- filterEtaMat(varCovMat, iivSpace[[which.min(objDfAll[[sortBy]]) & which.min(objDfAll$nParams)]])
    finalOFit <- fit$env$originalUi |> ini(finalVarCov)
    finalFit <- nlmixr(fit$env$originalUi, nlme::getData(fit), est = "focei", 
            control = nlmixr2est::foceiControl(mceta=10, covMethod = "", calcTables=FALSE))


    list(res = res, summary = objDfAll, finalFit = finalFit)
}

#' Filter covariance matrix
#' @param oMat covariance matrix 
#' @param filterStr single string with pattern of "-" between entities and "," for within entity correlation
#' @author Omar Elashkar
#' @noRd 
filterEtaMat <- function(oMat, filterStr, minIni=0.1){
        stopifnot(length(filterStr) == 1)
        resMat <- oMat
        resMat[] <- 0
        elem <- unlist(strsplit(filterStr, "-"))
        for(i in elem){
            elem2 <- unlist(strsplit(i, ","))
            x <- elem2[1]
            y <- ifelse(length(elem2) == 1, elem2[1], elem2[2])
            resMat[x,y] <- oMat[x, y]
            resMat[y, x] <- oMat[y, x]  # Ensure symmetry
        }
        resMat[resMat < minIni & resMat !=0] <- minIni

        stopifnot(isSymmetric(resMat))

        resMat
    }

#' Check if a model has any individual variability (eta) parameters
#' @return TRUE if any eta parameters are present, FALSE otherwise
#' @author Omar Elashkar
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
#' @author Omar Elashkar
#' @noRd
hasUnFixedEta <- function(ui){
    etadf <- ui$iniDf[!is.na(ui$iniDf$neta1) & !ui$iniDf$fix,]
    # thetadf <- ui$iniDf[is.na(ui$iniDf$neta1) & is.na(ui$iniDf$condition),]

    nrow(etadf) > 0
}


#' Generate all combinations of input parameters
#' @param params character vector of parameter names
#' @author Omar Elashkar
#' @return char vector of possible combinations separated by "-". Correlations are marked by ",".
#' @noRd 
iivCombn <- function(params) {
    # Generate all combinations of input parameters
    all_combinations <- unlist(lapply(1:length(params), function(i) {
    combn(params, i, FUN = function(x) paste(x, collapse = "-"), simplify = TRUE)
    }))

    # Generate results with/without all possible corr per all available etas
    all_results <- unlist(lapply(all_combinations, function(combo) {
        individual <- strsplit(combo, "-")[[1]]
        if (length(individual) == 1) {
            return(combo)
        }
        # Generate all pairwise correlations
        all_corr <- combn(individual, 2, FUN = function(x) paste(x, collapse = ","), simplify = TRUE)
        
        # Generate subsets of correlations
        corr_subsets <- unlist(lapply(0:length(all_corr), function(k) {
            combn(all_corr, k, FUN = function(x) paste(x, collapse = "-"), simplify = TRUE)
        }))
        
        # Combine the base combination with each correlation subset
        paste0(combo, "-", corr_subsets)
    }))

    all_results
}

#' Add Individual Random Effects and Fix them to Small Value
#' @param ui model with no individual random effects added
#' @return nlmixr2 fit with individual random effects added on all fixed effects and fixed.
#' @author Omar Elashkar
#' @export 
addFixedEtas <- function(ui){ 
    if(inherits(ui, "nlmixr2FitCore")){ui <- ui$ui}
    if(inherits(ui, "function")){ui <- ui()}
    
    if(hasAnyEta(ui)){stop("Model already has IIV parameters. Must provide a model without etas.")}
    if(!inherits(ui, "rxUi")){stop("`ui` must be a model or fit")}

    # ui$singleTheta 
    # ui$muRefDataFrame
    # ui$noMuEtas

    thetaNames <- ui$iniDf[is.na(ui$iniDf$neta1) & is.na(ui$ini$condition), "name"]
    # add random effects 
    newMod <- nlmixr2lib::addEta(ui, thetaNames)

    # fix all random effects to 10^-4
    iniDf <- newMod$iniDf
    iniDf$fix[!is.na(iniDf$neta1)] <- TRUE
    iniDf$est[!is.na(iniDf$neta1)] <- 1e-4
    ini(newMod) <- iniDf

    newMod

}
