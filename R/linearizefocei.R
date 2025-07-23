#' Extract derivatives from NLME model fitted with FOCE
#'
#' @param fit object fitted with nlmixr2 of class nlmixr2.focei
#'
#' @author Omar Elashkar
#' @noRd
getDeriv <- function(fit){
    # if(fit$method != "FOCE") stop("This method requires FOCE method")
    rxode2::assertRxUiMixedOnly(fit, " for the procedure routine 'linerize'", .var.name=fit$modelName)

    ui <- rxode2::assertRxUi(fit)
    innerModel <- ui$foceiModel

    eta <- fit$eta

    eta <- eta[,-1]
    eta <- setNames(eta, paste0("ETA[", seq_along(eta), "]"))

    theta <- fit$theta
    reps <- ceiling(nrow(eta) / nrow(theta))
    # theta <- as.data.frame(theta)
    theta <- setNames(fit$theta, paste0("THETA[", seq_along(fit$theta), "]")) |>
        t() |>
        as.data.frame()

    theta <- theta[rep(1, nrow(eta)),,drop=FALSE]
    params_df <- cbind(theta, eta)

    oData <- nlme::getData(fit)
    names(oData) <- toupper(names(oData))

    
    stopifnot(all(c("rx_pred_", "rx_r_") %in% innerModel$inner$lhs))
    stopifnot(sum(grepl("rx__sens_rx_pred__BY_ETA_\\d+___", innerModel$inner$lhs)) == ncol(eta))
    stopifnot(sum(grepl("rx__sens_rx_r__BY_ETA_\\d+___", innerModel$inner$lhs)) >= ncol(eta))

    derv <- rxode2::rxSolve(innerModel$innerOeta, oData, params=params_df,
        addDosing=FALSE,
        keep = c("DV", setdiff(names(oData),  c("DV", "ID", "TIME", "DVID", "ADDL", "EVID", "AMT", "CMT")))) # add covariates
    
    eta <- fit$eta[,-1]
    # OPRED
    derv <- renameCol(derv, "OPRED", "rx_pred_")

    # D_ETA
    for(i in seq_len(ncol(eta))){
        currEta <- colnames(eta)[i]
        derv <- renameCol(derv, paste0("D_", currEta), paste0("rx__sens_rx_pred__BY_ETA_", i, "___"))
    }

    #O_ETA
    for(i in seq_len(ncol(eta))){
        currEta <- colnames(eta)[i]
        derv <- renameCol(derv, paste0("O_", currEta), paste0("rx__ETA", i))
    }

    #O_EPS
    derv <- renameCol(derv, "O_ResVar", "rx_r_")

    # O_IRES 
    derv$O_IRES <- fit$IRES

    # D_EPS
    derv$D_ResVar <- 1

    # D_EPSETA
    for(i in seq_len(ncol(eta))){
        currEta <- colnames(eta)[i]
        derv <- renameCol(derv, paste0("D_VAR_",currEta), paste0("rx__sens_rx_r__BY_ETA_", i, "___"))
    }

    # add OCMT col
    if(length(unique(fit$predDf$cond)) > 1){
        # derv$OCMT <- as.integer(derv$CMT)

        mv <- rxode2::rxModelVars(innerModel$inner) # model variables
        cmtName <- c(mv$state, mv$stateExtra) # + cmt names (states)
        cmtDf <- data.frame(CMT=seq_along(cmtName), cond=cmtName)
        predDf <- fit$predDf
        predDf$OCMT <- predDf$cmt
        cmtDf <- merge(cmtDf, predDf)
        cmtDf <- cmtDf[,c("OCMT", "CMT", "dvid")]
        
        derv$rxRow <- seq_along(derv$id) # dummy sort
        derv <- merge(derv, cmtDf) 
        derv <- derv[order(derv$rxRow), ]
        derv$rxRow <- NULL
        derv$CMT <- NULL
 
    }

    derv
}

renameCol <- function(df, new, old){
    names(df)[names(df) == old] <- new
    df
}

#' Generate a Linearization Model From Previous Fit 
#' @param fit a nonlinear model.
#' @param focei boolean. If TRUE, use FOCEI linearization with individual and residual linearization. Default is TRUE.
#' @param derivFct boolean. If TRUE, use normalization derivative factors. Default is FALSE.
#' @author Omar Elashkar
#' @noRd
linModGen <- function(fit, focei = TRUE, derivFct = FALSE){
    rxode2::assertRxUiMixedOnly(fit, " for the procedure routine 'linerize'", .var.name=fit$modelName)

    # errSym <- rxode2:::.rxGetVarianceForErrorType(rxUiDecompress(fit), fit$predDf)
    fitui <- rxode2::assertRxUi(fit)
    fitui <- rxode2::rxUiDecompress(fitui)
    assign("derivFct", derivFct, envir = fitui)
    extractError <- fitui$linearizeError
    # fitui <- rxode2::rxUiCompress(fitui)
    rm(fitui)
    epStr <- fit$predDf$var
    errNames <- fit$iniDf$name[!is.na(fit$iniDf$err)]
    nlmod <- fit$finalUi
    etaDf <- fit$eta[,-1]

    modelStr <- list()

    # mu block
    for(i in seq_along(colnames(etaDf))){
        currEta <- colnames(etaDf)[i]
        modelStr$muRef[i] <- paste0("mu_", currEta , " = " , "theta.", currEta , " + ", currEta)
    }
    # linearize eta
    for(i in seq_along(colnames(etaDf))){
        # D_ETAn * (-OETAn + eta.n)
        currEta <- colnames(etaDf)[i]
        modelStr$baseEta[i] <- paste0("base_",currEta , "=", "D_", currEta, "*(", "-O_", currEta, "+ mu_", currEta, ")")
    }
    # BASE_TERMS = base1 + ... + basen
    modelStr$baseEta[i+1] <- paste("BASE_TERMS =", paste(paste0("base_", colnames(etaDf)), collapse = " + "))
    modelStr$baseEta[i+2] <- paste0("y = BASE_TERMS + OPRED")
    # modelStr$baseEta[i+3] <- paste0("eps = y - DV")

    # linearize residuals
    for(i in seq_along(colnames(etaDf))){
        currEta <- colnames(etaDf)[i]
        # ERRn = D_VAR_ETA_1_n*(-O_ETAn + eta.n)
        modelStr$baseEps[i] <- paste0("err_",currEta , "=", "D_VAR_", currEta, "*(", "-O_", currEta, "+", currEta, ")")
    }

    errSym <- extractError$rxR2

    # for single ep, this will be single line
    # rxR2
    for(ii in seq_along(errSym)){
        modelStr$baseEps[i+ii] <- errSym[ii]
    }

    modelStr$baseEps[i+ii+1] <- "r <- sqrt(rxR2)"
    modelStr$baseEps[i+ii+2] <- paste0("foceiLin <- ", ifelse(focei, 1, 0))
    modelStr$baseEps[i+ii+3] <- paste("BASE_ERROR = foceiLin * fct * (", paste(paste0("err_", colnames(etaDf)), 
        collapse = " + "), ")/(2*r) + r")
    modelStr$baseEps[i+ii+4] <- "rxR = BASE_ERROR"

    for(i in seq_along(extractError$err)){
        modelStr$basePred[i] <-   extractError$err[i]
    }

    # iniDf already captures final estimates
    iniDf <- nlmod$iniDf[nlmod$iniDf$name %in% errNames | !is.na(nlmod$iniDf$neta1), ]

    # n_theta == n_errors + n_eta
    iniDf$ntheta[!is.na(iniDf$ntheta)] <- seq_along(na.omit(iniDf$ntheta)) # remove thetas except for err
    # add theta.etaname 
    thetaetadf <- data.frame(
        ntheta = seq_along(colnames(etaDf)) + max(iniDf$ntheta, na.rm = TRUE),
        neta1 = NA_real_,
        neta2 = NA_real_,
        name = paste0("theta.", colnames(etaDf)),
        lower = 0,
        est = 0,
        upper = Inf,
        fix = TRUE, 
        label = NA_character_,
        backTransform = NA_character_,
        condition = NA_character_,
        err = NA_character_
    )

    iniDf  <- rbind(iniDf, thetaetadf)
    iniDf  <- iniDf[order(iniDf$ntheta), ]

    nlmod <- rxode2::rxUiDecompress(nlmod)
    assign("iniDf", iniDf, envir = nlmod)
    v <- c(modelStr$muRef, modelStr$baseEta, modelStr$baseEps, 
                        extractError$tipred, modelStr$basePred) 

    nlmod$lstExpr <- as.list(str2lang(paste0("{", paste(v, collapse="\n"), "}"))[-1])
    nlmod <- nlmod$fun
    nlmod <- nlmod()

    nlmod
}

#' Perform linearization of a model fitted using FOCEI
#' @param fit fit of nonlinear model.
#' @param mceta a numeric vector for mceta to try. See decription.
#' @param relTol relative deviation tolerance between original and linearized models objective functions. Used for switching if focei = NA. See details.
#' @param focei Default is NA for automatic switch from FOCEI to FOCE if failed. See details.
#' @param derivFct boolean. If TRUE, turn on derivatives for linearization. Default is FALSE.
#' @param plot boolean. Print plot of linearized vs original
#' 
#' @details
#' 
#' mceta vector will be iterated over to find the best linearization if linearization failed. 
#' Escalating to next mceta will depend on the relative deviation `relTol` of the original and linearized models objective functions.
#' 
#' If `focei` is set to `NA`, the function will first try to linearize using FOCEI. 
#' If the relative deviation between original and linearized models objective functions is greater than `relTol`, it will switch to FOCE where residual linearization is skipped.
#' If `focei` is set to `TRUE`, the function will use FOCEI linearization with individual and residual linearization.
#' If `focei` is set to `FALSE`, the function will use FOCE linearization with residual linearization skipped.
#' 
#' If `derivFct` is set to `TRUE`, the function will use derivatives for linearization. 
#' This might be usefull to try with FOCEI.
#' 
#' `plot` argument can only print ggplot figure with default settings. 
#' If a user wish to capture the plot, one might use `linearizePlot()` call. 
#' 
#' 
#' 
#' @author Omar Elashkar
#' @export 
linearize <- function(fit, mceta=c(-1, 10, 100, 1000), relTol=0.4, focei = NA, derivFct = FALSE, plot = FALSE){ 
    checkmate::assertIntegerish(mceta, lower = -1, upper = 2000,  unique = TRUE)
    checkmate::assertNumeric(relTol, lower=0, upper=0.6)
    checkmate::assertLogical(plot)
    checkmate::assertLogical(derivFct)
    checkmate::assertLogical(focei, max.len = 1, any.missing=TRUE)

    derv <- getDeriv(fit)
    linMod <- linModGen(fit, focei = ifelse(is.na(focei), TRUE, focei), derivFct = derivFct)

    evalFun <- function(){
        # check linearization feasbility
        fitL <- evalLinModel(fit, linMod, derv)
        fitL_map <- evalLinModel(fit, linMod, derv, 1000)

        
        lObj_exact <- fitL$objDf$OBJF
        lObj_map <- fitL_map$objDf$OBJF
        oObj <- fit$objDf$OBJF
        
        list(
            lObj_exact = lObj_exact,
            lObj_map = lObj_map,
            oObj = oObj,
            message = paste("Non-Linear OFV: ", oObj, 
            "\n Eval Linear OFV (Exact):", lObj_exact, 
            "\n Eval Linear OFV (MAP):", lObj_map)
        )
    }

    firstEval <- evalFun()

    if(isTRUE(all.equal(firstEval$oObj, firstEval$lObj_map, tolerance = relTol))){
        message("Linearization evaluation matched. Linearization might be feasible ...")
    } else{
        message("Linearization evaluation mismatched by deltaOFV > 2%. Linearization might be difficult")
        message("Switching to linearization around predictions only (Variance linearization skipped)")
        linMod <- linMod |> model(foceiLin <- 0) 
        secondEval <- evalFun()
    }


    # fit linearized model
    for(i in seq_along(mceta)){
        fitL <- nlmixr(linMod, derv, est="focei",
            control = nlmixr2est::foceiControl(etaMat = fit, mceta=mceta[i], covMethod = "", calcTables=FALSE))
        oObj <- fit$objDf$OBJF
        lObj <- fitL$objDf$OBJF

        relDev <- abs((oObj-lObj)/lObj)

        if(relDev <=relTol){
            break
        } else{
            message(paste("Non-Linear OFV: ", oObj, "Linear OFV:", lObj, "mceta:", mceta[i], "Relative OFV Dev", relDev*100, "%"))
            if(i != length(mceta)){
                message(paste("Using mceta = ", mceta[i], "provided inadequate linearization. Trying mceta = ", mceta[i+1], "..."))
            } else{
                warning("Linearization was inadequate for the given tolerance. Try increasing the tolerance, refine the model or fit with mceta")
                # TODO automatically switch to FOCE if not already
            }
        }
    }

    fitL <- nlmixr2est::addTable(fitL)
    nlme::getVarCov(fitL)

    if (exists("secondEval")) {
        finalEval <- secondEval
    } else {
        finalEval <- firstEval
    }
    m <- paste(
        "Linearization method: ", ifelse(is.na(focei), 
            "Auto", ifelse(focei, "FOCE (individual + residual)", "Variance Skipped")), "\n",
        "Linearization Relative Tolerance: ", relTol, "\n",
        ifelse(is.na(focei) & exists("secondEval"), 
            "Linearization method switched automatically to FOCE (residual linearization skipped)", ""),
        "\n", paste(finalEval$message, collapse = "\n"),
        "\nFitted Linear OFV: ", lObj,
        "\nRelative OFV Dev: ", round(relDev, 4) * 100, "%\n",
        "\nmceta: ", mceta[i],
        "\nAll Eta Estimated:", hasAllEta(fit)
    )
    message("Linearization Summary:")
    message(m)

    tmpEnv <- fitL$env
    assign("message", c(fitL$message, m), envir=tmpEnv)
    v <- c("nlmixr2Linearize", class(fitL))
    attr(v, ".foceiEnv") <- tmpEnv
    class(fitL) <- v
    
    if(plot){
        print(linearizePlot(fit, fitL))
    }
    
    fitL
}

#' Plot Original Versus Linear Models iObj and Etas
#' @param nl non-linear fitting object
#' @param l linear fitting object
#' @author Omar Elashkar
#' @return ggplot object
#' @export
linearizePlot <- function(nl, lin){
    # TODO assertion that both objects are siblings. Random effects and errors names same

    stopifnot(inherits(nl, "nlmixr2FitCore"))
    stopifnot(inherits(lin, "nlmixr2Linearize"))

    l <- function(x, descr) {
        x$factor <- descr
        x
    }
    originalIval <- l(nl$etaObf, "original")
    linearIval <- l(lin$etaObf, "linear")

    fig <- rbind(originalIval, linearIval) |> 
        tidyr::pivot_longer(cols = c(-c("ID", "factor")), names_to = "parameter", 
        values_to = "value") |> 
        tidyr::pivot_wider(names_from = "factor", values_from = "value") |> 
        ggplot2::ggplot(aes(x = .data[["original"]], y=.data[["linear"]])) +
        ggplot2::geom_smooth(se = FALSE, method = "lm") +
        ggplot2::geom_point() + 
        ggplot2::facet_wrap("parameter", scales = "free") 
    fig
}

#' Evaluate A Linear Model Without Estimation 
#' @param fit fit of nonlinear model.
#' @param linMod linear model to evaluate
#' @param derv data frame of derivatives
#' @param innerIter number of inner iterations to use. Default is 0. 
#' 
#' `evalLinModel()` evaluates a linear model without estimation. 
#' This is useful for checking the linearization feasibility under given non-linear model and its fit. 
#' `innerIter=0` is equivalent to NONMEM `$ESTIMATION MAXITER=0` and will not perform any inner iterations.
#' 
#' @author Omar Elashkar
#' @noRd 
evalLinModel <- function(fit, linMod, derv, innerIter = 0, covMethod = ""){
    fitL <- nlmixr(linMod, derv, est="focei", 
            control = nlmixr2est::foceiControl(etaMat = fit, mceta=-1, 
            covMethod = covMethod, calcTables=TRUE, 
            maxInnerIterations=innerIter, maxOuterIterations=0L))
    fitL
}

#' Check Linearization Match 
#' @param fit fit of nonlinear model.
#' @param linFit fit of linear model.
#' @param tol relative tolerance for matching. Default is 0.05.
#' @author Omar Elashkar
#' @noRd
isLinearizeMatch <- function(fit, linFit, tol = 0.05){
    
    est <-setNames(fit$iniDf$est, fit$iniDf$name)
    estLin <- setNames(linFit$iniDf$est, linFit$iniDf$name)
    est <- est[names(estLin)]
    
    deltaOfv <- all.equal(fit$objDf$OBJF, linFit$objDf$OBJF, tol = tol)
    deltaOmega <- all.equal(fit$omega, linFit$omega, tol = tol)
    deltaEta <- all.equal(fit$eta, linFit$eta, tol = tol)
    deltaErr <- all.equal(est, estLin, tol = tol)

    list(
        ofv = list(isTRUE(deltaOfv), deltaOfv),
        omega = list(isTRUE(deltaOmega), deltaOmega),
        eta = list(isTRUE(deltaEta), deltaEta),
        err = list(isTRUE(deltaErr), deltaErr)
    )
}


.printInnerExp <- function(mod){
    ui <- rxode2::assertRxUi(mod)
    innerModel <- ui$foceiModel
    summary(innerModel$inner)
}

#'Parse covariate information from expression and data
#'@param expr formula
#'@param oData dataframe
#'@param effect character
#'@return dataframe of covariates
#'@author Omar Elashkar
#'@noRd 
parseCovExpr <- function(expr, oData, effect){
    checkmate::assertFormula(expr)
    checkmate::assertChoice(effec, c("linear", "power", "exp", "hockyStick"))
    checkmate::assertDataFrame(oData)
    
    # ensure all functions are either + or / 
    if (expr[[1]] != as.name("~")){
        stop("Expression must be a formula using '~'")
        }

    currentCovDf <- covExprDf(expr)

    # assert all cov in df 
    if(!all(currentCovDf$covariate %in% names(oData))) {
        stop("Not all covariates are present in the data")
        } else {
        message("All covariates exists in the model data")
        }
    
    ## type is extracted from getData() ==> cont or cat only 
    currentCovDf$type <- sapply(currentCovDf$covariate, function(x) class(oData[[x]]))
    if(!(currentCovDf$type %in% c("numeric", "logical", "factor", "character"))){
        stop("Covariate type cannot be identified")
    }
    currentCovDf$type <- ifelse(currentCovDf$type == "numeric", "cont", "cat")

    ## each of normfactor must be either numeric or 1 if cont  
    # TODO support mean/median/NA
    currentCovDf$normfactor <- ifelse(is.na(currentCovDf$normfactor), 1, currentCovDf$normfactor)
    if(any(currentCovDf$normfactor != 1 &  currentCovDf$normfactor == "cat")){
        stop("Categorical covariates does not have normalization factor")
    }

    ## add col levels for cat covars
    currentCovDf$levels <- unlist(lapply(1:nrow(currentCovDf), function(x){
                                if(currentCovDf$type[x] == "cat"){
                                    covName <- currentCovDf$covariate[x]
                                    unique(oData[[covName]])
                                } else{
                                    NA
                                }
    }))

    ## add cols min and max for cont
    currentCovDf$min <- unlist(lapply(1:nrow(currentCovDf), function(x){
                                if(currentCovDf$type[x] == "cont"){
                                    covName <- currentCovDf$covariate[x]
                                    min(oData[[covName]])
                                } else{
                                    NA
                                }
    }))
    
    currentCovDf$max <- unlist(lapply(1:nrow(currentCovDf), function(x){
                                if(currentCovDf$type[x] == "cont"){
                                    covName <- currentCovDf$covariate[x]
                                    max(oData[[covName]])
                                } else{
                                    NA
                                }
    }))

    currentCovDf$expr <- unlist(lapply(1:nrow(currentCovDf), function(x){
                                if(currentCovDf$type[x] == "cont"){
                                    param <- currentCovDf$param[x]
                                    covName <- currentCovDf$covariate[x]
                                    normfactor <- currentCovDf$normfactor[x]
                                    covTheta <- paste0("theta", param, covName)
                                    if(effect == "power"){
                                        paste0(param, covName,  "= (" ,  covName, "/", normfactor, ")^",  covTheta)
                                    }
                                } else{
                                    # TODO cat
                                    NA
                                }
    }))

    currentCovDf$effect <- effect

    currentCovDf

}

#' Add Covariate to Model Fit (Generic)
#'
#' @param fit Model fit object
#' @param expr Expression eg. CL ~ WT/70 + AGE/80 + ... .
#' @param normDefault Default normalization for continuous covariates. "mean", "median" or NA. Default "median"
#' @param mutiply boolean. If TRUE, multiply the covariate by the parameter, if FALSE, add the covariate to the parameter.
#' @param effect character or list of characters of "linear", "piece_lin", "exp", "power". see details 
#' 
#' `effect` and `normaDefault` are only used if covariate is continuous.
#' `normaDefault` call will be skipped if what covpar expression is normalized by other value
#' 
#' @export
addCovariate <- function(fit, expr, effect) {
    UseMethod("addCovariate")
}

#'@export 
addCovariate.default <- function(fit, expr, effect) {
    stop("addCovariate is not supported for this object")
}

#'@export 
addCovariate.nlmixr2Linearize <- function(fit, expr, effect) {
    covParseDf <- parseCovExpr(expr, nlme::getData(fit), effect = effect)
    covParseDf$Deriv <- paste0("D_", currentCovDf$param)

    ## TODO check all param in model @mattfidler
    # if(!all(currentCovDf$param %in% fit$lhs)) {
    #     stop("get lhs")
    # }
    
    # get D_eta of parameter
}



covExprDf <- function(expr) {
    # Extract parameter
    param <- as.character(expr[[2]])

    # Extract covariate and normfactor parts
    cov_norm <- expr[[3]]
    
    # A helper function to handle nested calls inside the expression
    handle_part <- function(part, covariates, normfactors) {
        if (inherits(part, "call")) {
            if (as.character(part[[1]]) == "/") {
            covariates <- c(covariates, as.character(part[[2]]))
            normfactors <- c(normfactors, as.numeric(as.character(part[[3]])))
        } else if (as.character(part[[1]]) == "+") {
            results1 <- handle_part(part[[2]], covariates, normfactors)
            results2 <- handle_part(part[[3]], results1$covariates, results1$normfactors)
            covariates <- results2$covariates
            normfactors <- results2$normfactors
        }
        } else {
            covariates <- c(covariates, as.character(part))
            normfactors <- c(normfactors, NA)
        }
        return(list(covariates = covariates, normfactors = normfactors))
}

    results <- handle_part(cov_norm, covariates = c(), normfactors = c())

    # Create a data frame with param, covariates, and normfactors
    result <- data.frame(
        param = param,
        covariate = results$covariates,
        normfactor = results$normfactors,
        stringsAsFactors = FALSE
    )

    result
}


#' Automated Search for Inter-individual Variability
#' @fit fit `nlmixr2` object
#'@author Omar Elashkar
#'@export 
iivSearch <- function(fit){
    UseMethod("iivSearch")
}

#'@export 
iivSearch.default <- function(fit){
    stop("Currently supports only linearized models")
}


#'@export 
iivSearch.nlmixr2Linearize <- function(fit){
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
        fit <- nlmixr(newMod, nlme::getData(fit), est = "focei")
        summ <- fit$objDf[, c("OBJF", "AIC", "BIC")]
        summ$search <- x
        list(fit = fit, objDf = summ, search = x)
    })

    objDfAll <- do.call(rbind, lapply(res, function(x) x$objDf))

    list(res = res, summary = objDfAll) # TODO fit and return original model with correct 
}

#' Filter covariance matrix
#' @param oMat covariance matrix 
#' @param filterStr single string with pattern of "-" between entities and "," for within entity correlation

#'@author Omar Elashkar
#'@noRd 
filterEtaMat <- function(oMat, filterStr){
        stopifnot(length(filterStr) == 1)
        resMat <- oMat
        resMat[] <- 0
        elem <- unlist(strsplit(filterStr, "-"))
        for(i in elem){
            elem2 <- unlist(strsplit(i, ","))
            x <- elem2[1]
            y <- ifelse(length(elem2) == 1, elem2[1], elem2[2])
            resMat[x,y] <- oMat[x, y]
        }
        resMat
    }

#' Check if not all fixed effects has eta 
#' @author Omar Elashkar
#' @noRd
hasAllEta <- function(ui){
    etadf <- ui$iniDf[!is.na(ui$iniDf$neta1),]
    thetadf <- ui$iniDf[is.na(ui$iniDf$neta1) & is.na(ui$ini$condition),]
    # TODO remove iov 
    # TODO remove corr
    nrow(etadf) == nrow(thetadf)
}


#'Create a combination 
#'@params char vector of IIV names
#'@author Omar Elashkar
#'@return char vector of possible combinations separated by "-". Correlations are marked by ",".
#'@noRd 
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
