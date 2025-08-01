#' Extract derivatives from NLME model fitted with FOCE
#'
#' @param fit object fitted with nlmixr2 of class nlmixr2.focei
#'
#' @author Omar Elashkar
#' @noRd
getDeriv <- function(fit){
    # if(fit$method != "FOCE") stop("This method requires FOCE method")
    rxode2::assertRxUiMixedOnly(fit, " for the procedure routine 'linearize'", .var.name=fit$modelName)

    ui <- rxode2::assertRxUi(fit)
    innerModel <- ui$foceiModel

    eta <- fit$eta
    eta <- eta[,-1, drop=FALSE]
    eta <- setNames(eta, paste0("ETA[", seq_along(eta), "]"))

    theta <- fit$theta
    reps <- ceiling(nrow(eta) / nrow(theta))
    # theta <- as.data.frame(theta)
    theta <- setNames(fit$theta, paste0("THETA[", seq_along(fit$theta), "]")) %>% 
        t() %>% 
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
        keep = c("DV", setdiff(names(oData),  c("DV", "ID", "TIME", "DVID", 
                                                "ADDL", "EVID", "AMT", "CMT", "RATE", "OCC")))) # add covariates

    eta <- fit$eta[,-1, drop=FALSE]
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
    #if(length(unique(fit$predDf$cond)) > 1){
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
    #}


    derv
}

renameCol <- function(df, new, old){
    names(df)[names(df) == old] <- new
    df
}

#' Generate a Linearization Model From Previous Fit
#' @param ui rx model or fit object.
#' @param focei boolean. If TRUE, use FOCEI linearization with individual and residual linearization. Default is TRUE.
#' @param derivFct boolean. If TRUE, use normalization derivative factors. Default is FALSE.
#' @author Omar Elashkar
#' @export
linModGen <- function(ui, focei = TRUE, derivFct = FALSE){

    if(is.function(ui)){
        ui <- ui()
    }
    if(inherits(ui, "nlmixr2FitCore")){
        stopifnot(all(ui$iniDf == ui$ui$iniDf, na.rm = TRUE))
        ui <- ui$ui
    }
    nlmod <- ui
    
    rxode2::assertRxUiMixedOnly(ui, " for the procedure routine 'linearize'", .var.name=ui$modelName)

    # errSym <- rxode2:::.rxGetVarianceForErrorType(rxUiDecompress(ui), ui$predDf)
    uiEnv <- rxode2::assertRxUi(ui)
    uiEnv <- rxode2::rxUiDecompress(uiEnv)
    assign("derivFct", derivFct, envir = uiEnv)
    on.exit({
      rm("derivFct", envir = uiEnv)
    })
    extractError <- uiEnv$linearizeError
    # fitui <- rxode2::rxUiCompress(fitui)
    epStr <- ui$predDf$var
    errNames <- ui$iniDf$name[!is.na(ui$iniDf$err)]
    errEstim <- ui$iniDf$est[!is.na(ui$iniDf$err)]
    
    noAdderrNames <- ui$iniDf$name[!is.na(ui$iniDf$err) & ui$iniDf$err != "add"]
    noAdderrEstim <- ui$iniDf$est[!is.na(ui$iniDf$err) & ui$iniDf$err != "add"]
    etaNames <- ui$eta

    modelStr <- list()

    # mu block
    for(i in seq_along(etaNames)){
        currEta <- etaNames[i]
        modelStr$muRef[i] <- paste0("mu_", currEta , " = " , "theta.", currEta , " + ", currEta)
    }
    # linearize eta
    for(i in seq_along(etaNames)){
        # D_ETAn * (-OETAn + eta.n)
        currEta <- etaNames[i]
        modelStr$baseEta[i] <- paste0("base_",currEta , "=", "D_", currEta, "*(", "-O_", currEta, "+ mu_", currEta, ")")
    }
    # BASE_TERMS = base1 + ... + basen
    modelStr$baseEta[i+1] <- paste("BASE_TERMS =", paste(paste0("base_", etaNames), collapse = " + "))
    modelStr$baseEta[i+2] <- paste0("y = BASE_TERMS + OPRED")
    # modelStr$baseEta[i+3] <- paste0("eps = y - DV")

    # linearize residuals
    for(i in seq_along(etaNames)){
        currEta <- etaNames[i]
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
    modelStr$baseEps[i+ii+3] <- paste("BASE_ERROR = foceiLin * fct * (", paste(paste0("err_", etaNames),
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
    iniDf <- addThetaToIniDf(iniDf, paste0("theta.", etaNames), ini=0, fix = TRUE)
    if(derivFct){
        if(derivFct & length(noAdderrNames) < 1){stop("Having `derivFct= TRUE` is irrelevant for this error model")}
        iniDf <- addThetaToIniDf(iniDf, paste0(noAdderrNames, ".l"), ini=noAdderrEstim , fix = TRUE)
    }

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
linearize <- function(fit, mceta=c(-1, 10, 100, 1000), relTol=0.4, focei = NA,
    derivFct = FALSE, plot = FALSE){
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
        cli::cli_alert_info("Linearization evaluation matched. Linearization might be feasible ...")
    } else{
        cli::cli_alert_warning("Linearization evaluation mismatched by deltaOFV > 2%. Linearization might be difficult")
        cli::cli_alert_warning("Switching to linearization around predictions only (Variance linearization skipped)")
        linMod <- linMod %>%  model(foceiLin <- 0)
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
            cli::cli_alert_info(paste("Non-Linear OFV: ", oObj, "Linear OFV:", lObj, "mceta:", mceta[i], "Relative OFV Dev", relDev*100, "%"))
            if(i != length(mceta)){
                cli::cli_alert_info(paste("Using mceta = ", mceta[i], "provided inadequate linearization. Trying mceta = ", mceta[i+1], "..."))
            } else{
                if(is.na(focei) & !exists("secondEval")){ # final switch iteration (after estimation) 
                    cli::cli_alert_info("Switching to FOCE after full FOCEI estimation failed ...")
                    secondEval <- evalFun()
                    linMod <- linMod %>%  model(foceiLin <- 0)
                    fitL <- nlmixr(linMod, derv, est="focei",
                        control = nlmixr2est::foceiControl(etaMat = fit, mceta=10, covMethod = "", calcTables=FALSE))
                    oObj <- fit$objDf$OBJF
                    lObj <- fitL$objDf$OBJF

                    relDev <- abs((oObj-lObj)/lObj)
                    
                    if(relDev <=relTol){
                        cli::cli_alert_warning("FOCE Linearization was inadequate for the given tolerance. Try increasing the tolerance, refine the model or fit with mceta")
                    }

                } else{
                    cli::cli_alert_warning("Linearization was inadequate for the given tolerance. Try increasing the tolerance, refine the model or fit with mceta")

                }
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
        "\nLinearized Model Runtime:", sum(fit$time),
        "\nNon-Linearized Model Runtime:", sum(fitL$time)
    )
    cli::cli_alert_info("Linearization Summary:")
    cli::cli_alert_info(m)

    tmpEnv <- fitL$env
    assign("message", c(fitL$message, m), envir=tmpEnv)
    assign("originalUi", fit$finalUi, envir=tmpEnv)

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
#' @return ggplot object
#' @author Omar Elashkar
#' @export
linearizePlot <- function(nl, lin){
    stopifnot(inherits(nl, "nlmixr2FitCore"))
    stopifnot(inherits(lin, "nlmixr2FitCore"))
    stopifnot(inherits(lin, "nlmixr2Linearize"))

    stopifnot(all(nl$iniDf$name[!is.na(nl$iniDf$neta1)] == lin$iniDf$name[!is.na(lin$iniDf$neta1)]))

    l <- function(x, descr) {
        x$factor <- descr
        x
    }
    originalIval <- l(nl$etaObf, "original")
    linearIval <- l(lin$etaObf, "linear")

    fig <- rbind(originalIval, linearIval) %>% 
        tidyr::pivot_longer(cols = c(-c("ID", "factor")), names_to = "parameter",
        values_to = "value") %>% 
        tidyr::pivot_wider(names_from = "factor", values_from = "value") %>% 
        ggplot2::ggplot(aes(x = .data[["original"]], y=.data[["linear"]])) +
        ggplot2::geom_smooth(se = FALSE, method = "lm") +
        ggplot2::geom_point() +
        ggplot2::facet_wrap("parameter", scales = "free") + 
        ggplot2::labs(title = paste("Model OBJ:", nl$objDf$OBJF, "Linear OBJ:", lin$objDf$OBJF))
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
#' @return match list for OFV, omega, eta, and residual variance terms
#' @author Omar Elashkar
#' @export
isLinearizeMatch <- function(fit, linFit, tol = 0.05){

    est <-setNames(fit$iniDf$est, fit$iniDf$name)
    estLin <- setNames(linFit$iniDf$est, linFit$iniDf$name)
    # est <- est[names(estLin)]

    estErrNl <- est[!is.na(fit$iniDf$err)]
    estErrLin <- estLin[!is.na(linFit$iniDf$ntheta) &  linFit$iniDf$fix == FALSE]

    deltaOfv <- all.equal(fit$objDf$OBJF, linFit$objDf$OBJF, tol = tol)
    deltaOmega <- all.equal(fit$omega, linFit$omega, tol = tol)
    deltaEta <- all.equal(fit$eta, linFit$eta, tol = tol)
    deltaErr <- all.equal(estErrNl, estErrLin, tol = tol)
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

#' Parse Covariate Expression and Extract Data Type
#'@param expr expression to parse. Must have no operation but division and addition
#'@param oData dataframe
#'@param effect character of effect type. One of "linear", "power", "exp", "hockyStick"
#'@return data frame with columns: param, covariate, normfactor, type, levels, min, max, expr, effect
#'@author Omar Elashkar
#'@noRd
parseCovExpr <- function(expr, oData, effect){
    checkmate::assertFormula(expr)
    checkmate::assertDataFrame(oData)
    checkmate::assertChoice(effect, c("linear", "power", "exp", "hockyStick"))


    currentCovDf <- covExprDf(expr)

    # assert all cov in df
    if(!all(currentCovDf$covariate %in% names(oData))) {
        stop("Not all covariates are present in the data")
        } else {
        cli::cli_alert_info("All covariates exists in the model data")
        }

    ## type is extracted from getData() ==> cont or cat only
    currentCovDf$type <- sapply(currentCovDf$covariate, function(x) class(oData[[x]]))
    if(!all(currentCovDf$type %in% c("numeric", "logical", "factor", "character"))){
        stop("Covariate type cannot be identified")
    }
    currentCovDf$type <- ifelse(currentCovDf$type == "numeric", "cont", "cat")


    ## add col levels for cat covars
    refCalc <- lapply(1:nrow(currentCovDf), function(x){
                                if(currentCovDf$type[x] == "cat"){
                                    covName <- currentCovDf$covariate[x]
                                    freqTable <- prop.table(table(oData[[covName]]))
                                    refLevel <- names(freqTable)[which.max(freqTable)]
                                    refFreq <- freqTable[refLevel]
                                    list(refLevel=refLevel, refFreq=refFreq)
                                } else{
                                    list(refLevel=NA, refFreq=NA)
                                }
    })

    currentCovDf$refLevel <- lapply(refCalc, function(x) x$refLevel)
    currentCovDf$refFreq <- unlist(lapply(refCalc, function(x) x$refFreq))
    
    currentCovDf$levels <- lapply(1:nrow(currentCovDf), function(x){
                                if(currentCovDf$type[x] == "cat"){
                                    covName <- currentCovDf$covariate[x]
                                    list(unique(oData[[covName]]))
                                } else{
                                    NA
                                }
    })
    

    # cont covariaties
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


    currentCovDf$mean <- unlist(lapply(1:nrow(currentCovDf), function(x){
                                if(currentCovDf$type[x] == "cont"){
                                    covName <- currentCovDf$covariate[x]
                                    mean(oData[[covName]])
                                } else{
                                    NA
                                }
    }))


    currentCovDf$median <- unlist(lapply(1:nrow(currentCovDf), function(x){
                                if(currentCovDf$type[x] == "cont"){
                                    covName <- currentCovDf$covariate[x]
                                    median(oData[[covName]])
                                } else{
                                    NA
                                }
    }))



    ## each of normfactor must be either numeric parse, summary stat or 1 if cont
    currentCovDf$normFactor <- ifelse(is.na(currentCovDf$normFactor), "1", currentCovDf$normFactor)
    currentCovDf$normFactor <- ifelse(currentCovDf$normFactor == "median", currentCovDf$median, currentCovDf$normFactor)
    currentCovDf$normFactor <- ifelse(currentCovDf$normFactor == "mean", currentCovDf$mean, currentCovDf$normFactor)
    currentCovDf$normFactor <- as.numeric(currentCovDf$normFactor)
    if(any(currentCovDf$normFactor != 1 &  currentCovDf$normFactor == "cat")){
        stop("Categorical covariates does not have normalization factor")
    }


    currentCovDf$expr <- lapply(1:nrow(currentCovDf), function(x){
                                param <- currentCovDf$param[x]
                                covName <- currentCovDf$covariate[x]
                                normFactor <- currentCovDf$normFactor[x]
                                covTheta <- paste0("cov.", param, covName)
                                eqName <- paste0("eq.", covTheta)


                                if(currentCovDf$type[x] == "cont"){

                                    if(effect == "power"){
                                        xpr <- paste0(eqName,  "= (" ,  covName, "/", normFactor, ")^",  covTheta)
                                    }
                                    if(effect == "exp"){
                                        xpr <- paste0(eqName,  "= exp(", covTheta, "* (", covName, "-", normFactor, "))")
                                    }
                                    if(effect == "linear"){
                                        xpr <- paste0(eqName,  "= 1 + ", covTheta, "* (", covName, "-", normFactor, ")")
                                    }
                                    if(effect == "hockyStick"){
                                        # FIXME must have 2 covTheta
                                        xpr <- substitute({
                                            if (covName <= normFactor) eqName <- 1+ covTheta * (covName - normFactor)
                                            if (covName > normFactor) eqName <-  1+ covTheta * (covName - normFactor)
                                                },
                                            list(covName = as.name(covName),  
                                                eqName = as.name(eqName), 
                                                covTheta=as.name(covTheta), 
                                                normFactor = normFactor)
                                        )
                                        xpr
                                    }
                                } else{
                                    # FIXME Cat is not fully working model. Need to standarize if must be factor to start with to avoid manual work
                                    refLevel <- currentCovDf$refLevel[x]
                                    refFreq <- currentCovDf$refFreq[x]
                                    catLevels <- unlist(currentCovDf$levels[x])
                                    covIndicator <- paste0(param, covName, ".Ind") 
                                    
                                    if_exprs <- lapply(catLevels, function(lvl) {
                                        substitute(
                                            if (covName == L) covIndicator <- L,
                                                list(L = lvl, covName = as.name(covName), covIndicator = as.name(covIndicator))
                                                )
                                            
                                        })

                                    cov_eq <- substitute( paramcovName <- 1 + covTheta * (refFreq - covIndicator), 
                                        list(covTheta = as.name(covTheta), 
                                            refFreq = refFreq, 
                                            covIndicator = as.name(covIndicator)

                                        ))

                                    xpr <- c(if_exprs, list(cov_eq))
                                }
                                xpr
    })

    currentCovDf$effect <- effect

    currentCovDf

}

#' Add Covariate to Model Fit (Generic)
#'
#' @param fit Model fit object
#' @param expr Expression, or vector or list of expressions. eg. CL ~ WT/70 + AGE/80 + ... .
#' @param ref  Default normalization for continuous covariates. "mean", "median" or NA. Default "median"
#' @param effect character or list of characters of "linear", "piece_lin", "exp", "power". see details
#'
#' `effect` and `normaDefault` are only used if covariate is continuous.
#' `normaDefault` call will be skipped if what covpar expression is normalized by other value
#' 
#' @author Omar Elashkar
#' @export
addCovariate <- function(fit, expr, effect, ref ) {
    UseMethod("addCovariate")
}


#'@rdname addCovariate
#'@export
addCovariate.default <- function(fit, expr, effect="power", ref  = "median") {
    stop("addCovariate is not supported for this object")
}

#'@rdname addCovariate
#'@export
addCovariate.rxUi <- function(fit, expr, effect = "power", ref = "median"){
    if(is.null(oData)){
        stop("Use addData2Rx() first to get covariate adding from this model.")
    }
    
    stop("addCovariate is not supported for this object") # TODO covariate adding with normal models

}

#'@rdname addCovariate
#'@export
addCovariate.nlmixr2FitCore <- function(fit, expr, effect = "power", ref = "median"){
    ui <- fit$ui
    addCovariate(ui)
}

#'@rdname addCovariate
#'@export
addCovariate.nlmixr2Linearize <- function(fit, expr, effect="power", ref  = "median") {

    ui <- fit$ui
    if(is.list(expr)){
        covParseDf <- lapply(expr, function(x) parseCovExpr(x, oData = nlme::getData(fit), effect = effect))
        covParseDf <- do.call(rbind, covParseDf)
    } else{
        covParseDf <- parseCovExpr(expr, nlme::getData(fit), effect = effect)
    }

    if(!any(covParseDf$param %in% ui$eta)){
        cli::cli_alert_info("eta parameters: ", paste(ui$eta, collapse = ", "))
        stop("For linearized models, covariates are added to eta parameters only")
    }
    covParseDf$Deriv <- paste0("D_", covParseDf$param)
    
    covRef <- paste0("cov", covParseDf$param, " = ", 
        covParseDf$Deriv, "*1*(" , covParseDf$param, covParseDf$covariate , " - 1)") # TODO support more covariate effects
    covTermsPrev <- grep("covTerms = ", ui$lstChr, value = F)
    if(length(covTermsPrev) == 0){ # first time add
        covTermLine <- paste0("covTerms = ", 
            paste0("cov", covParseDf$param, collapse = "+"))
            
        v <- c(covParseDf$expr, covRef, covTermLine, ui$lstChr)
    } else { # second time
        covTermLineModel <- ui$lstChr
        # replace line 
        covTermLineModel[covTermsPrev] <- paste0(covTermLineModel[covTermsPrev], " + ", paste0("cov", covParseDf$param, collapse = "+"))
        v <- c(covParseDf$expr, covRef, covTermLineModel)
    }


    iniDf <- addThetaToIniDf(fit$iniDf, paste0("cov.", covParseDf$param, covParseDf$covariate), 1, fix = FALSE)
    newMod <- fit$ui
    newMod <- rxode2::rxUiDecompress(newMod)
    assign("iniDf", iniDf, envir = newMod)

    newMod$lstExpr <- as.list(str2lang(paste0("{", paste(v, collapse="\n"), "}"))[-1])
    newMod <- newMod$fun
    newMod <- newMod()

    newMod <- newMod %>% model(y = BASE_TERMS+ OPRED + covTerms)
    # rxUi
    
    linEnv <- rxUiDecompress(newMod)
    linEnv$ui <- newMod
    linEnv$originalUi <- fit$originalUi
    linEnv$origData <- nlme::getData(fit)
    linEnv$message <- fit$message
    linEnv <- rxUiCompress(linEnv)

    class(linEnv) <- c("nlmixr2Linearize", class(linEnv))

    linEnv
}


#' Parse covariate expression into a data frame
#' @param expr expression to parse
#' @return data frame with columns: param, covariate, normfactor
#' @author Omar Elashkar
#' @noRd
covExprDf <- function(expr) {
    # ensure all functions are either + or /
    if (expr[[1]] != as.name("~")){
        stop("Expression must be a formula using '~'")
        }
    
    # Extract parameter
    param <- as.character(expr[[2]])

    # Extract covariate and normfactor parts
    cov_norm <- expr[[3]]

    # A helper function to handle nested calls inside the expression
    handle_part <- function(part, covariates, normFactors) {
        if (inherits(part, "call")) {
            if (as.character(part[[1]]) == "/") {
            covariates <- c(covariates, as.character(part[[2]]))
            normFactors <- c(normFactors, as.character(part[[3]]))
        } else if (as.character(part[[1]]) == "+") {
            results1 <- handle_part(part[[2]], covariates, normFactors)
            results2 <- handle_part(part[[3]], results1$covariates, results1$normFactors)
            covariates <- results2$covariates
            normFactors <- results2$normFactors
        }
        } else {
            covariates <- c(covariates, as.character(part))
            normFactors <- c(normFactors, NA)
        }
        return(list(covariates = covariates, normFactors = normFactors))
}

    results <- handle_part(cov_norm, covariates = c(), normFactors = c())

    # Create a data frame with param, covariates, and normfactors
    result <- data.frame(
        param = param,
        covariate = results$covariates,
        normFactor = results$normFactors,
        stringsAsFactors = FALSE
    )

    result
}


#' Add a theta to initial parameter data frame
#' @param iniDf initial parameter data frame
#' @param thetaname vector of the thetas to add
#' @param ini initial value for the theta
#' @param fix boolean. If TRUE, the theta will be fixed. Default is FALSE.
#' @return updated initial parameter data frame with the new theta added
#' @author Omar Elashkar
#' @noRd
addThetaToIniDf <- function(iniDf, thetaname, ini, fix = TRUE){
    checkmate::assertCharacter(thetaname, any.missing = FALSE, min.len = 1)
    tryCatch({
        stopifnot(!any(thetaname %in% iniDf$name))
        }, error = function(e) {
            duplicated_names <- thetaname[thetaname %in% iniDf$name]
            stop(paste("Duplicated names found:", paste(duplicated_names, collapse = ", ")))
})


    newTheta <- data.frame(
        ntheta = seq_along(thetaname) + max(iniDf$ntheta, na.rm = TRUE),
        neta1 = NA_real_,
        neta2 = NA_real_,
        name = thetaname,
        lower = -Inf,
        est = ini,
        upper = Inf,
        fix = fix,
        label = NA_character_,
        backTransform = NA_character_,
        condition = NA_character_,
        err = NA_character_
    )
    iniDf <- rbind(iniDf, newTheta)

    iniDf[order(iniDf$ntheta), ]
}


copyEnv <- function(env) {
    list2env(as.list(env, all.names = TRUE), parent = parent.env(env))
}


getData.rxUi <- function(x){
    x$origData
}


#' Add Data to RxUi model
#' @param ui a model.
#' @param data a data frame to attach.
#' The data can then be accessed via `nlme::getData(ui)`
#' @return rxUi.
#' @author Omar Elashkar
#' @export
addData2Rx <- function(ui, data){
    
    if(is.function(ui)){
        ui <- ui()
    }

    checkmate::assertClass(ui, "RxUi")
    checkmate::assertDataframe(data)
    
    ui <- rxUiDecompress(ui)
    ui$origData <- data
    rxUiCompress(ui)
}


#' Covariate Finding Using Different Algorithms
#'@author Omar Elashkar
#'@noRd
# TODO 
covSearch <- function(ui, covSpace, method = "scm", pValBack=0.01, pValForward=0.05){
    UseMethod("covSearch")
}

