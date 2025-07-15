

#' Extract derivatives from NLME model fitted with FOCE
#'
#' @param fit object fitted with nlmixr2 of class nlmixr2.focei
#'
#' @author Omar Elashkar
#'
#' @noRd
getDeriv <- function(fit){
    if(fit$method != "FOCE") stop("This method requires FOCE method")
    rxode2::assertRxUiMixedOnly(fit, " for the procedure routine 'linerize'", .var.name=fit$modelName)

    ui <- rxode2::assertRxUi(fit)
    innerModel <- ui$foceiModel

    eta <- fit$eta

    eta <- eta[,-1]
    eta <- setNames(eta, paste0("ETA[", seq_along(eta), "]"))

    theta <- fit$theta
    reps <- ceiling(nrow(eta) / nrow(theta))
    theta <- setNames(fit$theta, paste0("THETA[", seq_along(fit$theta), "]")) |>
        t() |>
        as.data.frame()

    theta <- theta[rep(1, nrow(eta)),]
    params_df <- cbind(theta, eta)

    oData <- nlme::getData(fit)
    names(oData) <- toupper(names(oData))

    
    stopifnot(all(c("rx_pred_", "rx_r_") %in% innerModel$inner$lhs))
    stopifnot(sum(grepl("rx__sens_rx_pred__BY_ETA_\\d+___", innerModel$inner$lhs)) == ncol(eta))
    stopifnot(sum(grepl("rx__sens_rx_r__BY_ETA_\\d+___", innerModel$inner$lhs)) >= ncol(eta))

    derv <- rxode2::rxSolve(innerModel$innerOeta, oData, params=params_df,
        addDosing=FALSE,
        keep = c("DV"))
        # keep = c("DV", setdiff(names(oData),  c("ID", "TIME", "EVID", "AMT", "CMT"))))
    
    
    # OPRED
    derv <- renameCol(derv, "OPRED", "rx_pred_")

    # D_ETA
    for(i in seq_len(ncol(fit$eta))){
        derv <- renameCol(derv, paste0("D_ETA", i), paste0("rx__sens_rx_pred__BY_ETA_", i, "___"))
    }

    #O_ETA
    for(i in seq_len(ncol(fit$eta))){
        derv <- renameCol(derv, paste0("O_ETA", i), paste0("rx__ETA", i))
    }

    #O_EPS
    derv <- renameCol(derv, "O_ResVar", "rx_r_")

    # O_IRES 
    derv$O_IRES <- fit$IRES

    # D_EPS
    derv$D_ResVar <- 1

    # D_EPSETA
    for(i in seq_len(ncol(eta))){
        derv <- renameCol(derv, paste0("D_VAR_ETA_",1, "_", i), paste0("rx__sens_rx_r__BY_ETA_", i, "___"))
    }

    # add OCMT col
    if(length(unique(fit$predDf$cond)) > 1){
        # derv$OCMT <- as.integer(derv$CMT)

        mv <- rxode2::rxModelVars(innerModel$inner)
        cmtName <- c(mv$state, mv$stateExtra)
        cmtDf <- data.frame(CMT=seq_along(cmtName), cond=cmtName)
        predDf <- fit$predDf
        predDf$OCMT <- predDf$cmt
        cmtDf <- merge(cmtDf, predDf)
        cmtDf <- cmtDf[,c("OCMT", "CMT", "dvid")]

        derv <- merge(derv, cmtDf)
        derv$CMT <- NULL
 
    }

    derv
}

renameCol <- function(df, new, old){
    names(df)[names(df) == old] <- new
    df
}

linModGen <- function(fit, focei = TRUE){
    rxode2::assertRxUiMixedOnly(fit, " for the procedure routine 'linerize'", .var.name=fit$modelName)

    # errSym <- rxode2:::.rxGetVarianceForErrorType(rxUiDecompress(fit), fit$predDf)
    extractError <- fit$linearizeError
    epStr <- fit$predDf$var
    errNames <- fit$iniDf$name[!is.na(fit$iniDf$err)]
    nlmod <- fit$finalUi
    etaDf <- fit$eta[,-1]

    modelStr <- list()

    # linearize eta
    for(i in seq_along(colnames(etaDf))){
        # D_ETAn * (-OETAn + eta.n)
        modelStr$baseEta[i] <- paste0("base",i , "=", "D_ETA", i, "*(", "-O_ETA", i, "+", colnames(etaDf)[i], ")")
    }
    # BASE_TERMS = base1 + ... + basen
    modelStr$baseEta[i+1] <- paste("BASE_TERMS =", paste(paste0("base", seq(ncol(etaDf))), collapse = " + "))
    modelStr$baseEta[i+2] <- paste0("y = BASE_TERMS + OPRED")
    # modelStr$baseEta[i+3] <- paste0("eps = y - DV")

    # linearize residuals
    for(i in seq_along(colnames(etaDf))){
        # ERRn = D_VAR_ETA_1_n*(-O_ETAn + eta.n)
        modelStr$baseEps[i] <- paste0("err",i , "=", "D_VAR_ETA_1_", i, "*(", "-O_ETA", i, "+", colnames(etaDf)[i], ")")
    }

    errSym <- extractError$rxR2

    # for single ep, this will be single line
    # rxR2
    for(ii in seq_along(errSym)){
        modelStr$baseEps[i+ii] <- errSym[ii]
    }

    modelStr$baseEps[i+ii+1] <- "r <- sqrt(rxR2)"
    modelStr$baseEps[i+ii+2] <- paste0("foceiLin <- ", ifelse(focei, 1, 0))
    modelStr$baseEps[i+ii+3] <- paste("BASE_ERROR = foceiLin * (", paste(paste0("err", seq(ncol(etaDf))), collapse = " + "), ")/(2*r) + r")
    modelStr$baseEps[i+ii+4] <- "rxR = BASE_ERROR"

    for(i in seq_along(extractError$err)){
        modelStr$basePred[i] <-   extractError$err[i]
    }

    # iniDf already captures final estimates
    iniDf <- nlmod$iniDf[nlmod$iniDf$name %in% c(colnames(etaDf), errNames) , ]

    # n_theta == n_errors
    iniDf$ntheta[!is.na(iniDf$ntheta)] <- seq_along(na.omit(iniDf$ntheta))

    ini(nlmod) <- iniDf
    v <- c(modelStr$baseEta, modelStr$baseEps, 
                        extractError$tipred, modelStr$basePred) 
    # model(nlmod) <-  v

    nlmod <- rxode2::rxUiDecompress(nlmod)
    nlmod$lstExpr <- as.list(str2lang(paste0("{", paste(v, collapse="\n"), "}"))[-1])
    nlmod <- rxode2::rxUiCompress(nlmod)
    nlmod <- nlmod$fun
    nlmod <- nlmod()

    nlmod
}

#' Perform linearization of a model fitted using FOCEI
#' @param fit fit of nonlinear model.
#' @param mceta a numeric vector for mceta to try. See decription.
#' @param relTol relative deviation tolerance between original and linearized models objective functions. Used for switching if focei = NA. See details.
#' @param focei Default is NA for automatic switch from FOCEI to FOCE if failed. See details.
#' @param plot boolean. Print plot of linearized vs original
#' 
#' `plot` argument can only print ggplot figure with default settings. 
#' If a user wish to capture the plot, one might use `linearizePlot()` call. 
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
    linMod <- linModGen(fit, focei = ifelse(is.na(focei), TRUE, focei))

    evalFun <- function(){
        # check linearization feasbility
        fitL <- evalLinModel(fit, linMod, derv)
        fitL_map <- evalLinModel(fit, linMod, derv, 1000)
        
        list(
            lObj_exact = fitL$objDf$OBJF,
            lObj_map = fitL_map$objDf$OBJF,
            oObj = fit$objDf$OBJF,
            message = paste("Non-Linear OFV: ", oObj, 
            "Eval Linear OFV (Exact):", lObj_exact, 
            "Eval Linear OFV (MAP):", lObj_map)
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
            }
        }
    }

    if(plot){
        print(linearizePlot(fit, fitL))
    }

    fitL <- nlmixr2est::addTable(fitL)
    nlme::getVarCov(fitL)

    message("---Linearization Summary---")
    message("Relative Tolerance:", relTol)
    message("Linearization method:", ifelse(is.na(focei), "Auto", ifelse(focei, "FOCEI", "Variance Skipped")))
    message("Skipping variance during linearization:", ifelse(is.null(secondEval), FALSE, focei))
    message(paste("Non-Linear OFV: ", oObj))
    message(paste("Linear OFV (Exact):", ifelse(is.null(secondEval), firstEval$lObj_exact, secondEval$lObj_exact)))
    message(paste("Linear OFV (MAP):", ifelse(is.null(secondEval), firstEval$lObj_exact, secondEval$lObj_exact)))
    message(paste("mceta:", mceta[i]))
    message(paste("Relative OFV Dev", relDev*100, "%"))
    fitL
}

#' Plot original versus linear iobj and etas
#' @param nl non-linear fitting object
#' @param l linear fitting object
#' @author Omar Elashkar
#' @return ggplot object
#' @export
linearizePlot <- function(nl, lin){
    # TODO assertion that both objects are siblings. Random effects and errors are names same

    stopifnot(inherits(nl, "nlmixr2FitCore"))
    stopifnot(inherits(lin, "nlmixr2FitCore"))

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

evalLinModel <- function(fit, linMod, derv, innerInter = 0){
    fitL <- nlmixr(linMod, derv, est="focei", 
            control = nlmixr2est::foceiControl(etaMat = fit, mceta=-1, 
            covMethod = "", calcTables=FALSE, maxInnerIterations=innerInter, maxOuterIterations=0L))
    fitL
}


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
