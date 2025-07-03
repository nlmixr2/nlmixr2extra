

#' Extract derivatives from NLME model fitted with FOCE
#'
#' @param fit object fitted with nlmixr2 of class nlmixr2.focei
#'
#' @author Omar Elashkar
#'
#' @noRd
getDerv <- function(fit){
    if(fit$method != "FOCE") stop("This method requires FOCE method")
    stopifnot(inherits(fit, "nlmixr2.focei"))



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

    o_data <- nlme::getData(fit)
    names(o_data) <- toupper(names(o_data))

    stopifnot(all(c("rx_pred_", "rx_r_") %in% innerModel$inner$lhs))
    stopifnot(sum(grepl("rx__sens_rx_pred__BY_ETA_\\d+___", innerModel$inner$lhs)) == ncol(eta))
    stopifnot(sum(grepl("rx__sens_rx_r__BY_ETA_\\d+___", innerModel$inner$lhs)) >= ncol(eta))

    derv <- rxSolve(innerModel$innerOeta, o_data, params=params_df,
        addDosing=TRUE,
        keep = c("DV", setdiff(names(o_data),  c("ID", "TIME", "EVID", "AMT", "CMT"))))

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

    # D_EPS
    derv$D_ResVar <- 1

    ## TODO Add stops if error models different from what is going to be supported

    # D_EPSETA
    for(i in seq_len(ncol(eta))){
        derv <- renameCol(derv, paste0("D_VAR_ETA_",1, "_", i), paste0("rx__sens_rx_r__BY_ETA_", i, "___"))
    }

    derv
}

renameCol <- function(df, new, old){
    names(df)[names(df) == old] <- new
    df
}

linModGen <- function(fit){

    # if(nrow(fit$predDf) != 1){
    #     stop("Mutiple endpoints linearization is not supported")
    # }

    if(!(fit$predDf$errType %in% c("add", "prop", "add + prop"))){
        stop("Error model not supported")
    }

    errSym <- rxode2:::.rxGetVarianceForErrorType(rxUiDecompress(fit), fit$predDf)
    epStr <- fit$predDf$var
    errNames <- fit$iniDf$name[!is.na(fit$iniDf$err)]
    nlmod <- fit$finalUi
    eta_df <- fit$eta[,-1]

    modelstr <- list()

    # linearize eta
    for(i in seq_along(colnames(eta_df))){
        # D_ETAn * (-OETAn + eta.n)
        modelstr$baseEta[i] <- paste0("base",i , "=", "D_ETA", i, "*(", "-O_ETA", i, "+", colnames(eta_df)[i], ")")
    }
    # BASE_TERMS = base1 + ... + basen
    modelstr$baseEta[i+1] <- paste("BASE_TERMS =", paste(paste0("base", seq(ncol(eta_df))), collapse = " + "))
    modelstr$baseEta[i+2] <- paste0("F = BASE_TERMS + OPRED")

    # linearize residuals
    for(i in seq_along(colnames(eta_df))){
        # ERRn = D_VAR_ETA_1_n*(-O_ETAn + eta.n)
        modelstr$baseEps[i] <- paste0("err",i , "=", "D_VAR_ETA_1_", i, "*(", "-O_ETA", i, "+", colnames(eta_df)[i], ")")
    }


    errSym <- gsub("rx_pred_f_", "OPRED", deparse1(errSym))
    modelstr$baseEps[i+1] <- paste("BASE_ERROR = (", paste(paste0("err", seq(ncol(eta_df))), collapse = " + "), ")")
    modelstr$baseEps[i+1] <-  paste(modelstr$baseEps[i+1], "+", errSym)
    # modelstr$baseEps[i+2] <- paste0("R2 = BASE_ERROR +", errSym, ")")
    modelstr$baseEps[i+2] <- "R2 = BASE_ERROR"

    modelstr$basePred[1] <- paste0(epStr ,"  = F")
    modelstr$basePred[2] <- paste0(epStr, " ~ add(R2) + var()")

    # iniDf already captures final estimates
    inidf <- nlmod$iniDf[nlmod$iniDf$name %in% c(colnames(eta_df), errNames) , ]

    # n_theta == n_errors
    inidf$ntheta[!is.na(inidf$ntheta)] <- seq_along(na.omit(inidf$ntheta))

    # newMod <- modelExtract(nlmod, endpoint=FALSE)
    # model(nlmod) <- model(newMod)
    ini(nlmod) <- inidf
    model(nlmod) <-   c(modelstr$baseEta, modelstr$baseEps, modelstr$basePred)

    nlmod
}

#' Perform linearization of a model fitted using FOCEI
#' @author Omar Elashkar
#' @export
linearize <- function(fit, mceta=c(-1, 10, 100, 1000), relTol){  # TODO mceta

    derv <- getDerv(fit)
    linMod <- linModGen(fit)
    fitL <- nlmixr(linMod, derv, est="focei",
        control = nlmixr2::foceiControl(etaMat = as.matrix(fit$eta[-1])))

    oObj <- fit$objDf$OBJF
    lObj <- fitL$objDf$OBJF

    message(paste("Non-Linear OFV: ", oObj, "Linear OFV:", lObj))

    fitL
}
