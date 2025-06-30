#' Extract derivatives from NLME model fitted with FOCE
#' @param fit object fitted with nlmixr2 of class nlmixr2.focei
#' @author Omar Elashkar
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
        derv <- renameCol(derv, paste0("D_EPSETA_",1, "_", i), paste0("rx__sens_rx_r__BY_ETA_", i, "___"))
    }
    
    derv
}

renameCol <- function(df, new, old){
    names(df)[names(df) == old] <- new
    df
}
