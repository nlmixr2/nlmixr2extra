#' Factor for the additive expression
#'
#' @param env environment for the rxode2 model
#'
#' @param pred1 The pred1 endpoint that is being considered
#'
#' @return quote(fct<-0), since with additive error, the residuals do
#'   not depend on the predictions.
#' @noRd
#' @author Matthew L. Fidler
.rxGetFctForErrorAdd <- function(env, pred1) {
  str2lang("fct <- 0")
}

#' Get the variance factor for the proportional error
#'
#' @param env rxode2 environment
#' @param pred1 prediction endpoint that is being considered
#' @return factor for the proportional error
#' @noRd
#' @author Matthew L. Fidler
.rxGetFctForErrorProp <- function(env, pred1) {
  .f <- .rxGetVarianceForErrorPropOrPowF(env, pred1)
  if (!is.na(pred1$b)) {
    .p1 <- str2lang(pred1$b)
  } else {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err %in% c("prop", "propF", "propT") & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p1 <- str2lang(env$iniDf$name[.w])
      .plast <- str2lang(paste0(env$iniDf$name[.w], ".l"))
    } else {
      stop("cannot find proportional standard deviation", call.=FALSE)
    }
  }
  if (pred1$variance) {
    bquote(fct <- .(.p1)/(.plast))
  } else {
    bquote(fct <- (.(.p1))^2/((.plast)^2))
  }
}
#' Get the factor for the power error
#'
#' @param env environment of the rxode2 model
#' @param pred1 The pred1 endpoint that is being considered
#' @return A expression for the current factor for linearization
#' @noRd
#' @author Matthew L. Fidler
.rxGetFctForErrorPow <- function(env, pred1) {
  .f <- .rxGetVarianceForErrorPropOrPowF(env, pred1)
  .cnd <- pred1$cond
  if (!is.na(pred1$b)) {
    .p1 <- str2lang(pred1$b)
  } else {
    .w <- which(env$iniDf$err %in% c("pow", "powF", "powT") & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p1 <- str2lang(env$iniDf$name[.w])
      .plast <- str2lang(paste0(env$iniDf$name[.w], ".l"))
    } else {
      stop("cannot find power standard deviation", call.=FALSE)
    }
  }
  if (pred1$variance) {
    bquote(fct <- .(.p1)/(.plast))
  } else {
    bquote(fct <- (.(.p1))^2/((.plast)^2))
  }
}

#' Get the additive + proportional factor
#'
#' @param env environment of the rxode2 model
#' @param pred1 The pred1 endpoint that is being considered
#' @return A expression for the current factor for linearization
#' @noRd
#' @author Matthew L. Fidler
.rxGetFctForErrorAddProp <- function(env, pred1) {
  if (!is.na(pred1$a)) {
    .p1 <- str2lang(pred1$a)
    .p1last <- str2lang(paste0(pred1$a, ".l"))
  } else {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err %in% c("add", "lnorm", "probitNorm", "logitNorm") & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p1 <- str2lang(env$iniDf$name[.w])
      .p1last <- str2lang(paste0(env$iniDf$name[.w], ".l"))
    } else {
      stop("cannot find additive standard deviation", call.=FALSE)
    }
  }
  if (!is.na(pred1$b)) {
    .p2 <- str2lang(pred1$b)
    .p2last <- str2lang(paste0(pred1$b, ".l"))
  } else {
    .cnd <- pred1$cond
    .w <- which(env$iniDf$err %in% c("prop", "propT", "propF") & env$iniDf$condition == .cnd)
    if (length(.w) == 1L) {
      .p2 <- str2lang(env$iniDf$name[.w])
      .p2last <- str2lang(paste0(env$iniDf$name[.w], ".l"))
    } else {
      stop("cannot find proportional standard deviation", call.=FALSE)
    }
  }
  if (pred1$addProp == "default") {
    .addProp <- rxGetControl(env, "addProp", getOption("rxode2.addProp", "combined2"))
  } else {
    .addProp <- pred1$addProp
  }
  if (pred1$variance) {
    if (.addProp == "combined2") {
      bquote(fct <- .(.p2)/.(.p2last))
    } else {
      bquote(fct <- (.(.p2)*OPRED + sqrt(.(.p2)*(.p1)))/(.(.p2last)*OPRED + sqrt(.(.p2last)*(.p1last))))
    }
  } else {
    if (.addProp == "combined2") {
      bquote(fct <- .(p2)^2 / (.p2last)^2)
    } else {
      bquote(fct <- (.(.p2)^2*OPRED + .(.p2)*(.p1))/(.(.p2last)^2*OPRED + .(.p2last)*(.p1last)))
    }
  }
}

#' Get the factor adjustment for a specific error type
#'
#' @param env  environment of the rxode2 model
#' @param pred1 The pred1 endpoint that is being considered
#' @return A expression for the current factor for linearization
#' @noRd
#' @author Matthew L. Fidler
.rxGetFctForErrorType <- function (env, pred1) {
  if (env$derivFct) {
    switch(as.character(pred1$errType),
           add = .rxGetFctForErrorAdd(env, pred1),
           prop = .rxGetFctForErrorProp(env, pred1),
           pow = .rxGetFctForErrorPow(env, pred1),
           `add + prop` = .rxGetFctForErrorAddProp(env, pred1),
           `add + pow` = quote(fct <- 1))
  } else {
    quote(fct <- 1)
  }
}


#' Replace rx_pred_f and rx_pred_ with the right OPRED values
#'
#' @param expr expression to replace
#' @param env rxode2 ui model
#' @param pred1 rxode2 single endpoint value
#' @return expression with OPRED for linearized model
#' @noRd
#' @author Matthew L. Fidler
.replaceFwithOpred <- function(expr, env, pred1, y="OPRED") {
  if (is.name(expr)) {
    if (identical(expr, quote(rx_pred_f_))) {
      # These need to be back-transformed to the original scale.
      # Because of this it needs to use the transformation
      .yj <- as.double(pred1$transform) - 1
      if (.yj == 2) {
        # No transformation
        str2lang(y)
      } else if (.yj == 3) {
        # lognormal transformation
        str2lang(paste0("exp(", y, ")"))
      } else {
        # Otherwise use the inverse transformation
        # rxTBSi(DV, rx_lambda_, rx_yj_, rx_low_, rx_hi_)
        .lambda <- deparse1(rxode2::.rxGetLambdaFromPred1AndIni(env, pred1))
        .low <- deparse1(rxode2::.rxGetLowBoundaryPred1AndIni(env, pred1))
        .hi <- deparse1(rxode2::.rxGetHiBoundaryPred1AndIni(env, pred1))
        str2lang(sprintf("rxTBSi(%s, %s, %s, %s, %s)", y,
                         .lambda, .yj, .low, .hi))
      }
    } else if (identical(expr, quote(rx_pred_))) {
      # This is the transformed predictions
      str2lang("OPRED")
    } else {
      expr
    }
  } else if (is.call(expr)) {
    as.call(c(list(expr[[1]]),
              lapply(expr[-1], .replaceFwithOpred, env=env, pred1=pred1, y=y)))
  } else {
    expr
  }
}

#' This creates the linearization error lines from a rxui model
#'
#' @param line line to parse
#'
#' @return error lines for posology
#'
#' @export
#'
#' @keywords internal
#'
#' @author Matthew L. Fidler
linearizeErrorLines <- function(line) {
  UseMethod("linearizeErrorLines")
}

#' @rdname linearizeErrorLines
#' @export
linearizeErrorLines.norm <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  ret <- vector("list", 1)
  list(
    bquote(rxR2 <- .(.replaceFwithOpred(rxode2::.rxGetVarianceForErrorType(env, pred1), env=env, pred1=pred1, y="OPRED"))),
    .rxGetFctForErrorType(env, pred1))
}

#' @rdname linearizeErrorLines
#' @export
linearizeErrorLines.default  <- function(line) {
  stop("distribution not supported", call.=FALSE)
}


# This handles the errors for focei
linearizeErrorLinesObject <- function(x, line) {
  pred_df <- x$predDf
  if (line > nrow(pred_df)) {
    return(NULL)
  }
  pred_line <- pred_df[line, ]
  ret <- list(x, pred_line, line)
  class(ret) <- c(paste(pred_line$distribution), "linearizeErrorLines")
  ret
}

#' @rdname linearizeErrorLines
#' @export
linearizeErrorLines.rxUi <- function(line) {
  pred_df <- line$predDf
  lapply(seq_along(pred_df$cond), function(c) {
    mod <- linearizeErrorLinesObject(line, c)
    linearizeErrorLines(mod)
  })
}


rxUiGet.linearizeError <- function(x, ...) {
  .ui <- x[[1]]
  .errLines <- linearizeErrorLines(.ui)
  .predDf <- .ui$predDf
  .expr <-
  .tipred <- lapply(seq_along(.predDf$cmt),
                    function(i) {
                      .replaceFwithOpred(expr=str2lang("TIPRED <- rx_pred_f_"),
                                                  env=.ui,
                                                  pred1=.predDf[i, ],
                                                  y="y")
                    })
  if (all(vapply(seq_along(.tipred), function(i) {
    identical(.tipred[[i]], str2lang("TIPRED <- y"))
  }, logical(1)))) {
    .tipred <- "TIPRED <- y"
  } else {
    .tipred <- vapply(seq_along(.predDf$cmt),
                      function(i) {
                        paste(deparse(as.call(list(quote(`if`), as.call(list(quote(`==`),
                                                                             quote(OCMT), as.numeric(.predDf$cmt[i]))),
                                                   as.call(c(list(quote(`{`)),
                                                             .tipred[[i]]))))), collapse="\n")

                      }, character(1))

    .tipred <- strsplit(paste(.tipred, collapse="\n"), "\n")[[1]]
  }
  if(length(.predDf$cmt) > 1){
    .rxR2 <- vapply(seq_along(.predDf$cmt),
                    function(i) {
                      paste(deparse(as.call(list(quote(`if`), as.call(list(quote(`==`),
                                                                           quote(OCMT), as.numeric(.predDf$cmt[i]))),
                                                 as.call(c(list(quote(`{`)),
                                                           .errLines[[i]]))))), collapse="\n")

                    }, character(1))
  } else {
    .rxR2 <- vapply(seq_along(.errLines[[1]]),
                    function(i) {
                      paste(deparse(.errLines[[1]][[i]]), collapse="\n")
                    }, character(1))
  }

  .errModel <-
    vapply(seq_along(.predDf$cmt),
           function(i) {
             .pred1 <- .predDf[i, ]
             .lambda <- deparse1(rxode2::.rxGetLambdaFromPred1AndIni(.ui, .pred1))
             .low <- deparse1(rxode2::.rxGetLowBoundaryPred1AndIni(.ui, .pred1))
             .hi <- deparse1(rxode2::.rxGetHiBoundaryPred1AndIni(.ui, .pred1))
             .transform <- paste0(.predDf$transform[i])
             .first <- switch(.transform,
                              boxCox="add(rxR)",
                              yeoJohnson="add(rxR)",
                              untransformed="add(rxR)",
                              lnorm="lnorm(rxR)",
                              logit=sprintf("logitNorm(rxR, %s, %s)", .low, .hi),
                              `logit + yeoJohnson`=sprintf("logitNorm(rxR, %s, %s)", .low, .hi),

                              probit=sprintf("probitNorm(rxR, %s, %s)", .low, .hi),
                              `probit + yeoJohnson`=sprintf("probitNorm(rxR, %s, %s)", .low, .hi),
                              `logit + boxCox`=sprintf("logitNorm(rxR, %s, %s)", .low, .hi),
                              `probit + boxCox`=sprintf("probitNorm(rxR, %s, %s)", .low, .hi)
                              )

             .last <- switch(.transform,
                             boxCox=sprintf(" + boxCox(%s) + dv()", .lambda),
                             yeoJohnson=sprintf("+ yeoJohnson(%s) + dv()", .lambda),
                             untransformed="",
                             lnorm="+dv()",
                             logit="+dv()",
                             `logit + yeoJohnson`=sprintf("+yeoJohnson(%s)+dv()", .lambda),

                             probit="+dv()",
                             `probit + yeoJohnson`=sprintf("+yeoJohnson(%s)+dv()", .lambda),
                             `logit + boxCox`=sprintf(" + boxCox(%s) + dv()", .lambda),
                             `probit + boxCox`=sprintf(" + boxCox(%s) + dv()", .lambda)
                             )
             deparse1(str2lang(paste0("y", i, " ~ ",  .first, .last)))
           }, character(1), USE.NAMES=FALSE)

  .errModel <- c(vapply(seq_along(.predDf$cmt),
                      function(i) {
                        deparse1(str2lang(sprintf("y%s <- y", i)))
                      }, character(1), USE.NAMES=FALSE),
                 .errModel)

  list(rxR2=strsplit(paste(.rxR2, collapse="\n"), "\n")[[1]],
       tipred=.tipred,
       err=.errModel)
}
