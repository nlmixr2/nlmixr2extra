#' @importFrom Rcpp evalCpp
#' @importFrom stats setNames approxfun cov cov2cor qchisq qnorm quantile AIC median pchisq sd
#' @importFrom utils head
#' @useDynLib nlmixr2extra, .registration=TRUE
#' @noRd
.getUiFunFromIniAndModel <- function(ui, ini, model) {
  .ls <- ls(ui$meta, all.names=TRUE)
  .ret <- vector("list", length(.ls) + 3)
  .ret[[1]] <- quote(`{`)
  for (.i in seq_along(.ls)) {
    .ret[[.i + 1]] <- eval(parse(text=paste("quote(", .ls[.i], "<-", deparse1(ui$meta[[.ls[.i]]]), ")")))
  }
  .len <- length(.ls)
  .ret[[.len + 2]] <- ini
  .ret[[.len + 3]] <- model
  .retf <- function(){}
  body(.retf) <- as.call(.ret)
  .retf
}

#' Linearly re-parameterize the model to be less sensitive to rounding errors
#'
#' @param fit A nlmixr2 fit to be preconditioned
#' @param estType Once the fit has been linearly reparameterized,
#'   should a "full" estimation, "posthoc" estimation or simply a
#'   estimation of the covariance matrix "none" before the fit is
#'   updated
#' @param ntry number of tries before giving up on a pre-conditioned
#'   covariance estimate
#'
#' @return A nlmixr2 fit object that was preconditioned to stabilize
#'   the variance/covariance calculation
#'
#' @export
#'
#' @references Aoki Y, Nordgren R, Hooker AC. Preconditioning of
#'   Nonlinear Mixed Effects Models for Stabilisation of
#'   Variance-Covariance Matrix Computations. AAPS
#'   J. 2016;18(2):505-518. doi:10.1208/s12248-016-9866-5
#'
#' Format a number for inclusion in a generated model line
#'
#' Uses the full round-trip precision of a double so the preconditioning
#' coefficients are not degraded by the trip through text.
#'
#' @param x A single numeric value
#' @return A character string
#' @noRd
.preCondNum <- function(x) {
  sprintf("%.17g", x)
}

#' Build the preconditioned reparameterization lines
#'
#' Each original theta is rewritten as the linear combination of the
#' preconditioned thetas given by a row of `pre`, i.e.
#' \code{d0[i] = sum_j pre[i, j] * nlmixr2Pre_<d0[j]>}.
#'
#' Built by string assembly rather than symbolically: symengine cannot parse an
#' identifier containing a \code{.}, so a conventional name like \code{add.sd}
#' (as \code{nlmixr2Pre_add.sd}) raised "SymEngine exception: Parse error" and
#' made `preconditionFit()` unusable for most models (#124).
#'
#' @param pre The preconditioning matrix from `preCondInv()`
#' @param d0 Character vector of the original parameter names, in the row and
#'   column order of `pre`
#' @return A single string of newline-separated model lines
#' @noRd
.preCondModExtra <- function(pre, d0) {
  .d <- paste0("nlmixr2Pre_", d0)
  .lines <- vapply(seq_along(d0), function(.i) {
    .coef <- pre[.i, ]
    # drop exact zeros; symengine dropped them too, and they only bloat the model
    .keep <- which(.coef != 0)
    if (length(.keep) == 0L) {
      return(paste0(d0[.i], "=0"))
    }
    # parenthesise every coefficient so a negative one cannot produce `+-1.5`
    .terms <- paste0("(", vapply(.coef[.keep], .preCondNum, character(1)), ")*",
                     .d[.keep])
    paste0(d0[.i], "=", paste(.terms, collapse = "+"))
  }, character(1))
  paste(.lines, collapse = "\n")
}

#' Expand a preconditioning matrix to the full covariance parameter space
#'
#' `preCondInv()` is built from `fit$R`, which spans only the population
#' parameters, while the fit covariance also carries the omega elements.  The
#' reparameterization leaves those untouched, so the transform that takes the
#' preconditioned covariance back to the original scale is `pre` on the theta
#' block and the identity elsewhere; keeping it a single matrix means the
#' theta/omega cross-covariances transform correctly too.
#'
#' @param pre The preconditioning matrix from `preCondInv()`
#' @param covNames Row names of the preconditioned fit covariance
#' @param d0 Original parameter names, in the row/column order of `pre`
#' @return A square matrix the size of the covariance
#' @noRd
.preCondExpand <- function(pre, covNames, d0) {
  .idx <- match(paste0("nlmixr2Pre_", d0), covNames)
  if (anyNA(.idx)) {
    stop("could not find the preconditioned parameters (",
         paste(paste0("nlmixr2Pre_", d0)[is.na(.idx)], collapse = ", "),
         ") in the preconditioned fit covariance",
         call. = FALSE)
  }
  .ret <- diag(length(covNames))
  .ret[.idx, .idx] <- pre
  .ret
}

preconditionFit <- function(fit, estType = c("full", "posthoc", "none"),
                            ntry = 10L) {
  nlmixrWithTiming("covariance", {
    if (!exists("R", fit$env)) {
      stop("this assumes a covariance matrix with a R matrix",
           call. = FALSE
           )
    }
    .R <- fit$R
    .covMethod <- ""
    .i <- 1
    while (.i < ntry & .covMethod != "r,s") {
      .i <- .i + 1
      pre <- preCondInv(.R)
      d0 <- dimnames(fit$R)[[1]]
      modExtra <- .preCondModExtra(pre, d0)
      preInv <- solve(pre)

      .ini <- as.data.frame(fit$ui$iniDf)
      newEst <- setNames(as.vector(preInv %*% matrix(fit$theta[d0])), d0)
      for (v in d0) {
        .w <- which(.ini$name == v)
        .ini$lower[.w] <- -Inf
        .ini$upper[.w] <- Inf
        .ini$est[.w] <- newEst[v]
        .ini$name[.w] <- paste0("nlmixr2Pre_", v)
      }
      .ini <- as.expression(lotri::as.lotri(.ini))
      .ini[[1]] <- quote(`ini`)
      .newModel <- eval(parse(text = paste0("quote(model({", modExtra, "\n", fit$ui$fun.txt, "}))")))
      .newModel <- .getUiFunFromIniAndModel(fit$ui, .ini, .newModel)
      .newModel <- .newModel()
      .ctl <- fit$foceiControl
      estType <- match.arg(estType)
      if (estType == "none") {
        .ctl$maxInnerIterations <- 0
        .ctl$maxOuterIterations <- 0
        .ctl$boundTol <- 0
        .ctl$etaMat <- as.matrix(fit$eta[, -1])
        .ctl$calcTables <- FALSE
        .ctl$compress <- FALSE
      } else if (estType == "posthoc") {
        .ctl$maxOuterIterations <- 0
        .ctl$boundTol <- 0
        .ctl$calcTables <- FALSE
        .ctl$compress <- FALSE
      } else if (estType == "full") {
        .ctl$boundTol <- 0
        .ctl$calcTables <- FALSE
        .ctl$compress <- FALSE
      }
      .ctl$covMethod <- 1L
      ## FIXME compare objective functions
      newFit <- suppressWarnings(nlmixr2est::nlmixr2(.newModel, nlme::getData(fit), est = "focei", control = .ctl))
      .R <- newFit$R
      .covMethod <- newFit$covMethod
    }
    if (.covMethod != "r,s") {
      stop("preconditioning failed after ", ntry, "tries",
           call. = FALSE
           )
    }
    # `pre` only spans the parameters in fit$R; widen it to the covariance's own
    # parameter space so a model with random effects conforms (fit$R is thetas
    # only, fit$cov also carries the omega elements)
    .covNames <- dimnames(newFit$cov)[[1]]
    .A <- .preCondExpand(pre, .covNames, d0)
    cov <- .A %*% newFit$cov %*% t(.A)
    # back to the original parameter names, omega names untouched
    .nm <- sub("^nlmixr2Pre_", "", .covNames)
    dimnames(cov) <- list(.nm, .nm)
    assign("precondition", cov, envir = fit$env)
    .setCov(fit, covMethod = cov)
    assign("covMethod", "precondition", envir=fit$env)
  }, envir=fit)
  return(invisible(fit$env$precondition))
}
