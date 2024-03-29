#' Build formula for the brms from the fit
#'
#' @param fit compiled rxode2 nlmir2 model fit
#' @param covarsVec  character vector of covariates that need to be added
#' @param inicovarsVec covariates vector to add in initial model, default empty
#' @param eta the eta parameter name to construct formula
#' @return formula for the brms
#' @noRd
.buildbrmsFormula <- function(fit,covarsVec,eta,inicovarsVec=NULL){
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  } else {
    ui <- fit$finalUiEnv
  }

  checkmate::assert_character(covarsVec)

  formula <- brms::bf(stats::as.formula(paste0(eta, " ~ a + b")),
                      stats::as.formula(paste0("a ~ ", paste(c("0", covarsVec), collapse = " + "))),
                      stats::as.formula(paste0("b ~ ", paste(c("1", inicovarsVec), collapse = " + "))),
                      nl = TRUE)

  if (!inherits(formula, c("brmsformula" ,"bform"))) {
    stop("BRMS Formula construction fails ")
  } else {
    return(formula)
  }
}

#' Global scale for regularized horseshoe prior
#' @param fit compiled rxode2 nlmir2 model fit
#' @param covarsVec  character vector of covariates that need to be added
#' @param p0 expected number of covariate terms, default 2
#' @return Global shrinkage parameter
#' @noRd
.calTau0 <- function(fit,covarsVec,p0=2){
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  else {
    ui <- fit$finalUiEnv
  }

  checkmate::assert_character(covarsVec)
  D <- length(covarsVec)-1
  n <- nrow(nlme::getData(fit))
  # Global shrinkage parameter, equation
  tau0 <- p0/(D-p0) / sqrt(n)

  if (!is.finite(tau0)| tau0 < 0  ){tau0 <- 0.25}

  tau0
}

normal <- function(...){}
horseshoe <- function(...){}
#' @import utils
utils::globalVariables("tau0")

#' Build Horseshoe prior
#'
#' @param tauZero Global shrinkage parameter
#' @return Prior for the string
#' @noRd
.horseshoePrior <- function(tau0){
  # Check if tau0 is valid
  checkmate::assert_double(tau0)
  priorString <- c(brms::prior(horseshoe(df = 1, scale_global = tau0,df_global = 1), class ="b",nlpar = "a"),
                   brms::prior(normal(0, 10), class = "b", nlpar = "b"))
  # stan variable for parsing
  stanvars <- brms::stanvar(tau0, name='tau0')
  list(priorString, stanvars)
}

lasso <- function(...){}

#' Build Lasso prior
#'
#' @param dfr	Degrees of freedom of the chi-square prior of the inverse tuning parameter
#' @param scalep Scale of the lasso prior
#' @return Prior for the string
#' @noRd
.lassoPrior <- function(df = 1,  scale = 1){

  # Check if given parameters are valid
  checkmate::assert_double(df)
  checkmate::assert_double(scale)

  priorString <- c(brms::prior(lasso(df = 1, scale = 1), class ="b",nlpar = "a"),
                   brms::prior(normal(0, 10), class = "b", nlpar = "b"))

  # stan variable for parsing
  stanvars <- brms::stanvar(df, name='df')+brms::stanvar(scale, name='scale')
  list(priorString, stanvars)
}

#' BRMS model fits for the eta
#'
#' @param fit compiled rxode2 nlmir2 model fit
#' @param covarsVec  character vector of covariates that need to be added
#' @param inicovarsVec covariates vector to add in initial model, default empty
#' @param priorVar priorstring and stanvars list for the covariates
#'  other parameters passed to brm(): warmup = 1000, iter = 2000, chains = 4, cores = 4,
#'  control = list(adapt_delta = 0.99, max_treedepth = 15)
#' @return list of the fitted models
#' @noRd
.fitbrmsModel <- function(fit,covarsVec,inicovarsVec=NULL,priorVar,warmup = 1000, iter = 2000, chains = 4,cores = 2,
                          control = list(adapt_delta = 0.99, max_treedepth = 15),seed=1015){
  #Normalized covariate data
  data <- nlme::getData(fit)
  covData <- normalizedData(data,covarsVec)
  # Extract eta parameters
  etaData <- fit$eta
  etaVector <- colnames(etaData[grepl('eta', colnames(etaData))])
  # Extract Individual column
  uidCol <- .idColumn(data)
  # Make a combined data set of eta parameters and covariate parameters
  combData <- merge(covData,etaData,by=uidCol)

  brms_formulas <- list()
  # brms formulas for all length of the eta parameters

  brms_formulas <- lapply(etaVector, .buildbrmsFormula,fit=fit,covarsVec=covarsVec,inicovarsVec=inicovarsVec)

  # Run brms on all eta parameters
  brms_models <- list()

  brms_models <- suppressWarnings(lapply(brms_formulas,brms::brm,data = combData,family = stats::gaussian(),prior =priorVar[[1]],
                                         stanvars = priorVar[[2]],warmup = warmup, iter = iter, chains = chains,cores = cores,
                                         control = control,seed=seed))
  names(brms_models) <- etaVector
  brms_models
}

#' Create Summary data frame from the BRMS models
#'
#' @param modelList List of BRMS model fits
#' @return Summary data frame of all covariates
#' @noRd
.brmSummarydf <- function(all_models){
  # Check if the model list is named
  checkmate::assert_list(all_models,min.len = 1,names = "named")
  # Construct data frame of estimates by adding eta and covariate column
  dfsList <- a <- lapply(names(all_models),function(x) {
    data.frame(eta=x,covariate= gsub("a_|b_","",rownames(summary(all_models[[x]])$fixed)),
               summary(all_models[[x]])$fixed,row.names = NULL)})
  # Merge Estimates for all eta parameters
  summaryDf <- do.call("rbind", dfsList)
  summaryDf <- summaryDf[!(summaryDf$covariate=="Intercept"),]
  summaryDf
}

tau0 <- NULL

#' Create Horseshoe summary posterior estimates
#'
#' @param fit compiled rxode2 nlmir2 model fit
#' @param covarsVec  character vector of covariates that need to be added
#' @param ...   other parameters passed to brm(): warmup = 1000, iter = 2000, chains = 4, cores = 4,
#'  control = list(adapt_delta = 0.99, max_treedepth = 15)
#' @return Horse shoe Summary data frame of all covariates
#' @export
#' @author  Vishal Sarsani, Christian Bartels
#' @examples
#' \dontrun{
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45; label("Ka")
#'     tcl <- log(c(0, 2.7, 100)); label("Cl")
#'     tv <- 3.45; label("V")
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' d <- nlmixr2data::theo_sd
#' fit <- nlmixr2(one.cmt, d, est = "saem", control = list(print = 0))
#' covarsVec <- c("WT")
#'
#' # Horseshoe summary posterior estimates:
#'
#' #hsDf <- horseshoeSummardf(fit,covarsVec,cores=2)
#' #brms sometimes may throw a Error in sink(type = “output”)
#' #Issue Should be fixed by uninstalling and re-installing rstan
#' }
horseshoeSummardf <- function(fit,covarsVec,...){
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  checkmate::assert_character(covarsVec)
  # Global shrinkage prior estimate
  assignInMyNamespace("tau0", .calTau0 (fit,covarsVec,p0=2))
  # Get prior String
  priorString <- .horseshoePrior(tau0)
  # Fit BRMS models

  .horseshoeModels <-.fitbrmsModel(fit,covarsVec,priorVar = priorString,inicovarsVec=NULL,
                                   warmup = 1000, iter = 2000, chains = 4,cores = 2,
                                   control = list(adapt_delta = 0.99, max_treedepth = 15),seed=1015)

  # Extract Summary of models
  horseshoeSummary <-  .brmSummarydf(.horseshoeModels)
  horseshoeSummary
}

#' Create Lasso summary posterior estimates
#' @param fit compiled rxode2 nlmir2 model fit
#' @param covarsVec  character vector of covariates that need to be added
#' @param ...   other parameters passed to brm(): warmup = 1000, iter = 2000, chains = 4, cores = 4,
#'  control = list(adapt_delta = 0.99, max_treedepth = 15)
#' @return Horse shoe Summary data frame of all covariates
#' @export
#' @author  Vishal Sarsani, Christian Bartels
#'
#' @examples
#' \dontrun{
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45; label("Ka")
#'     tcl <- log(c(0, 2.7, 100)); label("Cl")
#'     tv <- 3.45; label("V")
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' d <- nlmixr2data::theo_sd
#' fit <- nlmixr2(one.cmt, d, est = "saem", control = list(print = 0))
#' covarsVec <- c("WT")
#'
#' # Horseshoe summary posterior estimates:
#'
#' #lassoDf <- lassoSummardf(fit,covarsVec,cores=2)
#' #brms sometimes may throw a Error in sink(type = “output”)
#' #Issue Should be fixed by uninstalling and re-installing rstan
#' }
lassoSummardf <- function(fit,covarsVec,...){
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit")
  }
  checkmate::assert_character(covarsVec)


  # Get prior String
  priorString <- .lassoPrior(df=1,scale=1)
  # Fit BRMS models

  .lassoModels <-.fitbrmsModel(fit,covarsVec,priorVar = priorString,inicovarsVec=NULL,
                               warmup = 1000, iter = 2000, chains = 4,cores = 2,
                               control = list(adapt_delta = 0.99, max_treedepth = 15),seed=1015)

  # Extract Summary of models
  lassoSummary <-  .brmSummarydf(.lassoModels)
  lassoSummary
}
