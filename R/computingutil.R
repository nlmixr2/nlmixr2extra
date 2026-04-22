# Function to calculate population standard deviation from group means
#' @param x vector of values to calculate standard deviation.
#' @return population standard deviation
#' @noRd
.sd.p <- function(x) {
  sd(x)*sqrt((length(x)-1)/length(x))
}

.bootstrapEnv <- new.env(parent=emptyenv())
.bootstrapEnv$nSampIndiv <- 0L

#' Function to return pop mean, pop std of a given covariate
#'
#' @param data given data frame
#' @param covariate the covariate that needs popmean,stdev
#'
#' @return list containing values of population mean, standard deviation
#' @noRd
.popMeanStd <- function(data, covariate) {

  checkmate::assertDataFrame(data,col.names = "named")
  checkmate::assertCharacter(covariate,len = 1,any.missing = FALSE )

  if (inherits(try(str2lang(covariate)),"try-error")) {
    stop("`varName` must be a valid R expression",call. = FALSE)
  }

  .new <- intersect(names(data), covariate)
  if (length(.new) == 0L) stop("covariate specified not in original dataset", call.=FALSE)

  #extract Individual ID from data frame

  uidCol <- .idColumn(data)
  # mean by groups (Individual)
  groupMeans <- with(data, ave(get(covariate),get(uidCol), FUN = function(x) mean(x, na.rm = TRUE)))
  # pop mean
  popMean <- mean(groupMeans, na.rm=TRUE)
  # pop std
  popStd <- .sd.p (groupMeans)

  .meanStd <-c(popMean,popStd)
  names(.meanStd) <- c("popmean","popstd")
  .meanStd
}

#' Function to return normalization column of a covariates
#'
#' @param data given a data frame
#' @param covariate the covariate that needs normalization
#'
#' @return data frame with normalized covariate
#' @noRd
.normalizeDf <- function(data, covariate,sub=TRUE) {
  checkmate::assertDataFrame(data,col.names = "named")
  checkmate::assertCharacter(covariate,len = 1,any.missing = FALSE )

  if (inherits(try(str2lang(covariate)),"try-error")) {
    stop("`varName` must be a valid R expression",call. = FALSE)
  }
  .new <- intersect(names(data), covariate)
  if (length(.new) == 0L) stop("covariate specified not in original dataset")

  if (is.factor(data[[covariate]])) {
    data
  } else {
    # Column name for the standardized covariate
    datColNames <- paste0("normalized_", covariate)
    # popMean
    .popMean = .popMeanStd(data,covariate)[[1]]
    # popStdev
    .popStd = .popMeanStd(data,covariate)[[2]]
    # add standardized covariate values to the data frame
    data[,datColNames] <- (data[,covariate]-.popMean)/(.popStd)
    data
  }
}

#' Function to return data of normalized covariates
#'
#' @param data a dataframe with covariates to normalize
#' @param covarsVec a list of covariate names (parameters) that need to be estimates
#' @param replace replace the original covariate data with normalized data for easier updated model.
#'
#' @return data frame with all normalized covariates
#' @author Vishal Sarsani
#' @export
#'
#' @examples
#'
#' d <- nlmixr2data::theo_sd
#' d$SEX <-0
#' d$SEX[d$ID<=6] <-1
#'
#' covarsVec <- c("WT")
#'
#' # Normalized covariate (replaced)
#' df1 <- normalizedData(d,covarsVec,replace=TRUE)
#'
#' # Normalized covariate (without replacement)
#' df2 <- normalizedData(d,covarsVec,replace=FALSE)
normalizedData <- function(data,covarsVec,replace=TRUE) {
  checkmate::assert_character(covarsVec)
  .normalizedDFs <- lapply(covarsVec,.normalizeDf,data=data)

  # final data frame of normalized covariates
  if(replace) {
    .dat <- Reduce(merge,.normalizedDFs)
    dropnormPrefix <- function(x) {
      colnames(x) <- gsub("normalized_", "", colnames(x))
      x
    }
    catCheck <- intersect(covarsVec,names(Filter(is.factor, data)))
    .dat <- cbind(.dat[ , !names(.dat) %in% covarsVec],subset(.dat,select=catCheck))
    .finalDf <- dropnormPrefix(.dat)
  } else {
    .finalDf <- Reduce(merge,.normalizedDFs)
  }
  .finalDf
}

#' Stratified cross-validation fold generator function, inspired from the caret
#'
#' @param data data frame used in the analysis
#' @param nfold number of k-fold cross validations. Default is 5
#' @param stratVar  Stratification Variable. Default is NULL and ID is used for CV
#' @return return data.frame with the fold column attached
#' @author Vishal Sarsani, caret
#' @export
#'
#' @examples
#' d <- nlmixr2data::theo_sd
#' d$SEX <-0
#' d$SEX[d$ID<=6] <-1
#'
#' covarsVec <- c("WT")
#'
#' # Stratified cross-validation data with CMT
#' df1 <- foldgen(d, nfold=5, stratVar="CMT")
#'
#' # Stratified cross-validation data with ID (individual)
#' df2 <- foldgen(d, nfold=5, stratVar=NULL)
foldgen <-  function(data,nfold=5,stratVar=NULL) {
  # check if data frame
  checkmate::assert_data_frame(data,min.cols = 7)
  # check if user want to stratify on a variable , if not default is on individual
  if(!is.null(stratVar)) {
    checkmate::assertCharacter(stratVar,len = 1,any.missing = FALSE)
    stratCheck <- intersect(names(data), stratVar)
    if(!is.null(stratCheck)) {
      y <- data[,stratCheck]
    } else {
      stop(paste0(stratVar, "not in the data to stratify"),
           call.=FALSE)
    }
  } else {
    # extract ID column from the data frame
    ID <- .idColumn(data)
    # Extract list of individuals
    y <- unique(data[,ID])
  }
  ## Group based on magnitudes and sample within groups
  if(is.numeric(y)) {
    ## Group the numeric data based on their magnitudes
    ## and sample within those groups.

    ## When the number of samples is low, we may have
    ## issues further slicing the numeric data into
    ## groups. The number of groups will depend on the
    ## ratio of the number of folds to the sample size.
    ## At most, we will use quantiles. If the sample
    ## is too small, we just do regular unstratified
    ## CV
    cuts <- floor(length(y)/nfold)
    if(cuts < 2) cuts <- 2
    if(cuts > 5) cuts <- 5
    y <- cut(y,
             unique(quantile(y,
                             probs = seq(0, 1, length = cuts),
                             na.rm=TRUE)),
             include.lowest = TRUE)
  }

  if (nfold < length(y)) {
    ## reset levels so that the possible levels and
    ## the levels in the vector are the same
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))

    ## For each class, balance the fold allocation as far
    ## as possible, then resample the remainder.
    ## The final assignment of folds is also randomized.
    for(i in 1:seq_along(numInClass)) {
      ## create a vector of integers from 1:k as many times as possible without
      ## going over the number of samples in the class. Note that if the number
      ## of samples in a class is less than k, nothing is producd here.
      seqVector <- rep(1:nfold, numInClass[i] %/% nfold)
      ## add enough random integers to get  length(seqVector) == numInClass[i]
      if(numInClass[i] %% nfold > 0) {
        seqVector <- c(seqVector, sample(1:nfold, numInClass[i] %% nfold))
      }
      ## shuffle the integers for fold assignment and assign to this classes's data
      foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
    }
  } else {
    foldVector <- seq(along = y)
  }

  out <- split(seq(along = y), foldVector)
  names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), sep = "")
  out <- foldVector

  if(!is.null(stratVar)) {
    out <- cbind(data,fold=out)
  } else {
    indv <- unique(data[,ID])
    out <- merge(data,cbind(ID=indv,fold=out))
  }
  out
}

#' Sample from uniform distribution by optim
#'
#' @param xvec A vector of min,max values . Ex:c(10,20)
#' @param N Desired number of values
#' @param medValue  Desired Median
#' @param floorT boolean indicating whether to round up
#'
#' @return Samples with approx desired median.
#' @author Vishal Sarsani
#' @export
#'
#' @examples
#' # Simulate 1000 creatine clearance values with median of 71.7 within range of c(6.7,140)
#' creatCl <- optimUnisampling(xvec=c(6.7,140), N=1000, medValue = 71.7, floorT=FALSE)
optimUnisampling <- function(xvec,N=1000,medValue,floorT=TRUE) {
  #Function to calculate distance between sampling median and desired
  fun <- function(xvec, N=1000) {
    xmin <- xvec[1]
    xmax <- xvec[2]
    if (floorT) {
      x <- floor(stats::runif(N, xmin, xmax))
    } else {
      x <- stats::runif(N, xmin, xmax)
    }
    xdist <- (median(x)-medValue)^2
    xdist
  }
  # Optimization
  xr <- stats::optim(xvec, fun)
  xrmin <- xr$par[[1]]
  xrmax <- xr$par[[2]]
  sampled <- stats::runif(N, min = xr$par[[1]], max = xr$par[[2]])
  if (xrmin==xvec[1] && xrmax==xvec[2] && floorT) {
    return(floor(sampled))
  }
  else if (xrmin==xvec[1] && xrmax==xvec[2]) {
    return(sampled)
  }
  optimUnisampling(xvec,N=1000,medValue)
}

#' Format confidence bounds for a variable into bracketed notation using string formatting
#'
#' @param var a list of values for the varaible
#' @param confLower the lower bounds for each of the values
#' @param confUpper the upper bounds for each of the values
#' @param sigdig the number of significant digits
#' @author Vipul Mann
#' @noRd
addConfboundsToVar <- function(var, confLower, confUpper, sigdig = 3) {
  res <- lapply(seq_along(var), function(idx) {
    paste0(
      signif(var[idx], sigdig),
      " (",
      signif(confLower[idx], sigdig),
      ", ",
      signif(confUpper[idx], sigdig),
      ")"
    )
  })
  unlist(res)
}

#' Write bootstrap output files and print a console summary
#'
#' Writes two files to \code{outputDir}:
#' \describe{
#'   \item{\code{bootstrap_results.csv}}{One row per successful replicate with
#'     replicate number and all fixed-effect parameter estimates (both the
#'     log/model-scale Estimate and the Back-transformed value).}
#'   \item{\code{bootstrap_summary.txt}}{Human-readable summary including the
#'     random seed used for resampling, original parameter estimates, and
#'     bootstrap confidence intervals.}
#' }
#' Also prints a concise summary to the console.
#'
#' @param fit nlmixr2 fit object (already augmented with bootstrap results)
#' @param bootSummary result of \code{getBootstrapSummary()}
#' @param outputDir path to the permanent output directory
#' @param fitName model name string
#' @param nboot total number of bootstrap replicates requested
#' @param ci confidence interval level (e.g. 0.95)
#' @return invisible \code{NULL}
#' @noRd
.writeBootstrapOutput <- function(fit, bootSummary, outputDir, fitName, nboot, ci) {

  if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

  # -- 1. Per-replicate results table -----------------------------------------
  # Read each modelsEnsemble_*.rds, extract replicate index from filename,
  # skip failed-replicate placeholders (length == 0).
  mod_files <- list.files(outputDir, pattern = "modelsEnsemble_[0-9]+\\.rds",
                          full.names = FALSE)
  if (length(mod_files) > 0L) {
    mod_idx_raw <- as.integer(gsub("modelsEnsemble_([0-9]+)\\.rds", "\\1", mod_files))
    ord <- order(mod_idx_raw)
    mod_files  <- mod_files[ord]
    mod_indices <- mod_idx_raw[ord]

    row_list <- lapply(seq_along(mod_files), function(i) {
      m <- tryCatch(readRDS(file.path(outputDir, mod_files[i])),
                    error = function(e) list())
      if (length(m) == 0L || is.null(m$parFixedDf)) return(NULL)
      par_df <- m$parFixedDf
      row <- data.frame(replicate = mod_indices[i], stringsAsFactors = FALSE)
      for (par in rownames(par_df)) {
        col_est <- gsub("[^A-Za-z0-9_]", ".", par)
        row[[paste0(col_est, ".Estimate")]]        <- par_df[par, "Estimate"]
        row[[paste0(col_est, ".BackTransformed")]] <- par_df[par, "Back-transformed"]
      }
      row
    })
    row_list <- Filter(Negate(is.null), row_list)

    if (length(row_list) > 0L) {
      results_df <- do.call(rbind, row_list)
      csv_path <- file.path(outputDir, "bootstrap_results.csv")
      utils::write.csv(results_df, csv_path, row.names = FALSE)
    } else {
      results_df <- NULL
      csv_path   <- NULL
    }
    n_ok <- length(row_list)
  } else {
    results_df <- NULL
    csv_path   <- NULL
    n_ok       <- 0L
  }

  # -- 2. Random seed ----------------------------------------------------------
  seed_file <- file.path(outputDir, "bootstrap_seed.rds")
  boot_seed <- if (file.exists(seed_file)) readRDS(seed_file) else NULL

  seed_kind <- if (!is.null(boot_seed)) RNGkind()[1] else "unknown"
  seed_preview <- if (!is.null(boot_seed) && length(boot_seed) >= 8L) {
    paste(boot_seed[seq_len(8L)], collapse = " ")
  } else if (!is.null(boot_seed)) {
    paste(boot_seed, collapse = " ")
  } else {
    "(not recorded)"
  }

  # -- 3. Parameter summary lines ---------------------------------------------
  par_names  <- rownames(fit$parFixedDf)
  orig_est   <- fit$parFixedDf$Estimate
  orig_bt    <- fit$parFixedDf$`Back-transformed`
  boot_mean  <- bootSummary$parFixedDf$mean[, 1]
  boot_se    <- bootSummary$parFixedDf$stdDev[, 1]
  boot_lo    <- bootSummary$parFixedDf$confLower[, 2]  # back-transformed scale
  boot_hi    <- bootSummary$parFixedDf$confUpper[, 2]

  par_summary_df <- data.frame(
    Parameter        = par_names,
    Orig.Estimate    = orig_est,
    Orig.BackTrans   = orig_bt,
    Boot.Mean        = boot_mean,
    Boot.SE          = boot_se,
    Boot.CI.Lower    = boot_lo,
    Boot.CI.Upper    = boot_hi,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  colnames(par_summary_df) <- c(
    "Parameter",
    "Orig.Estimate",
    "Orig.Back-transformed",
    "Boot.Mean",
    "Boot.SE",
    sprintf("Boot.CI.Lower (%.0f%%)", ci * 100),
    sprintf("Boot.CI.Upper (%.0f%%)", ci * 100)
  )

  # -- 4. Omega diagonal summary -----------------------------------------------
  omega_mean <- round(diag(bootSummary$omega$mean), 4)
  omega_lo   <- round(diag(bootSummary$omega$confLower), 4)
  omega_hi   <- round(diag(bootSummary$omega$confUpper), 4)
  omega_df   <- data.frame(
    Parameter    = names(omega_mean),
    Boot.Mean    = omega_mean,
    Boot.CI.Lower = omega_lo,
    Boot.CI.Upper = omega_hi,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  colnames(omega_df) <- c(
    "Parameter", "Boot.Mean",
    sprintf("Boot.CI.Lower (%.0f%%)", ci * 100),
    sprintf("Boot.CI.Upper (%.0f%%)", ci * 100)
  )

  # -- 5. Write bootstrap_summary.txt -----------------------------------------
  txt_path <- file.path(outputDir, "bootstrap_summary.txt")
  .w <- function(...) cat(..., file = txt_path, append = TRUE, sep = "")
  .ln <- function(...) .w(..., "\n")

  if (file.exists(txt_path)) file.remove(txt_path)

  .ln("Bootstrap Summary")
  .ln("=================")
  .ln()
  .ln(sprintf("Run date            : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  .ln(sprintf("Model               : %s", fitName))
  .ln(sprintf("Estimation method   : %s", tryCatch(fit$est, error = function(e) "unknown")))
  .ln(sprintf("Total replicates    : %d", nboot))
  .ln(sprintf("Successful          : %d", n_ok))
  .ln(sprintf("Failed              : %d", nboot - n_ok))
  .ln(sprintf("Confidence interval : %.0f%%", ci * 100))
  .ln(sprintf("Output directory    : %s", outputDir))
  .ln()
  .ln("Random seed (state before resampling)")
  .ln(sprintf("  RNG kind    : %s", seed_kind))
  .ln(sprintf("  Seed[1:8]   : %s", seed_preview))
  .ln("  (full RNG state saved in bootstrap_seed.rds)")
  .ln()
  .ln(strrep("-", 78))
  .ln("Fixed-effect parameter summary")
  .ln(strrep("-", 78))
  .ln(paste(utils::capture.output(print(par_summary_df, row.names = FALSE)),
           collapse = "\n"))
  .ln()
  .ln(strrep("-", 78))
  .ln(sprintf("Omega diagonal (BSV variances) - %.0f%% CI", ci * 100))
  .ln(strrep("-", 78))
  .ln(paste(utils::capture.output(print(omega_df, row.names = FALSE)),
           collapse = "\n"))

  # -- 6. Console summary -----------------------------------------------------
  cli::cli_rule(left = "Bootstrap output")
  cli::cli_inform(c("i" = "Output directory  : {outputDir}"))
  if (!is.null(csv_path)) {
    cli::cli_inform(c("i" = "Results table     : bootstrap_results.csv ({n_ok} replicates, {length(par_names)} parameters)"))
  } else {
    cli::cli_warn(c("!" = "Results table     : (no successful replicates \u2014 CSV not written)"))
  }
  cli::cli_inform(c("i" = "Summary file      : bootstrap_summary.txt"))
  cli::cli_rule(left = "Parameter summary")
  cli::cli_inform(c(
    "i" = "Fixed effects ({sprintf('%.0f', ci * 100)}% CI, back-transformed scale):"
  ))
  print(par_summary_df[, c("Parameter",
                            "Orig.Back-transformed",
                            sprintf("Boot.CI.Lower (%.0f%%)", ci * 100),
                            sprintf("Boot.CI.Upper (%.0f%%)", ci * 100))],
        row.names = FALSE)
  cli::cli_rule()

  invisible(list(seed = boot_seed, results = results_df))
}

#' Bootstrap nlmixr2 fit
#'
#' Bootstrap input dataset and rerun the model to get confidence bounds and aggregated parameters
#'
#' @param fit the nlmixr2 fit object
#' @param nboot an integer giving the number of bootstrapped models to
#'   be fit; default value is 200
#' @param nSampIndiv an integer specifying the number of samples in
#'   each bootstrapped sample; default is the number of unique
#'   subjects in the original dataset
#' @param stratVar Variable in the original dataset to stratify on;
#'   This is useful to distinguish between sparse and full sampling
#'   and other features you may wish to keep distinct in your
#'   bootstrap
#' @param pvalues a vector of pvalues indicating the probability of
#'   each subject to get selected; default value is NULL implying that
#'   probability of each subject is the same
#' @param restart a boolean that indicates if a previous session has
#'   to be restarted; default value is FALSE
#' @param fitName Name of fit to be saved (by default the variable
#'   name supplied to fit)
#' @param stdErrType This gives the standard error type for the
#'   updated standard errors; The current possibilities are: `"perc"`
#'   which gives the standard errors by percentiles (default), `"sd"`
#'   which gives the standard errors by the using the normal
#'   approximation of the mean with standard devaition, or `"se"`
#'   which uses the normal approximation with standard errors
#'   calculated with `nSampIndiv`
#' @param ci Confidence interval level to calculate.  Default is 0.95
#'   for a 95 percent confidence interval
#' @param plotHist A boolean indicating if a histogram plot to assess
#'   how well the bootstrap is doing.  By default this is turned off
#'   (`FALSE`)
#' @param outputDir character; path of the subdirectory used to store all
#'   bootstrap cache files and output tables.  When \code{NULL} (default), the
#'   name is derived automatically as \code{<fitName>_boot_<N>} where \code{N}
#'   is one greater than the number of existing \code{<fitName>_boot_*}
#'   directories in the current working directory (so repeated runs produce
#'   \code{fit_boot_1}, \code{fit_boot_2}, ...).  When \code{restart = FALSE}
#'   and a prior numbered directory exists, the most recently numbered directory
#'   is reused for resuming.  Providing a path explicitly overrides the
#'   automatic scheme.
#'
#' @param updateFit logical; if \code{TRUE} (default) the bootstrap results are
#'   stored on the original \code{fit} object as a single new \code{$bootstrap}
#'   list element without modifying any pre-existing slots.  The list contains:
#'   \code{summary} (augmented \code{parFixedDf}), \code{omegaSummary},
#'   \code{covMatrix}, \code{corMatrix}, \code{bias} (list with
#'   \code{parFixed} and \code{omega} elements), \code{nboot}, \code{ci},
#'   \code{outputDir}, and \code{timestamp}.  Set to \code{FALSE} to leave
#'   \code{fit} completely unchanged.
#'
#' @param returnType this describes the return type
#'
#' - `model` this is the default and returns an updated model boostrapped
#'   confidence intervals and standard errors updated in the estimates.
#'
#' - `fitList` this returns the list of the full fits from the bootstrap.
#'
#' - `modelList` this returns the model list (which abbreviates the
#'    parameters and messages) used for the bootstrap summary
#'    statistics.
#'
#' @param workers integer or \code{NULL}; number of parallel workers to use
#'   when fitting bootstrap replicates.  When \code{NULL} (default), the
#'   current \code{future} plan is used unchanged.  A positive integer
#'   temporarily switches to a \code{multisession} plan with that many workers
#'   for the duration of the bootstrap run.
#'
#' @author Vipul Mann, Matthew Fidler
#' @return A named list with elements:
#' \describe{
#'   \item{\code{seed}}{Integer vector: full RNG state before resampling (also
#'     saved as \file{bootstrap_seed.rds} in \code{outputDir}).}
#'   \item{\code{results}}{Data frame: one row per successful replicate with
#'     columns \code{replicate}, \code{<par>.Estimate}, and
#'     \code{<par>.BackTransformed} for every fixed-effect parameter (also
#'     written as \file{bootstrap_results.csv}).}
#'   \item{\code{summary}}{Data frame: original \code{parFixedDf} augmented
#'     with columns \code{Bootstrap Estimate}, \code{Bootstrap SE},
#'     \code{Bootstrap \%RSE}, \code{Bootstrap CI Lower},
#'     \code{Bootstrap CI Upper}, and \code{Bootstrap Back-transformed}.}
#'   \item{\code{model}}{Character: the \code{fitName} used for the run.}
#'   \item{\code{outputDir}}{Character: absolute path to the output directory
#'     containing all cache files and result tables.}
#'   \item{\code{timestamp}}{POSIXct: date-time when the function returned.}
#' }
#' When \code{updateFit = TRUE}, the same bootstrap results are also attached to
#' the original \code{fit} object as \code{fit$bootstrap}.
#' @export
#' @examples
#' \dontrun{
#'
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45; label("Ka")
#'     tcl <- 1; label("Cl")
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
#' fit <- nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "saem", control = list(print = 0))
#'
#' withr::with_tempdir({ # Run example in temp dir
#'
#' runBootstrap(fit, nboot = 5, restart = TRUE) # overwrites any of the existing data or model files
#' runBootstrap(fit, nboot = 7) # resumes fitting using the stored data and model files
#'
#' # Note this resumes because the total number of bootstrap samples is not 10
#'
#' runBootstrap(fit, nboot=10)
#'
#' # Note the boostrap standard error and variance/covariance matrix is retained.
#' # If you wish to switch back you can change the covariance matrix by
#'
#' nlmixr2est::setCov(fit, "linFim")
#'
#' # And change it back again
#'
#' nlmixr2est::setCov(fit, "boot10")
#'
#' # This change will affect any simulations with uncertainty in their parameters
#'
#' # You may also do a chi-square diagnostic plot check for the bootstrap with
#' bootplot(fit)
#' })
#' }
runBootstrap <- function(fit,
                         nboot = 200,
                         nSampIndiv,
                         stratVar,
                         stdErrType = c("perc", "sd", "se"),
                         ci = 0.95,
                         pvalues = NULL,
                         restart = FALSE,
                         plotHist = FALSE,
                         fitName = as.character(substitute(fit)),
                         outputDir = NULL,
                         updateFit = TRUE,
                         returnType=c("model", "fitList", "modelList"),
                         workers = NULL) {

  stdErrType <- match.arg(stdErrType)
  returnType <- match.arg(returnType)
  checkmate::assertNumeric(ci, lower=0, upper=1, len=1, any.missing=FALSE, null.ok = FALSE)

  if (missing(stratVar)) {
    performStrat <- FALSE
  }
  else {
    if (!(stratVar %in% colnames(nlme::getData(fit)))) {
      cli::cli_alert_danger("{stratVar} not in data")
      stop("aborting ...stratifying variable not in data", call. = FALSE)
    }
    performStrat <- TRUE
  }

  # -- Resolve permanent output directory (SCM-style numbered scheme) ----------
  if (is.null(outputDir)) {
    .fit_nm <- gsub("[^A-Za-z0-9_.]", "_", fitName)
    .boot_pattern <- paste0("^", .fit_nm, "_boot_[0-9]+$")
    .existing <- sort(grep(
      .boot_pattern,
      list.dirs(".", full.names = FALSE, recursive = FALSE),
      value = TRUE
    ))
    if (!restart && length(.existing) > 0L) {
      # Resume: use the most recently numbered directory
      .nums <- as.integer(sub(paste0("^", .fit_nm, "_boot_"), "", .existing))
      outputDir <- .existing[which.max(.nums)]
    } else {
      # Fresh run: increment
      outputDir <- paste0(.fit_nm, "_boot_", length(.existing) + 1L)
    }
  }
  # Keep outputDir as a relative path so that modelBootstrap()'s legacy
  # paste0("./", output_dir, ...) constructions remain valid.
  # normalizePath() would convert to absolute and produce "./C:\..." on Windows.

  if (performStrat) {
    resBootstrap <-
      modelBootstrap(
        fit,
        nboot = nboot,
        nSampIndiv = nSampIndiv,
        stratVar = stratVar,
        pvalues = pvalues,
        restart = restart,
        fitName = fitName,
        outputDir = outputDir,
        workers = workers
      ) # multiple models

    modelsList <- resBootstrap[[1]]
    fitList <- resBootstrap[[2]]
  }
  else {
    resBootstrap <-
      modelBootstrap(
        fit,
        nboot = nboot,
        nSampIndiv = nSampIndiv,
        pvalues = pvalues,
        restart = restart,
        fitName = fitName,
        outputDir = outputDir,
        workers = workers
      ) # multiple models
    modelsList <- resBootstrap[[1]]
    fitList <- resBootstrap[[2]]
  }

  if (returnType == "fitList") {
    return(fitList)
  } else if (returnType == "modelList") {
    return(modelsList)
  }


  bootSummary <-
    getBootstrapSummary(modelsList, ci=ci, stdErrType=stdErrType) # aggregate values/summary

  nrws <- nrow(bootSummary$parFixedDf$mean)

  # Build summary data frame: original parFixedDf + Bootstrap columns
  summaryDf <- fit$parFixedDf

  estEst   <- unname(bootSummary$parFixedDf$mean[seq_len(nrws), 1])
  seBoot   <- unname(bootSummary$parFixedDf$stdDev[seq_len(nrws), 1])
  cLowerBT <- unname(bootSummary$parFixedDf$confLower[seq_len(nrws), 2])
  cUpperBT <- unname(bootSummary$parFixedDf$confUpper[seq_len(nrws), 2])
  estBT    <- unname(bootSummary$parFixedDf$mean[seq_len(nrws), 2])

  summaryDf["Bootstrap Estimate"]        <- estEst
  summaryDf["Bootstrap SE"]              <- seBoot
  summaryDf["Bootstrap %RSE"]            <- seBoot / estEst * 100
  summaryDf["Bootstrap CI Lower"]        <- cLowerBT
  summaryDf["Bootstrap CI Upper"]        <- cUpperBT
  summaryDf["Bootstrap Back-transformed"] <- estBT

  # -- Write permanent output files and print console summary -----------------
  out <- .writeBootstrapOutput(
    fit        = fit,
    bootSummary = bootSummary,
    outputDir  = outputDir,
    fitName    = fitName,
    nboot      = nboot,
    ci         = ci
  )

  abs_output_dir <- normalizePath(outputDir, mustWork = FALSE)

  # -- Optionally attach results to the original fit object -------------------
  # A single new $bootstrap slot is added; no pre-existing slot is modified.
  if (updateFit) {
    bias_par   <- abs(
      data.frame(
        "Estimate"         = fit$parFixedDf$Estimate,
        "Back-transformed" = fit$parFixedDf$`Back-transformed`,
        check.names = FALSE
      ) - bootSummary$parFixedDf$mean
    )
    bias_omega <- abs(fit$omega - bootSummary$omega$mean)

    assign("bootstrap", list(
      summary      = summaryDf,
      omegaSummary = bootSummary$omega,
      covMatrix    = bootSummary$omega$covMatrix,
      corMatrix    = bootSummary$omega$corMatrix,
      bias         = list(parFixed = bias_par, omega = bias_omega),
      nboot        = nboot,
      ci           = ci,
      outputDir    = abs_output_dir,
      timestamp    = Sys.time()
    ), envir = fit$env)
  }

  list(
    seed      = out$seed,
    results   = out$results,
    summary   = summaryDf,
    model     = fitName,
    outputDir = abs_output_dir,
    timestamp = Sys.time()
  )
}

#' Perform bootstrap-sampling from a given dataframe
#'
#' @param data the original dataframe object to sample from for bootstrapping
#'
#' @param nsamp an integer specifying the number of samples in each
#'   bootstrapped sample; default is the number of unique subjects in
#'   the original dataset
#'
#' @param uid_colname a string representing the unique ID of each
#'   subject in the data; default values is 'ID'
#'
#' @param pvalues a vector of pvalues indicating the probability of
#'   each subject to get selected; default value is NULL implying that
#'   probability of each subject is the same
#'
#' @return returns a bootstrap sampled dataframe object
#' @author Vipul Mann, Matthew Fidler
#'
#' @examples
#' sampling(data)
#' sampling(data, 10)
#' @noRd
sampling <- function(data,
                     nsamp=NULL,
                     uid_colname,
                     pvalues = NULL,
                     performStrat = FALSE,
                     stratVar) {
  checkmate::assert_data_frame(data)
  if (is.null(nsamp)) {
    nsamp <- length(unique(data[, uid_colname]))
  }
  else {
    checkmate::assert_integerish(nsamp,
                                 len = 1,
                                 any.missing = FALSE,
                                 lower = 2
                                 )
  }

  if (performStrat && missing(stratVar)) {
    print("stratVar is required for stratifying")
    stop("aborting... stratVar not specified", call. = FALSE)
  }

  checkmate::assert_integerish(nsamp,
                               lower = 2,
                               len = 1,
                               any.missing = FALSE
                               )

  if (missing(uid_colname)) {
    # search the dataframe for a column name of 'ID'
    colNames <- colnames(data)
    colNamesLower <- tolower(colNames)
    if ("id" %in% colNames) {
      uid_colname <- colNames[which("id" %in% colNamesLower)]
    }
    else {
      uid_colname <- "ID"
    }
  }
  else {
    checkmate::assert_character(uid_colname)
  }

  if (performStrat) {
    stratLevels <-
      as.character(unique(data[, stratVar])) # char to access freq. values

    dataSubsets <- lapply(stratLevels, function(x) {
      data[data[, stratVar] == x, ]
    })

    names(dataSubsets) <- stratLevels

    tab <- table(data[stratVar])
    nTab <- sum(tab)

    sampledDataSubsets <- lapply(names(dataSubsets), function(x) {
      dat <- dataSubsets[[x]]

      uids <- unique(dat[, uid_colname])
      uids_samp <- sample(
        list(uids),
        size = ceiling(nsamp * unname(tab[x]) / nTab),
        replace = TRUE,
        prob = pvalues
      )

      sampled_df <-
        data.frame(dat)[0, ] # initialize an empty dataframe with the same col names

      # populate dataframe based on sampled uids
      # new_id = 1
      .env <- environment()
      .env$new_id <- 1
      do.call(rbind, lapply(uids_samp, function(u) {
        data_slice <- dat[dat[, uid_colname] == u, ]
        start <- NROW(sampled_df) + 1
        end <- start + NROW(data_slice) - 1

        data_slice[uid_colname] <-
          .env$new_id # assign a new ID to the sliced dataframe
        .env$new_id <- .env$new_id + 1
        data_slice
      }))
    })
    do.call("rbind", sampledDataSubsets)
  }

  else {
    uids <- unique(data[, uid_colname])
    uids_samp <-
      sample(
        uids,
        size = nsamp,
        replace = TRUE,
        prob = pvalues
      )

    sampled_df <-
      data.frame(data)[0, ] # initialize an empty dataframe with the same col names

    # populate dataframe based on sampled uids
    # new_id = 1
    .env <- environment()
    .env$new_id <- 1

    do.call(rbind, lapply(uids_samp, function(u) {
      data_slice <- data[data[, uid_colname] == u, ]
      start <- NROW(sampled_df) + 1
      end <- start + NROW(data_slice) - 1

      data_slice[uid_colname] <-
        .env$new_id # assign a new ID to the sliced dataframe
      .env$new_id <- .env$new_id + 1
      data_slice
    }))
  }
}

#' Fitting multiple bootstrapped models without aggregaion; called by the function runBootstrap()
#'
#' @param fit the nlmixr2 fit object
#' @param nboot an integer giving the number of bootstrapped models to be fit; default value is 100
#' @param nSampIndiv an integer specifying the number of samples in each bootstrapped sample; default is the number of unique subjects in the original dataset
#' @param pvalues a vector of pvalues indicating the probability of each subject to get selected; default value is NULL implying that probability of each subject is the same
#' @param restart a boolean that indicates if a previous session has to be restarted; default value is FALSE
#'
#' @return a list of lists containing the different attributed of the fit object for each of the bootstrapped models
#' @author Vipul Mann, Matthew Fidler
#' @examples
#' modelBootstrap(fit)
#' modelBootstrap(fit, 5)
#' modelBootstrap(fit, 5, 20)
#' @noRd
modelBootstrap <- function(fit,
                           nboot = 100,
                           nSampIndiv=NULL,
                           stratVar,
                           pvalues = NULL,
                           restart = FALSE,
                           fitName = "fit",
                           outputDir = NULL,
                           workers = NULL) {
  nlmixr2est::assertNlmixrFit(fit)
  if (missing(stratVar)) {
    performStrat <- FALSE
    stratVar <- NULL
  } else {
    performStrat <- TRUE
  }

  data <- nlme::getData(fit)

  .w <- tolower(names(data)) == "id"
  uidCol <- names(data)[.w]

  checkmate::assert_integerish(nboot,
                               len = 1,
                               any.missing = FALSE,
                               lower = 1
                               )

  if (missing(nSampIndiv)) {
    nSampIndiv <- length(unique(data[, uidCol]))
  } else {
    checkmate::assert_integerish(
      nSampIndiv,
      len = 1,
      any.missing = FALSE,
      lower = 2
    )
  }

  # infer the ID column from data
  colNames <- names(data)
  colNamesLower <- tolower(colNames)
  if ("id" %in% colNamesLower) {
    uid_colname <- colNames[which("id" %in% colNamesLower)]
  } else {
    stop("cannot find the 'ID' column! aborting ...", call. = FALSE)
  }

  ui <- fit$finalUiEnv
  fitMeth <- getFitMethod(fit)

  bootData <- vector(mode = "list", length = nboot)

  if (is.null(fit$bootstrapMd5)) {
    bootstrapMd5 <- fit$md5
    assign("bootstrapMd5", bootstrapMd5, envir = fit$env)
  }

  output_dir <- if (!is.null(outputDir)) {
    outputDir
  } else {
    paste0("nlmixr2BootstrapCache_", fitName, "_", fit$bootstrapMd5)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  } else if (dir.exists(output_dir) && restart == TRUE) {
    unlink(output_dir, recursive = TRUE, force = TRUE) # unlink any of the previous directories
    dir.create(output_dir) # create a fresh directory
  }

  fnameBootDataPattern <-
    paste0("boot_data_", "[0-9]+", ".rds",
           sep = ""
           )
  fileExists <-
    list.files(paste0("./", output_dir), pattern = fnameBootDataPattern)

  if (length(fileExists) == 0) {
    restart <- TRUE
  }

  if (!restart) {
    # read saved bootData from boot_data files on disk
    if (length(fileExists) > 0) {
      cli::cli_alert_success("resuming bootstrap data sampling using data at {paste0('./', output_dir)}")

      bootData <- lapply(fileExists, function(x) {
        readRDS(paste0("./", output_dir, "/", x, sep = ""))
      })

      startCtr <- length(bootData) + 1
    } else {
      cli::cli_alert_danger(
        cli::col_red(
          "need the data files at {.file {paste0(getwd(), '/', output_dir)}} to resume"
        )
      )
      stop("aborting...resume file missing", call. = FALSE)
    }
  } else {
    startCtr <- 1
  }

  # Save random seed before any new sampling so it can be reported later.
  # Only write on the first run (or after restart); preserve existing seed on resume.
  seed_file <- file.path(output_dir, "bootstrap_seed.rds")
  if (!file.exists(seed_file)) {
    saveRDS(.Random.seed, seed_file)
  }

  # Generate additional samples (if nboot>startCtr)
  if (nboot >= startCtr) {
    for (mod_idx in startCtr:nboot) {
      bootData[[mod_idx]] <- sampling(
        data,
        nsamp = nSampIndiv,
        uid_colname = uidCol,
        pvalues = pvalues,
        performStrat = performStrat,
        stratVar = stratVar
      )

      # save bootData in curr directory: read the file using readRDS()
      attr(bootData, "randomSeed") <- .Random.seed
      saveRDS(bootData[[mod_idx]],
              file = paste0(
                "./",
                output_dir,
                "/boot_data_",
                mod_idx,
                ".rds"))
    }
  }

  # check if number of samples in stored file is the same as the required number of samples
  fileExists <-
    list.files(paste0("./", output_dir), pattern = fnameBootDataPattern)
  bootData <- lapply(fileExists, function(x) {
    readRDS(paste0("./", output_dir, "/", x, sep = ""))
  })

  currBootData <- length(bootData)

  # Fitting models to bootData now
  fnameModelsEnsemblePattern <-
    paste0("modelsEnsemble_", "[0-9]+",
           ".rds",
           sep = "")
  modFileExists <-
    list.files(paste0("./", output_dir), pattern = fnameModelsEnsemblePattern)

  fnameFitEnsemblePattern <-
    paste0("fitEnsemble_", "[0-9]+",
           ".rds",
           sep = "")
  fitFileExists <- list.files(paste0("./", output_dir), pattern = fnameFitEnsemblePattern)

  # Determine which replicate indices are already complete
  completed_indices <- integer(0)
  if (!restart) {
    if (length(modFileExists) > 0 &&
          (length(fileExists) > 0)) {
      cli::cli_alert_success(
        "resuming bootstrap model fitting using data and models stored at {paste0(getwd(), '/', output_dir)}"
      )
      completed_indices <- as.integer(
        gsub("fitEnsemble_([0-9]+)\\.rds", "\\1", fitFileExists)
      )
      curr_num_models <- length(completed_indices)
      if (curr_num_models > nboot) {
        cli::cli_alert_danger(
          cli::col_red(
            "the model file already has {curr_num_models} models when max models is {nboot}; using only the first {nboot} model(s)"
          )
        )
        # Load and return immediately, capped at nboot
        modelsEnsembleLoaded <- lapply(modFileExists, function(x) {
          readRDS(paste0("./", output_dir, "/", x, sep = ""))
        })
        fitEnsembleLoaded <- lapply(fitFileExists, function(x) {
          readRDS(paste0("./", output_dir, "/", x, sep = ""))
        })
        modelsEnsembleLoaded <- Filter(function(x) length(x) > 0L,
                                       modelsEnsembleLoaded[seq_len(nboot)])
        fitEnsembleLoaded    <- Filter(function(x) !identical(x, NA),
                                       fitEnsembleLoaded[seq_len(nboot)])
        return(list(modelsEnsembleLoaded, fitEnsembleLoaded))
      } else if (curr_num_models == nboot) {
        cli::cli_alert_success(
          "all {nboot} bootstrap models already complete, loading from {paste0(getwd(), '/', output_dir)}"
        )
        modelsEnsembleLoaded <- lapply(modFileExists, function(x) {
          readRDS(paste0("./", output_dir, "/", x, sep = ""))
        })
        fitEnsembleLoaded <- lapply(fitFileExists, function(x) {
          readRDS(paste0("./", output_dir, "/", x, sep = ""))
        })
        modelsEnsembleLoaded <- Filter(function(x) length(x) > 0L, modelsEnsembleLoaded)
        fitEnsembleLoaded    <- Filter(function(x) !identical(x, NA), fitEnsembleLoaded)
        return(list(modelsEnsembleLoaded, fitEnsembleLoaded))
      } else {
        cli::cli_inform("estimating {nboot - curr_num_models} additional model(s) ...")
      }
    } else {
      cli::cli_alert_danger(
        cli::col_red(
          "need both the data and the model files at: {paste0(getwd(), '/', output_dir)} to resume"
        )
      )
      stop(
        "aborting...data and model files missing at: {paste0(getwd(), '/', output_dir)}",
        call. = FALSE
      )
    }
  }

  pending_indices <- setdiff(seq_len(nboot), completed_indices)

  # -- Startup banner ---------------------------------------------------------
  worker_desc <- if (is.null(workers)) {
    "sequential (using current future::plan())"
  } else if (identical(workers, "auto")) {
    if (requireNamespace("future", quietly = TRUE))
      paste0("parallel, auto (", future::availableCores(omit = 1L), " workers)")
    else
      "sequential (future not available)"
  } else if (identical(as.integer(workers), 1L)) {
    "sequential (workers = 1)"
  } else {
    paste0("parallel (", workers, " workers)")
  }
  resume_note <- if (length(completed_indices) > 0L)
    paste0(" (resuming; ", length(completed_indices), " already done)")
  else
    ""

  cli::cli_rule(left = "Bootstrap")
  cli::cli_inform(c(
    "i" = "Model             : {fitName}",
    "i" = "Estimation method : {fitMeth}",
    "i" = "Total replicates  : {nboot}{resume_note}",
    "i" = "Pending replicates: {length(pending_indices)}",
    "i" = "Subjects / sample : {nSampIndiv}",
    "i" = "Stratification    : {if (performStrat) paste0('yes (', stratVar, ')') else 'no'}",
    "i" = "Execution         : {worker_desc}",
    "i" = "Cache directory   : {output_dir}"
  ))
  cli::cli_rule()

  # get control settings for the 'fit' object and save computation effort by not computing the tables
  .ctl <- setQuietFastControl(fit$control)

  .withWorkerPlan(workers, {
    .plap(
      pending_indices,
      function(orig_idx) {
        boot_data <- bootData[[orig_idx]]
        n_subj <- length(unique(boot_data[[uidCol]]))
        cli::cli_rule(left = paste0("Replicate ", orig_idx, "/", nboot))
        cli::cli_inform(c(
          ">" = "Replicate {orig_idx}/{nboot} | method: {fitMeth} | {nrow(boot_data)} rows, {n_subj} subjects"
        ))

        fit_result <- tryCatch(
          {
            fit_r <- suppressWarnings(nlmixr2(ui,
                                              boot_data,
                                              est = fitMeth,
                                              control = .ctl))

            multiple_fits <- list(
              omega = fit_r$omega,
              parFixedDf = fit_r$parFixedDf[, c("Estimate", "Back-transformed")],
              message = fit_r$message,
              warnings = fit_r$warnings)

            list(.failed = FALSE, fit = fit_r, models = multiple_fits)
          },
          error = function(error_message) {
            message("error fitting the model")
            message(error_message)
            message("storing the models as NA ...")
            list(.failed = TRUE,
                 .reason = conditionMessage(error_message),
                 .idx = orig_idx)
          })

        if (!isTRUE(fit_result$.failed)) {
          saveRDS(
            fit_result$models,
            file = paste0(
              "./",
              output_dir,
              "/modelsEnsemble_",
              orig_idx,
              ".rds"))
          saveRDS(
            fit_result$fit,
            file = paste0(
              "./",
              output_dir,
              "/fitEnsemble_",
              orig_idx,
              ".rds"
            )
          )
        } else {
          # Save NA placeholder so the index is treated as complete on resume
          saveRDS(
            list(),
            file = paste0(
              "./",
              output_dir,
              "/modelsEnsemble_",
              orig_idx,
              ".rds"))
          saveRDS(
            NA,
            file = paste0(
              "./",
              output_dir,
              "/fitEnsemble_",
              orig_idx,
              ".rds"
            )
          )
        }

        invisible(NULL)
      },
      .label = function(orig_idx) {
        paste0("replicate ", orig_idx, "/", nboot)
      }
    )
  })

  # Report any failures (fits stored as NA)
  fitFileExists <- list.files(paste0("./", output_dir), pattern = fnameFitEnsemblePattern)
  n_failed <- sum(vapply(fitFileExists, function(f) {
    r <- tryCatch(readRDS(paste0("./", output_dir, "/", f)),
                  error = function(e) NULL)
    identical(r, NA)
  }, logical(1)))
  # Completion summary
  n_ok <- nboot - n_failed
  cli::cli_rule(left = "Bootstrap complete")
  if (n_failed == 0L) {
    cli::cli_inform(c(
      "v" = "All {nboot} replicates fitted successfully.",
      "i" = "Results saved to: {output_dir}"
    ))
  } else {
    cli::cli_warn(c(
      "!" = "{n_failed} of {nboot} bootstrap replicate(s) failed to fit.",
      "i" = "Failed replicates are stored as {.code NA} in the output."
    ))
    cli::cli_inform(c(
      "i" = "{n_ok}/{nboot} replicates succeeded.",
      "i" = "Results saved to: {output_dir}"
    ))
  }
  cli::cli_rule()

  modFileExists <-
    list.files(paste0("./", output_dir), pattern = fnameModelsEnsemblePattern)

  modelsEnsemble <- lapply(modFileExists, function(x) {
    readRDS(paste0("./", output_dir, "/", x, sep = ""))
  })
  # Remove empty-list placeholders saved for failed replicates so that
  # getBootstrapSummary (which uses names(fitList[[1]])) sees a valid entry first.
  modelsEnsemble <- Filter(function(x) length(x) > 0L, modelsEnsemble)

  fitFileExists <- list.files(paste0("./", output_dir), pattern = fnameFitEnsemblePattern)
  fitEnsemble <- lapply(fitFileExists, function(x) {
    readRDS(paste0("./", output_dir, "/", x, sep = ""))
  })
  # Remove NA placeholders for failed replicates from fitEnsemble as well.
  fitEnsemble <- Filter(function(x) !identical(x, NA), fitEnsemble)

  list(modelsEnsemble, fitEnsemble)
}

#' Get the nlmixr2 method used for fitting the model
#'
#' @param fit the nlmixr2 fit object
#'
#' @return returns a string representing the method used by nlmixr2 for fitting the given model
#'
#' @author Vipul Mann, Matthew Fidler
#' @examples
#' getFitMethod(fit)
#' @noRd
getFitMethod <- function(fit) {
  if (!(inherits(fit, "nlmixr2FitCore"))) {
    stop("'fit' needs to be a nlmixr2 fit", call. = FALSE)
  }
  fit$est
}

#' Extract all the relevant variables from a set of bootstrapped models
#'
#' @param fitlist a list of lists containing information on the multiple bootstrapped models; similar to the output of modelsBootstrap() function
#' @param id a character representing the variable of interest: OBJF, AIC, omega, parFixedDf, method, message, warnings
#'
#' @return returns a vector or list across of the variable of interest from all the fits/bootstrapped models
#'
#' @author Vipul Mann, Matthew Fidler
#' @examples
#' extractVars(fitlist, 1) # returns a vector of OBJF values
#' extractVars(fitlist, 4) # returns a list of dataframes containing parFixedDf values
#' @noRd
extractVars <- function(fitlist, id = "method") {
  if (id == "method") {
    # no lapply for 'method'
    unlist(unname(fitlist[[1]][id]))
  }
  else {
    # if id not equal to 'method'
    res <- lapply(fitlist, function(x) {
      x[[id]]
    })


    if (!(id == "omega" ||
            id == "parFixedDf")) {
      # check if all message strings are empty
      if (id == "message") {
        prev <- TRUE
        for (i in length(res)) {
          status <- (res[[i]] == "") && prev
          prev <- status
        }
        if (status == TRUE) {
          ""
        }
        else {
          # if non-empty 'message'
          unlist(res)
        }
      }

      else {
        # if id does not equal 'message'
        unlist(res)
      }
    }
    else {
      # if id equals 'omega' or 'parFixedDf
      res
    }
  }
}

#' Summarize the bootstrapped fits/models
#'
#' @param fitList a list of lists containing information on the multiple bootstrapped models; similar to the output of modelsBootstrap() function
#' @return returns aggregated quantities (mean, median, standard deviation, and variance) as a list for all the quantities
#' @author Vipul Mann, Matthew Fidler
#' @inheritParams runBootstrap
#' @examples
#' getBootstrapSummary(fitlist)
#' @noRd
getBootstrapSummary <- function(fitList,
                                ci = 0.95,
                                stdErrType = "perc") {
  checkmate::assertNumeric(ci, len=1, lower=0, upper=1, any.missing=FALSE, null.ok=FALSE)

  quantLevels <-
    c(0.5, (1 - ci)/2, 1 - (1 - ci)/2) # median, (1-ci)/2, 1-(1-ci)/2

  varIds <-
    names(fitList[[1]]) # number of different variables present in fitlist
  summaryList <- lapply(varIds, function(id) {
    # if (!(id %in% c("omega", "parFixedDf", "method", "message", "warnings"))) {
    #   varVec <- extractVars(fitList, id)
    #   mn <- mean(varVec)
    #   median <- median(varVec)
    #   sd <- sd(varVec)
    #
    #   c(
    #     mean = mn,
    #     median = median,
    #     stdDev = sd
    #   )
    # }
    if (id == "omega") {
      # omega estimates
      omegaMatlist <- extractVars(fitList, id)
      varVec <- simplify2array(omegaMatlist)
      mn <- apply(varVec, 1:2, mean, na.rm=TRUE)
      sd <- apply(varVec, 1:2, sd, na.rm=TRUE)

      quants <- apply(varVec, 1:2, function(x) {
        unname(quantile(x, quantLevels, na.rm=TRUE))
      })
      median <- quants[1, , ]
      confLower <- quants[2, , ]
      confUpper <- quants[3, , ]

      if (stdErrType != "perc") {
        confLower <- mn + qnorm(quantLevels[[2]]) * sd
        confUpper <- mn + qnorm(quantLevels[[3]]) * sd
      }

      # computing the covariance and correlation matrices
      # =======================================================
      parFixedOmegaBootVec <- list()
      parFixedlist <- extractVars(fitList, id = "parFixedDf")
      parFixedlistVec <- lapply(parFixedlist, function(x) {
        ret <- x$Estimate
        if (is.null(names(ret))) {
          ret <- stats::setNames(ret, rownames(x))
        }
        ret
      })
      parFixedlistVec <- do.call("rbind", parFixedlistVec)

      omgVecBoot <- list()
      omegaIdx <- seq_along(omegaMatlist)

      omgVecBoot <- lapply(omegaIdx, function(idx) {
        omgMat <- omegaMatlist[[idx]]
        omgVec <- omgMat[lower.tri(omgMat, TRUE)]
        omgVecBoot[[idx]] <- omgVec
      })
      omgVecBoot <- do.call("rbind", omgVecBoot)

      idxName <- 1
      namesList <- list()
      for (nam1 in colnames(omegaMatlist[[1]])) {
        for (nam2 in colnames(omegaMatlist[[1]])) {
          if (nam1 == nam2) {
            if (!(nam1 %in% namesList)) {
              namesList[idxName] <- nam1
              idxName <- idxName + 1
            }
          } else {
            nam <- paste0("(", nam1, ",", nam2, ")")
            namRev <- paste0("(", nam2, ",", nam1, ")")
            if (!(nam %in% namesList | namRev %in% namesList)) {
              namesList[idxName] <- nam
              idxName <- idxName + 1
            }
          }
        }
      }
      colnames(omgVecBoot) <- namesList

      .w <- which(vapply(namesList, function(x) {
        !all(omgVecBoot[, x] == 0)
      }, logical(1), USE.NAMES=FALSE))
      omgVecBoot <- omgVecBoot[, .w]


      parFixedOmegaCombined <- cbind(parFixedlistVec, omgVecBoot)

      covMatrix <- cov(parFixedOmegaCombined)
      w <- which(diag(covMatrix) == 0)
      if (length(w) > 0) {
        d <- dim(covMatrix)[1]
        corMatrix <- matrix(rep(0,d * d), d, d)
        corMatrix[-w, -w] <- cov2cor(covMatrix[-w, -w])
      } else {
        corMatrix <- cov2cor(covMatrix)
      }
      diag(corMatrix) <- sqrt(diag(covMatrix))
      dimnames(corMatrix) <- dimnames(covMatrix)
      lst <- list(
        mean = mn,
        median = median,
        stdDev = sd,
        confLower = confLower,
        confUpper = confUpper,
        covMatrixCombined = covMatrix,
        corMatrixCombined = corMatrix
      )
    }
    else if (id == "parFixedDf") {
      # parameter estimates (dataframe)
      varVec <- extractVars(fitList, id)
      mn <-
        apply(simplify2array(lapply(varVec, as.matrix)), 1:2, mean, na.rm = TRUE)
      sd <-
        apply(simplify2array(lapply(varVec, as.matrix)), 1:2, sd, na.rm = TRUE)

      quants <-
        apply(simplify2array(lapply(varVec, as.matrix)), 1:2, function(x) {
          unname(quantile(x, quantLevels, na.rm = TRUE))
        })

      median <- quants[1, , ]
      confLower <- quants[2, , ]
      confUpper <- quants[3, , ]

      if (stdErrType != "perc") {
        confLower <- mn + qnorm(quantLevels[[2]]) * sd
        confUpper <- mn + qnorm(quantLevels[[3]]) * sd
      }

      lst <- list(
        mean = mn,
        median = median,
        stdDev = sd,
        confLower = confLower,
        confUpper = confUpper
      )
    }

    else {
      # if id equals method, message, or warning
      extractVars(fitList, id)
    }
  })

  names(summaryList) <- varIds
  summaryList$nboot <- length(fitList)
  summaryList$warnings <- unique(summaryList$warnings)
  summaryList$message <- unique(summaryList$message)
  class(summaryList) <- "nlmixr2BoostrapSummary"
  summaryList
}

#' @export
print.nlmixr2BoostrapSummary <- function(x, ..., sigdig = NULL) {
  if (is.null(sigdig)) {
    if (any(names(x) == "sigdig")) {
      sigdig <- x$sigdig
    } else {
      sigdig <- 3
    }
  }

  # objf <- x$objf
  # aic <- x$aic
  message <- x$message
  warnings <- x$warnings

  omega <- x$omega
  parFixedDf <- x$parFixedDf

  nboot <- x$nboot
  cli::cli_h1(
    cli::col_red(
      "Summary of the bootstrap models (nboot: {nboot})"
    )
  )
  cli::cli_li(cli::col_magenta(
    cli::style_bold(
      "Omega matrices: mean, median, standard deviation, and confidence bousnds"
    ),
    cli::col_yellow(" (summary$omega)")
  ))

  lapply(seq_along(omega), function(x) {
    cli::cli_text(cli::col_green(paste0("$", names(omega)[x])))
    print(signif(omega[[x]], sigdig))
  })

  cli::cli_li(cli::col_magenta(
    cli::style_bold(
      "Estimated parameters: mean, median, standard deviation, and confidence bounds"
    ),
    cli::col_yellow(" (summary$parFixedDf)")
  ))
  lapply(seq_along(parFixedDf), function(x) {
    cli::cli_text(cli::col_yellow(paste0("$", names(parFixedDf)[x])))
    print(signif(parFixedDf[[x]], sigdig))
  })

  cli::cli_li(cli::cli_text(
    cli::bg_yellow(cli::style_bold("Messages")),
    cli::col_yellow(" (summary$message)")
  ))
  print(message)

  cli::cli_li(cli::cli_text(
    cli::bg_red(cli::style_bold(cli::col_white("Warnings"))),
    cli::col_yellow(" (summary$warnings)")
  ))
  print(warnings)

  cli::cli_h1("end")
  invisible(x)
}

#' Assign a set of variables to the nlmixr2 fit environment
#'
#' @param namedVars a named list of variables that need to be assigned to the given environment
#' @param fitobject the nlmixr2 fit object that contains its environment information
#' @noRd
#'
assignToEnv <- function(namedVars, fitobject) {
  if (!inherits(fitobject, "nlmixr2FitCore")) {
    stop("'fit' needs to be a nlmixr2 fit", call. = FALSE)
  }

  if (is.null(names(namedVars))) {
    stop("'namedVars needs to be a named list", call. = FALSE)
  }

  if (length(namedVars) != length(names(namedVars))) {
    stop("'namedVars does not have all the elements named", call. = FALSE)
  }
  env <- fitobject$env
  lapply(names(namedVars), function(x) {
    assign(x, namedVars[[x]], envir = env)
  })
}

#' @title Produce delta objective function for boostrap
#'
#' @param x fit object
#' @param ... other parameters
#' @return Fit traceplot or nothing.
#' @author Vipul Mann,  Matthew L. Fidler
#' @references
#' R Niebecker,  MO Karlsson. (2013)
#' *Are datasets for NLME models large enough for a bootstrap to provide reliable parameter uncertainty distributions?*
#' PAGE 2013.
#' <https://www.page-meeting.org/?abstract=2899>
#' @export
bootplot <- function(x, ...) {
  UseMethod("bootplot")
}

#' @rdname bootplot
#' @export
#' @importFrom ggplot2 .data
bootplot.nlmixr2FitCore <- function(x, ...) {
  .fitName <- as.character(substitute(x))
  if (inherits(x, "nlmixr2FitCore")) {
    if (exists("bootSummary", x$env) & (!exists(".bootPlotData", x$env))) {
      runBootstrap(x, x$bootSummary$nboot, plotHist = TRUE, fitName = .fitName)
    }
    if (exists(".bootPlotData", x$env)) {
      if (x$bootSummary$nboot != x$env$.bootPlotData$deltaN) {
        runBootstrap(x, x$bootSummary$nboot, plotHist = TRUE, fitName = .fitName)
      }
      .chisq <- x$env$.bootPlotData$chisq
      .dfD <- x$env$.bootPlotData$dfD
      .deltaN <- x$env$.bootPlotData$deltaN
      .df2 <- x$env$.bootPlotData$df2
      .plot <- ggplot2::ggplot(.chisq, ggplot2::aes(.data$quantiles, .data$deltaofv, color = .data$Distribution)) +
        ggplot2::geom_line() +
        ggplot2::ylab("\u0394 objective function") +
        ggplot2::geom_text(data = .dfD, ggplot2::aes(label = .data$label), hjust = 0) +
        ggplot2::xlab("Distribution quantiles") +
        ggplot2::scale_color_manual(values = c("red", "blue")) +
        rxode2::rxTheme() +
        ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")

      if (requireNamespace("ggtext", quietly = TRUE)) {
        .plot <- .plot +
          ggplot2::theme(
            plot.title = ggtext::element_markdown(),
            legend.position = "none"
          ) +
          ggplot2::labs(
            title = paste0(
              'Bootstrap <span style="color:blue; opacity: 0.2;">\u0394 objective function (', .deltaN,
              " models, df\u2248", .df2, ')</span> vs <span style="color:red; opacity: 0.2;">reference \u03C7\u00B2(df=',
              length(x$ini$est), ")</style>"
            ),
            caption = "\u0394 objective function curve should be on or below the reference distribution curve"
          )
      } else {
        .plot <- ggplot2::labs(
          title = paste0("Distribution of \u0394 objective function values for ", .deltaN, " df=", .df2, " models"),
          caption = "\u0394 objective function curve should be on or below the reference distribution curve"
        )
      }
      .plot
    } else {
      stop("this nlmixr2 object does not include boostrap distribution statics for comparison",
           call. = FALSE
           )
    }
  } else {
    stop("this is not a nlmixr2 object",
         call. = FALSE
         )
  }
}

#' Bootstrap nlmixr2 fit
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `bootstrapFit()` has been renamed to [runBootstrap()]. This alias is provided
#' for backward compatibility and will be removed in a future release.
#'
#' @inheritParams runBootstrap
#' @inherit runBootstrap return
#' @export
#' @keywords internal
bootstrapFit <- function(...) {
  .Deprecated("runBootstrap", package = "nlmixr2extra",
              msg = paste0("'bootstrapFit()' has been renamed to 'runBootstrap()'. ",
                           "Please update your code."))
  runBootstrap(...)
}
