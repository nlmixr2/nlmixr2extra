
#' Built-in covariate shape expression builders for SCM
#'
#' A named list of functions \code{function(col, center, level)} that return the
#' covariate expression string embedded inside the model body.  Each function is
#' multiplied by the covariate theta: \code{cov_theta * covExpr}.
#'
#' \describe{
#'   \item{power}{  \code{log(col/center)} — power/allometric scaling;
#'     on the exponentiated scale this gives
#'     \code{baseline * (col/center)^theta}.}
#'   \item{lin}{    \code{(col - center)} — linear on the log scale
#'     (exponential in natural scale);
#'     \code{baseline * exp(theta * (col - center))}.}
#'   \item{log}{    \code{log(col)} — uncentered log transform.}
#'   \item{identity}{\code{col} — raw covariate on the log scale.}
#'   \item{cat}{    \code{col_level} — pre-computed 0/1 indicator column
#'     (e.g. \code{sex_male}).  The column must already exist in the dataset,
#'     typically created by \code{addCatCovariates()}.  The theta multiplies
#'     this indicator directly on the log scale, giving mu-referencing
#'     compatible code: \code{log(par) = theta_pop + eta + theta_cat * indicator}.}
#' }
#'
#' Users may pass additional shapes via the \code{customShapes} argument of
#' \code{runSCM()}.  Custom shape functions must accept named arguments
#' \code{col}, \code{center}, and \code{level}.
#'
#' @noRd
.SCM_SHAPES <- list(
  power    = function(col, center, level = NULL) paste0("log(", col, "/", center, ")"),
  lin      = function(col, center, level = NULL) paste0("(", col, " - ", center, ")"),
  log      = function(col, center, level = NULL) paste0("log(", col, ")"),
  identity = function(col, center, level = NULL) col,
  cat      = function(col, center = NULL, level) {
    # Return the name of the pre-computed 0/1 indicator column (e.g. "sex_male").
    # The column must exist in the dataset before model fitting, typically
    # created by addCatCovariates().  This keeps the model body mu-referencing
    # compatible with log-linear parameterisation.
    paste0(col, "_", level)
  }
)

#' Parse an init specification into a canonical list(est, lower, upper)
#'
#' Accepts either a scalar numeric (lower/upper default to -5/5) or a named
#' list with elements \code{est} (or \code{init}), \code{lower}, and
#' \code{upper}.  \code{NULL} is treated as the null-effect default
#' (est=0.1, lower=-5, upper=5).
#' @param spec scalar, named list, or NULL
#' @return list(est, lower, upper)
#' @noRd
.parseInitSpec <- function(spec) {
  if (is.null(spec))
    return(list(est = 0.1, lower = -5, upper = 5))
  if (is.numeric(spec) && length(spec) == 1L)
    return(list(est = spec, lower = -5, upper = 5))
  if (is.list(spec)) {
    est   <- if (!is.null(spec[["est"]])) spec[["est"]] else
      if (!is.null(spec[["init"]])) spec[["init"]] else 0.1
    lower <- if (!is.null(spec[["lower"]])) spec[["lower"]] else -5
    upper <- if (!is.null(spec[["upper"]])) spec[["upper"]] else  5
    return(list(est = est, lower = lower, upper = upper))
  }
  list(est = 0.1, lower = -5, upper = 5)
}

#' Build the covariate expression string for a given shape name
#'
#' @param shape  character; key into \code{.SCM_SHAPES} or \code{customShapes}
#' @param col    character; data column name (e.g. \code{"wt"}, \code{"sex"})
#' @param center numeric or NULL; centering value for continuous shapes
#' @param level  character or NULL; level for categorical shapes
#' @param customShapes named list of additional shape builder functions; merged
#'   with \code{.SCM_SHAPES} (custom entries take precedence)
#' @return character expression string
#' @noRd
.scmCovExpr <- function(shape, col, center = NULL, level = NULL, customShapes = NULL) {
  all_shapes <- c(customShapes, .SCM_SHAPES)  # custom takes precedence
  fn <- all_shapes[[shape]]
  if (is.null(fn)) {
    stop("Unknown SCM shape '", shape, "'. Available: ",
         paste(names(all_shapes), collapse = ", "), call. = FALSE)
  }
  fn(col = col, center = center, level = level)
}

#' Rebuild a UI by applying a set of covariate relationships to a clean base UI
#'
#' Starting from \code{base_ui} (always the original model without any
#' covariates), applies each row of \code{pairs_df} in sequence via
#' \code{.builduiCovariate()}.  This avoids the tainted-UI problem that arises
#' when trying to add a second covariate to a model that already has one.
#'
#' @param base_ui  the clean base rxode2 UI (no covariates added)
#' @param pairs_df data frame with columns \code{var}, \code{covar}, and
#'   optionally \code{covExpr} and \code{init}
#' @return updated UI with all covariate relationships applied, or
#'   \code{base_ui} unchanged when \code{pairs_df} is \code{NULL} or empty
#' @noRd
.rebuildUiFromPairs <- function(base_ui, pairs_df) {
  if (is.null(pairs_df) || nrow(pairs_df) == 0L) return(base_ui)

  base_ui <- rxode2::rxUiDecompress(base_ui)

  # Single-pass expansion from the clean base model.  Calling .builduiCovariate()
  # in a loop fails because the first call produces a tainted UI (tv moves from
  # pureMuRef to taintMuRef), so the second call can no longer find tv in
  # muRefDataFrame and produces a malformed "+".  Here we build the full covDf
  # for ALL pairs at once and process the base model body in one pass.
  cov_rows <- lapply(seq_len(nrow(pairs_df)), function(i) {
    var_i   <- pairs_df$var[i]
    covar_i <- pairs_df$covar[i]
    expr_i  <- {
      cv <- if ("covExpr" %in% names(pairs_df)) pairs_df$covExpr[i] else NA_character_
      if (is.null(cv) || is.na(cv)) covar_i else cv
    }
    theta_i <- .getThetaName(base_ui, var_i)
    data.frame(
      theta              = theta_i,
      covariate          = expr_i,
      covariateParameter = paste0("cov_", covar_i, "_", var_i),
      stringsAsFactors   = FALSE
    )
  })
  .covDf <- do.call(rbind, cov_rows)

  .murefDf <- base_ui$muRefDataFrame
  .split   <- base_ui$getSplitMuModel
  .pars    <- c(names(.split$pureMuRef), names(.split$taintMuRef))
  .model   <- nlmixr2est::.saemDropMuRefFromModel(base_ui)

  .model_list <- lapply(.model, .expandRefMu,
                        murefDf = .murefDf, covDf = .covDf, pars = .pars)
  .model_list <- .injectCovBlocks(.model_list, .covDf, .murefDf)

  .newModel <- eval(parse(text = paste0(
    "quote(model({",
    paste0(as.character(.model_list), collapse = "\n"),
    "}))"
  )))

  .iniDf       <- base_ui$iniDf
  nthetaLength <- length(which(!is.na(.iniDf$ntheta)))
  new_ini_rows <- lapply(seq_len(nrow(.covDf)), function(i) {
    init_i <- if ("init" %in% names(pairs_df)) pairs_df$init[i] else NA_real_
    if (is.null(init_i) || is.na(init_i)) init_i <- 0.1
    lower_i <- if ("lower" %in% names(pairs_df)) pairs_df$lower[i] else NA_real_
    if (is.null(lower_i) || is.na(lower_i)) lower_i <- -5
    upper_i <- if ("upper" %in% names(pairs_df)) pairs_df$upper[i] else NA_real_
    if (is.null(upper_i) || is.na(upper_i)) upper_i <- 5
    data.frame(
      ntheta        = as.integer(nthetaLength + i),
      neta1         = NA_character_,
      neta2         = NA_character_,
      name          = .covDf$covariateParameter[i],
      lower         = lower_i,
      est           = init_i,
      upper         = upper_i,
      fix           = FALSE,
      label         = NA_character_,
      backTransform = NA_character_,
      condition     = NA_character_,
      err           = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  .ini <- rbind(.iniDf, do.call(rbind, new_ini_rows))
  .ini <- as.expression(lotri::as.lotri(.ini))
  .ini[[1]] <- quote(`ini`)

  rxode2::rxUiDecompress(.getUiFunFromIniAndModel(base_ui, .ini, .newModel)())
}

#' Add covariate
#'
#'
#' @param ui compiled rxode2 nlmir2 model or fit
#' @param varName  the variable name to which the given covariate is to be added
#' @param covariate the covariate that needs string to be constructed
#' @param add  boolean indicating if the covariate needs to be added or removed.
#' @param covExpr optional pre-built covariate expression string.  When
#'   \code{NULL} (default) the expression is derived from the covariate name.
#' @author Matthew Fidler, Vishal Sarsani
#' @export

addorremoveCovariate <- function(ui, varName, covariate, add = TRUE, covExpr = NULL) {

  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$ui
  }
  ui <- rxode2::assertRxUi(ui)


  checkmate::assertCharacter(varName, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(covariate, len = 1, any.missing = FALSE)

  if (inherits(try(str2lang(varName)), "try-error")) {
    stop("`varName` must be a valid R expression", call. = FALSE)
  }
  if (inherits(try(str2lang(covariate)), "try-error")) {
    stop("`varName` must be a valid R expression", call. = FALSE)
  }
  .pop <- .getThetaName(ui, varName = varName)
  .cov <- paste0("cov_", covariate, "_", varName)
  # covExpr is the expression embedded in the model (may differ from the
  # covariate name used for theta naming, e.g. "log(wt/70)" or
  # 'ifelse(sex=="male",1,0)').  Defaults to the raw covariate name.
  .expr <- if (!is.null(covExpr)) covExpr else covariate

  if (add) {
    .covdf <- rbind(ui$muRefCovariateDataFrame,
                    data.frame(theta = .pop, covariate = .expr, covariateParameter = .cov))
  } else {
    .covdf <- ui$muRefCovariateDataFrame[ui$muRefCovariateDataFrame$covariateParameter != .cov, ]
  }
  .split <- ui$getSplitMuModel
  .pars <- c(names(.split$pureMuRef), names(.split$taintMuRef))
  .model <- nlmixr2est::.saemDropMuRefFromModel(ui)
  .model_list <- lapply(.model, .expandRefMu, murefDf = ui$muRefDataFrame, covDf = .covdf, pars = .pars)
  if (nrow(.covdf) > 0L) {
    .model_list <- .injectCovBlocks(.model_list, .covdf, ui$muRefDataFrame)
  }
  .model_list
}




#' Given model expression, Expand population expression
#'
#' @param x model expression
#' @param murefDf MU referencing data frame
#' @param covDf covariate referencing data frame
#' @param pars theta parameters
#'
#' @noRd
#' @author Matthew Fidler, Vishal Sarsani
#' @noRd
.expandRefMu <- function(x, murefDf, covDf, pars) {
  if (is.name(x)) {
    currparam <- as.character(x)
    if (currparam %in% pars) {
      return(str2lang(.expandPopExpr(currparam, murefDf, covDf)))
    }
  } else if (is.call(x)) {
    return(as.call(c(list(x[[1]]), lapply(x[-1], .expandRefMu, murefDf = murefDf, covDf = covDf, pars = pars))))
  }
  x
}



#' Derive the PK parameter name from the population theta name
#'
#' @param popParam character; population theta name (e.g. \code{"tv"})
#' @param murefDf  the muRefDataFrame from the UI (columns \code{theta}, \code{eta})
#' @return character; PK parameter name (e.g. \code{"v"})
#' @noRd
.covVarNameFromTheta <- function(popParam, murefDf) {
  eta <- murefDf$eta[murefDf$theta == popParam]
  if (length(eta) == 0L) return(popParam)
  sub("^eta\\.", "", eta[[1L]])
}

#' Derive the intermediate variable name from a covariateParameter string
#'
#' Strips the \code{"cov_"} prefix so that e.g. \code{"cov_wt_power_v"} becomes
#' \code{"wt_power_v"} and \code{"cov_sex_male_v"} becomes \code{"sex_male_v"}.
#'
#' @param covariateParameter character; theta name from \code{muRefCovariateDataFrame}
#' @return character; intermediate variable name
#' @noRd
.covIntermediateName <- function(covariateParameter) {
  sub("^cov_", "", covariateParameter)
}

#' Build the covariate block lines for a single parameter
#'
#' Generates one intermediate assignment per covariate and a final
#' \code{cov_<var>} aggregate line.  Example output for weight (power) and
#' sex on \code{v}:
#' \preformatted{
#'   wt_power_v = log(wt/70.5) * cov_wt_power_v
#'   sex_male_v = sex_male * cov_sex_male_v
#'   cov_v = wt_power_v + sex_male_v
#' }
#'
#' @param popParam character; population theta name (e.g. \code{"tv"})
#' @param murefDf  muRefDataFrame
#' @param covDf    full covariate reference data frame
#' @param factor   optional scalar multiplied into each covariate theta (LASSO)
#' @return character vector of lines, or \code{character(0)} when no covariates
#' @noRd
.buildCovBlock <- function(popParam, murefDf, covDf, factor = NULL) {
  varName <- .covVarNameFromTheta(popParam, murefDf)
  rows    <- covDf[covDf$theta == popParam, , drop = FALSE]
  if (nrow(rows) == 0L) return(character(0L))

  int_names <- .covIntermediateName(rows$covariateParameter)

  int_lines <- if (!is.null(factor)) {
    paste0(int_names, " = ", rows$covariate, " * ", rows$covariateParameter, " * ", factor)
  } else {
    paste0(int_names, " = ", rows$covariate, " * ", rows$covariateParameter)
  }

  agg_line <- paste0("cov_", varName, " = ", paste(int_names, collapse = " + "))
  c(int_lines, agg_line)
}

#' Inject covariate block lines into the model expression list
#'
#' For each PK parameter that has entries in \code{covDf}, inserts the
#' intermediate variable assignments and the \code{cov_<var>} aggregate line
#' immediately before that parameter's assignment line.
#'
#' @param model_list list of R language objects (model body statements)
#' @param covDf      covariate reference data frame
#' @param murefDf    muRefDataFrame
#' @param factor     optional LASSO constraint factor
#' @return modified \code{model_list}
#' @noRd
.injectCovBlocks <- function(model_list, covDf, murefDf, factor = NULL) {
  params_with_covs <- unique(covDf$theta)
  strs <- vapply(model_list,
                 function(x) paste(deparse(x, width.cutoff = 500L), collapse = " "),
                 character(1L))

  # Collect (insertion_index, block_exprs) pairs; insert in reverse order so
  # earlier insertions do not shift the indices of later ones.
  insertions <- vector("list", length(params_with_covs))

  for (k in seq_along(params_with_covs)) {
    popParam <- params_with_covs[[k]]
    varName  <- .covVarNameFromTheta(popParam, murefDf)
    cov_var  <- paste0("cov_", varName)

    block <- .buildCovBlock(popParam, murefDf, covDf, factor)
    if (length(block) == 0L) next

    idx <- which(grepl(cov_var, strs, fixed = TRUE))
    if (length(idx) == 0L) next
    idx <- idx[[1L]]

    block_exprs <- lapply(block, function(ln) {
      tryCatch(str2lang(ln), error = function(e) parse(text = ln)[[1L]])
    })

    insertions[[k]] <- list(idx = idx, exprs = block_exprs)
  }

  insertions <- Filter(Negate(is.null), insertions)
  if (length(insertions) == 0L) return(model_list)

  # Sort descending so we insert from the bottom up
  ord <- order(vapply(insertions, `[[`, integer(1L), "idx"), decreasing = TRUE)
  insertions <- insertions[ord]

  for (ins in insertions) {
    idx  <- ins$idx
    exprs <- ins$exprs
    model_list <- c(
      model_list[seq_len(idx - 1L)],
      exprs,
      model_list[idx:length(model_list)]
    )
  }

  model_list
}

#' Expand population expression given the mu reference and covariate reference Data frames
#'
#' When covariates are present the reference is replaced with
#' \code{theta+eta+cov_<var>}; the detailed covariate expressions are
#' generated as separate lines by \code{.buildCovBlock()} /
#' \code{.injectCovBlocks()} in \code{addorremoveCovariate()}.
#'
#' @param popParam  population parameter variable
#' @param murefDf MU referencing data frame
#' @param covDf covariate referencing data frame
#' @param factor unused; kept for signature compatibility
#'
#' @return expanded expression string
#' @author Matthew Fidler,Vishal Sarsani
#' @noRd
.expandPopExpr <- function(popParam, murefDf, covDf, factor = NULL) {
  .par1 <- murefDf[murefDf$theta == popParam, ]
  .res  <- paste0(.par1$theta, "+", .par1$eta)
  .w    <- which(covDf$theta == popParam)
  if (length(.w) > 0L) {
    varName <- .covVarNameFromTheta(popParam, murefDf)
    .res <- paste0(.res, "+cov_", varName)
  }
  .res
}

#' get the population parameter from variable name
#'
#' @param ui compiled rxode2 nlmir2 model or fi
#' @param varName the variable name to which the given covariate is to be added
#'
#' @return population parameter variable
#'
#' @author Matthew Fidler
#' @noRd
.getThetaName <- function(ui, varName) {
  .split <- ui$getSplitMuModel
  if (varName %in% names(.split$pureMuRef)) {
    return(varName)
  }
  .w <- which(.split$pureMuRef == varName)
  if (length(.w) == 1) {
    return(names(.split$pureMuRef)[.w])
  }
  if (varName %in% names(.split$taintMuRef)) {
    return(varName)
  }
  # taintMuRef values are PK param names (e.g. "v") whose theta has a covariate;
  # the original check only looked at names (theta names), so we also check values
  .w2 <- which(.split$taintMuRef == varName)
  if (length(.w2) == 1) {
    return(names(.split$taintMuRef)[.w2])
  }
  # Fallback: use muRefDataFrame which maps theta names to eta names (e.g. tv -> eta.v).
  # This works even after a covariate has been added (tv moves from pureMuRef to
  # taintMuRef and taintMuRef values may not be plain PK param name strings).
  .mdf <- ui$muRefDataFrame
  .eta_name <- paste0("eta.", varName)
  .w3 <- which(.mdf$eta == .eta_name)
  if (length(.w3) >= 1) {
    return(.mdf$theta[.w3[1]])
  }
  stop("'", varName, "'", "has not been found in the model ui", call. = FALSE)
}


#' Given a data frame extract column corresponding to  Individual
#'
#' @param data given data frame
#'
#' @return column name of individual
#' @noRd
#' @author  Vishal Sarsani
.idColumn <- function(data) {
  #check if it is a dataframe
  checkmate::assertDataFrame(data, col.names = "named")
  # Extract individual ID from column names
  colNames <- colnames(data)
  colNamesLower <- tolower(colNames)
  if ("id" %in% colNamesLower) {
    uidCol <- colNames[which("id" %in% colNamesLower)]
  } else {
    uidCol <- "ID"
  }
  uidCol
}

#' Build ui from the covariate
#'
#' @import lotri
#' @param ui compiled rxode2 nlmir2 model or fit
#' @param varName  the variable name to which the given covariate is to be added
#' @param covariate the covariate name used for theta naming (e.g. \code{"wt"},
#'   \code{"sex_male"})
#' @param add boolean indicating if the covariate needs to be added or removed
#' @param type \code{"continuous"} (default) or \code{"categorical"}.  Controls
#'   how the model expression is generated.
#' @param center for continuous covariates, the centering value (typically the
#'   population median).  When supplied the model expression becomes
#'   \code{cov_theta * log(raw_col / center)}.
#' @param raw_col the data column name to use in the model expression.  Defaults
#'   to \code{covariate}.  For categorical covariates this is the original
#'   factor/character column (e.g. \code{"sex"}), while \code{covariate} carries
#'   the level-suffixed name used for the theta (e.g. \code{"sex_male"}).
#' @param level for categorical covariates, the level being tested (e.g.
#'   \code{"male"}).  The model expression becomes
#'   \code{cov_theta * ifelse(raw_col == "level", 1, 0)}.
#' @param covExpr optional pre-computed covariate expression string.  When
#'   supplied, the \code{type}/\code{center}/\code{level} shape logic is
#'   skipped and this string is used verbatim as the expression multiplied by
#'   the covariate theta.  Intended for use by \code{.fitCandidatePairs()} after
#'   \code{.expandShapes()} has already determined the expression.
#' @param init numeric starting estimate for the covariate theta; default
#'   \code{0} (no covariate effect).  Set via \code{.expandShapes()} from the
#'   \code{inits} argument of \code{runSCM()}.
#' @return ui with added or removed covariate
#' @noRd
#' @author  Vishal Sarsani
.builduiCovariate <- function(ui, varName, covariate, add = TRUE,
                              type = "continuous", center = NULL,
                              raw_col = NULL, level = NULL,
                              covExpr = NULL, init = 0.1,
                              lower = -5, upper = 5) {
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$finalUiEnv
  }
  ui <- rxode2::assertRxUi(ui)
  ui <- rxode2::rxUiDecompress(ui)

  checkmate::assertCharacter(varName,   len = 1, any.missing = FALSE)
  checkmate::assertCharacter(covariate, len = 1, any.missing = FALSE)

  if (inherits(try(str2lang(varName)),   "try-error")) {
    stop("`varName` must be a valid R expression",   call. = FALSE)
  }
  if (inherits(try(str2lang(covariate)), "try-error")) {
    stop("`covariate` must be a valid R expression", call. = FALSE)
  }

  # The data column used in the model expression (may differ from covariate for
  # categorical: covariate = "sex_male", raw_col = "sex", level = "male")
  rc <- if (!is.null(raw_col)) raw_col else covariate

  # Build the covariate expression that appears in the model body.
  # For continuous (centred) covariates the expression is log(col / center).
  # For categorical covariates the expression is the pre-computed 0/1 indicator
  # column name, which must already exist in the data.
  # For plain (uncentred) the raw column name is used as the expression.
  # When covExpr is supplied by the caller (e.g. from .expandShapes()), use it
  # directly rather than recomputing from type/center/level.
  if (is.null(covExpr)) {
    covExpr <- if (type == "categorical" && !is.null(level)) {
      paste0(rc, "_", level)
    } else if (type == "continuous" && !is.null(center)) {
      paste0("log(", rc, "/", center, ")")
    } else {
      rc
    }
  }

  covName <- paste0("cov_", covariate, "_", varName)

  # add covariate
  if (add) {
    lst <- addorremoveCovariate(ui, varName, covariate, add = TRUE, covExpr = covExpr)
    .newModel <- eval(parse(text = paste0("quote(model({", paste0(as.character(lst), collapse = "\n"), "}))")))
    nthetaLength <- length(which(!is.na(ui$iniDf$ntheta)))
    .ini <- ui$iniDf
    .ini <- rbind(.ini, data.frame(ntheta = as.integer(nthetaLength + 1), neta1 = NA_character_, neta2 = NA_character_,
                                   name = covName, lower = lower, est = init, upper = upper, fix = FALSE, label = NA_character_,
                                   backTransform = NA_character_, condition = NA_character_, err = NA_character_))
  } else {
    # remove covariate
    lst <- addorremoveCovariate(ui, varName, covariate, add = FALSE)
    .newModel <- eval(parse(text = paste0("quote(model({", paste0(as.character(lst), collapse = "\n"), "}))")))
    .ini <- ui$iniDf[ui$iniDf$name != covName, ]
  }

  # build ui
  .ini <- as.expression(lotri::as.lotri(.ini))
  .ini[[1]] <- quote(`ini`)
  rxode2::rxUiDecompress(.getUiFunFromIniAndModel(ui, .ini, .newModel)())
}

#' Build covInfo list from varsVec and covarsVec
#'
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec character vector of covariates that need to be added
#' @return covInfo list of covariate info
#' @author  Vishal Sarsani
#' @export
buildcovInfo <- function(varsVec, covarsVec) {
  checkmate::assert_character(varsVec, min.len = 1)
  checkmate::assert_character(covarsVec, min.len = 1)
  possiblePerms <- expand.grid(varsVec, covarsVec)
  possiblePerms <-
    list(
      as.character(possiblePerms[[1]]),
      as.character(possiblePerms[[2]])
    )
  names(possiblePerms) <- c("vars", "covars")
  covInfo <- list() # reversivle listVarName!!
  for (item in Map(list, possiblePerms$vars, possiblePerms$covars)) {
    listVarName <- paste0(item[[2]], item[[1]])
    covInfo[[listVarName]] <- list(varName = item[[1]], covariate = item[[2]])
  }
  covInfo
}



#' Build updated from the covariate and variable vector list
#'
#'
#' @param ui compiled rxode2 nlmir2 model or fit
#' @param varsVec character vector of variables that need to be added
#' @param covarsVec character vector of covariates that need to be added
#' @param add boolean indicating if the covariate needs to be added or removed
#' @param indep a boolean indicating if the covariates should be added independently, or
#'  sequentially (append to the previous model). only applicable to adding covariate
#' @return updated ui with added covariates
#' @author  Vishal Sarsani
#' @export
buildupatedUI <- function(ui, varsVec, covarsVec, add = TRUE, indep = FALSE) {
  if (inherits(ui, "nlmixr2FitCore")) {
    ui <- ui$finalUiEnv
  }
  ui <- rxode2::assertRxUi(ui)
  ui <- rxode2::rxUiDecompress(ui)
  # construct covInfo
  covInfo <- buildcovInfo(varsVec, covarsVec)
  # check if the covInfo is a list
  checkmate::assert_list(covInfo)
  covSearchRes <- list()
  covsAdded <- list() # to keep track of covariates added and store in a file
  if (add) {
    # Add covariates one after other
    if (indep) {
      ui_base <- ui  # snapshot before the loop so each candidate starts fresh
      for (i in seq_along(covInfo)) {
        x <- covInfo[[i]]
        covName <- paste0("cov_", x$covariate, "_", x$varName)
        ui_candidate <- tryCatch(
          {
            res <- do.call(.builduiCovariate, c(ui_base, x))
            res # to return 'ui'
          },
          error = function(error_message) {
            message("error  while  ADDING covariate")
            message(error_message)
            message("skipping this covariate")
            res # return NA otherwise (instead of NULL)
          }
        )
        covSearchRes[[i]] <- list(ui_candidate, c(x$covariate, x$varName), covName)[[1]]
      }
      covSearchRes
    } else {
      ## Add all at once
      covsAddedIdx <- 1
      for (x in covInfo) {
        covName <- paste0("cov_", x$covariate, "_", x$varName)
        if (length(covsAdded) == 0) {
          covsAdded[[covsAddedIdx]] <- c(x$covariate, x$varName)
          ui <- tryCatch(
            {
              res <- do.call(.builduiCovariate, c(ui, x))
              res # to return 'ui'
            },
            error = function(error_message) {
              message("error  while SIMULTANEOUSLY ADDING covariates")
              message(error_message)
              message("skipping this covariate")
              res # return NA otherwise (instead of NULL)
            }
          )
          covSearchRes[[covsAddedIdx]] <- list(ui, covsAdded[[covsAddedIdx]], covName)
          covsAddedIdx <- covsAddedIdx + 1
        } else {
          covsAdded[[covsAddedIdx]] <-
            c(covsAdded[[covsAddedIdx - 1]], x$covariate, x$varName)
          ui <- tryCatch(
            {
              res <- do.call(.builduiCovariate, c(ui, x))
              res # to return 'ui'
            },
            error = function(error_message) {
              message("error  while SIMULTANEOUSLY ADDING covariates")
              message(error_message)
              message("skipping this covariate")
              res # return NA otherwise (instead of NULL)
            }
          )
          covSearchRes[[covsAddedIdx]] <- list(ui, covsAdded[[covsAddedIdx]], covName)
          covsAddedIdx <- covsAddedIdx + 1
        }
        covSearchRes
      }
      covSearchRes[length(covSearchRes)][[1]][[1]]
    }
  } else {
    #remove covariate one after other
    for (i in seq_along(covInfo)) {
      x <- covInfo[[i]]
      covName <- paste0("cov_", x$covariate, "_", x$varName)
      ui <- tryCatch(
        {
          res <- do.call(.builduiCovariate, c(ui, x, add = FALSE))
          res # to return 'ui'
        },
        error = function(error_message) {
          message("error  while  ADDING covariate")
          message(error_message)
          message("skipping this covariate")
          res # return NA otherwise (instead of NULL)
        }
      )
      covSearchRes[[i]] <- list(ui, c(x$covariate, x$varName), covName)[[1]]
    }
    covSearchRes
  }
}


#' Make dummy variable cols and updated covarsVec
#'
#' @param data data frame used in the analysis
#' @param covarsVec character vector of covariates that need to be
#'   added
#' @param catcovarsVec character vector of categorical covariates that
#'   need to be added
#' @return return updated Data along with the updated covarsVec
#' @author Vishal Sarsani
#' @export
addCatCovariates <- function(data, covarsVec, catcovarsVec) {
  # check for valid inputs
  checkmate::assert_data_frame(data)
  checkmate::assert_character(covarsVec)
  checkmate::assert_character(catcovarsVec)
  #create new catcovarsvec
  newcatvars <- character(0)
  for (col in catcovarsVec) {
    if (is.factor(data[[col]])) {
      uniqueVals <- levels(data[[col]])
      if (any(is.na(data[[col]]))) {
        uniqueVals <- c(uniqueVals, NA)
      }
      # Else by order values appear.
    } else {
      uniqueVals <- unique(data[[col]])
    }
    uniqueVals <- as.character(uniqueVals)
    # Remove NA values and first dummy
    uniqueVals <- uniqueVals[!is.na(uniqueVals)][-1]
    for (uniqueValue in uniqueVals) {
      colname <- paste0(col, "_", uniqueValue)
      data[, colname] <- as.integer(as.character(data[[col]]) == uniqueValue)
      newcatvars <- c(newcatvars, colname)
    }
  }
  # Remove original categorical variables
  updatedData <- data[, !(names(data) %in% catcovarsVec)]
  #Update entire covarsvec with added categorical variables
  updcovarsVec <- c(covarsVec, newcatvars)
  list(updatedData, updcovarsVec)
}
