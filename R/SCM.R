#' Stepwise Covariate Model-selection (SCM) method
#'
#' @param fit an nlmixr2 'fit' object
#' @param varsVec character vector of PK parameters to which covariates may be
#'   added (e.g. \code{c("ka", "cl", "v")}).  Ignored when \code{pairsVec} is
#'   provided.
#' @param covarsVec character vector of continuous covariate names to test.
#'   Ignored when \code{pairsVec} is provided.
#' @param catvarsVec character vector of categorical covariate names.  Dummy
#'   columns are created via \code{addCatCovariates} and appended to
#'   \code{covarsVec} before building the candidate set.
#' @param pairsVec optional data frame with columns \code{var} and \code{covar}
#'   (or a list of \code{list(var=, covar=)} items) that explicitly enumerates
#'   the parameter–covariate pairs to test.  When supplied, the cartesian
#'   product of \code{varsVec} x \code{covarsVec} is bypassed.  Categorical
#'   dummy columns must already be present in the dataset; use the expanded
#'   covariate name (e.g. \code{covar = "sex_male"}).
#' @param pVal a named list with names \code{fwd} and \code{bck} for the
#'   forward and backward p-value thresholds; default
#'   \code{list(fwd = 0.05, bck = 0.01)}
#' @param inits named list mapping shape names to initial theta estimates for
#'   the covariate parameter, e.g.
#'   \code{list(power = 0.1, lin = 0.01, cat = 0.01)}.  Applied globally to
#'   all pairs; a per-pair scalar \code{init} element in a \code{pairsVec} list
#'   item takes precedence.  Shapes not listed default to \code{0}.  Suggested
#'   values: \code{0.1} for \code{"power"} (allometric exponent near null),
#'   \code{0.01} for other shapes.
#' @param missingToken sentinel value used to identify missing covariate
#'   observations in addition to \code{NA}.  Default \code{NA} (only true
#'   \code{NA} values are treated as missing).  Set to a numeric value (e.g.
#'   \code{-99}) or a character string (e.g. \code{"."}) to also flag that
#'   sentinel as missing.  When missing observations are detected for a
#'   covariate, the generated model expression is wrapped in an
#'   \code{ifelse()} guard so that the covariate contribution is set to the
#'   population-typical value: \code{0} for continuous covariates (equivalent
#'   to imputing the median, since \code{log(median/median) = 0}) or the mode
#'   of the observed levels for categorical covariates.
#' @param data data frame containing all subjects and columns required by the
#'   search (covariates, dose/observation records, etc.).  When \code{NULL}
#'   (default), \code{nlme::getData(fit)} is used, which may omit covariate
#'   columns that are not referenced by the base model.  Pass the data set
#'   explicitly whenever the base model does not reference a covariate column
#'   that will be used in the search (e.g. a categorical covariate such as
#'   \code{sex}).  No column transformations should be applied before passing
#'   the data: centering and other shape transformations are generated inside
#'   the model body by the SCM machinery.
#' @param outputDir character; path of the subdirectory used to store all
#'   fitted model \code{.rds} files.  When \code{NULL} (default), the name is
#'   derived automatically as \code{<fit_name>_scm_<N>} where \code{N} is one
#'   greater than the number of existing \code{<fit_name>_scm_*} directories in
#'   the current working directory (so repeated runs produce
#'   \code{fit_scm_1}, \code{fit_scm_2}, …).  If the resolved directory
#'   already exists, the search stops with an error; set
#'   \code{restart = TRUE} to back it up and start fresh instead.  Ignored
#'   when \code{saveModels = FALSE}.
#' @param saveModels logical; if \code{TRUE} (default) every fitted candidate
#'   model accepted during forward or backward search is written to
#'   \code{outputDir} as an \code{.rds} file.  Set to \code{FALSE} to run the
#'   search without writing any files.
#' @param verbose logical; if \code{TRUE}, print additional detail at each
#'   search step: the full candidate table (with covariate expression and
#'   starting estimate) before fitting, the complete results table for every
#'   candidate after fitting (sorted by p-value), and—when a covariate is
#'   accepted—the accepted model's fixed-effect estimates and diagonal omega
#'   variances.  Default \code{FALSE} (only the summary line and accepted-
#'   model best row are printed, as at present).
#' @param includedRelations optional specification of covariate–parameter
#'   relationships that must be present in the model at the start of backward
#'   elimination (analogous to PsN's \code{[included_relations]}).  Can be a
#'   data frame with columns \code{var} and \code{covar} (and optionally
#'   \code{shape} or \code{shapes}), or a list of
#'   \code{list(var=, covar=, shapes=)} items in the same format as
#'   \code{pairsVec}.  Relations already in the forward-search result are left
#'   untouched; missing ones are added to the model and a single re-fit is
#'   performed before backward elimination begins.  Included relations are
#'   eligible for removal during backward elimination like any other covariate.
#'   Ignored when \code{searchType = "forward"}.
#' @param shapes character vector of shape names to test for each continuous
#'   covariate–parameter pair.  Built-in shapes are \code{"power"} (default;
#'   \code{log(cov/median)}), \code{"lin"} (\code{cov - median}), \code{"log"}
#'   (\code{log(cov)}), and \code{"identity"} (raw covariate).  Per-pair
#'   overrides can be specified via the \code{shapes} element of individual
#'   \code{pairsVec} list items.  Categorical covariates always use the
#'   \code{"cat"} shape regardless of this setting.
#' @param customShapes named list of additional shape builder functions of the
#'   form \code{function(col, center, level) -> character}.  These are merged
#'   with the built-in shapes and take precedence over them.
#' @param searchType one of \code{'scm'}, \code{'forward'}, or
#'   \code{'backward'}; default \code{'scm'}
#' @param restart logical; if \code{TRUE} the existing cache is backed up and
#'   the search restarts from scratch; default \code{FALSE}
#'
#' @return A list with elements \code{summaryTable} (combined forward and
#'   backward results), \code{resFwd} (list of final fit and step table from
#'   forward search), and \code{resBck} (same for backward search).
#'
#' @export
#' @author Vipul Mann, Matthew Fidler, Vishal Sarsani
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
#' fit <- nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "saem", control = list(print = 0))
#' rxode2::.rxWithWd(tempdir(), {
#'
#'   auto1 <- covarSearchAuto(fit, varsVec = c("ka", "cl"), covarsVec = c("WT"))
#'
#'   # Specify exact pairs instead of all combinations:
#'   auto2 <- covarSearchAuto(fit,
#'     pairsVec = list(list(var = "cl", covar = "WT"), list(var = "ka", covar = "WT")))
#' })
#' }
covarSearchAuto <- function(fit,
                            varsVec           = NULL,
                            covarsVec         = NULL,
                            pVal              = list(fwd = 0.05, bck = 0.01),
                            catvarsVec        = NULL,
                            pairsVec          = NULL,
                            shapes            = "power",
                            customShapes      = NULL,
                            inits             = list(),
                            missingToken      = NA,
                            includedRelations = NULL,
                            data              = NULL,
                            outputDir         = NULL,
                            saveModels        = TRUE,
                            verbose           = FALSE,
                            searchType        = c("scm", "forward", "backward"),
                            restart           = FALSE) {
  if (!is.numeric(AIC(fit))) {
    cli::cli_alert_danger("the 'fit' object needs to have an objective function value associated with it")
    cli::cli_alert_info("try computing 'AIC(fitobject)' in console to compute and store the corresponding OBJF value")
    stop("aborting...objf value not associated with the current 'fit' object", call. = FALSE)
  }

  if (is.null(data)) data <- nlme::getData(fit)
  # When using varsVec/covarsVec, merge catvarsVec into covarsVec so that
  # buildPairs creates (var, catcovar) base rows that .enrichPairs then
  # expands into per-level rows.  pairsVec users must include cat covariates
  # explicitly (either pre-expanded with level info, or as raw column names
  # paired with catvarsVec for auto-expansion).
  if (!is.null(varsVec) && is.null(pairsVec)) {
    covarsVec <- c(covarsVec, catvarsVec)
  }
  pairs <- buildPairs(varsVec = varsVec, covarsVec = covarsVec, pairsVec = pairsVec)
  pairs <- .enrichPairs(pairs, data, catvarsVec = catvarsVec, missingToken = missingToken)
  pairs <- .expandShapes(pairs, shapes = shapes, customShapes = customShapes,
                         inits = inits)

  # Process includedRelations through the same pipeline so it carries covExpr
  included_pairs <- NULL
  if (!is.null(includedRelations)) {
    included_pairs <- buildPairs(pairsVec = includedRelations)
    included_pairs <- .enrichPairs(included_pairs, data, catvarsVec = catvarsVec,
                                   missingToken = missingToken)
    included_pairs <- .expandShapes(included_pairs,
                                    shapes = shapes, customShapes = customShapes,
                                    inits = inits)
  }

  searchType <- match.arg(searchType)

  if (!all(names(pVal) %in% c("fwd", "bck"))) {
    stop("pVal must be a list with names 'fwd' and 'bck'")
  }

  # ── Resolve and prepare output directory ──────────────────────────────────
  if (saveModels) {
    if (is.null(outputDir)) {
      # Auto-name: <fit_name>_scm_<N> where N = existing count + 1
      model_name <- gsub("[^A-Za-z0-9_.]", "_", as.character(substitute(fit)))
      n_existing <- sum(grepl(
        paste0("^", model_name, "_scm_[0-9]+$"),
        list.dirs(".", full.names = FALSE, recursive = FALSE)
      ))
      outputDir <- paste0(model_name, "_scm_", n_existing + 1L)
    }

    if (!dir.exists(outputDir)) {
      dir.create(outputDir)
    } else if (isTRUE(restart)) {
      backupDir <- paste0(outputDir, "_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
      file.rename(outputDir, backupDir)
      cli::cli_alert_info("Prior cache moved to {backupDir}")
      dir.create(outputDir)
    } else {
      stop(
        "Output directory '", outputDir, "' already exists. ",
        "Use restart = TRUE to back it up and start fresh, or set ",
        "saveModels = FALSE to run without saving.",
        call. = FALSE
      )
    }
    cli::cli_alert_info("Saving models to: {outputDir}")
  } else {
    # saveModels = FALSE: no directory needed; use a placeholder for forwardSearch/backwardSearch
    if (is.null(outputDir)) outputDir <- tempdir()
  }

  if (searchType == "scm") {
    resFwd <- forwardSearch(pairs, fit, data, pVal = pVal$fwd,
                            outputDir = outputDir, restart = restart,
                            saveModels = saveModels, verbose = verbose)
    resBck <- backwardSearch(pairs, fitorig = fit, fitupdated = resFwd[[1]],
                             data = data, pVal = pVal$bck, reFitCovars = FALSE,
                             context_pairs = resFwd[[3]],
                             includedRelations = included_pairs,
                             outputDir = outputDir, restart = restart,
                             saveModels = saveModels, verbose = verbose)
    summaryTable <- Reduce(rbind, list(resFwd[[2]], resBck[[2]]))
    return(list(summaryTable = summaryTable, resFwd = resFwd, resBck = resBck))

  } else if (searchType == "forward") {
    resFwd <- forwardSearch(pairs, fit, data, pVal = pVal$fwd,
                            outputDir = outputDir, restart = restart,
                            saveModels = saveModels, verbose = verbose)
    return(list(summaryTable = resFwd[[2]], resFwd = resFwd, resBck = NULL))

  } else {
    resBck <- backwardSearch(pairs, fitorig = fit, data = data,
                             pVal = pVal$bck, reFitCovars = TRUE,
                             includedRelations = included_pairs,
                             outputDir = outputDir, restart = restart,
                             saveModels = saveModels, verbose = verbose)
    return(list(summaryTable = resBck[[2]], resFwd = NULL, resBck = resBck))
  }
}

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Build a pairs data frame from varsVec/covarsVec or an explicit pairsVec
#'
#' @param varsVec character vector of PK parameter names
#' @param covarsVec character vector of covariate names
#' @param pairsVec optional data frame with columns \code{var} and \code{covar},
#'   or a list of \code{list(var=, covar=)} items
#' @return data frame with columns \code{var} and \code{covar}, one row per pair
#' @noRd
buildPairs <- function(varsVec = NULL, covarsVec = NULL, pairsVec = NULL) {
  if (!is.null(pairsVec)) {
    if (is.data.frame(pairsVec)) {
      stopifnot(all(c("var", "covar") %in% names(pairsVec)))
      keep <- intersect(c("var", "covar", "shapes", "inits"), names(pairsVec))
      df   <- pairsVec[, keep, drop = FALSE]
      # If caller used a 'shape' character column instead of a 'shapes' list-col,
      # convert it so .expandShapes() can dispatch per-row shapes correctly.
      if ("shape" %in% names(pairsVec) && !"shapes" %in% names(df)) {
        df$shapes <- as.list(pairsVec$shape)
      }
      return(df)
    }
    # Build column vectors manually to support optional list-columns
    n           <- length(pairsVec)
    var_vec     <- character(n)
    covar_vec   <- character(n)
    shapes_list <- vector("list", n)
    inits_list  <- vector("list", n)
    has_shapes  <- FALSE
    has_inits   <- FALSE
    for (k in seq_len(n)) {
      p            <- pairsVec[[k]]
      var_vec[k]   <- if (is.list(p)) p$var   else p[[1]]
      covar_vec[k] <- if (is.list(p)) p$covar else p[[2]]
      if (is.list(p) && !is.null(p$shapes)) {
        shapes_list[[k]] <- p$shapes
        has_shapes        <- TRUE
      }
      if (is.list(p) && !is.null(p$inits)) {
        inits_list[[k]] <- p$inits   # named list, e.g. list(power=0.75, lin=0.01)
        has_inits       <- TRUE
      }
    }
    df <- data.frame(var = var_vec, covar = covar_vec, stringsAsFactors = FALSE)
    if (has_shapes) df$shapes <- shapes_list
    if (has_inits)  df$inits  <- inits_list
    return(df)
  }
  if (is.null(varsVec) || is.null(covarsVec)) {
    stop("provide either pairsVec or both varsVec and covarsVec")
  }
  expand.grid(var = varsVec, covar = covarsVec, stringsAsFactors = FALSE)
}

#' Enrich a pairs data frame with covariate type, centering values, and levels
#'
#' For continuous covariates the median of the covariate column is computed from
#' \code{data} and stored in \code{center}; the model expression will be
#' \code{cov_theta * log(covariate / center)}.  For categorical covariates each
#' row is expanded into one row per non-reference level and the model expression
#' will be \code{cov_theta * ifelse(raw_col == "level", 1, 0)}.
#'
#' When \code{missingToken} is supplied (or the covariate column contains
#' \code{NA}s), three extra columns are added to support wrapping the
#' covariate expression in an \code{ifelse()} guard inside
#' \code{.expandShapes()}:
#' \describe{
#'   \item{has_missing}{logical — whether any observations are missing}
#'   \item{missing_check}{character — the \code{ifelse} condition string,
#'     e.g. \code{"is.na(wt)"} or \code{"is.na(wt) | wt == -99"}}
#'   \item{missing_fill}{integer — replacement value (0 for continuous;
#'     0 or 1 for categorical based on whether the observed mode equals
#'     the indicator level)}
#' }
#'
#' @param pairs data frame with at least \code{var} and \code{covar} columns
#' @param data  data frame used for fitting (source of medians and factor levels)
#' @param catvarsVec character vector of categorical covariate column names;
#'   rows whose \code{covar} matches an entry are expanded per level
#' @param missingToken sentinel value treated as missing in addition to
#'   \code{NA}; default \code{NA} (only true \code{NA}s detected)
#' @return enriched data frame with additional columns \code{type},
#'   \code{center}, \code{raw_col}, \code{level}, and (when missing values
#'   are present) \code{has_missing}, \code{missing_check}, \code{missing_fill}
#' @noRd
.enrichPairs <- function(pairs, data, catvarsVec = NULL, missingToken = NA) {
  if (!"type"    %in% names(pairs)) pairs$type    <- "continuous"
  if (!"center"  %in% names(pairs)) pairs$center  <- NA_real_
  if (!"raw_col" %in% names(pairs)) pairs$raw_col <- pairs$covar
  if (!"level"   %in% names(pairs)) pairs$level   <- NA_character_

  # Expand categorical covariate rows into one row per non-reference level
  if (!is.null(catvarsVec) && length(catvarsVec) > 0) {
    keep  <- list()
    extra <- list()
    for (i in seq_len(nrow(pairs))) {
      rc <- pairs$covar[i]
      if (rc %in% catvarsVec) {
        col_data <- data[[rc]]
        lvls <- if (is.factor(col_data)) {
          levels(col_data)
        } else {
          sort(unique(as.character(col_data[!is.na(col_data)])))
        }
        lvls <- lvls[-1]  # drop first (reference) level
        for (lev in lvls) {
          row          <- pairs[i, , drop = FALSE]
          row$covar    <- paste0(rc, "_", lev)
          row$type     <- "categorical"
          row$center   <- NA_real_
          row$raw_col  <- rc
          row$level    <- lev
          extra[[length(extra) + 1]] <- row
        }
      } else {
        keep[[length(keep) + 1]] <- pairs[i, , drop = FALSE]
      }
    }
    pairs <- do.call(rbind, c(keep, extra))
    rownames(pairs) <- NULL
  }

  # Compute medians for continuous covariates that don't already have a center
  for (i in seq_len(nrow(pairs))) {
    if (pairs$type[i] == "continuous" && is.na(pairs$center[i])) {
      col <- pairs$raw_col[i]
      if (col %in% names(data)) {
        pairs$center[i] <- median(data[[col]], na.rm = TRUE)
      }
    }
  }

  # ── Missing value detection ────────────────────────────────────────────────
  # For each row, check whether the covariate has NA or missingToken values.
  # When found, record three columns consumed by .expandShapes() to wrap the
  # covariate expression in an ifelse() guard:
  #
  #   has_missing   – logical flag
  #   missing_check – condition string for the ifelse() (e.g. "is.na(wt)")
  #   missing_fill  – replacement value when condition is TRUE:
  #                     continuous  → 0 (imputing to median → log(m/m) = 0)
  #                     categorical → 1 if mode == level, else 0
  pairs$has_missing   <- FALSE
  pairs$missing_fill  <- 0L
  pairs$missing_check <- NA_character_

  for (i in seq_len(nrow(pairs))) {
    rc <- pairs$raw_col[i]
    if (!rc %in% names(data)) next
    col_data <- data[[rc]]

    # Identify missing observations (NA and/or sentinel token)
    is_missing <- is.na(col_data)
    if (!is.na(missingToken)) {
      is_missing <- is_missing | (!is.na(col_data) & col_data == missingToken)
    }
    if (!any(is_missing)) next

    pairs$has_missing[i] <- TRUE

    # Build the ifelse() condition string
    check <- paste0("is.na(", rc, ")")
    if (!is.na(missingToken)) {
      # Quote the token if the data column is character/factor or the token
      # itself is a string; leave unquoted for numeric tokens on numeric columns
      quote_tok <- is.character(missingToken) ||
                   inherits(col_data, c("character", "factor"))
      if (quote_tok) {
        check <- paste0(check, ' | ', rc, ' == "', missingToken, '"')
      } else {
        check <- paste0(check, " | ", rc, " == ", missingToken)
      }
    }
    pairs$missing_check[i] <- check

    # Fill value: continuous → 0 (median imputation on log scale)
    #             categorical → 1 if mode of valid values equals this level
    if (pairs$type[i] == "continuous") {
      pairs$missing_fill[i] <- 0L
    } else {
      valid_vals <- as.character(col_data[!is_missing])
      if (length(valid_vals) == 0L) {
        pairs$missing_fill[i] <- 0L
      } else {
        mode_val <- names(sort(table(valid_vals), decreasing = TRUE))[[1L]]
        pairs$missing_fill[i] <- as.integer(mode_val == pairs$level[i])
      }
    }
  }

  pairs
}

#' Expand a pairs data frame into one row per (pair, shape)
#'
#' For continuous covariates each shape in \code{shapes} (or the per-row
#' \code{shapes} list-column when present) produces one output row.  The
#' \code{covar} column is suffixed with the shape name so that theta identifiers
#' remain unique (e.g. \code{"wt"} with shapes \code{c("power","lin")} yields
#' rows with \code{covar = "wt_power"} and \code{covar = "wt_lin"}).  Categorical
#' rows are not expanded; they always use the \code{"cat"} shape and their
#' \code{covar} name is unchanged.
#'
#' The \code{covExpr} column is pre-computed by \code{.scmCovExpr()} so that
#' downstream code can pass it directly to \code{.builduiCovariate()} without
#' repeating shape logic.
#'
#' @param pairs   enriched pairs data frame (output of \code{.enrichPairs()})
#' @param shapes  character vector of shape names applied globally to continuous
#'   pairs with no per-row override; default \code{"power"}
#' @param customShapes named list of additional shape builder functions; merged
#'   with \code{.SCM_SHAPES} (custom entries take precedence)
#' @param inits   named list mapping shape names to initial theta estimates,
#'   e.g. \code{list(power = 0.1, lin = 0.01, cat = 0.01)}.  Used as the
#'   global default; a per-pair \code{inits} list-column (set in
#'   \code{buildPairs()} from \code{pairsVec} items) takes precedence on a
#'   shape-by-shape basis, allowing different starting values per shape even
#'   within the same relationship.  Shapes not covered by either source
#'   default to \code{0}.
#' @return expanded pairs data frame with additional columns \code{shape},
#'   \code{covExpr}, and \code{init}; the \code{shapes} and \code{inits}
#'   list-columns (if present) are removed
#' @noRd
.expandShapes <- function(pairs, shapes = "power", customShapes = NULL,
                          inits = list()) {
  has_per_row   <- "shapes"      %in% names(pairs)
  has_per_inits <- "inits"       %in% names(pairs)
  has_miss_col  <- "has_missing" %in% names(pairs)
  result_rows <- vector("list", nrow(pairs) * length(shapes))
  out_idx     <- 0L

  for (i in seq_len(nrow(pairs))) {
    row      <- pairs[i, , drop = FALSE]
    cov_type <- row$type
    raw_col  <- row$raw_col
    center   <- if (is.na(row$center)) NULL else row$center
    level    <- if (is.na(row$level))  NULL else row$level
    base_covar <- row$covar

    # Per-pair inits: named list keyed by shape (e.g. list(power=0.75, lin=0.01))
    pair_inits <- if (has_per_inits) row$inits[[1]] else list()

    # Missing-value guard (populated by .enrichPairs())
    has_miss   <- has_miss_col && isTRUE(row$has_missing)
    miss_check <- if (has_miss) row$missing_check else NA_character_
    miss_fill  <- if (has_miss) row$missing_fill  else 0L

    if (cov_type == "categorical") {
      expr <- .scmCovExpr("cat", raw_col, level = level,
                          customShapes = customShapes)
      if (has_miss && !is.na(miss_check)) {
        expr <- paste0("ifelse(", miss_check, ", ", miss_fill, ", ", expr, ")")
      }
      row$shape   <- "cat"
      row$covExpr <- expr
      row$init    <- if (!is.null(pair_inits[["cat"]])) pair_inits[["cat"]] else
                       if (!is.null(inits[["cat"]])) inits[["cat"]] else 0
      if (has_per_row)   row$shapes <- NULL
      if (has_per_inits) row$inits  <- NULL
      out_idx <- out_idx + 1L
      result_rows[[out_idx]] <- row
    } else {
      # Determine which shapes to test for this continuous row
      row_shapes <- if (has_per_row && length(row$shapes[[1]]) > 0) {
        row$shapes[[1]]
      } else {
        shapes
      }
      for (sh in row_shapes) {
        r    <- row
        expr <- .scmCovExpr(sh, raw_col, center = center,
                            customShapes = customShapes)
        if (has_miss && !is.na(miss_check)) {
          expr <- paste0("ifelse(", miss_check, ", ", miss_fill, ", ", expr, ")")
        }
        r$shape   <- sh
        r$covExpr <- expr
        r$covar   <- paste0(base_covar, "_", sh)
        r$init    <- if (!is.null(pair_inits[[sh]])) pair_inits[[sh]] else
                       if (!is.null(inits[[sh]])) inits[[sh]] else 0
        if (has_per_row)   r$shapes <- NULL
        if (has_per_inits) r$inits  <- NULL
        out_idx <- out_idx + 1L
        result_rows[[out_idx]] <- r
      }
    }
  }

  result <- do.call(rbind, result_rows[seq_len(out_idx)])
  rownames(result) <- NULL
  result
}

#' Filter a pairs data frame to only those pairs present in a fitted model
#'
#' @param pairs data frame with columns \code{var} and \code{covar}
#' @param fit nlmixr2 fit object
#' @return subset of \code{pairs} whose \code{cov_<covar>_<var>} theta exists in \code{fit}
#' @noRd
.pairsInModel <- function(pairs, fit) {
  mrd <- tryCatch(fit$finalUiEnv$muRefCovariateDataFrame, error = function(e) NULL)
  if (!is.null(mrd) && nrow(mrd) > 0 && "covariateParameter" %in% names(mrd)) {
    present <- mrd$covariateParameter
  } else {
    ini_names <- tryCatch(fit$finalUiEnv$iniDf$name, error = function(e) character(0))
    present   <- ini_names[grepl("^cov_", ini_names)]
  }
  expected <- paste0("cov_", pairs$covar, "_", pairs$var)
  pairs[expected %in% present, , drop = FALSE]
}

#' Compute % reduction in omega variance for the eta associated with a variable
#'
#' @param fit_without fit without the covariate (higher omega expected)
#' @param fit_with    fit with the covariate
#' @param nam_var     variable name (e.g. \code{"cl"})
#' @return numeric % reduction, or \code{NA_real_} if the eta cannot be resolved
#' @noRd
.bsvReduction <- function(fit_without, fit_with, nam_var) {
  tryCatch({
    theta_name <- .getThetaName(fit_without$finalUiEnv, nam_var)
    mrd  <- fit_without$muRefDataFrame
    eta  <- mrd$eta[mrd$theta == theta_name]
    if (length(eta) == 0) return(NA_real_)
    eta  <- eta[1]
    om_without <- diag(fit_without$omega)[eta]
    om_with    <- diag(fit_with$omega)[eta]
    if (is.na(om_without) || om_without == 0) return(NA_real_)
    (om_without - om_with) / om_without * 100
  }, error = function(e) NA_real_)
}

#' Fit all candidate models for one SCM step
#'
#' Builds a candidate UI for every row of \code{pairs} (adding or removing the
#' covariate–parameter relationship), fits each one, and assembles a per-row
#' result list.  Rows whose candidate cannot be built or fitted are silently
#' dropped.
#'
#' @param pairs    data frame with columns \code{var} and \code{covar}
#' @param ui       current model UI (rxode2 environment)
#' @param fit      current reference fit (nlmixr2FitCore)
#' @param data     data frame passed to \code{nlmixr2}
#' @param pVal     p-value threshold (used only to compute the \code{qchisqr} label)
#' @param stepIdx  integer step counter stored in the stats table
#' @param add      logical; \code{TRUE} = forward (add covariate),
#'                 \code{FALSE} = backward (remove covariate)
#' @return list of \code{list(var, covar, fit, stats)}, one per successful candidate
#' @noRd
.fitCandidatePairs <- function(pairs, base_ui, context_pairs = NULL, fit, data, pVal, stepIdx, add) {
  Filter(Negate(is.null), lapply(seq_len(nrow(pairs)), function(i) {
    nam_var   <- pairs$var[i]
    nam_covar <- pairs$covar[i]
    covNames  <- paste0("cov_", nam_covar, "_", nam_var)

    # Extract pre-computed covExpr, shape, and init (set by .expandShapes())
    cov_expr  <- if ("covExpr" %in% names(pairs)) pairs$covExpr[i] else NULL
    cov_shape <- if ("shape"   %in% names(pairs)) pairs$shape[i]   else NA_character_
    cov_init  <- if ("init"    %in% names(pairs) && !is.na(pairs$init[i]))
                   pairs$init[i] else 0

    ui_cand <- tryCatch({
      if (add) {
        # Forward: rebuild context model from base, then add candidate
        ctx_ui <- .rebuildUiFromPairs(base_ui, context_pairs)
        .builduiCovariate(ctx_ui, nam_var, nam_covar, add = TRUE,
                          covExpr = cov_expr, init = cov_init)
      } else {
        # Backward: rebuild from base with all context pairs except this candidate
        cp_minus <- if (!is.null(context_pairs) && nrow(context_pairs) > 0) {
          context_pairs[
            !(context_pairs$var == nam_var & context_pairs$covar == nam_covar), ,
            drop = FALSE
          ]
        } else {
          NULL
        }
        .rebuildUiFromPairs(base_ui, cp_minus)
      }
    }, error = function(e) {
      message("error building candidate for ",
              nam_covar, "~", nam_var, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(ui_cand)) return(NULL)

    x <- tryCatch(
      suppressWarnings(nlmixr2(ui_cand, data, fit$est)),
      error = function(e) {
        message("error fitting candidate for ", nam_covar, "~", nam_var, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(x)) return(NULL)

    # dObjf > 0 means the model change improved the fit:
    #   forward:  candidate has lower OFV than reference  → fit$objf - x$objf
    #   backward: removal raises OFV relative to reference → x$objf - fit$objf
    dObjf <- if (add) fit$objf - x$objf else x$objf - fit$objf
    dof   <- if (add) {
      length(x$finalUiEnv$ini$est) - length(fit$finalUiEnv$ini$est)
    } else {
      length(fit$finalUiEnv$ini$est) - length(x$finalUiEnv$ini$est)
    }
    pchisqr <- if (dObjf > 0) 1 - pchisq(dObjf, df = dof) else 1

    # covarEffect: from the candidate fit when adding, from the full model when removing
    covarEffect  <- if (add) x$parFixedDf[covNames, "Estimate"] else fit$parFixedDf[covNames, "Estimate"]
    bsvReduction <- if (add) {
      .bsvReduction(fit_without = fit, fit_with = x,   nam_var = nam_var)
    } else {
      .bsvReduction(fit_without = x,   fit_with = fit, nam_var = nam_var)
    }

    stats <- data.frame(
      step         = stepIdx,
      covar        = nam_covar,
      var          = nam_var,
      shape        = cov_shape,
      objf         = x$objf,
      deltObjf     = dObjf,
      AIC          = x$AIC,
      BIC          = x$BIC,
      numParams    = length(x$finalUiEnv$ini$est),
      qchisqr      = qchisq(1 - pVal, dof),
      pchisqr      = pchisqr,
      included     = if (add) "no" else "",
      searchType   = if (add) "forward" else "backward",
      covNames     = covNames,
      covarEffect  = covarEffect,
      bsvReduction = bsvReduction,
      stringsAsFactors = FALSE
    )

    list(var = nam_var, covar = nam_covar, fit = x, stats = stats)
  }))
}

# ── Search functions ───────────────────────────────────────────────────────────

#' Forward covariate search
#'
#' @param pairs     data frame with columns \code{var} and \code{covar}
#' @param fit       base nlmixr2 fit object
#' @param data      data frame for fitting
#' @param pVal      p-value threshold for inclusion
#' @param outputDir cache directory path
#' @param restart   logical; \code{FALSE} = resume from cache (only when
#'   \code{saveModels = TRUE})
#' @param saveModels logical; if \code{FALSE} all \code{saveRDS}/\code{readRDS}
#'   operations are skipped
#' @param verbose logical; if \code{TRUE} print candidate table before each
#'   step, full results table after fitting, and parameter/omega detail when a
#'   covariate is accepted
#' @return \code{list(final_fit, resTableComplete)}
#' @noRd
#' @author Vipul Mann, Matthew Fidler, Vishal Sarsani
forwardSearch <- function(pairs, fit, data, pVal = 0.05, outputDir, restart = FALSE,
                          saveModels = TRUE, verbose = FALSE) {
  if (!inherits(fit, "nlmixr2FitCore")) stop("'fit' needs to be a nlmixr2 fit")
  if (saveModels && missing(outputDir)) stop("please specify outputDir for forward search")

  base_ui        <- fit$finalUiEnv  # fixed — never updated
  orig_pairs     <- pairs           # preserved for resume reconstruction
  accepted_pairs <- NULL
  resTableComplete <- NULL
  stepIdx <- 1

  # When not saving there is no cache to read; treat as a clean start
  if (!saveModels) restart <- TRUE

  if (!restart) {
    fileExistsTab  <- list.files(paste0("./", outputDir),
      pattern = "forward_step_[0-9]+_table_[a-zA-Z0-9_]+\\.rds")
    fileExistsFit  <- list.files(paste0("./", outputDir),
      pattern = "forward_step_[0-9]+_fit_[a-zA-Z0-9_]+\\.rds")
    fileExistsCT   <- list.files(paste0("./", outputDir),
      pattern = "forward_step_[0-9]+_completetable_[a-zA-Z0-9_]+\\.rds")

    if (length(fileExistsTab) == 0) {
      restart <- TRUE
    } else {
      resumeRows  <- lapply(fileExistsTab, function(f) readRDS(paste0(outputDir, "/", f)))
      resumeTable <- data.table::rbindlist(resumeRows)
      fit            <- readRDS(paste0(outputDir, "/", fileExistsFit[[length(fileExistsFit)]]))
      accepted_pairs <- .pairsInModel(orig_pairs, fit)
      resTableComplete <- readRDS(paste0(outputDir, "/", fileExistsCT[[length(fileExistsCT)]]))
      stepIdx <- unlist(resumeTable[nrow(resumeTable), ]$step) + 1
      for (i in seq_len(nrow(resumeTable))) {
        pairs <- pairs[!(pairs$var == resumeTable$var[i] & pairs$covar == resumeTable$covar[i]), ]
      }
      cli::cli_alert_success("loaded forward search data from disk, resuming ...")
    }
  }

  cli::cli_h1("starting forward search...")

  while (nrow(pairs) > 0) {
    # ── verbose: current model code + candidates before fitting ───────────
    if (verbose) {
      cli::cli_h2("Forward step {stepIdx}: current model")
      tryCatch(
        print(fit$finalUiEnv),
        error = function(e) cli::cli_alert_warning(
          "Model code unavailable: {conditionMessage(e)}")
      )
      cli::cli_h2("Forward step {stepIdx}: {nrow(pairs)} candidate(s) to test")
      v_cols <- intersect(
        c("var", "covar", "shape", "init", "covExpr"),
        names(pairs)
      )
      print(pairs[, v_cols, drop = FALSE])
    }

    results <- .fitCandidatePairs(pairs, base_ui, accepted_pairs, fit, data, pVal, stepIdx, add = TRUE)
    if (length(results) == 0) break

    resTable <- do.call(rbind, lapply(results, `[[`, "stats"))
    bestIdx  <- which.min(resTable$pchisqr)
    bestRow  <- resTable[bestIdx, , drop = FALSE]

    # ── verbose: show all candidate results ───────────────────────────────
    if (verbose) {
      cli::cli_h2("Forward step {stepIdx}: all candidate results")
      v_cols <- intersect(
        c("covar", "var", "shape", "objf", "deltObjf",
          "pchisqr", "AIC", "BIC", "bsvReduction", "covarEffect"),
        names(resTable)
      )
      print(resTable[order(resTable$pchisqr), v_cols, drop = FALSE])
    }

    if (bestRow$pchisqr <= pVal) {
      resTable[bestIdx, "included"] <- "yes"
      bestRow[, "included"] <- "yes"

      cli::cli_h1("best model at step {stepIdx}:")
      print(bestRow)

      fit <- results[[bestIdx]]$fit

      # ── verbose: updated model code + parameter estimates + omega ────────
      if (verbose) {
        cli::cli_h2("Accepted model — updated model code:")
        tryCatch(
          print(fit$finalUiEnv),
          error = function(e) cli::cli_alert_warning(
            "Model code unavailable: {conditionMessage(e)}")
        )
        cli::cli_h2("Accepted model — fixed effects:")
        print(fit$parFixedDf)
        cli::cli_h2("Accepted model — omega diagonal (BSV variances):")
        print(round(diag(fit$omega), 4))
      }

      acc_var   <- results[[bestIdx]]$var
      acc_covar <- results[[bestIdx]]$covar

      # Track accepted pair for rebuild-from-base at next forward step
      acc_row <- pairs[pairs$var == acc_var & pairs$covar == acc_covar, , drop = FALSE]
      if (nrow(acc_row) > 0L) {
        accepted_pairs <- rbind(accepted_pairs, acc_row[1L, , drop = FALSE])
      }

      # Remove the accepted pair AND all other shapes of the same raw covariate
      # on the same parameter — only one parameterisation of a given
      # covariate–parameter relationship may enter the model.
      if ("raw_col" %in% names(pairs)) {
        acc_raw <- pairs$raw_col[pairs$var == acc_var & pairs$covar == acc_covar]
        acc_raw <- if (length(acc_raw) > 0L) acc_raw[1L] else NA_character_
        if (!is.na(acc_raw)) {
          dropped <- pairs$covar[pairs$var == acc_var & pairs$raw_col == acc_raw &
                                   pairs$covar != acc_covar]
          if (length(dropped) > 0L)
            cli::cli_alert_info(
              "Dropping alternative shape(s) for {acc_raw}~{acc_var}: {paste(dropped, collapse=', ')}")
          pairs <- pairs[!(pairs$var == acc_var & pairs$raw_col == acc_raw), ]
        } else {
          pairs <- pairs[!(pairs$var == acc_var & pairs$covar == acc_covar), ]
        }
      } else {
        pairs <- pairs[!(pairs$var == acc_var & pairs$covar == acc_covar), ]
      }

      resTableComplete <- rbind(resTableComplete, resTable)
      if (saveModels) {
        key <- paste0(acc_covar, "_", acc_var)
        saveRDS(fit,              paste0(outputDir, "/forward_step_", stepIdx, "_fit_",           key, ".rds"))
        saveRDS(bestRow,          paste0(outputDir, "/forward_step_", stepIdx, "_table_",         key, ".rds"))
        saveRDS(resTableComplete, paste0(outputDir, "/forward_step_", stepIdx, "_completetable_", key, ".rds"))
      }

      cli::cli_h2("accepted {acc_covar}~{acc_var}")
      stepIdx <- stepIdx + 1
    } else {
      cli::cli_h1("OFV did not improve, exiting forward search ...")
      resTableComplete <- rbind(resTableComplete, resTable)
      break
    }
  }

  cli::cli_h2(cli::col_red("forward search complete"))
  list(fit, resTableComplete, accepted_pairs)
}

#' Backward covariate search
#'
#' @param pairs       data frame with columns \code{var} and \code{covar}
#' @param fitorig     the original base fit object
#' @param fitupdated  the fit after forward search (if any)
#' @param data        data frame for fitting
#' @param pVal        p-value threshold for retaining a covariate
#' @param reFitCovars logical; if \code{TRUE} all pairs are added simultaneously
#'   before backward elimination (used for standalone backward search)
#' @param includedRelations enriched pairs data frame (output of the pipeline
#'   in \code{covarSearchAuto}) of relationships that must be present in the
#'   model before backward elimination begins.  Any that are missing are added
#'   and the model is re-fit once.  \code{NULL} = no forced inclusions.
#' @param outputDir   cache directory path
#' @param restart     logical; \code{FALSE} = resume from cache (only when
#'   \code{saveModels = TRUE})
#' @param saveModels  logical; if \code{FALSE} all \code{saveRDS}/\code{readRDS}
#'   operations are skipped
#' @param verbose     logical; if \code{TRUE} print candidate table before each
#'   step, full results table after fitting, and parameter/omega detail when a
#'   covariate is retained
#' @return \code{list(final_fit, resTableComplete)}
#' @noRd
#' @author Vipul Mann, Matthew Fidler, Vishal Sarsani
backwardSearch <- function(pairs, fitorig, fitupdated, data,
                           pVal = 0.01, reFitCovars = FALSE,
                           context_pairs = NULL,
                           includedRelations = NULL,
                           outputDir, restart = FALSE,
                           saveModels = TRUE, verbose = FALSE) {
  if (!inherits(fitorig, "nlmixr2FitCore")) stop("'fitorig' needs to be a nlmixr2 fit")
  if (saveModels && missing(outputDir)) stop("please specify outputDir for backward search")

  base_ui <- fitorig$finalUiEnv  # fixed — never updated
  resTableComplete <- NULL
  stepIdx <- 1

  # Determine starting model
  if (!missing(fitupdated)) {
    fwd_thetas  <- fitupdated$iniDf$name[!is.na(fitupdated$iniDf$ntheta)]
    orig_thetas <- fitorig$iniDf$name[!is.na(fitorig$iniDf$ntheta)]
    if (all(fwd_thetas %in% orig_thetas)) {
      # No forward covariates accepted. Skip backward unless includedRelations
      # forces something into the model first.
      if (is.null(includedRelations) || nrow(includedRelations) == 0) {
        cli::cli_alert_warning("no covariates added in the forward search, skipping backward search")
        return(list(fitorig, NULL))
      }
      cli::cli_alert_info("no forward covariates accepted; proceeding for includedRelations")
    }
    fit <- fitupdated
  } else {
    fit <- fitorig
  }

  if (reFitCovars) {
    xmod <- buildupatedUI(base_ui, unique(pairs$var), unique(pairs$covar), indep = FALSE, add = TRUE)
    fit  <- suppressWarnings(nlmixr2(xmod, data, fit$est))
    context_pairs <- .pairsInModel(pairs, fit)
  }

  # ── includedRelations ────────────────────────────────────────────────────────
  # Ensure all specified relations are in the model.  Relations already present
  # (from forward search) are left untouched.  Any that are missing are added
  # to context_pairs and the model is rebuilt from base and re-fit once before
  # the backward loop begins.  Included relations are then unioned into `pairs`
  # so that .pairsInModel() can find them and they are eligible for removal.
  if (!is.null(includedRelations) && nrow(includedRelations) > 0) {
    # Derive context_pairs now if not yet set (needed to rebuild below)
    if (is.null(context_pairs)) {
      context_pairs <- .pairsInModel(pairs, fit)
    }

    cur_thetas <- fit$iniDf$name[!is.na(fit$iniDf$ntheta)]
    inc_keys   <- paste0("cov_", includedRelations$covar, "_", includedRelations$var)
    missing_rel <- includedRelations[!inc_keys %in% cur_thetas, , drop = FALSE]

    if (nrow(missing_rel) > 0) {
      cli::cli_alert_info(
        "Adding {nrow(missing_rel)} includedRelation(s) not in current model: \\
        {paste(missing_rel$covar, missing_rel$var, sep='~', collapse=', ')}")

      # Merge missing relations into context_pairs and rebuild from base
      new_context <- rbind(context_pairs, missing_rel)
      cp_key      <- paste0(new_context$var, "_", new_context$covar)
      context_pairs <- new_context[!duplicated(cp_key), , drop = FALSE]

      xmod <- .rebuildUiFromPairs(base_ui, context_pairs)
      fit  <- suppressWarnings(nlmixr2(xmod, data, fit$est))
    }

    # Union includedRelations into the candidate pairs so .pairsInModel()
    # picks them up.  De-duplicate on (covar, var) key.
    pair_keys <- paste0(pairs$covar, "_", pairs$var)
    new_rows  <- includedRelations[!inc_keys %in% pair_keys, , drop = FALSE]
    if (nrow(new_rows) > 0) {
      pairs <- rbind(pairs, new_rows)
    }
  }

  # Restrict to pairs that are actually in the current model
  pairs <- .pairsInModel(pairs, fit)
  if (nrow(pairs) == 0) {
    cli::cli_alert_warning("no candidate pairs found in the model, skipping backward search")
    return(list(fit, NULL))
  }

  # Derive context_pairs from the model if not yet set
  if (is.null(context_pairs)) {
    context_pairs <- pairs
  }

  # When not saving there is no cache to read; treat as a clean start
  if (!saveModels) restart <- TRUE

  if (!restart) {
    fileExistsTab <- list.files(paste0("./", outputDir),
      pattern = "backward_step_[0-9]+_table_[a-zA-Z0-9_]+\\.rds")
    fileExistsFit <- list.files(paste0("./", outputDir),
      pattern = "backward_step_[0-9]+_fit_[a-zA-Z0-9_]+\\.rds")
    fileExistsCT  <- list.files(paste0("./", outputDir),
      pattern = "backward_step_[0-9]+_completetable_[a-zA-Z0-9_]+\\.rds")

    if (length(fileExistsTab) == 0) {
      restart <- TRUE
    } else {
      resumeRows  <- lapply(fileExistsTab, function(f) readRDS(paste0(outputDir, "/", f)))
      resumeTable <- data.table::rbindlist(resumeRows)
      fit           <- readRDS(paste0(outputDir, "/", fileExistsFit[[length(fileExistsFit)]]))
      context_pairs <- .pairsInModel(pairs, fit)
      resTableComplete <- readRDS(paste0(outputDir, "/", fileExistsCT[[length(fileExistsCT)]]))
      stepIdx <- unlist(resumeTable[nrow(resumeTable), ]$step) + 1
      for (i in seq_len(nrow(resumeTable))) {
        pairs <- pairs[!(pairs$var == resumeTable$var[i] & pairs$covar == resumeTable$covar[i]), ]
      }
      cli::cli_alert_success("loaded backward search data from disk, resuming ...")
    }
  }

  cli::cli_h1("starting backward search...")

  while (nrow(pairs) > 0) {
    # ── verbose: current model code + candidates before fitting ───────────
    if (verbose) {
      cli::cli_h2("Backward step {stepIdx}: current model")
      tryCatch(
        print(fit$finalUiEnv),
        error = function(e) cli::cli_alert_warning(
          "Model code unavailable: {conditionMessage(e)}")
      )
      cli::cli_h2("Backward step {stepIdx}: {nrow(pairs)} candidate(s) to test")
      v_cols <- intersect(
        c("var", "covar", "shape", "init", "covExpr"),
        names(pairs)
      )
      print(pairs[, v_cols, drop = FALSE])
    }

    results <- .fitCandidatePairs(pairs, base_ui, context_pairs, fit, data, pVal, stepIdx, add = FALSE)
    if (length(results) == 0) break

    resTable <- do.call(rbind, lapply(results, `[[`, "stats"))
    bestIdx  <- which.min(resTable$pchisqr)
    bestRow  <- resTable[bestIdx, , drop = FALSE]

    # ── verbose: show all candidate results ───────────────────────────────
    if (verbose) {
      cli::cli_h2("Backward step {stepIdx}: all candidate results")
      v_cols <- intersect(
        c("covar", "var", "shape", "objf", "deltObjf",
          "pchisqr", "AIC", "BIC", "bsvReduction", "covarEffect"),
        names(resTable)
      )
      print(resTable[order(resTable$pchisqr), v_cols, drop = FALSE])
    }

    if (bestRow$pchisqr <= pVal) {
      resTable[bestIdx, "included"] <- "yes"
      bestRow[, "included"] <- "yes"

      cli::cli_h1("best model at step {stepIdx}:")
      print(bestRow)

      # fit is now the reduced model (with the removed covariate dropped)
      fit <- results[[bestIdx]]$fit

      # ── verbose: updated model code + parameter estimates + omega ────────
      if (verbose) {
        cli::cli_h2("Retained model — updated model code:")
        tryCatch(
          print(fit$finalUiEnv),
          error = function(e) cli::cli_alert_warning(
            "Model code unavailable: {conditionMessage(e)}")
        )
        cli::cli_h2("Retained model — fixed effects:")
        print(fit$parFixedDf)
        cli::cli_h2("Retained model — omega diagonal (BSV variances):")
        print(round(diag(fit$omega), 4))
      }

      ret_var   <- results[[bestIdx]]$var
      ret_covar <- results[[bestIdx]]$covar
      pairs     <- pairs[!(pairs$var == ret_var & pairs$covar == ret_covar), ]

      # Update context_pairs to reflect the covariate that was removed
      context_pairs <- context_pairs[
        !(context_pairs$var == ret_var & context_pairs$covar == ret_covar), ,
        drop = FALSE
      ]

      resTableComplete <- rbind(resTableComplete, resTable)
      if (saveModels) {
        key <- paste0(ret_covar, "_", ret_var)
        saveRDS(fit,              paste0(outputDir, "/backward_step_", stepIdx, "_fit_",           key, ".rds"))
        saveRDS(bestRow,          paste0(outputDir, "/backward_step_", stepIdx, "_table_",         key, ".rds"))
        saveRDS(resTableComplete, paste0(outputDir, "/backward_step_", stepIdx, "_completetable_", key, ".rds"))
      }

      cli::cli_h2("retaining {ret_covar}~{ret_var}")
      stepIdx <- stepIdx + 1
    } else {
      cli::cli_h1("OFV did not improve, exiting backward search ...")
      resTableComplete <- rbind(resTableComplete, resTable)
      break
    }
  }

  cli::cli_h2(cli::col_red("backward search complete"))
  list(fit, resTableComplete)
}
