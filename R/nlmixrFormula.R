#' A simple formula-based interface for nlmixr2
#'
#' @details
#'
#' The formula is given with different notation than typical formulas.
#' The formula notation is inspired by and similar to
#' \code{lme4::nlmer()}.  It is a 3-part formula:
#' \code{dependentVariable~predictorEquation~randomEffects}.
#'
#' The \code{dependentVariable} is any variable in the dataset.  It
#' may not include any math; for example, \code{log(DV)} is not
#' allowed.
#'
#' The \code{predictorEquation} is any valid math, and it will be used
#' directly in the nlmixr2 model.
#'
#' The \code{randomEffects} are one or more random effect parameters
#' defined by putting the parameter in parentheses and putting a
#' vertical bar and the grouping parameter.  Only one grouping
#' parameter is allowed for all random effects.  An example would be
#' \code{(slope|ID)} to estimate a random effect parameter named
#' "slope" for each "ID" in the data.
#'
#' @param object The formula defining the model (see details)
#'
#' @param data The data to fit
#'
#' @param start A named list of starting estimates.  The names define
#'   the parameters in the model.  If a single parameter estimate is
#'   desired, it can be given here.  If a parameter estimate per
#'   factor level is desired, either a single starting estimate can be
#'   given across all factor levels or one estimate may be given per
#'   factor level.  (Specify the factors with the \code{param}
#'   argument.)
#'
#' @param param A formula or list of two-sided formulas giving the
#'   model used for parameters.  If a parameter is a simple fixed
#'   effect, only, then it should not be included here.  The
#'   right-hand side of a \code{param} formula names one or more
#'   columns in \code{data}: factor columns get a separate fixed effect
#'   per level; numeric columns get a linear slope. Multiple covariates
#'   on the same parameter can be combined with \code{+}, e.g.
#'   \code{b ~ z + w}.
#'
#' @param paramLink Optional named character vector mapping parameter
#'   names to a link function. \code{"identity"} (the default) emits
#'   \code{<param> <- <linear combination>}; \code{"log"} wraps the
#'   linear combination in \code{exp()} so the parameter is
#'   strictly-positive on the natural scale. Only parameters that have
#'   a \code{param} entry may appear here.
#'
#' @inheritDotParams nlmixr2est::nlmixr
#'
#' @param residualModel The residual model formula to use as a
#'   one-sided formula. The default is \code{~ add(addSd)}; richer
#'   models such as \code{~ add(addSd) + prop(propSd)} are supported as
#'   long as the corresponding sigma parameter names appear in
#'   \code{start}.
#'
#' @return The model fit from \code{nlmixr2est::nlmixr2()}
#'
#' @examples
#' \dontrun{
#' nlmixrFormula(
#'   height ~ (Asym+AsymRe)+(R0-(Asym+AsymRe))*exp(-exp(lrc)*age) ~ (AsymRe|Seed),
#'   data = Loblolly,
#'   start = list(Asym = 103, R0 = -8.5, lrc = -3.3, addSd=1),
#'   est="focei"
#' )
#' }
#' @export
nlmixrFormula <- function(object, data, start, param=NULL, paramLink=NULL, ..., residualModel=~add(addSd)) {
  if (missing(data)) {
    data <- NULL
  }
  built <- .nlmixrFormulaBuild(
    object=object, data=data, start=start, param=param,
    paramLink=paramLink, residualModel=residualModel
  )
  if (is.null(built$data)) {
    nlmixr2est::nlmixr(object=built$modelFun, ...)
  } else {
    nlmixr2est::nlmixr(object=built$modelFun, data=built$data, ...)
  }
}

#' Assemble the synthetic nlmixr2 model function for \code{nlmixrFormula}
#'
#' Performs the formula parse, data prep, parameter expansion, and ini/model
#' assembly without invoking \code{nlmixr2est::nlmixr}. Extracted so integration
#' tests can inspect the assembled model body without paying the cost of a fit
#' or simulation.
#'
#' @inheritParams nlmixrFormula
#' @return A list with elements \code{modelFun} (the synthetic
#'   \code{function() { ini(); model() }}), \code{data} (the data after
#'   \code{TIME}/\code{DV}/\code{ID} rename), \code{ini} (the assembled ini
#'   body), and \code{model} (the assembled model body).
#' @noRd
.nlmixrFormulaBuild <- function(object, data, start, param=NULL, paramLink=NULL, residualModel=~add(addSd)) {
  parsedFormula <- .nlmixrFormulaParser(object)

  # Setup the random effects
  ranefGroup <- NULL
  for (currentRanef in parsedFormula$ranef) {
    if (is.null(ranefGroup)) {
      ranefGroup <- currentRanef$ranefGroup
    } else if (!identical(currentRanef$ranefGroup, ranefGroup)) {
      # TODO: support multiple grouping variables (e.g. IOV (slope|id/occ));
      # see follow-up issue.
      stop(
        "Only one random-effect grouping variable is supported; got both '",
        deparse(ranefGroup), "' and '", deparse(currentRanef$ranefGroup), "'."
      )
    }
  }
  data <- .nlmixrFormulaDataPrep(data, dvName=parsedFormula$DV, idName=ranefGroup)
  startAll <- .nlmixrFormulaExpandStartParam(start=start, param=param, paramLink=paramLink, data=data)
  iniFixed <- .nlmixrFormulaSetupIniFixed(start=startAll)
  iniComplete <-
    .nlmixrFormulaSetupIniRandom(
      ranefDefinition=parsedFormula$ranef,
      base=iniFixed
    )
  modelComplete <-
    .nlmixrFormulaSetupModel(
      start=startAll,
      predictor=parsedFormula$predictor,
      residualModel=residualModel
    )

  modelFun <- function() {
    ini()
    model()
  }
  body(modelFun)[[2]][[2]] <- iniComplete
  body(modelFun)[[3]][[2]] <- modelComplete
  list(modelFun=modelFun, data=data, ini=iniComplete, model=modelComplete)
}

#' Perform any required data modifications for the nlmixrFormula interface
#'
#' @inheritParams nlmixr2est::nlmixr2
#' @param dvName,idName The name of the DV and ID columns for the dataset,
#'   respectively
#' @return A data frame modified, as needed for nlmixrFormula; if data is NULL,
#'   return NULL
.nlmixrFormulaDataPrep <- function(data, dvName, idName) {
  if (is.null(data)) return(NULL)
  # Add "TIME", if needed to the data
  if (!("TIME" %in% names(data))) {
    data$TIME <- seq_len(nrow(data))
  }
  # Setup DV
  data <- .renameOrOverwrite(data, newName="DV", oldName=dvName)
  # Setup ID
  if (length(idName) != 0) {
    data <- .renameOrOverwrite(data, newName="ID", oldName=idName)
  }
  data
}

#' Rename a column in a dataset
#'
#' If \code{newName} already exists in \code{data} and would be overwritten with
#' a different source column, this is treated as a hard error so that user data
#' is never silently clobbered.
#'
#' @param data The dataset to modify
#' @param newName,oldName The new and old column names
#' @return data with \code{data[[newName]] <- data[[charOld]]}
#' @keywords Internal
#' @examples
#' \dontrun{
#' .renameOrOverwrite(data.frame(A=1), newName="B", oldName="A")
#' }
.renameOrOverwrite <- function(data, newName, oldName) {
  charOld <- as.character(oldName)
  checkmate::assertChoice(charOld, choices=names(data))
  if (charOld != newName) {
    if (newName %in% names(data)) {
      stop(
        "Cannot rename column '", charOld, "' to '", newName,
        "': '", newName, "' already exists in the data. ",
        "Rename or remove the existing column first."
      )
    }
    data[[newName]] <- data[[charOld]]
  }
  data
}

#' Parse the formula to extract the dependent variable, predictor, and random
#' effects
#'
#' @param object The formula to parse
#' @return A list with names of "DV", "predictor", and "ranef" which are each
#'   part of the formula broken down into calls
#' @keywords Internal
#' @examples
#' \dontrun{
#' .nlmixrFormulaParser(a~b+c~(c|id))
#' }
.nlmixrFormulaParser <- function(object) {
  checkmate::assertClass(object, classes = "formula")
  # Confirm that the formula looks like it should

  # TODO: consider reformulas::splitForm for richer (v1+v2|g) and || syntax;
  # see follow-up issue.

  # Confirm that it is a two-sided formula
  if (length(object) == 2) {
    stop("formula must be two-sided")
  }

  # break the formula up into its component parts
  checkmate::assertTRUE(length(object) == 3)
  checkmate::assertTRUE(identical(as.name("~"), object[[1]]))
  # Split out the left and right-hand sides of the formula before determining if
  # there are random effects
  lhsFirst <- object[[2]]
  rhsFirst <- object[[3]]
  if (length(lhsFirst) == 3 && identical(as.name("~"), lhsFirst[[1]])) {
    # the model has random effects
    dvPart <- lhsFirst[[2]]
    predictorPart <- lhsFirst[[3]]
    ranefPart <- .nlmixrFormulaParserRanef(rhsFirst)
  } else {
    # the model does not have random effects
    dvPart <- lhsFirst
    predictorPart <- rhsFirst
    ranefPart <- NULL
  }

  # We cannot yet handle DV that are not a NAME (no manipulations are done
  # within the data)
  if (!is.name(dvPart)) {
    stop("formula left-hand-side must be a single variable, not ", deparse(dvPart))
  }

  list(
    DV=dvPart,
    predictor=list(predictorPart),
    ranef=ranefPart
  )
}

#' Parse the random effects part of a formula
#'
#' @param object The formula to parse
#' @return An unnamed list with one element per random effect.  The elements
#'   each have names of "ranefVar", "ranefGroup", and "start" indicating the
#'   variable with the random effect, the grouping variable for the random
#'   effect, and the starting value for the standard deviation of the random
#'   effect.
#' @keywords Internal
#' @examples
#' \dontrun{
#' .nlmixrFormulaParserRanef(str2lang("(c|id)+(d|id2)"))
#' }
.nlmixrFormulaParserRanef <- function(object) {
  if (!is.call(object)) {
    stop("Random-effect expression must be a call, not ", deparse(object))
  }
  if (object[[1]] == as.name("+")) {
    # multiple random effects being parsed, separate them
    ret <-
      append(
        .nlmixrFormulaParserRanef(object[[2]]),
        .nlmixrFormulaParserRanef(object[[3]])
      )
  } else if (object[[1]] == as.name("(")) {
    checkmate::assertTRUE(length(object) == 2)
    # work inside the parentheses
    ret <- .nlmixrFormulaParserRanef(object[[2]])
  } else if (object[[1]] == as.name("|")) {
    checkmate::assertTRUE(length(object) == 3)
    if (!is.name(object[[2]])) {
      stop("Random-effect variable must be a single name, not ", deparse(object[[2]]))
    }
    if (!is.name(object[[3]])) {
      stop("Random-effect grouping variable must be a single name, not ", deparse(object[[3]]))
    }
    ret <-
      list(
        list(
          ranefVar=object[[2]],
          ranefGroup=object[[3]],
          # TODO: allow start[[ranefVar]] to set this from the user; see
          # follow-up issue.
          start=1
        )
      )
  } else {
    stop("Invalid random effect: ", deparse(object))
  }
  ret
}

#' Expand the \code{param} argument into a named list of covariate columns
#'
#' @param param NULL, a two-sided formula, or a list of two-sided formulas.
#'   The right-hand side of each formula may be a single column name or
#'   multiple column names joined by \code{+}.
#' @return A named list of character vectors mapping parameter name to one or
#'   more covariate column names. Entries with the same name are concatenated.
#'   Returns an empty named list when \code{param} is NULL.
#' @noRd
.paramExpand <- function(param) {
  if (is.null(param)) {
    return(stats::setNames(list(), character()))
  }
  if (is.list(param)) {
    out <- stats::setNames(list(), character())
    for (entry in param) {
      out <- .paramExpandMerge(out, .paramExpand(entry))
    }
    return(out)
  }
  if (inherits(param, "formula")) {
    checkmate::assertTRUE(length(param) == 3)
    varnames <- all.vars(param[[2]])
    covCols <- .paramExpandRhs(param[[3]])
    out <- stats::setNames(
      replicate(length(varnames), covCols, simplify = FALSE),
      varnames
    )
    return(out)
  }
  stop("Cannot expand param") # nocov
}

#' Walk a \code{param} formula RHS and return its covariate column names
#'
#' Accepts a bare \code{name} or any \code{+}-joined tree of names. Anything
#' else (e.g. an interaction term \code{a*b}) is rejected.
#' @noRd
.paramExpandRhs <- function(rhs) {
  if (is.name(rhs)) {
    return(as.character(rhs))
  }
  if (is.call(rhs) && identical(rhs[[1]], as.name("+"))) {
    return(c(.paramExpandRhs(rhs[[2]]), .paramExpandRhs(rhs[[3]])))
  }
  stop(
    "Invalid right-hand side in `param`: ", deparse(rhs),
    ". Only single column names or `+`-joined names are supported."
  )
}

#' Merge two named lists of character vectors, concatenating duplicate names
#' @noRd
.paramExpandMerge <- function(a, b) {
  for (nm in names(b)) {
    a[[nm]] <- c(a[[nm]], b[[nm]])
  }
  a
}

#' Expand parameters to include their covariate representations, if applicable.
#'
#' @param start the starting values for the model
#' @param startName The base name for the parameter
#' @param startValue The initial value for the base parameter
#' @param param The parameter in the model
#' @param paramLink Named character vector giving the link function for each
#'   parameter that has a \code{param} entry. Recognised values: \code{"identity"}
#'   (default) emits \code{<param> <- <linear combination>}; \code{"log"} wraps
#'   the linear combination in \code{exp()}.
#' @param data The dataset
#' @return A list with the ini and model parts needed to use the parameter
#' @keywords Internal
.nlmixrFormulaExpandStartParam <- function(start, param, paramLink=NULL, data) {
  paramExpanded <- .paramExpand(param)
  checkmate::assertSubset(names(paramExpanded), choices=names(start))
  if (!is.null(paramLink)) {
    checkmate::assertCharacter(paramLink, names="named")
    checkmate::assertSubset(paramLink, choices=c("identity", "log"))
    checkmate::assertSubset(names(paramLink), choices=names(start))
  }
  ret <- list()
  for (currentStart in names(start)) {
    currentLink <- if (currentStart %in% names(paramLink)) {
      paramLink[[currentStart]]
    } else {
      NULL
    }
    ret <-
      append(
        ret,
        list(.nlmixrFormulaExpandStartParamSingle(
          startName=currentStart,
          startValue=start[[currentStart]],
          param=paramExpanded[[currentStart]],
          link=currentLink,
          data=data
        ))
      )
  }
  ret
}

#' @describeIn dot-nlmixrFormulaExpandStartParam Expand a single parameter in a
#'   model using dataset covariates, if applicable
#' @param link Optional link function for this parameter (\code{"identity"} or
#'   \code{"log"}). Defaults to \code{"identity"}.
.nlmixrFormulaExpandStartParamSingle <- function(startName, startValue, param=NULL, link=NULL, data=NULL) {
  checkmate::assertString(startName)
  checkmate::assertNumeric(startValue)
  if (is.null(link) || identical(link, NA_character_)) {
    link <- "identity"
  }
  checkmate::assertChoice(link, choices=c("identity", "log"))

  if (length(param) == 0 || (length(param) == 1 && is.na(param))) {
    # Scalar fixed effect with no covariate model
    checkmate::assertNumeric(startValue, len=1)
    if (!identical(link, "identity")) {
      stop(
        "paramLink for '", startName, "' is '", link,
        "' but the parameter has no covariate model in `param`."
      )
    }
    return(list(
      ini=list(str2lang(paste(startName, "<-", startValue))),
      model=list(NULL)
    ))
  }

  if (is.null(data)) {
    stop("data must be given when parameters are not single fixed effects")
  }
  checkmate::assertCharacter(param, min.len=1)
  checkmate::assertSubset(param, choices=names(data))

  # Build one fragment per covariate column. A factor fragment emits its own
  # baseline level so we *omit* the global pop.<startName> intercept when any
  # factor is present; a sole continuous covariate carries its own intercept.
  #
  # startValue is consumed as a queue: each factor takes n_levels values, each
  # continuous takes one slope. When only one covariate exists we accept the
  # familiar scalar shorthands (factor: length 1 = baseline-only; continuous:
  # length 1 = intercept-only with slope=0; length 2 = c(intercept, slope)).
  hasFactor <- any(vapply(param, function(col) is.factor(data[[col]]), logical(1)))
  singleCovariate <- length(param) == 1

  # Pre-validate covariate types (also generates the per-covariate consume size)
  consumeSizes <- integer(length(param))
  for (idx in seq_along(param)) {
    col <- param[idx]
    if (any(is.na(data[[col]]))) {
      stop("NA found in data column: ", col)
    }
    if (is.factor(data[[col]])) {
      consumeSizes[idx] <- length(levels(data[[col]]))
    } else if (is.numeric(data[[col]])) {
      consumeSizes[idx] <- 1L
    } else if (is.character(data[[col]])) {
      stop(
        "Column '", col, "' in `data` is character; convert it to a factor ",
        "(e.g. data$", col, " <- factor(data$", col, ")) before passing to ",
        "nlmixrFormula()."
      )
    } else {
      stop(
        "Unsupported column type for parameter covariate '", col,
        "': ", paste(class(data[[col]]), collapse="/")
      )
    }
  }

  if (singleCovariate) {
    # Defer length validation to the single-covariate helper so the existing
    # scalar shorthands keep working.
    sliceList <- list(startValue)
  } else {
    expectedLen <- sum(consumeSizes) + if (hasFactor) 0L else 1L
    if (length(startValue) != expectedLen) {
      stop(
        "For parameter '", startName, "' with covariates ",
        paste0("'", param, "'", collapse=", "),
        ", start must have length ", expectedLen,
        " (one per factor level plus one slope per continuous covariate",
        if (!hasFactor) " plus one intercept" else "",
        "); got ", length(startValue), "."
      )
    }
    # Slice up startValue: factors take their level count, continuous take 1.
    sliceList <- vector("list", length(param))
    pos <- 1L
    for (idx in seq_along(param)) {
      n <- consumeSizes[idx]
      sliceList[[idx]] <- startValue[seq.int(pos, length.out = n)]
      pos <- pos + n
    }
  }

  iniLines <- list()
  rhsFragments <- character()
  for (idx in seq_along(param)) {
    col <- param[idx]
    slice <- sliceList[[idx]]
    if (is.factor(data[[col]])) {
      part <- .nlmixrFormulaExpandStartParamFactor(startName, slice, col, data)
    } else {
      # Continuous covariate
      if (singleCovariate) {
        part <- .nlmixrFormulaExpandStartParamContinuous(
          startName, slice, col, includeIntercept = TRUE
        )
      } else {
        # In the mixed case a factor already supplies the intercept; pass only
        # the slope to the continuous helper.
        part <- .nlmixrFormulaExpandStartParamContinuous(
          startName, slice, col, includeIntercept = FALSE
        )
      }
    }
    iniLines <- c(iniLines, part$ini)
    rhsFragments <- c(rhsFragments, part$rhs)
  }

  rhs <- paste(rhsFragments, collapse=" + ")
  if (identical(link, "log")) {
    rhs <- paste0("exp(", rhs, ")")
  }
  modelLine <- str2lang(sprintf("%s <- %s", startName, rhs))

  list(ini=iniLines, model=list(modelLine))
}

# Provide start and ini lines for one factor covariate on a parameter. Returns
# a list with `ini` (language objects to insert into the ini block) and `rhs`
# (a string fragment to be combined into the model line by the dispatcher).
.nlmixrFormulaExpandStartParamFactor <- function(startName, startValue, param, data) {
  paramLabel <- levels(data[[param]])
  checkmate::assertNumeric(startValue)
  if (!length(startValue) %in% c(1, length(paramLabel))) {
    stop(
      "For factor covariate '", param, "' on parameter '", startName,
      "', length(start[['", startName, "']]) must be 1 or ",
      length(paramLabel), " (number of factor levels); got ",
      length(startValue), "."
    )
  }
  if (length(startValue) == 1 && !is.ordered(data[[param]])) {
    message("ordering the parameters by factor frequency: ", startName, " with parameter ", param)
    # decreasing=TRUE gives a stable sort: most-frequent first, with original
    # factor level order preserved on ties
    paramLabel <- names(sort(summary(data[[param]]), decreasing = TRUE))
  }
  iniLines <- list()
  modelString <- character()
  for (idx in seq_along(paramLabel)) {
    currentParamName <- make.names(paste(startName, param, paramLabel[idx]))
    currentStartValue <-
      if (length(startValue) == 1 && idx == 1) {
        startValue
      } else if (length(startValue) == 1) {
        0
      } else {
        startValue[idx]
      }
    iniLines <- c(iniLines, list(str2lang(paste(currentParamName, "<-", currentStartValue))))
    # Setup for mu-referencing where the base value is estimated for all levels
    # and other values are estimated as differences from the base
    modelStringCurrent <-
      if (idx == 1) {
        currentParamName
      } else {
        sprintf('%s * (%s == "%s")', currentParamName, param, paramLabel[idx])
      }
    modelString <- c(modelString, modelStringCurrent)
  }
  list(
    ini=iniLines,
    rhs=paste(modelString, collapse=" + ")
  )
}

#' Provide start and ini lines for one continuous covariate on a parameter
#'
#' Emits a slope parameter named \code{cov_<param>_<startName>} (matching the
#' naming used elsewhere in nlmixr2extra; see \code{R/parsingutil.R:30}) and,
#' when \code{includeIntercept = TRUE}, an intercept named
#' \code{pop.<startName>}. The intercept is suppressed when a factor covariate
#' on the same parameter already provides one, to avoid an unidentifiable
#' duplicate intercept.
#'
#' \code{startValue} of length 1 is taken as the intercept (slope defaults to
#' 0). Length 2 is taken as \code{c(intercept, slope)}.
#'
#' @noRd
.nlmixrFormulaExpandStartParamContinuous <- function(startName, startValue, param, includeIntercept) {
  checkmate::assertNumeric(startValue)
  checkmate::assertString(param)
  checkmate::assertFlag(includeIntercept)

  if (includeIntercept) {
    if (!length(startValue) %in% c(1, 2)) {
      stop(
        "For continuous covariate '", param, "' on parameter '", startName,
        "', length(start[['", startName, "']]) must be 1 (intercept only) ",
        "or 2 (intercept, slope); got ", length(startValue), "."
      )
    }
    intercept <- startValue[1]
    slope <- if (length(startValue) == 2) startValue[2] else 0
  } else {
    # The intercept is supplied by another covariate (typically a factor) on
    # the same parameter; this helper only contributes the slope term.
    if (length(startValue) != 1) {
      stop(
        "When combined with another covariate on parameter '", startName,
        "', continuous covariate '", param, "' must contribute a single ",
        "slope value; got ", length(startValue), "."
      )
    }
    slope <- startValue[1]
  }

  covName <- paste0("cov_", param, "_", startName)
  iniLines <- list(str2lang(paste(covName, "<-", slope)))
  rhs <- sprintf("%s * %s", covName, param)

  if (includeIntercept) {
    popName <- paste0("pop.", startName)
    iniLines <- c(
      list(str2lang(paste(popName, "<-", intercept))),
      iniLines
    )
    rhs <- paste(popName, "+", rhs)
  }

  list(ini=iniLines, rhs=rhs)
}

#' Setup the ini() part of the model for fixed effects
#'
#' @param start The starting estimates for the model
#' @param base The initial basis for the ini definition
#' @return The inside of the ini() part of the model
#' @keywords Internal
.nlmixrFormulaSetupIniFixed <- function(start, base=str2lang("{}")) {
  checkmate::assertList(start, min.len = 1)
  for (idxOuter in seq_along(start)) {
    for (idxInner in seq_along(start[[idxOuter]]$ini)) {
      base[[length(base) + 1]] <- start[[idxOuter]]$ini[[idxInner]]
    }
  }
  base
}

#' @describeIn dot-nlmixrFormulaSetupIniFixed Setup the ini() part of the model for fixed effects
#' @param ranefDefinition The random effect definitions
.nlmixrFormulaSetupIniRandom <- function(ranefDefinition, base=str2lang("{}")) {
  if (!is.null(ranefDefinition)) {
    checkmate::assertList(ranefDefinition)
    for (currentRanef in ranefDefinition) {
      base[[length(base) + 1]] <-
        str2lang(paste(currentRanef$ranefVar, "~", currentRanef$start))
    }
  }
  base
}

#' Setup the model() part of the model
#'
#' Each parameter's per-level model line (if any) was already built upstream by
#' \code{\link{.nlmixrFormulaExpandStartParamSingle}}; this function only
#' concatenates those lines with the predictor and residual lines.
#'
#' @param start The starting estimates (used when fixed effects need
#'   more definition like for factors)
#' @param predictor The predictor from the formula
#' @param residualModel The residual model definition
#' @param predictorVar The variable in the data for the predictor (the left hand
#'   side of the equation)
#' @return the interior of the model()
#' @noRd
#' @author William Denney
.nlmixrFormulaSetupModel <- function(start, predictor, residualModel, predictorVar="value") {
  checkmate::assertClass(residualModel, classes = "formula")
  if (length(residualModel) != 2) {
    stop("residualModel must be a one-sided formula (e.g. ~ add(addSd))")
  }

  predLine <- str2lang(paste(predictorVar, " <- ."))
  predLine[[3]] <- predictor[[1]]
  residualLine <- str2lang(paste(predictorVar, "~."))
  residualLine[[3]] <- residualModel[[2]]

  base <- str2lang("{}")
  # Add any lines needed from the 'start'
  for (currentStart in start) {
    if (!is.null(currentStart$model[[1]])) {
      base[[length(base) + 1]] <- currentStart$model[[1]]
    }
  }
  base[[length(base) + 1]] <- predLine
  base[[length(base) + 1]] <- residualLine
  base
}

#' @export
nlmixr2.formula <- function(object, data=NULL, est = NULL, control = NULL,
                            table = nlmixr2est::tableControl(), ...,
                            start=NULL, param=NULL, paramLink=NULL,
                            residualModel=~add(addSd),
                            save = NULL, envir = parent.frame()) {
  .lst <- c(list(object=object,
                 data=data,
                 start=start,
                 param=param,
                 paramLink=paramLink,
                 est=est,
                 control=control,
                 table=table),
            list(...),
            list(save=save,
                 envir=envir,
                 residualModel=residualModel))
  do.call("nlmixrFormula", .lst)
}
