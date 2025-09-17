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
#'   effect, only, then it should not be included here.  If a
#'   parameter should have a separate estimate per level of a factor,
#'   give that as the two-sided formula here.
#'
#' @inheritDotParams nlmixr2est::nlmixr
#'
#' @param residualModel The residual model formula to use as a
#'   one-sided formula.
#'
#' @return The model fit from \code{nlmixr2est::nlmixr2()}
#'
#' @examples
#' nlmixrFormula(
#'   height ~ (Asym+AsymRe)+(R0-(Asym+AsymRe))*exp(-exp(lrc)*age) ~ (AsymRe|Seed),
#'   data = Loblolly,
#'   start = list(Asym = 103, R0 = -8.5, lrc = -3.3, addSd=1),
#'   est="focei"
#' )
#' @export
nlmixrFormula <- function(object, data, start, param=NULL, ..., residualModel=~add(addSd)) {
  parsedFormula <- .nlmixrFormulaParser(object)

  # Setup the random effects
  ranefGroup <- NULL
  for (currentRanef in parsedFormula$ranef) {
    if (is.null(ranefGroup)) {
      ranefGroup <- currentRanef$ranefGroup
    } else {
      # only one random effect grouping variable is allowed
      stopifnot(currentRanef$ranefGroup == ranefGroup)
    }
  }
  if (missing(data)) {
    data <- NULL
  }
  data <- .nlmixrFormulaDataPrep(data, dvName=parsedFormula$DV, idName=ranefGroup)
  startAll <- .nlmixrFormulaExpandStartParam(start=start, param=param, data=data)
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

  # Build up the model function
  modelOut <- function() {
    ini()
    model()
  }
  body(modelOut)[[2]][[2]] <- iniComplete
  body(modelOut)[[3]][[2]] <- modelComplete
  if (is.null(data)) {
    nlmixr2est::nlmixr(object=modelOut, ...)
  } else {
    nlmixr2est::nlmixr(object=modelOut, data=data, ...)
  }
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

#' Rename a column in a dataset, optionally overwriting it if the column does
#' not exist
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
  stopifnot(charOld %in% names(data))
  if (charOld != newName) {
    if (newName %in% names(data)) {
      warning(newName, " is in the data and will be overwritten with ", oldName)
    }
    data[[newName]] <- data[[charOld]]
  }
  stopifnot(newName %in% names(data))
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
  stopifnot(inherits(object, "formula"))
  # Confirm that the formula looks like it should

  # Confirm that it is a two-sided formula
  if (length(object) == 2) {
    stop("formula must be two-sided")
  }

  # break the formula up into its component parts
  stopifnot(length(object) == 3)
  stopifnot(identical(as.name("~"), object[[1]]))
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
  stopifnot(is.call(object))
  if (object[[1]] == as.name("+")) {
    # multiple random effects being parsed, separate them
    ret <-
      append(
        .nlmixrFormulaParserRanef(object[[2]]),
        .nlmixrFormulaParserRanef(object[[3]])
      )
  } else if (object[[1]] == as.name("(")) {
    stopifnot(length(object) == 2)
    # work inside the parentheses
    ret <- .nlmixrFormulaParserRanef(object[[2]])
  } else if (object[[1]] == as.name("|")) {
    stopifnot(is.call(object))
    stopifnot(length(object) == 3)
    stopifnot(is.name(object[[2]]))
    stopifnot(is.name(object[[3]]))
    ret <-
      list(
        list(
          ranefVar=object[[2]],
          ranefGroup=object[[3]],
          # Give the user control of start in the future
          start=1
        )
      )
  } else {
    stop("Invalid random effect: ", deparse(object))
  }
  ret
}

.paramExpand <- function(param) {
  if (is.null(param)) {
    character()
  } else if (is.list(param)) {
    unlist(lapply(X = param, FUN = .paramExpand))
  } else if (inherits(param, "formula")) {
    stopifnot(length(param) == 3)
    stopifnot(is.name(param[[3]]))
    varnames <- all.vars(param[[2]])
    setNames(rep(as.character(param[[3]]), length(varnames)), nm=varnames)
  } else {
    stop("Cannot expand param") # nocov
  }
}

#' Expand parameters to include their factor representations, if applicable.
#'
#' @param start the starting values for the model
#' @param startName The base name for the parameter
#' @param startValue The initial value for the base parameter
#' @param param The parameter in the model
#' @param data The dataset
#' @return A list with the ini and model parts needed to use the parameter
#' @keywords Internal
.nlmixrFormulaExpandStartParam <- function(start, param, data) {
  paramExpanded <- .paramExpand(param)
  stopifnot(all(names(paramExpanded) %in% names(start)))
  ret <- list()
  for (currentStart in names(start)) {
    ret <-
      append(
        ret,
        list(.nlmixrFormulaExpandStartParamSingle(
          startName=currentStart,
          startValue=start[[currentStart]],
          param=paramExpanded[currentStart], # this will be NA if it does not exist
          data=data
        ))
      )
  }
  ret
}

#' @describeIn dot-nlmixrFormulaExpandStartParam Expand a single parameter in a model using dataset factors, if applicable
.nlmixrFormulaExpandStartParamSingle <- function(startName, startValue, param, data) {
  stopifnot(is.character(startName))
  stopifnot(is.numeric(startValue))
  if (is.na(param)) {
    stopifnot(length(startValue) == 1)
    ret <-
      list(
        ini=list(str2lang(paste(startName, "<-", startValue))),
        model=list(NULL)
      )
  } else if (is.null(data)) {
    stop("data must be given when parameters are not single fixed effects")
  } else {
    stopifnot(is.character(param))
    stopifnot(length(param) == 1)
    stopifnot(param %in% names(data))
    if (any(is.na(data[[param]]))) {
      stop("NA found in data column: ", param)
    }
    if (is.factor(data[[param]])) {
      .nlmixrFormulaExpandStartParamFactor(startName, startValue, param, data)
    } else {
      stop("Can only handle factors for fixed effect grouping levels")
    }
  }
  ret
}

# Provide start and model setup for factor parameters
.nlmixrFormulaExpandStartParamFactor <- function(startName, startValue, param, data) {
  paramLabel <- levels(data[[param]])
  stopifnot(length(startValue) %in% c(1, length(paramLabel)))
  if (length(startValue) == 1 && !is.ordered(data[[param]])) {
    message("ordering the parameters by factor frequency: ", startName, " with parameter ", param)
    paramLabel <- names(rev(sort(summary(data[[param]]))))
  }
  ret <- list(ini=list(), model=list())
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
    ret$ini <-
      append(
        ret$ini,
        .nlmixrFormulaExpandStartParamSingle(startName=currentParamName, startValue=currentStartValue, param=NA)$ini
      )
    # Setup for mu-referencing where the base value is estimated for all levels
    # and other values are estimated as differences from the base
    modelStringCurrent <-
      if (idx == 1) {
        currentParamName
      } else {
        sprintf("%s*(%s == '%s')", currentParamName, param, paramLabel[idx])
      }
    modelString <- c(modelString, modelStringCurrent)
  }
  ret$model <-
    list(str2lang(sprintf(
      "%s <- %s", startName, paste(modelString, collapse="+")
    )))
  ret
}

#' Setup the ini() part of the model for fixed effects
#'
#' @param start The starting estimates for the model
#' @param base The initial basis for the ini definition
#' @return The inside of the ini() part of the model
#' @keywords Internal
.nlmixrFormulaSetupIniFixed <- function(start, base=str2lang("{}")) {
  stopifnot(length(start) > 0)
  stopifnot(is.list(start))
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
    stopifnot(is.list(ranefDefinition))
    for (currentRanef in ranefDefinition) {
      base[[length(base) + 1]] <-
        str2lang(paste(currentRanef$ranefVar, "~", currentRanef$start))
    }
  }
  base
}

#' Setup the model() part of the model
#'
#' @param start The starting estimates (used when fixed effects need more
#'   definition like for factors)
#' @param predictor The predictor from the formula
#' @param residualModel The residual model definition
#' @param predictorVar The variable in the data for the predictor
#' @param data The data used in the model
#' @return the interior of the model()
#' @keywords Internal
.nlmixrFormulaSetupModel <- function(start, predictor, residualModel, predictorVar="value", data) {
  stopifnot(inherits(residualModel, "formula"))
  # a one-sided formula
  stopifnot(length(residualModel) == 2)

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
