#' A simple formula-based interface for nlmixr2
#' 
#' @details
#' The formula is given with different notation than typical formulas.  The
#' formula notation is inspired by and similar to \code{lme4::nlmer()}.  It is a
#' 3-part formula: \code{dependentVariable~predictorEquation~randomEffects}.
#' 
#' The \code{dependentVariable} is any variable in the dataset.  It may not
#' include any math; for example, \code{log(DV)} is not allowed.
#' 
#' The \code{predictorEquation} is any valid math, and it will be used directly
#' in the nlmixr2 model.
#' 
#' The \code{randomEffects} are one or more random effect parameters defined by
#' putting the parameter in parentheses and putting a vertical bar and the
#' grouping parameter.  Only one grouping parameter is allowed for all random
#' effects.  An example would be \code{(slope|ID)} to estimate a random effect
#' parameter named "slope" for each "ID" in the data.
#' 
#' @param object The formula defining the model (see details)
#' @param data The data to fit
#' @param start A named list of starting estimates.  The names define the
#'   parameters in the model.  If a single parameter estimate is desired, it can
#'   be given here.  If a parameter estimate per factor level is desired, either
#'   a single starting estimate can be given across all factor levels or one
#'   estimate may be given per factor level.  (Specify the factors with the
#'   \code{param} argument.)
#' @param param A formula or list of two-sided formulas giving the model used
#'   for parameters.  If a parameter is a simple fixed effect, only, then it
#'   should not be included here.  If a parameter should have a separate
#'   estimate per level of a factor, give that as the two-sided formula here.
#' @inheritDotParams nlmixr2::nlmixr
#' @param residualModel The residual model formula to use as a one-sided formula.
#' @return The model fit from \code{nlmixr()}
#' @examples
#' nlmixrFormula(
#'   height ~ (Asym+AsymRe)+(R0-(Asym+AsymRe))*exp(-exp(lrc)*age) ~ (AsymRe|Seed),
#'   data = Loblolly,
#'   start = list(Asym = 103, R0 = -8.5, lrc = -3.3, addErr=1),
#'   est="focei"
#' )
#' @export
nlmixrFormula <- function(object, data, start, param=NULL, ..., residualModel=~add(addErr)) {
  parsedFormula <- .nlmixrFormulaParser(object)
  # Add "TIME", if needed to the data
  if (!("TIME" %in% names(data))) {
    data$TIME <- seq_len(nrow(data))
  }
  # Setup DV
  data <- .renameOrOverwrite(data, newName="DV", oldName=parsedFormula$DV)
  
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
  data <- .renameOrOverwrite(data, newName="ID", oldName=ranefGroup)
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
      residualModel=residualModel,
      param=param
    )
  
  # Build up the model function  
  modelOut <- function() {
    ini()
    model()
  }
  body(modelOut)[[2]][[2]] <- iniComplete
  body(modelOut)[[3]][[2]] <- modelComplete
  nlmixr2est::nlmixr(object=modelOut, data=data, ...)
}

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

.nlmixrFormulaParser <- function(object) {
  stopifnot(inherits(object, "formula"))
  # Confirm that the formula looks like it should and break it up into its
  # component parts
  stopifnot(length(object) == 3)
  stopifnot(identical(as.name("~"), object[[1]]))
  ranefPart <- object[[3]]
  mainFormulaPart <- object[[2]]
  stopifnot(length(mainFormulaPart) == 3)
  stopifnot(identical(as.name("~"), mainFormulaPart[[1]]))
  dvPart <- mainFormulaPart[[2]]
  predictorPart <- mainFormulaPart[[3]]
  
  # We cannot yet handle DV that are not a NAME (no manipulations are done
  # within the data)
  stopifnot(is.name(dvPart))
  
  list(
    DV=dvPart,
    predictor=list(predictorPart),
    ranef=.nlmixrFormulaParserRanef(ranefPart)
  )
}

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
    innerRanefSpec <- object[[2]]
    stopifnot(is.call(innerRanefSpec))
    stopifnot(length(innerRanefSpec) == 3)
    stopifnot(innerRanefSpec[[1]] == as.name("|"))
    stopifnot(is.name(innerRanefSpec[[2]]))
    stopifnot(is.name(innerRanefSpec[[3]]))
    ret <-
      list(
        list(
          ranefVar=innerRanefSpec[[2]],
          ranefGroup=innerRanefSpec[[3]],
          # Give the user control of start in the future
          start=1
        )
      )
  } else {
    stop("Invalid random effect") # nocov
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
  } else {
    stopifnot(is.character(param))
    stopifnot(length(param) == 1)
    stopifnot(param %in% names(data))
    if (any(is.na(data[[param]]))) {
      stop("NA found in data column:", param)
    }
    if (is.factor(data[[param]])) {
      paramLabel <- levels(data[[param]])
      stopifnot(length(startValue) %in% c(1, length(paramLabel)))
      ret <- list(ini=list(), model=list())
      modelString <- character()
      for (idx in seq_along(paramLabel)) {
        currentParamName <- make.names(paste(startName, param, paramLabel[idx]))
        currentStartValue <-
          if (length(startValue) == 1) {
            startValue
          } else {
            startValue[idx]
          }
        ret$ini <-
          append(
            ret$ini,
            .nlmixrFormulaExpandStartParamSingle(startName=currentParamName, startValue=currentStartValue, param=NA)$ini
          )
        modelString <-
          c(modelString, sprintf("%s*(%s == '%s')", currentParamName, param, paramLabel[idx]))
      }
      ret$model <-
        list(str2lang(sprintf(
          "%s <- %s", startName, paste(modelString, collapse="+")
        )))
    } else {
      stop("Can only handle factors for fixed effect grouping levels")
    }
  }
  ret
}

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

.nlmixrFormulaSetupIniRandom <- function(ranefDefinition, base=str2lang("{}")) {
  stopifnot(is.list(ranefDefinition))
  for (currentRanef in ranefDefinition) {
    base[[length(base) + 1]] <-
      str2lang(paste(currentRanef$ranefVar, "~", currentRanef$start))
  }
  base
}

.nlmixrFormulaSetupModel <- function(start, predictor, residualModel, predictorVar="value", param, data) {
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
