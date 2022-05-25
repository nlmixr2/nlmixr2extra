#' A simple formula-based interface for nlmixr2
#' 
#' @details
#' The formula is given with different notation than typical formulas.  The
#' formula notation is inspired by and similar to \code{lme4::nlmer()}.  It is a
#' 3-part formula: \code{dependent_variable~predictor_equation~random_effects}.
#' 
#' The \code{dependent_variable} is any variable in the dataset.  It may not
#' include any math; for example, \code{log(DV)} is not allowed.
#' 
#' The \code{predictor_equation} is any valid math, and it will be used directly
#' in the nlmixr2 model.
#' 
#' The \code{random_effects} are one or more random effect parameters defined by
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
#'   estimate may be given per factor level.  (Specify the factors with the \code{param} argument.)
#' @param ... Passed to \code{nlmixr()}
#' @param residual_model The residual model formula to use.
#' @return The model fit from \code{nlmixr()}
#' @examples
#' nlmixr_formula(
#'   height ~ (Asym+Asym_re)+(R0-(Asym+Asym_re))*exp(-exp(lrc)*age) ~ (Asym_re|Seed),
#'   data = Loblolly,
#'   start = list(Asym = 103, R0 = -8.5, lrc = -3.3, add_err=1),
#'   est="focei"
#' )
#' @export
nlmixr_formula <- function(object, data, start, param=NULL, ..., residual_model=~add(add_err)) {
  parsed_formula <- nlmixr_formula_parser(object)
  # Add "TIME", if needed to the data
  if (!("TIME" %in% names(data))) {
    data$TIME <- seq_len(nrow(data))
  }
  # Setup DV
  data <- rename_or_overwrite(data, new_name="DV", old_name=parsed_formula$DV)
  
  # Setup the random effects
  ranef_group <- NULL
  for (current_ranef in parsed_formula$ranef) {
    if (is.null(ranef_group)) {
      ranef_group <- current_ranef$ranef_group
    } else {
      # only one random effect grouping variable is allowed
      stopifnot(current_ranef$ranef_group == ranef_group)
    }
  }
  data <- rename_or_overwrite(data, new_name="ID", old_name=ranef_group)
  start_all <- nlmixr_formula_expand_start_param(start=start, param=param, data=data)
  ini_fixed <- nlmixr_formula_setup_ini_fixed(start=start_all)
  ini_complete <-
    nlmixr_formula_setup_ini_random(
      ranef_definition=parsed_formula$ranef,
      base=ini_fixed
    )
  model_complete <-
    nlmixr_formula_setup_model(
      start=start_all,
      predictor=parsed_formula$predictor,
      residual_model=residual_model,
      param=param
    )
  
  # Build up the model function  
  model_out <- function() {
    ini()
    model()
  }
  body(model_out)[[2]][[2]] <- ini_complete
  body(model_out)[[3]][[2]] <- model_complete
  nlmixr2est::nlmixr(object=model_out, data=data, ...)
}

rename_or_overwrite <- function(data, new_name, old_name) {
  char_old <- as.character(old_name)
  stopifnot(char_old %in% names(data))
  if (char_old != new_name) {
    if (new_name %in% names(data)) {
      warning(new_name, " is in the data and will be overwritten with ", old_name)
    }
    data[[new_name]] <- data[[char_old]]
  }
  stopifnot(new_name %in% names(data))
  data
}

nlmixr_formula_parser <- function(object) {
  stopifnot(inherits(object, "formula"))
  # Confirm that the formula looks like it should and break it up into its
  # component parts
  stopifnot(length(object) == 3)
  stopifnot(identical(as.name("~"), object[[1]]))
  ranef_part <- object[[3]]
  main_formula_part <- object[[2]]
  stopifnot(length(main_formula_part) == 3)
  stopifnot(identical(as.name("~"), main_formula_part[[1]]))
  dv_part <- main_formula_part[[2]]
  predictor_part <- main_formula_part[[3]]
  
  # We cannot yet handle DV that are not a NAME (no manipulations are done
  # within the data)
  stopifnot(is.name(dv_part))
  
  list(
    DV=dv_part,
    predictor=list(predictor_part),
    ranef=nlmixr_formula_parser_ranef(ranef_part)
  )
}

nlmixr_formula_parser_ranef <- function(object) {
  stopifnot(is.call(object))
  if (object[[1]] == as.name("+")) {
    # multiple random effects being parsed, separate them
    ret <-
      append(
        nlmixr_formula_parser_ranef(object[[2]]),
        nlmixr_formula_parser_ranef(object[[3]])
      )
  } else if (object[[1]] == as.name("(")) {
    stopifnot(length(object) == 2)
    inner_ranef_spec <- object[[2]]
    stopifnot(is.call(inner_ranef_spec))
    stopifnot(length(inner_ranef_spec) == 3)
    stopifnot(inner_ranef_spec[[1]] == as.name("|"))
    stopifnot(is.name(inner_ranef_spec[[2]]))
    stopifnot(is.name(inner_ranef_spec[[3]]))
    ret <-
      list(
        list(
          ranef_var=inner_ranef_spec[[2]],
          ranef_group=inner_ranef_spec[[3]],
          # Give the user control of start in the future
          start=1
        )
      )
  } else {
    stop("Invalid random effect") # nocov
  }
  ret
}

param_expand <- function(param) {
  if (is.null(param)) {
    character()
  } else if (is.list(param)) {
    unlist(lapply(param, param_expand))
  } else if (inherits(param, "formula")) {
    stopifnot(length(param) == 3)
    stopifnot(is.name(param[[3]]))
    varnames <- all.vars(param[[2]])
    setNames(rep(as.character(param[[3]]), length(varnames)), nm=varnames)
  } else {
    stop("Cannot expand param") # nocov
  }
}

nlmixr_formula_expand_start_param_single <- function(start_name, start_value, param, data) {
  stopifnot(is.character(start_name))
  stopifnot(is.numeric(start_value))
  if (is.na(param)) {
    stopifnot(length(start_value) == 1)
    ret <-
      list(
        ini=list(str2lang(paste(start_name, "<-", start_value))),
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
      param_label <- levels(data[[param]])
      stopifnot(length(start_value) %in% c(1, length(param_label)))
      ret <- list(ini=list(), model=list())
      model_string <- character()
      for (idx in seq_along(param_label)) {
        current_param_name <- make.names(paste(start_name, param, param_label[idx]))
        current_start_value <-
          if (length(start_value) == 1) {
            start_value
          } else {
            start_value[idx]
          }
        ret$ini <-
          append(
            ret$ini,
            nlmixr_formula_expand_start_param_single(start_name=current_param_name, start_value=current_start_value, param=NA)$ini
          )
        model_string <-
          c(model_string, sprintf("%s*(%s == '%s')", current_param_name, param, param_label[idx]))
      }
      ret$model <-
        list(str2lang(sprintf(
          "%s <- %s", start_name, paste(model_string, collapse="+")
        )))
    } else {
      stop("Can only handle factors for fixed effect grouping levels")
    }
  }
  ret
}

nlmixr_formula_expand_start_param <- function(start, param, data) {
  param_expanded <- param_expand(param)
  stopifnot(all(names(param_expanded) %in% names(start)))
  ret <- list()
  for (current_start in names(start)) {
    ret <-
      append(
        ret,
        list(nlmixr_formula_expand_start_param_single(
          start_name=current_start,
          start_value=start[[current_start]],
          param=param_expanded[current_start], # this will be NA if it does not exist
          data=data
        ))
      )
  }
  ret
}

nlmixr_formula_setup_ini_fixed <- function(start, base=str2lang("{}")) {
  stopifnot(length(start) > 0)
  stopifnot(is.list(start))
  for (idx_outer in seq_along(start)) {
    for (idx_inner in seq_along(start[[idx_outer]]$ini)) {
      base[[length(base) + 1]] <- start[[idx_outer]]$ini[[idx_inner]]
    }
  }
  base
}

nlmixr_formula_setup_ini_random <- function(ranef_definition, base=str2lang("{}")) {
  stopifnot(is.list(ranef_definition))
  for (current_ranef in ranef_definition) {
    base[[length(base) + 1]] <-
      str2lang(paste(current_ranef$ranef_var, "~", current_ranef$start))
  }
  base
}

nlmixr_formula_setup_model <- function(start, predictor, residual_model, predictor_var="value", param, data) {
  stopifnot(inherits(residual_model, "formula"))
  # a one-sided formula
  stopifnot(length(residual_model) == 2)
  
  pred_line <- str2lang(paste(predictor_var, " <- ."))
  pred_line[[3]] <- predictor[[1]]
  residual_line <- str2lang(paste(predictor_var, "~."))
  residual_line[[3]] <- residual_model[[2]]
  
  base <- str2lang("{}")
  # Add any lines needed from the 'start'
  for (current_start in start) {
    if (!is.null(current_start$model[[1]])) {
      base[[length(base) + 1]] <- current_start$model[[1]]
    }
  }
  base[[length(base) + 1]] <- pred_line
  base[[length(base) + 1]] <- residual_line
  base
}
