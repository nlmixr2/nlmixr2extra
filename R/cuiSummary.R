#' Endpoint summaries for a simulated clinical utility index
#'
#' Each summary reduces a simulation to one number per regimen: the mean
#' over virtual subjects of a per-subject value.
#'
#' - `cuiProb(x)`: mean of a simulated probability (must be in `[0, 1]`).
#' - `cuiMean(x)`: mean of a simulated value at a time, or with `window`
#'   the mean over subjects of each subject's time-weighted average (for
#'   example Cavg).
#' - `cuiAuc(x, window)`: mean over subjects of the integral of `x` over
#'   `window`.
#' - `cuiThreshold(cond)`: fraction of subjects for which `cond` is `TRUE`.
#'
#' `x` and `cond` are a column of the model's simulation (bare name or
#' string), or an expression of columns and constants, such as
#' `bmdchg > -1`.  A string matching a column exactly is that column (for
#' nonsyntactic names), otherwise it is parsed as an expression.  Names
#' that are model columns win over variables of the same name in the
#' calling environment (with a warning); write `.env$cut` for the variable
#' or `.data$cut` for the column.  Other names must be single constants in
#' the calling environment; they are replaced by their values when [cui()]
#' runs, so the result keeps no environment.
#'
#' Window integrals are computed exactly by adding an integrator state to
#' the simulation model, so `x` must be valid in the model.
#'
#' @param x Column or expression, as described above.
#' @param cond Logical expression of columns and constants.
#' @param ... Must be empty; `time` and `window` are given by name.
#' @param time Single time (model time units) to evaluate at; `NULL` uses
#'   each assessment time given to [cui()].
#' @param window Time interval `c(from, to)` for averages and integrals.
#' @param .x,.cond A string or language object instead of `x` / `cond`,
#'   for programmatic use.
#' @return A `cuiSummary` object for [cuiEndpoint()].
#' @family Clinical utility
#' @export
#' @examples
#' cuiProb(presp)
#' cuiThreshold(bmdchg > -1, time = 2016)
#' cuiMean(cp, window = c(0, 168))
#' cut <- -1
#' cuiThreshold(bmdchg > .env$cut)
cuiProb <- function(x, ..., time = NULL, .x = NULL) {
  .cuiNoDots(...)
  .cuiSummaryNew("prob", .cuiCaptureExpr(substitute(x), missing(x), .x, parent.frame(), "x"),
                 parent.frame(), time = time)
}

#' @rdname cuiProb
#' @export
cuiMean <- function(x, ..., time = NULL, window = NULL, .x = NULL) {
  .cuiNoDots(...)
  if (!is.null(time) && !is.null(window)) {
    stop("give 'time' or 'window', not both", call. = FALSE)
  }
  .cuiSummaryNew("mean", .cuiCaptureExpr(substitute(x), missing(x), .x, parent.frame(), "x"),
                 parent.frame(), time = time, window = window)
}

#' @rdname cuiProb
#' @export
cuiAuc <- function(x, ..., window, .x = NULL) {
  .cuiNoDots(...)
  if (missing(window)) stop("'window' is required", call. = FALSE)
  .cuiSummaryNew("auc", .cuiCaptureExpr(substitute(x), missing(x), .x, parent.frame(), "x"),
                 parent.frame(), window = window)
}

#' @rdname cuiProb
#' @param dose For `cuiResponse()`: the column of a fitted regression model
#'   (see [cuiModelPredict()]) that takes the regimen dose, as a bare name or
#'   a string, or missing for a model driven only by covariates and chained
#'   exposures; `cuiResponse()` is the mean predicted response of such a
#'   model over the population.
#' @param .dose The dose column as a string, for programmatic use.
#' @export
cuiResponse <- function(dose, ..., .dose = NULL) {
  .cuiNoDots(...)
  .doseVar <- if (missing(dose) && is.null(.dose)) {
    NULL
  } else {
    .cuiCaptureName(substitute(dose), missing(dose), .dose, "dose")
  }
  .s <- .cuiSummaryNew("response", NULL, NULL)
  .s$doseVar <- .doseVar
  .s
}

#' @rdname cuiProb
#' @export
cuiThreshold <- function(cond, ..., time = NULL, .cond = NULL) {
  .cuiNoDots(...)
  .cuiSummaryNew("threshold",
                 .cuiCaptureExpr(substitute(cond), missing(cond), .cond, parent.frame(), "cond"),
                 parent.frame(), time = time)
}

.cuiNoDots <- function(...) {
  if (...length() > 0) {
    stop("unexpected argument(s); give 'time' and 'window' by name", call. = FALSE)
  }
}

.cuiSummaryNew <- function(type, expr, env, time = NULL, window = NULL) {
  if (!is.null(time)) checkmate::assertNumber(time, finite = TRUE, .var.name = "time")
  if (!is.null(window)) {
    checkmate::assertNumeric(window, len = 2, finite = TRUE, any.missing = FALSE,
                             .var.name = "window")
    if (window[1] >= window[2]) stop("'window' must be c(from, to) with from < to", call. = FALSE)
  }
  structure(list(type = type, expr = expr, env = env, time = time, window = window),
            class = "cuiSummary")
}

#' @export
format.cuiSummary <- function(x, ...) {
  if (x$type == "response") {
    return(if (is.null(x$doseVar)) "response()" else sprintf("response(%s = dose)", x$doseVar))
  }
  .e <- if (is.character(x$expr)) x$expr else paste(deparse(x$expr), collapse = " ")
  .when <- if (!is.null(x$window)) {
    sprintf(" over [%s, %s]", format(x$window[1]), format(x$window[2]))
  } else if (!is.null(x$time)) {
    sprintf(" at %s", format(x$time))
  } else {
    ""
  }
  sprintf("%s(%s)%s", x$type, .e, .when)
}

#' @export
print.cuiSummary <- function(x, ...) {
  cat("cuiSummary: ", format(x), "\n", sep = "")
  invisible(x)
}

#' Check a value substituted for a name is a single constant
#' @noRd
.cuiCheckConst <- function(val, nm, ep = NULL) {
  if (is.function(val) || !is.atomic(val) || length(val) != 1L ||
    !(is.numeric(val) || is.character(val) || is.logical(val))) {
    stop(sprintf("'%s'%s must be a model column or a single numeric, character or logical constant",
                 nm, if (is.null(ep)) "" else sprintf(" in endpoint '%s'", ep)), call. = FALSE)
  }
  val
}

#' Apply a function to the arguments of a call, skipping empty arguments
#' @noRd
.cuiMapArgs <- function(e, f) {
  if (length(e) < 2L) return(e)
  for (.i in 2:length(e)) {
    if (identical(e[[.i]], quote(expr = ))) next
    .new <- f(e[[.i]])
    if (!identical(.new, e[[.i]])) e[[.i]] <- .new
  }
  e
}

.cuiIsDollar <- function(e, lhs) {
  is.call(e) && identical(e[[1]], as.name("$")) && identical(e[[2]], as.name(lhs))
}

#' Capture a summary expression given as a name, call or string
#'
#' Strings are kept as they are until the model columns are known;
#' `.env$name` references are replaced by their values now.
#' @noRd
.cuiCaptureExpr <- function(expr, isMissing, se, env, what) {
  if (!is.null(se)) {
    if (!isMissing) {
      stop(sprintf("supply only one of '%s' and '.%s'", what, what), call. = FALSE)
    }
    expr <- se
  } else if (isMissing) {
    stop(sprintf("'%s' is required", what), call. = FALSE)
  }
  if (is.character(expr)) {
    checkmate::assertString(expr, min.chars = 1, .var.name = what)
    return(expr)
  }
  if (!is.symbol(expr) && !is.call(expr)) {
    stop(sprintf("'%s' must be a column name, an expression or a string", what), call. = FALSE)
  }
  .cuiFreezeEnv(expr, env)
}

.cuiFreezeEnv <- function(e, env) {
  if (.cuiIsDollar(e, ".env")) {
    .nm <- as.character(e[[3]])
    return(.cuiCheckConst(get0(.nm, envir = env, inherits = TRUE), .nm))
  }
  if (is.call(e)) e <- .cuiMapArgs(e, function(a) .cuiFreezeEnv(a, env))
  e
}

#' Resolve a summary against the model's columns, dropping its environment
#'
#' @param s A `cuiSummary`.
#' @param cols Columns of the model's simulation.
#' @param ep Endpoint name, for messages.
#' @return The summary with a literal expression and `env = NULL`.
#' @noRd
.cuiResolveSummary <- function(s, cols, ep) {
  .e <- s$expr
  if (is.character(.e)) {
    .e <- if (.e %in% cols) as.name(.e) else str2lang(.e)
    .e <- .cuiFreezeEnv(.e, s$env)
  }
  if (s$type %in% c("prob", "mean", "auc") && is.symbol(.e) && !(as.character(.e) %in% cols)) {
    stop(sprintf("'%s' in endpoint '%s' is not a model column; columns: %s",
                 as.character(.e), ep, paste(cols, collapse = ", ")), call. = FALSE)
  }
  s$expr <- .cuiResolveSyms(.e, cols, s$env, ep)
  s$env <- NULL
  s
}

.cuiResolveSyms <- function(e, cols, env, ep) {
  if (is.symbol(e)) {
    .nm <- as.character(e)
    .val <- if (is.null(env)) NULL else get0(.nm, envir = env, inherits = TRUE)
    if (.nm %in% cols) {
      if (!is.null(.val) && !is.function(.val)) {
        warning(sprintf(paste0("'%s' in endpoint '%s' is a model column and also a variable; ",
                               "using the column (write .env$%s for the variable)"),
                        .nm, ep, .nm), call. = FALSE)
      }
      return(e)
    }
    if (is.null(.val)) {
      stop(sprintf("'%s' in endpoint '%s' is not a model column or a constant; columns: %s",
                   .nm, ep, paste(cols, collapse = ", ")), call. = FALSE)
    }
    return(.cuiCheckConst(.val, .nm, ep))
  }
  if (.cuiIsDollar(e, ".data")) {
    .nm <- as.character(e[[3]])
    if (!(.nm %in% cols)) {
      stop(sprintf("'.data$%s' in endpoint '%s' is not a model column", .nm, ep), call. = FALSE)
    }
    return(as.name(.nm))
  }
  if (is.call(e)) e <- .cuiMapArgs(e, function(a) .cuiResolveSyms(a, cols, env, ep))
  e
}
