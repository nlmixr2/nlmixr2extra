#' Use fitted regression models as clinical utility index endpoints
#'
#' Fitted models other than nlmixr2 fits (for example `glm()`, `lm()`,
#' `lme4::glmer()`, `glmmTMB::glmmTMB()` or `brms::brm()`) can be CUI
#' endpoints, in [cuiDoseResponse()] (`cuiBinary(tox, model = fit)`) and in
#' [cui()] (`cuiEndpoint(eff, fit, cuiResponse(dose = DOSE))`).  They are
#' used through four generics, so any model class works once it has
#' methods:
#'
#' - `cuiModelPredict(fit, newdata, draw)`: the predicted mean response per
#'   row of `newdata` (a probability for a binary endpoint), at the
#'   estimates when `draw` is `NULL`, otherwise at one parameter draw from
#'   `cuiModelDraws()`.  Predictions are population-level (random effects
#'   at zero).
#' - `cuiModelDraws(fit, n)`: `n` parameter draws, used for
#'   `uncertainty = "cov"` in [cui()]: coefficients from a multivariate
#'   normal on the estimates and their covariance, or posterior draws for
#'   `brmsfit`.  Called inside a seeded stream.
#' - `cuiModelRefit(fit, data)`: the model refit to new data, used by the
#'   bootstrap of [cuiDoseResponse()]; the default calls
#'   `update(fit, data = data)`.
#' - `cuiModelVars(fit)`: the variables the predictions need in `newdata`.
#'
#' Methods are provided for `lm` (including `glm`), `merMod` (lme4),
#' `glmmTMB` (conditional model without zero inflation for draws) and
#' `brmsfit`.  The defaults use `predict(type = "response")`,
#' `update()` and the model terms, and have no parameter draws.
#'
#' @param fit A fitted model.
#' @param newdata Data frame of predictors.
#' @param draw `NULL` for the estimates, or one element of
#'   `cuiModelDraws()`.
#' @param n Number of parameter draws.
#' @param data Data to refit to.
#' @param ... For methods.
#' @return `cuiModelPredict()`: numeric vector, one value per row of
#'   `newdata`; `cuiModelDraws()`: a list of `n` draws;
#'   `cuiModelRefit()`: the refit model; `cuiModelVars()`: character vector.
#' @family Clinical utility
#' @export
#' @examples
#' d <- data.frame(dose = rep(c(0, 10, 20, 40), each = 25))
#' set.seed(1)
#' d$eff <- rbinom(nrow(d), 1, plogis(-2 + 0.08 * d$dose))
#' fit <- glm(eff ~ dose, binomial, data = d)
#' cuiModelVars(fit)
#' cuiModelPredict(fit, data.frame(dose = c(0, 40)))
#' dr <- cuiModelDraws(fit, 2)
#' cuiModelPredict(fit, data.frame(dose = c(0, 40)), draw = dr[[1]])
cuiModelPredict <- function(fit, newdata, draw = NULL, ...) {
  UseMethod("cuiModelPredict")
}

#' @rdname cuiModelPredict
#' @export
cuiModelDraws <- function(fit, n, ...) {
  UseMethod("cuiModelDraws")
}

#' @rdname cuiModelPredict
#' @export
cuiModelRefit <- function(fit, data, ...) {
  UseMethod("cuiModelRefit")
}

#' @rdname cuiModelPredict
#' @export
cuiModelVars <- function(fit, ...) {
  UseMethod("cuiModelVars")
}

#' @rdname cuiModelPredict
#' @export
cuiModelPredict.default <- function(fit, newdata, draw = NULL, ...) {
  if (!is.null(draw)) .cuiNoDraws(fit)
  as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
}

#' @rdname cuiModelPredict
#' @export
cuiModelDraws.default <- function(fit, n, ...) {
  .cuiNoDraws(fit)
}

#' @rdname cuiModelPredict
#' @export
cuiModelRefit.default <- function(fit, data, ...) {
  .cuiUpdateIn(fit, "data", data)
}

#' Refit a model to new data, evaluating its call where its formula was made
#'
#' So fits built inside functions (with local family, weights or formula
#' variables) can be refit.
#' @noRd
.cuiUpdateIn <- function(fit, arg, data, ...) {
  .call <- stats::update(fit, ..., evaluate = FALSE)
  .call[[arg]] <- quote(.cuiRefitData)
  .f <- tryCatch(stats::formula(fit), error = function(e) NULL)
  .parent <- if (is.null(.f) || is.null(environment(.f))) globalenv() else environment(.f)
  .env <- new.env(parent = .parent)
  assign(".cuiRefitData", data, envir = .env)
  eval(.call, .env)
}

#' @rdname cuiModelPredict
#' @export
cuiModelVars.default <- function(fit, ...) {
  # an offset given as an argument is not part of the terms
  unique(c(all.vars(stats::delete.response(stats::terms(fit))),
           if (is.call(fit$call)) all.vars(fit$call$offset)))
}

.cuiNoDraws <- function(fit) {
  stop(sprintf(paste0("no parameter draws for a model of class '%s'; write cuiModelDraws() ",
                      "and cuiModelPredict() methods for it"), class(fit)[1]), call. = FALSE)
}

#' Multivariate normal draws of fixed effects
#' @noRd
.cuiCoefDraws <- function(beta, vc, n) {
  .ok <- !is.na(beta)
  beta <- beta[.ok]
  vc <- as.matrix(vc)[names(beta), names(beta), drop = FALSE]
  .r <- chol(vc)
  lapply(seq_len(n), function(i) beta + drop(stats::rnorm(length(beta)) %*% .r))
}

#' Factor levels of a model frame, for predicting on data missing levels
#' @noRd
.cuiFrameLevels <- function(frame) {
  .f <- Filter(function(x) is.factor(x) || is.character(x), frame)
  lapply(.f, function(x) levels(factor(x)))
}

#' Inverse link of a fixed-effects design times a coefficient draw
#' @noRd
.cuiLinkPredict <- function(tt, newdata, beta, linkinv, xlev = NULL, contrasts = NULL) {
  tt <- stats::delete.response(tt)
  if (!is.null(attr(tt, "offset"))) {
    stop("models with an offset are not supported for parameter draws", call. = FALSE)
  }
  # only factors of these terms (not, say, random-effect grouping variables)
  if (!is.null(xlev)) xlev <- xlev[intersect(names(xlev), all.vars(tt))]
  .mf <- stats::model.frame(tt, newdata, xlev = xlev, na.action = stats::na.pass)
  .x <- stats::model.matrix(tt, .mf, contrasts.arg = contrasts)
  .miss <- setdiff(names(beta), colnames(.x))
  if (length(.miss) > 0) {
    stop(sprintf("design matrix lacks coefficient(s) %s", paste(.miss, collapse = ", ")),
         call. = FALSE)
  }
  as.numeric(linkinv(drop(.x[, names(beta), drop = FALSE] %*% beta)))
}

#' @rdname cuiModelPredict
#' @export
cuiModelPredict.lm <- function(fit, newdata, draw = NULL, ...) {
  if (is.null(draw)) {
    return(as.numeric(stats::predict(fit, newdata = newdata, type = "response")))
  }
  if (is.call(fit$call) && !is.null(fit$call$offset)) {
    stop("models with an offset are not supported for parameter draws", call. = FALSE)
  }
  .cuiLinkPredict(stats::terms(fit), newdata, draw, stats::family(fit)$linkinv,
                  xlev = fit$xlevels, contrasts = fit$contrasts)
}

#' @rdname cuiModelPredict
#' @export
cuiModelDraws.lm <- function(fit, n, ...) {
  .cuiCoefDraws(stats::coef(fit), stats::vcov(fit), n)
}

#' @rdname cuiModelPredict
#' @export
cuiModelPredict.merMod <- function(fit, newdata, draw = NULL, ...) {
  if (is.null(draw)) {
    return(as.numeric(stats::predict(fit, newdata = newdata, type = "response", re.form = NA)))
  }
  .cuiLinkPredict(stats::terms(fit, fixed.only = TRUE), newdata, draw,
                  stats::family(fit)$linkinv, xlev = .cuiFrameLevels(fit@frame),
                  contrasts = attr(lme4::getME(fit, "X"), "contrasts"))
}

#' @rdname cuiModelPredict
#' @export
cuiModelDraws.merMod <- function(fit, n, ...) {
  if (!requireNamespace("lme4", quietly = TRUE)) stop("needs the 'lme4' package", call. = FALSE)
  .cuiCoefDraws(lme4::fixef(fit), stats::vcov(fit), n)
}

#' @rdname cuiModelPredict
#' @export
cuiModelVars.merMod <- function(fit, ...) {
  all.vars(stats::delete.response(stats::terms(fit, fixed.only = TRUE)))
}

#' @rdname cuiModelPredict
#' @export
cuiModelPredict.glmmTMB <- function(fit, newdata, draw = NULL, ...) {
  if (is.null(draw)) {
    return(as.numeric(stats::predict(fit, newdata = newdata, type = "response", re.form = NA)))
  }
  # these terms keep predvars, so poly() and scale() use the training basis
  .cuiLinkPredict(fit$modelInfo$terms$cond$fixed, newdata, draw, stats::family(fit)$linkinv,
                  xlev = .cuiFrameLevels(fit$frame),
                  contrasts = attr(glmmTMB::getME(fit, "X"), "contrasts"))
}

#' @rdname cuiModelPredict
#' @export
cuiModelDraws.glmmTMB <- function(fit, n, ...) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) stop("needs the 'glmmTMB' package", call. = FALSE)
  .zi <- fit$modelInfo$allForm$ziformula
  if (!is.null(.zi) && !identical(deparse(.zi), "~0")) {
    stop("parameter draws for glmmTMB need a model without zero inflation", call. = FALSE)
  }
  .cuiCoefDraws(glmmTMB::fixef(fit)$cond, stats::vcov(fit)$cond, n)
}

#' @rdname cuiModelPredict
#' @export
cuiModelVars.glmmTMB <- function(fit, ...) {
  all.vars(stats::delete.response(stats::terms(stats::formula(fit, fixed.only = TRUE))))
}

#' @rdname cuiModelPredict
#' @export
cuiModelPredict.brmsfit <- function(fit, newdata, draw = NULL, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("needs the 'brms' package", call. = FALSE)
  if (is.null(draw)) {
    return(colMeans(brms::posterior_epred(fit, newdata = newdata, re_formula = NA)))
  }
  as.numeric(brms::posterior_epred(fit, newdata = newdata, re_formula = NA, draw_ids = draw))
}

#' @rdname cuiModelPredict
#' @export
cuiModelDraws.brmsfit <- function(fit, n, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) stop("needs the 'brms' package", call. = FALSE)
  .nd <- brms::ndraws(fit)
  as.list(sample.int(.nd, n, replace = n > .nd))
}

#' @rdname cuiModelPredict
#' @export
cuiModelRefit.brmsfit <- function(fit, data, ...) {
  stats::update(fit, newdata = data, refresh = 0)
}

#' @rdname cuiModelPredict
#' @export
cuiModelVars.brmsfit <- function(fit, ...) {
  all.vars(stats::delete.response(stats::terms(.cuiNoBars(stats::formula(fit)$formula))))
}

#' Drop random-effect terms such as (1 | id) from a formula
#' @noRd
.cuiNoBars <- function(f) {
  .bar <- function(e) {
    if (is.call(e) && identical(e[[1]], as.name("("))) return(.bar(e[[2]]))
    is.call(e) && (identical(e[[1]], as.name("|")) || identical(e[[1]], as.name("||")))
  }
  .drop <- function(e) {
    if (!is.call(e) || !identical(e[[1]], as.name("+"))) return(if (.bar(e)) NULL else e)
    .l <- .drop(e[[2]])
    .r <- if (length(e) > 2) .drop(e[[3]]) else NULL
    if (is.null(.l)) return(.r)
    if (is.null(.r)) return(.l)
    call("+", .l, .r)
  }
  .rhs <- .drop(f[[length(f)]])
  if (is.null(.rhs)) .rhs <- 1
  f[[length(f)]] <- .rhs
  f
}
