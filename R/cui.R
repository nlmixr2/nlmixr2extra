#' Capture a name given as a symbol or a string
#'
#' @param expr The result of `substitute()` on the argument.
#' @param isMissing Whether the argument was missing.
#' @param se The standard-evaluation twin (a string or `NULL`).
#' @param what Argument name, for error messages.
#' @return A single string.
#' @noRd
.cuiCaptureName <- function(expr, isMissing, se, what) {
  if (!is.null(se)) {
    if (!isMissing) {
      stop(sprintf("supply only one of '%s' and '.%s'", what, what), call. = FALSE)
    }
    checkmate::assertString(se, min.chars = 1, .var.name = paste0(".", what))
    return(se)
  }
  if (isMissing) {
    stop(sprintf("'%s' is required", what), call. = FALSE)
  }
  if (is.symbol(expr)) {
    return(as.character(expr))
  }
  if (checkmate::testString(expr, min.chars = 1)) {
    return(expr)
  }
  stop(sprintf("'%s' must be a name or a single string", what), call. = FALSE)
}

#' Define a clinical utility index endpoint
#'
#' @param name Endpoint name, as a bare name (`resp`) or a string
#'   (`"resp"`).
#' @param fit The fit the endpoint is simulated from (used by the
#'   simulation-based CUI; `NULL` when the endpoint values are supplied
#'   directly to [cuiFromValues()]).
#' @param summary How the endpoint is summarized from a simulation (used
#'   by the simulation-based CUI).
#' @param utility A utility from [cuiIdentity()], [cuiLinear()] or
#'   [cuiStep()] mapping the endpoint to `[0, 1]`.
#' @param direction `"higher"` when larger utilities are better, `"lower"`
#'   to use `1 - utility` (for example a probability of an adverse event).
#' @param ... Not used.
#' @param .name Endpoint name as a string, for programmatic use instead of
#'   `name`.
#' @return A `cuiEndpoint` object.
#' @family Clinical utility
#' @export
#' @examples
#' cuiEndpoint(tox, direction = "lower")
#' cuiEndpoint("eff")
#' nm <- "tol"
#' cuiEndpoint(.name = nm, utility = cuiLinear(0, 50))
cuiEndpoint <- function(name, fit = NULL, summary = NULL, utility = cuiIdentity(),
                        direction = c("higher", "lower"), ..., .name = NULL) {
  .name <- .cuiCaptureName(substitute(name), missing(name), .name, "name")
  if (!inherits(utility, "cuiUtility")) {
    stop("'utility' must come from cuiIdentity(), cuiLinear() or cuiStep()", call. = FALSE)
  }
  direction <- match.arg(direction)
  structure(list(name = .name, fit = fit, summary = summary,
                 utility = utility, direction = direction),
            class = "cuiEndpoint")
}

#' @export
print.cuiEndpoint <- function(x, ...) {
  cat(sprintf("cuiEndpoint '%s': %s utility, %s is better\n",
              x$name, format(x$utility), x$direction))
  invisible(x)
}

#' Candidate regimens for a clinical utility index
#'
#' @param dose Dose of each regimen.
#' @param ii Dosing interval (model time units), or `NA`.
#' @param duration Treatment duration (model time units), or `NA`.
#' @param units Dose units, or `NA`.
#' @param label Unique regimen labels; by default built from `dose` and
#'   `ii`.
#' @return A `cuiRegimens` data frame with columns `label`, `dose`, `ii`,
#'   `duration` and `units`.
#' @family Clinical utility
#' @export
#' @examples
#' cuiRegimens(dose = c(150, 200), ii = c(24, 12), units = "mg")
cuiRegimens <- function(dose, ii = NA_real_, duration = NA_real_,
                        units = NA_character_, label = NULL) {
  checkmate::assertNumeric(dose, min.len = 1)
  .ret <- data.frame(dose = dose, ii = as.numeric(ii), duration = as.numeric(duration),
                     units = as.character(units), stringsAsFactors = FALSE)
  if (is.null(label)) {
    label <- .cuiNumLabel(.ret$dose)
    .hasIi <- !is.na(.ret$ii)
    label[.hasIi] <- paste0(label[.hasIi], " q", .cuiNumLabel(.ret$ii[.hasIi]))
  }
  label <- as.character(label)
  if (length(label) != nrow(.ret)) {
    stop("'label' must have one entry per regimen", call. = FALSE)
  }
  if (anyDuplicated(label) != 0L || anyNA(label)) {
    stop("regimen labels must be unique and not missing", call. = FALSE)
  }
  .ret <- cbind(data.frame(label = label, stringsAsFactors = FALSE), .ret)
  class(.ret) <- c("cuiRegimens", "data.frame")
  .ret
}

#' Label numbers one at a time, so labels match however a vector is formatted
#'
#' @param x Numeric vector.
#' @return Character labels (`10`, `2.5`, `1e+05`).
#' @noRd
.cuiNumLabel <- function(x) {
  vapply(x, function(.x) format(.x, digits = 15, trim = TRUE), character(1), USE.NAMES = FALSE)
}

#' Normalize the eligibility specification to a full grid
#'
#' @param eligible `NULL` (all eligible), a character vector of eligible
#'   regimen labels, a logical vector named by regimen, or a data frame
#'   listing eligible `(regimen, assessment)` pairs.
#' @param regimens Regimen labels.
#' @param assessments Assessment labels.
#' @return Data frame `(assessment, regimen, eligible)`.
#' @noRd
.cuiEligible <- function(eligible, regimens, assessments) {
  .grid <- expand.grid(regimen = regimens, assessment = assessments,
                       stringsAsFactors = FALSE)[, c("assessment", "regimen")]
  if (is.null(eligible)) {
    .grid$eligible <- TRUE
  } else if (is.data.frame(eligible)) {
    if (!all(c("regimen", "assessment") %in% names(eligible))) {
      stop("an 'eligible' data frame needs 'regimen' and 'assessment' columns", call. = FALSE)
    }
    if ("eligible" %in% names(eligible)) {
      if (!checkmate::testLogical(eligible[["eligible"]], any.missing = FALSE)) {
        stop("the 'eligible' column of an 'eligible' data frame must be TRUE/FALSE", call. = FALSE)
      }
      eligible <- eligible[eligible[["eligible"]], , drop = FALSE]
    }
    .reg <- as.character(eligible[["regimen"]])
    .ass <- as.character(eligible[["assessment"]])
    .bad <- setdiff(.reg, regimens)
    if (length(.bad) > 0) {
      stop(sprintf("unknown regimen(s) in 'eligible': %s", paste(.bad, collapse = ", ")), call. = FALSE)
    }
    .bad <- setdiff(.ass, assessments)
    if (length(.bad) > 0) {
      stop(sprintf("unknown assessment(s) in 'eligible': %s", paste(.bad, collapse = ", ")),
           call. = FALSE)
    }
    .key <- paste(.ass, .reg, sep = "\r")
    .grid$eligible <- paste(.grid$assessment, .grid$regimen, sep = "\r") %in% .key
  } else if (is.logical(eligible)) {
    if (is.null(names(eligible)) || !setequal(names(eligible), regimens) || anyNA(eligible)) {
      stop("a logical 'eligible' must be named with every regimen", call. = FALSE)
    }
    .grid$eligible <- unname(eligible[.grid$regimen])
  } else if (is.character(eligible)) {
    .bad <- setdiff(eligible, regimens)
    if (length(.bad) > 0) {
      stop(sprintf("unknown regimen(s) in 'eligible': %s", paste(.bad, collapse = ", ")), call. = FALSE)
    }
    .grid$eligible <- .grid$regimen %in% eligible
  } else {
    stop("'eligible' must be NULL, a character or named logical vector, or a data frame", call. = FALSE)
  }
  .none <- vapply(assessments, function(a) !any(.grid$eligible[.grid$assessment == a]), logical(1))
  if (any(.none)) {
    stop(sprintf("no eligible regimen at assessment(s): %s",
                 paste(assessments[.none], collapse = ", ")), call. = FALSE)
  }
  .grid
}

#' Clinical utility index from endpoint values
#'
#' Builds a clinical utility index (CUI) from endpoint summaries that are
#' already computed, for example observed response rates per dose or
#' model predictions.  For each assessment, regimen and draw, every
#' endpoint is mapped to a utility and the utilities are combined with the
#' normalized weights.
#'
#' `draw = 0` holds the point estimate.  Draws `1, 2, ...` (from a
#' bootstrap or parameter uncertainty) give percentile intervals and the
#' selection frequency: the fraction of draws in which a regimen has the
#' highest CUI among the eligible regimens, with ties within `tieTol`
#' sharing the draw equally.
#'
#' @param values Data frame with columns `endpoint`, `regimen` and
#'   `value`, and optionally `assessment` (default `"all"`) and `draw`
#'   (default `0`).  Every assessment needs exactly one value per endpoint,
#'   regimen and draw, and must include draw 0.
#' @param endpoints List of [cuiEndpoint()] objects; by default every
#'   endpoint in `values` with [cuiIdentity()] and `direction = "higher"`.
#' @param regimens A [cuiRegimens()] table whose labels match
#'   `values$regimen`; by default built from the regimens in `values`.
#' @param weights Named nonnegative endpoint weights (any scale; they are
#'   normalized to sum to one).  `NULL` weights endpoints equally.
#' @param combine `"additive"` (weighted mean of utilities) or
#'   `"geometric"` (weighted geometric mean over positive-weight endpoints;
#'   a zero utility gives a zero CUI).
#' @param eligible Regimens the recommendation may choose from: `NULL`
#'   (all), a character vector of regimen labels, a logical vector named by
#'   regimen, or a data frame listing eligible `regimen` and `assessment`
#'   pairs (with an optional logical `eligible` column, as in the
#'   `eligible` element of the result).
#' @param level Confidence level of the percentile intervals.
#' @param tieTol CUI differences up to `tieTol` count as ties.
#' @return An `nlmixr2Cui` object.
#' @family Clinical utility
#' @export
#' @examples
#' # CUI-MET example 1 (Zhang et al. 2025): modeled probabilities per dose
#' p <- data.frame(dose = 1:5,
#'                 tox = c(0.038, 0.058, 0.106, 0.232, 0.533),
#'                 eff = c(0.022, 0.119, 0.343, 0.568, 0.681),
#'                 tol = c(0.147, 0.206, 0.280, 0.368, 0.466))
#' values <- data.frame(regimen = rep(p$dose, 3),
#'                      endpoint = rep(c("tox", "eff", "tol"), each = 5),
#'                      value = c(p$tox, p$eff, p$tol))
#' res <- cuiFromValues(values,
#'   endpoints = list(cuiEndpoint(tox, direction = "lower"),
#'                    cuiEndpoint(eff), cuiEndpoint(tol)),
#'   weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
#' res
#' update(res, weights = c(tox = 0.3, eff = 0.2, tol = 0.5))
cuiFromValues <- function(values, endpoints = NULL, regimens = NULL, weights = NULL,
                          combine = c("additive", "geometric"), eligible = NULL,
                          level = 0.95, tieTol = sqrt(.Machine$double.eps)) {
  .cuiNew(values = values, endpoints = endpoints, regimens = regimens,
          weights = weights, combine = match.arg(combine), eligible = eligible,
          level = level, tieTol = tieTol, draws = NULL,
          settings = list(source = "values"))
}

#' Validate inputs and build an `nlmixr2Cui` object
#'
#' @param draws Optional draw manifest `(draw, status, source, message)`;
#'   failed draws must have no values.
#' @param settings List of provenance entries.
#' @param batch Optional draw-0 endpoint values per subject batch
#'   `(assessment, regimen, endpoint, batch, value)` for the Monte Carlo
#'   standard error.
#' @inheritParams cuiFromValues
#' @return An `nlmixr2Cui` object.
#' @noRd
.cuiNew <- function(values, endpoints, regimens, weights, combine, eligible,
                    level, tieTol, draws, settings, batch = NULL) {
  checkmate::assertDataFrame(values, min.rows = 1)
  if (!all(c("endpoint", "regimen", "value") %in% names(values))) {
    stop("'values' needs columns 'endpoint', 'regimen' and 'value'", call. = FALSE)
  }
  checkmate::assertNumber(level, lower = 0, upper = 1)
  checkmate::assertNumber(tieTol, lower = 0, finite = TRUE)
  # [[ ]] rather than $, which partial-matches (a 'draws' column is not 'draw')
  .regIn <- values[["regimen"]]
  .assIn <- values[["assessment"]]
  .drawIn <- values[["draw"]]
  if (!is.numeric(values[["value"]])) {
    stop("'values$value' must be numeric", call. = FALSE)
  }
  if (!is.null(.drawIn) && !checkmate::testIntegerish(.drawIn)) {
    stop("'values$draw' must be whole numbers", call. = FALSE)
  }
  .v <- data.frame(
    assessment = if (is.null(.assIn)) "all" else as.character(.assIn),
    regimen = if (is.numeric(.regIn)) .cuiNumLabel(.regIn) else as.character(.regIn),
    draw = if (is.null(.drawIn)) 0L else as.integer(.drawIn),
    endpoint = as.character(values[["endpoint"]]),
    value = values[["value"]],
    stringsAsFactors = FALSE
  )
  if (anyNA(.v)) {
    stop("'values' must not contain missing values", call. = FALSE)
  }
  if (any(.v$draw < 0)) {
    stop("'draw' must be 0 (point estimate) or a positive integer", call. = FALSE)
  }
  # endpoints
  if (is.null(endpoints)) {
    endpoints <- lapply(unique(.v$endpoint), function(e) cuiEndpoint(.name = e))
  }
  if (inherits(endpoints, "cuiEndpoint")) endpoints <- list(endpoints)
  if (!all(vapply(endpoints, inherits, logical(1), what = "cuiEndpoint"))) {
    stop("'endpoints' must be a list of cuiEndpoint() objects", call. = FALSE)
  }
  .epNames <- vapply(endpoints, `[[`, character(1), "name")
  if (anyDuplicated(.epNames) != 0L) {
    stop("endpoint names must be unique", call. = FALSE)
  }
  if (!setequal(.epNames, unique(.v$endpoint))) {
    stop(sprintf("endpoints in 'values' (%s) do not match 'endpoints' (%s)",
                 paste(unique(.v$endpoint), collapse = ", "),
                 paste(.epNames, collapse = ", ")), call. = FALSE)
  }
  # regimens
  if (is.null(regimens)) {
    if (is.numeric(.regIn)) {
      .d <- sort(unique(.regIn))
      regimens <- cuiRegimens(dose = .d, label = .cuiNumLabel(.d))
    } else {
      .l <- unique(.v$regimen)
      regimens <- cuiRegimens(dose = rep(NA_real_, length(.l)), label = .l)
    }
  }
  if (!inherits(regimens, "cuiRegimens")) {
    stop("'regimens' must come from cuiRegimens()", call. = FALSE)
  }
  .bad <- setdiff(unique(.v$regimen), regimens$label)
  if (length(.bad) > 0) {
    stop(sprintf("regimen(s) in 'values' not in 'regimens': %s", paste(.bad, collapse = ", ")),
         call. = FALSE)
  }
  .assess <- unique(.v$assessment)
  if (is.factor(.assIn)) {
    .assess <- intersect(levels(.assIn), .assess)
  }
  # complete grid: one value per endpoint x regimen x draw in every assessment
  if (anyDuplicated(.v[, c("assessment", "regimen", "draw", "endpoint")]) != 0L) {
    stop("'values' has more than one value for an endpoint/regimen/assessment/draw", call. = FALSE)
  }
  .drawSet <- sort(unique(.v$draw))
  if (!(0L %in% .drawSet)) {
    stop("'values' must include draw 0 (the point estimate)", call. = FALSE)
  }
  .nExpect <- length(.epNames) * nrow(regimens) * length(.drawSet)
  for (.a in .assess) {
    .va <- .v[.v$assessment == .a, ]
    if (nrow(.va) != .nExpect || !setequal(unique(.va$draw), .drawSet)) {
      stop(sprintf(paste0("'values' is incomplete at assessment '%s': need every endpoint ",
                          "x regimen x draw (%d values), found %d"),
                   .a, .nExpect, nrow(.va)), call. = FALSE)
    }
  }
  # draw manifest
  if (is.null(draws)) {
    draws <- data.frame(draw = .drawSet, status = "ok", source = "values",
                        message = "", stringsAsFactors = FALSE)
  }
  .failed <- draws$draw[draws$status != "ok"]
  if (any(.failed %in% .drawSet)) {
    stop("failed draws must not have values", call. = FALSE)
  }
  .obj <- structure(list(
    schema = 1L,
    endpoints = lapply(endpoints, function(e) {
      list(name = e$name, utility = e$utility, direction = e$direction,
           summary = if (inherits(e$summary, "cuiSummary")) format(e$summary))
    }),
    regimens = regimens,
    assessments = .assess,
    values = .v,
    draws = draws,
    weights = weights,
    combine = combine,
    eligibleSpec = eligible,
    level = level,
    tieTol = tieTol,
    settings = c(settings, list(version = as.character(utils::packageVersion("nlmixr2extra")))),
    batch = batch
  ), class = "nlmixr2Cui")
  .cuiCompute(.obj)
}

#' Compute utilities, CUI, intervals and selection for an `nlmixr2Cui`
#'
#' Uses only the stored endpoint values, so it never simulates or refits.
#'
#' @param obj An `nlmixr2Cui` object.
#' @return `obj` with `weightsNorm`, `eligible`, `utility`, `cui` and
#'   `summary` filled in.
#' @noRd
.cuiCompute <- function(obj) {
  .epNames <- vapply(obj$endpoints, `[[`, character(1), "name")
  .w <- .cuiNormWeights(obj$weights, .epNames)
  .el <- .cuiEligible(obj$eligibleSpec, obj$regimens$label, obj$assessments)
  .v <- obj$values
  .u <- numeric(nrow(.v))
  for (.ep in obj$endpoints) {
    .sel <- .v$endpoint == .ep$name
    .u[.sel] <- .cuiApplyUtility(.ep$utility, .ep$direction, .v$value[.sel], .ep$name)
  }
  # order so each block of K rows is one (assessment, draw, regimen)
  .ord <- order(match(.v$assessment, obj$assessments), .v$draw,
                match(.v$regimen, obj$regimens$label), match(.v$endpoint, .epNames))
  .v <- .v[.ord, ]
  .v$utility <- .u[.ord]
  rownames(.v) <- NULL
  .k <- length(.epNames)
  .m <- matrix(.v$utility, ncol = .k, byrow = TRUE, dimnames = list(NULL, .epNames))
  .cui <- .v[seq(1, nrow(.v), by = .k), c("assessment", "regimen", "draw")]
  rownames(.cui) <- NULL
  .cui$cui <- .cuiCombine(.m, .w, obj$combine)
  .cui$eligible <- .el$eligible[match(paste(.cui$assessment, .cui$regimen, sep = "\r"),
                                      paste(.el$assessment, .el$regimen, sep = "\r"))]
  # selection among eligible regimens in each uncertainty draw
  .tol <- obj$tieTol
  .nOk <- sum(obj$draws$status == "ok" & obj$draws$draw > 0)
  .d <- .cui[.cui$draw > 0 & .cui$eligible, ]
  .share <- numeric(0)
  if (nrow(.d) > 0) {
    .grp <- paste(.d$assessment, .d$draw, sep = "\r")
    .mx <- ave(.d$cui, .grp, FUN = max)
    .tied <- as.numeric(.d$cui >= .mx - .tol)
    .d$share <- .tied / ave(.tied, .grp, FUN = sum)
  }
  .alpha <- (1 - obj$level) / 2
  .sum <- .el
  .sum <- .sum[order(match(.sum$assessment, obj$assessments),
                     match(.sum$regimen, obj$regimens$label)), ]
  .key <- paste(.sum$assessment, .sum$regimen, sep = "\r")
  .c0 <- .cui[.cui$draw == 0, ]
  .sum$est <- .c0$cui[match(.key, paste(.c0$assessment, .c0$regimen, sep = "\r"))]
  .sum$mean <- NA_real_
  .sum$lower <- NA_real_
  .sum$upper <- NA_real_
  .sum$selFreq <- NA_real_
  if (.nOk > 0) {
    .cd <- .cui[.cui$draw > 0, ]
    .cdKey <- paste(.cd$assessment, .cd$regimen, sep = "\r")
    .sum$mean <- vapply(.key, function(k) mean(.cd$cui[.cdKey == k]), numeric(1), USE.NAMES = FALSE)
    .sum$lower <- vapply(.key, function(k) unname(stats::quantile(.cd$cui[.cdKey == k], .alpha)),
                         numeric(1), USE.NAMES = FALSE)
    .sum$upper <- vapply(.key, function(k) unname(stats::quantile(.cd$cui[.cdKey == k], 1 - .alpha)),
                         numeric(1), USE.NAMES = FALSE)
    .dKey <- paste(.d$assessment, .d$regimen, sep = "\r")
    .sum$selFreq <- vapply(.key, function(k) sum(.d$share[.dKey == k]) / .nOk,
                           numeric(1), USE.NAMES = FALSE)
    .sum$selFreq[!.sum$eligible] <- NA_real_
  }
  # recommended (best point estimate) and most frequently selected regimens
  .sum$best <- FALSE
  .sum$mostSelected <- NA
  for (.a in obj$assessments) {
    .s <- .sum$assessment == .a & .sum$eligible
    .mx <- max(.sum$est[.s])
    .sum$best[.s] <- .sum$est[.s] >= .mx - .tol
    if (.nOk > 0) {
      .sum$mostSelected[.sum$assessment == .a] <- FALSE
      .sum$mostSelected[.s] <- .sum$selFreq[.s] >= max(.sum$selFreq[.s]) - .tol
    }
  }
  rownames(.sum) <- NULL
  obj$weightsNorm <- .w
  obj$eligible <- .el
  obj$utility <- .v
  obj$cui <- .cui
  obj$summary <- .sum
  obj$mcse <- .cuiMcse(obj, .w, .sum)
  obj
}

#' Monte Carlo standard error of the CUI difference of the two best regimens
#'
#' Uses the draw-0 endpoint values per subject batch; batches hold the same
#' subjects for every regimen, so the difference is paired.
#'
#' @return Data frame `(assessment, first, second, diff, se)` or `NULL`.
#' @noRd
.cuiMcse <- function(obj, w, sum) {
  .b <- obj$batch
  if (is.null(.b) || nrow(.b) == 0) return(NULL)
  .epNames <- names(w)
  .u <- numeric(nrow(.b))
  for (.ep in obj$endpoints) {
    .sel <- .b$endpoint == .ep$name
    .u[.sel] <- .cuiApplyUtility(.ep$utility, .ep$direction, .b$value[.sel], .ep$name)
  }
  .ret <- lapply(obj$assessments, function(.a) {
    .s <- sum[sum$assessment == .a & sum$eligible, ]
    if (nrow(.s) < 2) return(NULL)
    .top <- .s$regimen[order(-.s$est)][1:2]
    .cuiB <- lapply(.top, function(.r) {
      .sel <- .b$assessment == .a & .b$regimen == .r
      .m <- tapply(.u[.sel], list(.b$batch[.sel], .b$endpoint[.sel]), identity)
      .cuiCombine(.m[, .epNames, drop = FALSE], w, obj$combine)
    })
    .d <- .cuiB[[1]] - .cuiB[[2]]
    data.frame(assessment = .a, first = .top[1], second = .top[2], diff = mean(.d),
               se = stats::sd(.d) / sqrt(length(.d)), stringsAsFactors = FALSE)
  })
  .ret <- do.call(rbind, .ret)
  if (is.null(.ret)) return(NULL)
  rownames(.ret) <- NULL
  .ret
}

#' Re-weight or re-restrict a clinical utility index
#'
#' Recomputes utilities, the CUI, intervals and selection frequencies from
#' the stored endpoint values; nothing is simulated or refit.
#'
#' @param object An `nlmixr2Cui` object.
#' @param weights,combine,eligible,level,tieTol As in [cuiFromValues()];
#'   omitted arguments keep their current values (`weights = NULL` means
#'   equal weights).
#' @param ... Not used.
#' @return The updated `nlmixr2Cui` object.
#' @family Clinical utility
#' @export
update.nlmixr2Cui <- function(object, weights, combine, eligible, level, tieTol, ...) {
  if (!missing(weights)) object["weights"] <- list(weights)
  if (!missing(combine)) object$combine <- match.arg(combine, c("additive", "geometric"))
  if (!missing(eligible)) object["eligibleSpec"] <- list(eligible)
  if (!missing(level)) {
    checkmate::assertNumber(level, lower = 0, upper = 1)
    object$level <- level
  }
  if (!missing(tieTol)) {
    checkmate::assertNumber(tieTol, lower = 0, finite = TRUE)
    object$tieTol <- tieTol
  }
  .cuiCompute(object)
}

#' @export
print.nlmixr2Cui <- function(x, ..., n = 10L) {
  .nOk <- sum(x$draws$status == "ok" & x$draws$draw > 0)
  .nFail <- sum(x$draws$status != "ok")
  cat("Clinical utility index (", x$combine, ")\n", sep = "")
  if (!is.null(x$settings$description)) cat(x$settings$description, "\n", sep = "")
  .ep <- vapply(x$endpoints, function(e) {
    sprintf("%s%s [%s%s] %.3g", e$name,
            if (is.null(e$summary)) "" else paste0(" = ", e$summary), format(e$utility),
            ifelse(e$direction == "lower", ", lower is better", ""),
            x$weightsNorm[[e$name]])
  }, character(1))
  cat("Endpoints [utility] weight:\n", paste0("  ", .ep, collapse = "\n"), "\n", sep = "")
  if (.nOk > 0 || .nFail > 0) {
    cat(sprintf("Uncertainty draws: %d ok, %d failed; %g%% percentile intervals\n",
                .nOk, .nFail, 100 * x$level))
  }
  for (.a in x$assessments) {
    .s <- x$summary[x$summary$assessment == .a, ]
    if (length(x$assessments) > 1) cat("\nAssessment: ", .a, "\n", sep = "")
    .show <- .s
    if (nrow(.show) > n) {
      .show <- .show[order(-.show$est), ][seq_len(n), ]
      cat(sprintf("Top %d of %d regimens:\n", n, nrow(.s)))
    }
    .tab <- data.frame(regimen = .show$regimen, cui = round(.show$est, 3),
                       stringsAsFactors = FALSE)
    if (.nOk > 0) {
      .tab$lower <- round(.show$lower, 3)
      .tab$upper <- round(.show$upper, 3)
      .tab$selFreq <- round(.show$selFreq, 3)
    }
    if (!all(.show$eligible)) .tab$eligible <- .show$eligible
    print(.tab, row.names = FALSE)
    cat("Recommended (highest CUI): ", paste(.s$regimen[.s$best], collapse = ", "), "\n", sep = "")
    .mc <- x$mcse[!is.null(x$mcse) & x$mcse$assessment == .a, , drop = FALSE]
    if (!is.null(.mc) && nrow(.mc) == 1) {
      cat(sprintf("Monte Carlo SE of CUI(%s) - CUI(%s): %.2g\n", .mc$first, .mc$second, .mc$se))
    }
    if (.nOk > 0) {
      cat("Most frequently selected: ",
          paste(.s$regimen[!is.na(.s$mostSelected) & .s$mostSelected], collapse = ", "),
          "\n", sep = "")
    }
  }
  invisible(x)
}

#' @param x An `nlmixr2Cui` object.
#' @param row.names,optional Not used.
#' @param type `"summary"` (one row per assessment and regimen),
#'   `"cui"` (CUI per draw), `"utility"` (endpoint values and utilities
#'   per draw) or `"observed"` (observed rates, from [cuiDoseResponse()]).
#' @rdname update.nlmixr2Cui
#' @export
as.data.frame.nlmixr2Cui <- function(x, row.names = NULL, optional = FALSE, ...,
                                     type = c("summary", "cui", "utility", "observed")) {
  type <- match.arg(type)
  if (type == "observed" && is.null(x$observed)) {
    stop("no observed data; only cuiDoseResponse() results have them", call. = FALSE)
  }
  x[[type]]
}
