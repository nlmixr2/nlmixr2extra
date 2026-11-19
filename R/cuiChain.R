#' Exposure from a PK model for exposure-response endpoints
#'
#' Simulates a PK model for each regimen and virtual subject of [cui()] and
#' passes one value per subject (for example Cavg) to the exposure-response
#' models that use it as a covariate.  Exposure-response predictions are
#' made per subject before averaging over the population, and every
#' downstream endpoint sees the same PK simulation.
#'
#' Exposures are subject-level scalars; a time-varying exposure needs a
#' joint PK/PD model instead.
#'
#' @param fit The PK model: an `nlmixr2` fit or an `rxUi`.
#' @param ... Named summaries, one per exposure covariate, each with an
#'   explicit `time` or `window`, for example
#'   `CAVG = cuiMean(cp, window = c(0, 168))`.
#' @return A `cuiExposure` object for the `exposure` argument of [cui()].
#' @family Clinical utility
#' @export
#' @examples
#' \dontrun{
#' cuiExposure(fitPK, CAVG = cuiMean(cp, window = c(0, 168)))
#' }
cuiExposure <- function(fit, ...) {
  if (!inherits(fit, "nlmixr2FitCore") && !inherits(fit, "rxUi")) {
    stop("'fit' must be an nlmixr2 fit or an rxUi", call. = FALSE)
  }
  .sums <- list(...)
  if (length(.sums) == 0 || is.null(names(.sums)) || any(names(.sums) == "") ||
    anyDuplicated(names(.sums)) != 0L) {
    stop("give exposures as uniquely named summaries, for example CAVG = cuiMean(cp, window = c(0, 24))",
         call. = FALSE)
  }
  for (.n in names(.sums)) {
    .s <- .sums[[.n]]
    if (!inherits(.s, "cuiSummary")) {
      stop(sprintf("exposure '%s' must be a cuiMean(), cuiAuc(), cuiProb() or cuiThreshold() summary", .n),
           call. = FALSE)
    }
    if (is.null(.s$time) && is.null(.s$window)) {
      stop(sprintf("exposure '%s' needs an explicit 'time' or 'window'", .n), call. = FALSE)
    }
  }
  structure(list(fit = fit, summaries = .sums), class = "cuiExposure")
}

#' @export
print.cuiExposure <- function(x, ...) {
  cat("cuiExposure:\n")
  for (.n in names(x$summaries)) cat("  ", .n, " = ", format(x$summaries[[.n]]), "\n", sep = "")
  invisible(x)
}

#' Bootstrap parameter uncertainty for a simulated clinical utility index
#'
#' Uses bootstrap refits, for example from
#' `bootstrapFit(fit, returnType = "fitList")`, as the parameter draws of
#' [cui()]: draw `b` simulates with every parameter (thetas, omega) of
#' bootstrap fit `b`.  Each ensemble is matched to the endpoint or exposure
#' fit it was bootstrapped from.  When several ensembles were refit on the
#' same resamples (same original subjects), draws pair the fits of the same
#' resample; otherwise their draws are independent.
#'
#' @param ... One or more bootstrap ensembles (lists of fits; failed refits
#'   are `NA`).
#' @param minSuccess Minimum fraction of each ensemble's refits that must
#'   have succeeded.
#' @return A `cuiBootstrap` object for the `uncertainty` argument of
#'   [cui()].
#' @family Clinical utility
#' @export
#' @examples
#' \dontrun{
#' boot <- bootstrapFit(fit, nboot = 200, returnType = "fitList")
#' res <- cui(endpoints, regimens, population, uncertainty = cuiBootstrap(boot))
#' }
cuiBootstrap <- function(..., minSuccess = 0.9) {
  checkmate::assertNumber(minSuccess, lower = 0, upper = 1)
  .ens <- list(...)
  if (length(.ens) == 0) stop("give at least one bootstrap ensemble", call. = FALSE)
  .ens <- lapply(seq_along(.ens), function(.i) {
    .e <- .ens[[.i]]
    if (inherits(.e, "nlmixr2FitCore")) {
      stop("pass the bootstrap fits, from bootstrapFit(..., returnType = \"fitList\")", call. = FALSE)
    }
    if (!is.list(.e) || length(.e) == 0) {
      stop("each bootstrap ensemble must be a non-empty list of fits", call. = FALSE)
    }
    .items <- lapply(.e, .cuiBootItem)
    .ok <- vapply(.items, `[[`, logical(1), "ok")
    if (mean(.ok) < minSuccess) {
      stop(sprintf("only %d of %d refits of bootstrap ensemble %d succeeded (minSuccess = %g)",
                   sum(.ok), length(.ok), .i, minSuccess), call. = FALSE)
    }
    .md5 <- unique(stats::na.omit(vapply(.items[.ok], `[[`, character(1), "parentMd5")))
    if (length(.md5) > 1L) {
      stop(sprintf("bootstrap ensemble %d mixes refits of different fits", .i), call. = FALSE)
    }
    .items <- .items[.ok]
    .idx <- vapply(.items, `[[`, numeric(1), "index")
    if (!anyNA(.idx)) .items <- .items[order(.idx)]
    list(parentMd5 = if (length(.md5) == 1L) .md5 else NA_character_, items = .items,
         nFit = length(.ok), nFailed = sum(!.ok))
  })
  structure(list(ensembles = .ens, minSuccess = minSuccess), class = "cuiBootstrap")
}

#' @export
print.cuiBootstrap <- function(x, ...) {
  cat("cuiBootstrap:", length(x$ensembles), "ensemble(s)\n")
  for (.e in x$ensembles) {
    cat(sprintf("  %d of %d refits succeeded\n", length(.e$items), .e$nFit))
  }
  invisible(x)
}

#' Parameters and resample identity of one bootstrap refit
#' @noRd
.cuiBootItem <- function(f) {
  if (!inherits(f, "nlmixr2FitCore")) {
    return(list(ok = FALSE, parentMd5 = NA_character_, index = NA_real_))
  }
  .info <- f$bootstrapResample
  list(ok = TRUE, theta = f$theta, omega = f$omega,
       parentMd5 = if (is.null(.info$parentMd5)) NA_character_ else as.character(.info$parentMd5),
       index = if (is.null(.info$index)) NA_real_ else as.numeric(.info$index),
       resample = if (is.null(.info$ids)) NA_character_ else digest::digest(.info$ids))
}

#' Pair the refits of several bootstrap ensembles into draws
#'
#' Ensembles refit on the same resamples are paired by resample; otherwise
#' the b-th refit of each ensemble forms draw b (independent draws).
#'
#' @param ens List of ensembles from `cuiBootstrap()`.
#' @return List with `draws` (matrix, one row per draw, one column per
#'   ensemble, of item indices) and `paired` (logical).
#' @noRd
.cuiBootPair <- function(ens) {
  .res <- lapply(ens, function(e) vapply(e$items, `[[`, character(1), "resample"))
  if (length(ens) > 1L && !anyNA(unlist(.res))) {
    .common <- Reduce(intersect, .res)
    if (length(.common) > 0) {
      .m <- vapply(.res, function(r) match(.common, r), integer(length(.common)))
      return(list(draws = matrix(.m, nrow = length(.common)), paired = TRUE))
    }
  }
  .n <- min(vapply(ens, function(e) length(e$items), integer(1)))
  list(draws = matrix(rep(seq_len(.n), times = length(ens)), nrow = .n), paired = FALSE)
}
