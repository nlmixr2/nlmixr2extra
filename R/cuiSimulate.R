#' Target population for a simulated clinical utility index
#'
#' Takes one row of covariates per subject from a data set (or the data a
#' fit was estimated on) and optionally resamples `n` subjects with
#' replacement, keeping each subject's covariates together.  Every
#' endpoint of [cui()] is simulated in this one population.
#'
#' @param data A data frame or an `nlmixr2` fit (its estimation data).
#' @param n Number of subjects to resample; `NULL` keeps every subject once.
#' @param covs Covariate columns to keep (first record of each subject).
#' @param id Subject id column (matched without regard to case).
#' @param seed Seed for the resampling, or `NULL` for the current stream.
#' @return A `cuiPopulation` data frame with `id` and the covariates.
#' @family Clinical utility
#' @export
#' @examples
#' pop <- cuiPopulation(nlmixr2data::theo_sd, n = 50, covs = "WT", seed = 1)
#' head(pop)
cuiPopulation <- function(data, n = NULL, covs = NULL, id = "ID", seed = NULL) {
  if (inherits(data, "nlmixr2FitCore")) data <- data$origData
  checkmate::assertDataFrame(data, min.rows = 1)
  checkmate::assertCharacter(covs, null.ok = TRUE, any.missing = FALSE)
  .idCol <- names(data)[tolower(names(data)) == tolower(id)]
  if (length(.idCol) != 1L) {
    stop(sprintf("id column '%s' not found in 'data'", id), call. = FALSE)
  }
  .miss <- setdiff(covs, names(data))
  if (length(.miss) > 0) {
    stop(sprintf("covariate(s) not in 'data': %s", paste(.miss, collapse = ", ")), call. = FALSE)
  }
  .first <- data[!duplicated(data[[.idCol]]), c(.idCol, covs), drop = FALSE]
  if (anyNA(.first[, covs, drop = FALSE])) {
    stop("covariates must not be missing on the first record of a subject", call. = FALSE)
  }
  .rows <- seq_len(nrow(.first))
  if (!is.null(n)) {
    checkmate::assertCount(n, positive = TRUE)
    .rows <- .cuiWithSeed(seed, .rows[sample.int(length(.rows), n, replace = TRUE)])
  }
  .ret <- data.frame(id = seq_along(.rows), .first[.rows, covs, drop = FALSE],
                     stringsAsFactors = FALSE)
  rownames(.ret) <- NULL
  class(.ret) <- c("cuiPopulation", "data.frame")
  .ret
}

#' Derive a reproducible seed for one random stream
#'
#' Streams are keyed (for example by draw and model) so a draw gets the same
#' random numbers however the run is split up.
#' @noRd
.cuiStreamSeed <- function(seed, ...) {
  strtoi(substr(digest::digest(list(seed, ...), algo = "xxhash32"), 1, 7), 16L)
}


#' Build the per-model simulation pieces for a group of endpoints
#'
#' @param fit The group's fit (an `nlmixr2` fit or an `rxUi`).
#' @param eps Endpoints sharing this fit.
#' @param exposures Named exposure summaries computed from this fit.
#' @param group Group number (keys its random streams).
#' @return List with the simulation model, resolved summaries and parameters.
#' @noRd
.cuiModelGroup <- function(fit, eps, group, exposures = list()) {
  .isFit <- inherits(fit, "nlmixr2FitCore")
  .ui <- if (.isFit) fit$finalUi else rxode2::assertRxUi(fit)
  .ini <- .ui$iniDf
  .theta <- setNames(.ini$est[!is.na(.ini$ntheta)], .ini$name[!is.na(.ini$ntheta)])
  .struct <- .ini$name[!is.na(.ini$ntheta) & is.na(.ini$err) & !.ini$fix]
  .omega <- .ui$omega
  if (!is.null(.omega) && !is.matrix(.omega)) {
    stop("multi-level random effects (for example IOV) are not supported by cui()", call. = FALSE)
  }
  .mv <- rxode2::rxModelVars(.ui)
  .cols <- unique(c(.mv$lhs, .mv$state, .mv$params, "time"))
  .lines <- list()
  .resolve <- function(s, nm) {
    .s <- .cuiResolveSummary(s, .cols, nm)
    if (!is.null(.s$window)) {
      .state <- paste0("cuiInt", length(.lines) + 1L)
      while (.state %in% .cols) .state <- paste0(.state, "x")
      .lines[[length(.lines) + 1L]] <<- bquote(d / dt(.(as.name(.state))) <- .(.s$expr))
      .s$state <- .state
    }
    .s
  }
  .sums <- lapply(eps, function(.ep) .resolve(.ep$summary, .ep$name))
  names(.sums) <- vapply(eps, `[[`, character(1), "name")
  .expo <- lapply(names(exposures), function(.n) .resolve(exposures[[.n]], .n))
  names(.expo) <- names(exposures)
  .sim <- .ui
  for (.l in .lines) {
    .sim <- suppressMessages(eval(bquote(rxode2::model(.(.sim), .(.l), append = TRUE))))
  }
  list(kind = "nlmixr2", group = group, fit = fit, isFit = .isFit, ui = .sim, theta = .theta,
       struct = .struct,
       lower = setNames(.ini$lower, .ini$name), upper = setNames(.ini$upper, .ini$name),
       omega = .omega, covs = .ui$allCovs, summaries = .sums, exposures = .expo,
       # bootstrapMd5 is only set in memory; a reloaded fit still has md5
       md5 = if (.isFit) {
         if (is.null(fit$bootstrapMd5)) fit$md5 else fit$bootstrapMd5
       })
}

#' Model group for endpoints from a fitted regression model
#'
#' Predictions come from the [cuiModelPredict()] generics; the regimen dose
#' goes into each summary's dose column.
#' @noRd
.cuiFittedGroup <- function(fit, eps, group) {
  .sums <- lapply(eps, `[[`, "summary")
  names(.sums) <- vapply(eps, `[[`, character(1), "name")
  .vars <- cuiModelVars(fit)
  for (.n in names(.sums)) {
    if (!is.null(.sums[[.n]]$doseVar) && !(.sums[[.n]]$doseVar %in% .vars)) {
      stop(sprintf("the model of endpoint '%s' does not use the dose column '%s'", .n,
                   .sums[[.n]]$doseVar), call. = FALSE)
    }
  }
  list(kind = "fitted", group = group, fit = fit, isFit = FALSE, ui = NULL, theta = NULL,
       struct = NULL, omega = NULL, md5 = NULL, vars = .vars,
       # variables every endpoint needs besides its own dose column
       covs = unique(unlist(lapply(.sums, function(s) setdiff(.vars, s$doseVar)))),
       summaries = .sums, exposures = list())
}

#' Cholesky factor of omega, with zero-variance etas held at zero
#' @noRd
.cuiOmegaFactor <- function(omega) {
  .pos <- diag(omega) > 0
  .r <- matrix(0, nrow(omega), ncol(omega))
  if (any(.pos)) .r[.pos, .pos] <- chol(omega[.pos, .pos, drop = FALSE])
  .r
}

#' Draw structural thetas from the fit's covariance, respecting bounds
#' @noRd
.cuiDrawTheta <- function(mg, covChol, seed, draw) {
  .cuiWithSeed(.cuiStreamSeed(seed, "theta", mg$group, draw), {
    .est <- mg$theta[mg$struct]
    for (.try in 1:100) {
      .t <- .est + drop(stats::rnorm(length(.est)) %*% covChol)
      if (all(.t >= mg$lower[mg$struct] & .t <= mg$upper[mg$struct])) break
      if (.try == 100) stop("could not draw thetas within their bounds in 100 tries", call. = FALSE)
    }
    .ret <- mg$theta
    .ret[mg$struct] <- .t
    .ret
  })
}

#' Event records for every regimen x subject
#' @noRd
.cuiEvents <- function(regimens, nSub, obsTimes, cmt, dose = TRUE) {
  .maxT <- max(obsTimes)
  .ev <- lapply(seq_len(nrow(regimens)), function(.r) {
    .ids <- (.r - 1L) * nSub + seq_len(nSub)
    .reg <- regimens[.r, ]
    .obs <- data.frame(id = rep(.ids, each = length(obsTimes)),
                       time = rep(obsTimes, times = nSub),
                       amt = 0, evid = 0L, ii = 0, addl = 0L)
    if (!dose || .reg$dose == 0) return(.obs)
    .ii <- if (is.na(.reg$ii)) 0 else .reg$ii
    .addl <- 0L
    if (.ii > 0) {
      .addl <- if (is.na(.reg$duration)) {
        as.integer(floor(.maxT / .ii))
      } else {
        max(0L, as.integer(ceiling(.reg$duration / .ii)) - 1L)
      }
    }
    rbind(data.frame(id = .ids, time = 0, amt = .reg$dose, evid = 1L, ii = .ii, addl = .addl),
          .obs)
  })
  .ev <- do.call(rbind, .ev)
  if (!is.null(cmt)) .ev$cmt <- cmt
  .ev[order(.ev$id, .ev$time), ]
}

#' Simulate one draw of one model group for all regimens
#'
#' @param extra Optional data frame with one row per regimen x subject
#'   (regimen-major), such as chained exposures.
#' @return Data frame of simulated rows with `.reg` and `.sub` columns.
#' @noRd
.cuiSolveDraw <- function(mg, theta, eta, cov, regimens, obsTimes, cmt, extra = NULL) {
  .nSub <- nrow(eta)
  .nReg <- nrow(regimens)
  .p <- as.data.frame(matrix(rep(theta, each = .nSub * .nReg), nrow = .nSub * .nReg,
                             dimnames = list(NULL, names(theta))))
  if (ncol(eta) > 0) {
    .p <- cbind(.p, eta[rep(seq_len(.nSub), times = .nReg), , drop = FALSE])
  }
  .popCovs <- setdiff(mg$covs, names(extra))
  if (length(.popCovs) > 0) {
    .p <- cbind(.p, cov[rep(seq_len(.nSub), times = .nReg), .popCovs, drop = FALSE])
  }
  if (!is.null(extra)) .p <- cbind(.p, extra)
  rownames(.p) <- NULL
  # rxode2 would quietly fill in (or sample) anything missing here
  .missing <- setdiff(rxode2::rxModelVars(mg$ui)$params, names(.p))
  if (length(.missing) > 0) {
    stop(sprintf("internal error: no value for model parameter(s) %s",
                 paste(.missing, collapse = ", ")), call. = FALSE)
  }
  # a model driven by a chained exposure gets the exposure, not the doses
  .ev <- .cuiEvents(regimens, .nSub, obsTimes, cmt, dose = is.null(extra))
  .s <- suppressMessages(rxode2::rxSolve(mg$ui, params = .p, events = .ev, nStud = 1,
                                         thetaMat = NULL, simVariability = FALSE, sigma = NA,
                                         returnType = "data.frame"))
  # rxode2 drops the id column when there is a single subject
  .id <- if (is.null(.s$id)) rep(1L, nrow(.s)) else as.integer(.s$id)
  .s$.reg <- (.id - 1L) %/% .nSub + 1L
  .s$.sub <- (.id - 1L) %% .nSub + 1L
  .s
}

#' Per-subject values of one summary at one time, ordered by regimen then subject
#' @noRd
.cuiSubjectValues <- function(s, sim, time, nReg, nSub, ep) {
  .at <- function(t) {
    .rows <- sim[abs(sim$time - t) <= 1e-8 * max(1, abs(t)), , drop = FALSE]
    if (nrow(.rows) != nReg * nSub) {
      stop(sprintf("endpoint '%s': expected one record per subject at time %s, found %d of %d",
                   ep, format(t), nrow(.rows), nReg * nSub), call. = FALSE)
    }
    .rows[order(.rows$.reg, .rows$.sub), , drop = FALSE]
  }
  if (!is.null(s$window)) {
    .v <- .at(s$window[2])[[s$state]] - .at(s$window[1])[[s$state]]
    if (s$type == "mean") .v <- .v / (s$window[2] - s$window[1])
    return(.v)
  }
  .rows <- .at(time)
  .v <- eval(s$expr, .rows, globalenv())
  if (length(.v) == 1L) .v <- rep(.v, nrow(.rows))
  if (length(.v) != nrow(.rows)) {
    stop(sprintf("endpoint '%s' must give one value per subject", ep), call. = FALSE)
  }
  if (s$type == "threshold") {
    if (!is.logical(.v) || anyNA(.v)) {
      stop(sprintf("endpoint '%s': cuiThreshold() needs a TRUE/FALSE condition", ep), call. = FALSE)
    }
    return(as.numeric(.v))
  }
  if (!is.numeric(.v) || anyNA(.v)) {
    stop(sprintf("endpoint '%s' must give numeric values without NA", ep), call. = FALSE)
  }
  if (s$type == "prob" && any(.v < 0 | .v > 1)) {
    stop(sprintf("endpoint '%s': cuiProb() values must be in [0, 1]", ep), call. = FALSE)
  }
  .v
}

#' Evaluate one summary at one time for every regimen (and batch)
#'
#' @return List with `value` (per regimen) and `batch` (regimen x batch).
#' @noRd
.cuiEvalSummary <- function(s, sim, time, nReg, nSub, batch, ep) {
  .cuiAggregate(.cuiSubjectValues(s, sim, time, nReg, nSub, ep), nReg, nSub, batch)
}

#' Mean per regimen (and per batch) of per-subject values, regimen-major
#' @noRd
.cuiAggregate <- function(v, nReg, nSub, batch) {
  .m <- matrix(v, nrow = nSub, ncol = nReg)
  .ret <- list(value = colMeans(.m))
  if (!is.null(batch)) {
    .ret$batch <- rowsum(.m, batch, reorder = TRUE) / as.vector(table(batch))
  }
  .ret
}

#' Group endpoint and exposure fits; endpoints sharing a fit share a group
#'
#' @return List of model groups, providers of exposures first, and the
#'   group index of every endpoint.
#' @noRd
.cuiGroups <- function(endpoints, exposure) {
  .fits <- c(lapply(endpoints, `[[`, "fit"), lapply(exposure, `[[`, "fit"))
  .key <- lapply(.fits, function(f) if (inherits(f, "nlmixr2FitCore")) f$env else f)
  .isSim <- vapply(.fits, function(f) inherits(f, "nlmixr2FitCore") || inherits(f, "rxUi"),
                   logical(1))
  if (length(exposure) > 0 && !all(.isSim[length(endpoints) + seq_along(exposure)])) {
    stop("exposures come from nlmixr2 fits or rxUi models", call. = FALSE)
  }
  .grp <- integer(length(.fits))
  for (.i in seq_along(.fits)) {
    .grp[.i] <- .i
    for (.j in seq_len(.i - 1L)) {
      if (identical(.key[[.i]], .key[[.j]])) {
        .grp[.i] <- .grp[.j]
        break
      }
    }
  }
  .nEp <- length(endpoints)
  .epGrp <- .grp[seq_len(.nEp)]
  .exGrp <- .grp[.nEp + seq_along(exposure)]
  .ids <- unique(.grp)
  .groups <- lapply(seq_along(.ids), function(.k) {
    .g <- .ids[.k]
    .exps <- list()
    for (.x in which(.exGrp == .g)) .exps <- c(.exps, exposure[[.x]]$summaries)
    if (anyDuplicated(names(.exps)) != 0L) {
      stop("exposure names must be unique", call. = FALSE)
    }
    if (!.isSim[.g]) return(.cuiFittedGroup(.fits[[.g]], endpoints[.epGrp == .g], .k))
    .cuiModelGroup(.fits[[.g]], endpoints[.epGrp == .g], .k, .exps)
  })
  .exNames <- unlist(lapply(.groups, function(mg) names(mg$exposures)))
  if (anyDuplicated(.exNames) != 0L) stop("exposure names must be unique", call. = FALSE)
  for (mg in .groups) {
    if (length(mg$exposures) > 0 && any(mg$covs %in% .exNames)) {
      stop("an exposure model cannot itself use a chained exposure", call. = FALSE)
    }
  }
  .provider <- vapply(.groups, function(mg) length(mg$exposures) > 0, logical(1))
  .ord <- c(which(.provider), which(!.provider))
  list(groups = .groups[.ord], epGroup = match(match(.epGrp, .ids), .ord), exposures = .exNames)
}

#' Simulate and summarize one draw for every model group
#'
#' @param params Per-group list of `theta` and `omega` for this draw.
#' @param rows Population rows of this draw's subjects.
#' @param etaDraw Draw key of the standardized random effects.
#' @return List of per-endpoint value (and draw-0 batch) tables.
#' @noRd
.cuiRunDraw <- function(d, groups, params, rows, etaDraw, population, regimens, times,
                        obsTimes, cmt, bt, seed) {
  .n <- length(rows)
  .nReg <- nrow(regimens)
  .extra <- NULL
  .out <- list()
  for (.k in seq_along(groups)) {
    .mg <- groups[[.k]]
    # params follow the order of groups; .mg$group only keys random streams
    .p <- params[[.k]]
    if (identical(.mg$kind, "fitted")) {
      .out <- c(.out, .cuiFittedDraw(.mg, .p$draw, rows, population, regimens, .extra, times,
                                     d, bt))
      next
    }
    .eta <- if (is.null(.p$omega) || nrow(.p$omega) == 0) {
      matrix(numeric(0), nrow = .n, ncol = 0)
    } else {
      .z <- .cuiWithSeed(.cuiStreamSeed(seed, "eta", .mg$group, etaDraw),
                         matrix(stats::rnorm(.n * nrow(.p$omega)), nrow = .n))
      .e <- .z %*% .cuiOmegaFactor(.p$omega)
      colnames(.e) <- rownames(.p$omega)
      .e
    }
    .need <- intersect(.mg$covs, names(.extra))
    .sim <- .cuiSolveDraw(.mg, .p$theta, .eta, population[rows, , drop = FALSE], regimens,
                          obsTimes, cmt, if (length(.need) > 0) .extra[, .need, drop = FALSE])
    for (.x in names(.mg$exposures)) {
      .s <- .mg$exposures[[.x]]
      .v <- .cuiSubjectValues(.s, .sim, .s$time, .nReg, .n, .x)
      .extra <- if (is.null(.extra)) setNames(data.frame(.v), .x) else cbind(.extra, setNames(data.frame(.v), .x))
    }
    for (.ep in names(.mg$summaries)) {
      .s <- .mg$summaries[[.ep]]
      for (.a in names(times)) {
        .t <- if (!is.null(.s$time)) .s$time else times[[.a]]
        .r <- .cuiEvalSummary(.s, .sim, .t, .nReg, .n, bt, .ep)
        .out[[length(.out) + 1L]] <- list(
          value = data.frame(assessment = .a, regimen = regimens$label, draw = d,
                             endpoint = .ep, value = .r$value, stringsAsFactors = FALSE),
          batch = if (!is.null(.r$batch)) {
            data.frame(assessment = .a, regimen = rep(regimens$label, each = nrow(.r$batch)),
                       endpoint = .ep, batch = rep(seq_len(nrow(.r$batch)), .nReg),
                       value = as.vector(.r$batch), stringsAsFactors = FALSE)
          })
      }
    }
  }
  .out
}

#' Predict and summarize the endpoints of a fitted-model group for one draw
#' @noRd
.cuiFittedDraw <- function(mg, draw, rows, population, regimens, extra, times, d, bt) {
  .n <- length(rows)
  .nReg <- nrow(regimens)
  .out <- list()
  for (.ep in names(mg$summaries)) {
    .s <- mg$summaries[[.ep]]
    .nd <- data.frame(row.names = seq_len(.n * .nReg))
    for (.v in setdiff(mg$vars, .s$doseVar)) {
      .nd[[.v]] <- if (.v %in% names(extra)) {
        extra[[.v]]
      } else {
        rep(population[[.v]][rows], times = .nReg)
      }
    }
    if (!is.null(.s$doseVar)) .nd[[.s$doseVar]] <- rep(regimens$dose, each = .n)
    .v <- cuiModelPredict(mg$fit, .nd, draw)
    if (!is.numeric(.v) || length(.v) != nrow(.nd) || anyNA(.v)) {
      stop(sprintf("the model of endpoint '%s' must predict one number per subject", .ep),
           call. = FALSE)
    }
    .r <- .cuiAggregate(.v, .nReg, .n, bt)
    for (.a in names(times)) {
      .out[[length(.out) + 1L]] <- list(
        value = data.frame(assessment = .a, regimen = regimens$label, draw = d,
                           endpoint = .ep, value = .r$value, stringsAsFactors = FALSE),
        batch = if (!is.null(.r$batch)) {
          data.frame(assessment = .a, regimen = rep(regimens$label, each = nrow(.r$batch)),
                     endpoint = .ep, batch = rep(seq_len(nrow(.r$batch)), .nReg),
                     value = as.vector(.r$batch), stringsAsFactors = FALSE)
        })
    }
  }
  .out
}

#' Set up the draw cache directory
#'
#' `NA` disables caching, `NULL` derives `nlmixr2CuiCache_<md5>` from `key`.
#' @noRd
.cuiCacheSetup <- function(cacheDir, restart, key) {
  if (length(cacheDir) == 1L && is.na(cacheDir)) return(NULL)
  if (is.null(cacheDir)) cacheDir <- paste0("nlmixr2CuiCache_", digest::digest(key))
  checkmate::assertString(cacheDir, min.chars = 1, .var.name = "cacheDir")
  if (dir.exists(cacheDir) && restart) unlink(cacheDir, recursive = TRUE, force = TRUE)
  if (!dir.exists(cacheDir)) dir.create(cacheDir, recursive = TRUE)
  .keyFile <- file.path(cacheDir, "key.rds")
  .k <- digest::digest(key)
  if (file.exists(.keyFile)) {
    if (!identical(readRDS(.keyFile), .k)) {
      stop(sprintf("cache '%s' belongs to a different simulation; use restart = TRUE or another cacheDir",
                   cacheDir), call. = FALSE)
    }
  } else {
    saveRDS(.k, .keyFile)
  }
  cacheDir
}

.cuiCacheFile <- function(dir, d) file.path(dir, sprintf("draw_%05d.rds", d))

.cuiCacheRead <- function(dir, d, check) {
  if (is.null(dir)) return(NULL)
  .f <- .cuiCacheFile(dir, d)
  if (!file.exists(.f)) return(NULL)
  .ret <- try(readRDS(.f), silent = TRUE)
  # a draw whose parameters changed (for example new bootstrap refits) is redone
  if (inherits(.ret, "try-error") || !identical(.ret$check, check)) NULL else .ret
}

.cuiCacheWrite <- function(dir, d, value) {
  if (is.null(dir)) return(invisible(NULL))
  # write then rename, so an interrupted write never leaves a partial draw
  .tmp <- tempfile("draw", tmpdir = dir, fileext = ".tmp")
  saveRDS(value, .tmp)
  file.rename(.tmp, .cuiCacheFile(dir, d))
  invisible(NULL)
}

#' Simulated clinical utility index from nlmixr2 models
#'
#' Simulates each endpoint's model for every candidate regimen in one
#' target population, summarizes the endpoints per regimen, and combines
#' them into a clinical utility index (CUI) as in [cuiFromValues()].
#'
#' The point estimate (draw 0) uses the fits' estimates.  Parameter
#' uncertainty gives intervals for the expected CUI and selection
#' frequencies:
#' - `uncertainty = "cov"` samples the structural thetas from each fit's
#'   covariance (residual-error and fixed parameters, and omega, stay at
#'   their estimates);
#' - `uncertainty = cuiBootstrap(...)` uses every parameter of bootstrap
#'   refits, pairing models refit on the same resamples.
#'
#' Every regimen and every draw sees the same virtual subjects (covariates
#' and standardized random effects), so regimen comparisons are not blurred
#' by simulation noise.  With `mode = "trial"` (as in Winzenborg et al.
#' 2021), parameters stay at their estimates and each draw is a new trial of
#' `nSub` subjects resampled from the population, so the intervals describe
#' trial-to-trial variation.
#'
#' Exposure-response endpoints can take a per-subject exposure from a PK
#' model with `exposure = cuiExposure(...)`.
#'
#' Simulations use model predictions (no residual error).  The Monte Carlo
#' precision of the comparison between the two best regimens is estimated
#' from `nBatch` batches of subjects; a warning suggests a larger
#' population when its standard error exceeds `mcTol`.
#'
#' With a `seed`, every draw is cached (by default under
#' `nlmixr2CuiCache_<md5>` in the working directory), so a rerun loads the
#' simulations, an interrupted run resumes and a larger `nDraw` only
#' simulates the new draws.  The result itself keeps the simulated endpoint
#' values: `update()` re-weights, [cuiSensitivity()] and `plot()` work from
#' them without simulating again.
#'
#' @param endpoints List of [cuiEndpoint()] objects, each with a `fit`
#'   (an `nlmixr2` fit, or an `rxUi` with `uncertainty = "none"`) and a
#'   `summary` ([cuiProb()], [cuiMean()], [cuiAuc()], [cuiThreshold()]).
#'   Endpoints sharing a fit are simulated together.
#' @param regimens A [cuiRegimens()] table; each regimen doses `dose` at
#'   time 0, repeated every `ii` until `duration` (or the last observation).
#' @param population A [cuiPopulation()], or a number of subjects for a
#'   model without covariates.
#' @param times Named assessment times (model time units), for example
#'   `c(m3 = 2016, m6 = 4032)` in hours; summaries without their own time
#'   are evaluated at each.  `NULL` uses the regimens' common `duration`.
#' @param uncertainty `"none"`, `"cov"` (parameter uncertainty from the
#'   fits' covariance) or a [cuiBootstrap()] object.
#' @param nDraw Number of uncertainty draws or simulated trials (with a
#'   bootstrap, 0 uses every available draw).
#' @param mode `"expectation"` (expected endpoints in the population) or
#'   `"trial"` (finite simulated trials).
#' @param nSub Subjects per simulated trial in `mode = "trial"`; defaults
#'   to the population size.
#' @param exposure A [cuiExposure()] or a list of them.  Models that use a
#'   chained exposure receive it as a covariate and no dose records.
#' @param cmt Dosing compartment, or `NULL` for the model default.
#' @param nBatch Number of subject batches for the Monte Carlo standard
#'   error.
#' @param mcTol Warn when the Monte Carlo standard error of the CUI
#'   difference between the two best regimens exceeds this.
#' @param seed Seed, or `NULL` to take one from the current random stream.
#' @param minSuccess Minimum fraction of draws that must simulate.
#' @param cacheDir Directory caching each draw.  `NULL` (default) uses
#'   `nlmixr2CuiCache_<md5>` from the simulation's inputs when a `seed` is
#'   given (without a seed the draws cannot be reproduced, so nothing is
#'   cached); `NA` caches nothing.
#' @param restart Discard the cache and simulate every draw again.
#' @inheritParams cuiFromValues
#' @return An `nlmixr2Cui` object.
#' @references Winzenborg I, Soliman AM, Shebley M (2021).  A Personalized
#'   Medicine Approach Using Clinical Utility Index and Exposure-Response
#'   Modeling Informed by Patient Preferences Data.  CPT Pharmacometrics
#'   Syst Pharmacol 10:40-47.
#' @family Clinical utility
#' @export
#' @examples
#' \dontrun{
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45; tcl <- 1; tv <- 3.45
#'     eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
#'     cp <- linCmt()
#'     pae <- expit(-6 + 0.4 * cp)
#'     peff <- expit(-2 + 0.8 * cp)
#'     cp ~ add(add.sd)
#'   })
#' }
#' fit <- nlmixr2(one.cmt, nlmixr2data::theo_sd, "focei")
#' res <- cui(list(cuiEndpoint(eff, fit, cuiProb(peff, time = 4)),
#'                 cuiEndpoint(ae, fit, cuiProb(pae, time = 2), direction = "lower")),
#'            regimens = cuiRegimens(dose = seq(0, 600, by = 100)),
#'            population = 500, uncertainty = "cov", nDraw = 100, seed = 1)
#' res
#' }
cui <- function(endpoints, regimens, population, times = NULL, weights = NULL,
                combine = c("additive", "geometric"), eligible = NULL,
                uncertainty = "none", nDraw = 0L, mode = c("expectation", "trial"),
                nSub = NULL, exposure = NULL, cmt = NULL, nBatch = 10L, mcTol = 0.01,
                seed = NULL, minSuccess = 0.9, cacheDir = NULL, restart = FALSE,
                level = 0.95, tieTol = sqrt(.Machine$double.eps)) {
  # a cache is only useful for draws that can be reproduced
  if (is.null(seed) && is.null(cacheDir)) cacheDir <- NA
  .boot <- NULL
  if (inherits(uncertainty, "cuiBootstrap")) {
    .boot <- uncertainty
    uncertainty <- "bootstrap"
  } else {
    uncertainty <- match.arg(uncertainty, c("none", "cov"))
  }
  mode <- match.arg(mode)
  combine <- match.arg(combine)
  # one base seed from the current stream; every keyed stream derives from it
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  checkmate::assertCount(nDraw)
  checkmate::assertCount(nBatch)
  checkmate::assertNumber(mcTol, lower = 0)
  checkmate::assertNumber(minSuccess, lower = 0, upper = 1)
  checkmate::assertFlag(restart)
  if (inherits(endpoints, "cuiEndpoint")) endpoints <- list(endpoints)
  if (!all(vapply(endpoints, function(e) {
    if (!inherits(e, "cuiEndpoint") || !inherits(e$summary, "cuiSummary") || is.null(e$fit)) {
      return(FALSE)
    }
    .sim <- inherits(e$fit, "nlmixr2FitCore") || inherits(e$fit, "rxUi")
    .sim == (e$summary$type != "response")
  }, logical(1)))) {
    stop("each endpoint needs cuiEndpoint(name, fit, summary): an nlmixr2 fit or rxUi with a ",
         "cuiProb()/cuiMean()/cuiAuc()/cuiThreshold() summary, or another fitted model ",
         "with cuiResponse()", call. = FALSE)
  }
  if (inherits(exposure, "cuiExposure")) exposure <- list(exposure)
  if (!all(vapply(exposure, inherits, logical(1), what = "cuiExposure"))) {
    stop("'exposure' must be a cuiExposure() or a list of them", call. = FALSE)
  }
  if (!inherits(regimens, "cuiRegimens")) stop("'regimens' must come from cuiRegimens()", call. = FALSE)
  if (anyNA(regimens$dose) || any(regimens$dose < 0)) {
    stop("regimen doses must be nonnegative numbers", call. = FALSE)
  }
  if (checkmate::testCount(population, positive = TRUE)) {
    population <- structure(data.frame(id = seq_len(population)),
                            class = c("cuiPopulation", "data.frame"))
  }
  if (!inherits(population, "cuiPopulation")) {
    stop("'population' must come from cuiPopulation() or be a number of subjects", call. = FALSE)
  }
  # modes and draws
  if (mode == "trial") {
    if (uncertainty != "none") stop("mode = \"trial\" uses uncertainty = \"none\"", call. = FALSE)
    if (nDraw < 1) stop("mode = \"trial\" needs nDraw >= 1 simulated trials", call. = FALSE)
    if (is.null(nSub)) nSub <- nrow(population)
    checkmate::assertCount(nSub, positive = TRUE)
  } else if (uncertainty == "none" && nDraw > 0) {
    stop("nDraw > 0 needs uncertainty = \"cov\", a cuiBootstrap() or mode = \"trial\"", call. = FALSE)
  } else if (uncertainty == "cov" && nDraw < 1) {
    stop("uncertainty = \"cov\" needs nDraw >= 1", call. = FALSE)
  }
  # assessments
  if (is.null(times)) {
    .dur <- unique(regimens$duration)
    .needsTime <- vapply(endpoints, function(e) {
      e$summary$type != "response" && is.null(e$summary$time) && is.null(e$summary$window)
    }, logical(1))
    if (any(.needsTime) && (length(.dur) != 1L || is.na(.dur))) {
      stop("give 'times', or a time/window in every summary (regimens have no common duration)",
           call. = FALSE)
    }
    times <- c(all = if (length(.dur) == 1L && !is.na(.dur)) .dur else NA_real_)
  } else {
    checkmate::assertNumeric(times, finite = TRUE, any.missing = FALSE, min.len = 1)
    if (is.null(names(times))) names(times) <- .cuiNumLabel(times)
    if (anyDuplicated(names(times)) != 0L) stop("assessment names must be unique", call. = FALSE)
  }
  # model groups, providers of chained exposures first
  .g <- .cuiGroups(endpoints, exposure)
  .groups <- .g$groups
  .clash <- intersect(.g$exposures, names(population))
  if (length(.clash) > 0) {
    stop(sprintf("exposure name(s) also in 'population': %s", paste(.clash, collapse = ", ")),
         call. = FALSE)
  }
  .miss <- setdiff(unlist(lapply(.groups, `[[`, "covs")), c(names(population), .g$exposures))
  if (length(.miss) > 0) {
    stop(sprintf("model covariate(s) not in 'population': %s", paste(.miss, collapse = ", ")),
         call. = FALSE)
  }
  # parameter draws per group
  .fitDraws <- lapply(.groups, function(mg) {
    if (uncertainty != "cov" || !identical(mg$kind, "fitted")) return(NULL)
    .cuiWithSeed(.cuiStreamSeed(seed, "model", mg$group), cuiModelDraws(mg$fit, nDraw))
  })
  .covChol <- lapply(.groups, function(mg) {
    if (uncertainty != "cov" || identical(mg$kind, "fitted")) return(NULL)
    if (!mg$isFit || is.null(mg$fit$cov)) {
      stop("uncertainty = \"cov\" needs fits with a covariance matrix", call. = FALSE)
    }
    .cov <- mg$fit$cov
    if (!all(mg$struct %in% rownames(.cov))) {
      stop("the fit covariance does not cover every estimated structural theta", call. = FALSE)
    }
    tryCatch(chol(.cov[mg$struct, mg$struct, drop = FALSE]),
             error = function(e) stop("the fit covariance is not positive definite", call. = FALSE))
  })
  .bootDraws <- NULL
  .bootGroup <- rep(NA_integer_, length(.groups))
  .bootNote <- ""
  if (uncertainty == "bootstrap") {
    .ens <- .boot$ensembles
    .gMd5 <- lapply(.groups, `[[`, "md5")
    for (.i in seq_along(.ens)) {
      .m <- which(vapply(.gMd5, function(m) identical(m, .ens[[.i]]$parentMd5), logical(1)))
      if (length(.m) == 0 && is.na(.ens[[.i]]$parentMd5) && length(.ens) == 1L &&
        length(.groups) == 1L && identical(.groups[[1]]$kind, "nlmixr2")) {
        .m <- 1L
      }
      if (length(.m) != 1L) {
        stop(sprintf(paste0("bootstrap ensemble %d does not belong to an endpoint or exposure fit ",
                            "(refit it with bootstrapFit() on that fit)"), .i), call. = FALSE)
      }
      if (!is.na(.bootGroup[.m])) stop("two bootstrap ensembles for the same fit", call. = FALSE)
      .it <- .ens[[.i]]$items[[1]]
      if (!setequal(names(.it$theta), names(.groups[[.m]]$theta)) ||
        !setequal(rownames(.it$omega), rownames(.groups[[.m]]$omega))) {
        stop(sprintf("bootstrap ensemble %d has different parameters than its fit", .i), call. = FALSE)
      }
      .bootGroup[.m] <- .i
    }
    if (anyNA(.bootGroup)) {
      warning("models without a bootstrap ensemble keep their estimated parameters", call. = FALSE)
    }
    .pair <- .cuiBootPair(.ens)
    .bootDraws <- .pair$draws
    if (nDraw == 0) nDraw <- nrow(.bootDraws)
    if (nDraw > nrow(.bootDraws)) {
      stop(sprintf("nDraw = %d but only %d bootstrap draws are available", nDraw, nrow(.bootDraws)),
           call. = FALSE)
    }
    .bootNote <- if (length(.ens) > 1L) {
      if (.pair$paired) "; ensembles paired by resample" else "; independent ensembles"
    } else {
      ""
    }
  }
  .params <- function(d) {
    lapply(seq_along(.groups), function(.k) {
      .mg <- .groups[[.k]]
      if (identical(.mg$kind, "fitted")) {
        return(list(draw = if (d > 0 && uncertainty == "cov") .fitDraws[[.k]][[d]]))
      }
      if (d > 0 && uncertainty == "cov") {
        return(list(theta = .cuiDrawTheta(.mg, .covChol[[.k]], seed, d), omega = .mg$omega))
      }
      if (d > 0 && uncertainty == "bootstrap" && !is.na(.bootGroup[.k])) {
        .it <- .boot$ensembles[[.bootGroup[.k]]]$items[[.bootDraws[d, .bootGroup[.k]]]]
        return(list(theta = .it$theta[names(.mg$theta)],
                    omega = .it$omega[rownames(.mg$omega), rownames(.mg$omega), drop = FALSE]))
      }
      list(theta = .mg$theta, omega = .mg$omega)
    })
  }
  # observation times: assessments, fixed summary times and window ends
  .obsTimes <- times[!is.na(times)]
  for (.mg in .groups) for (.s in c(.mg$summaries, .mg$exposures)) {
    .obsTimes <- c(.obsTimes, .s$time, .s$window)
  }
  .obsTimes <- sort(unique(unname(.obsTimes)))
  # cache keyed on everything that changes what a draw is; nDraw is left out
  # so a larger nDraw resumes
  .cache <- .cuiCacheSetup(cacheDir, restart, list(
    lapply(.groups, function(mg) {
      .sums <- lapply(c(mg$summaries, mg$exposures),
                      function(s) s[c("type", "expr", "time", "window", "doseVar")])
      if (identical(mg$kind, "fitted")) return(list(digest::digest(mg$fit), .sums))
      list(mg$ui$iniDf, mg$ui$lstExpr, .sums)
    }),
    as.data.frame(population), as.data.frame(regimens), times, mode, nSub, uncertainty,
    .covChol, cmt, nBatch, seed,
    as.character(utils::packageVersion("rxode2")),
    as.character(utils::packageVersion("nlmixr2extra"))))
  # draws
  .nPop <- nrow(population)
  .vals <- list()
  .batch <- NULL
  .status <- rep("ok", nDraw + 1L)
  .msg <- rep("", nDraw + 1L)
  for (.d in 0:nDraw) {
    .pd <- .params(.d)
    .check <- digest::digest(.pd)
    .res <- .cuiCacheRead(.cache, .d, .check)
    if (is.null(.res)) {
      .trial <- mode == "trial" && .d > 0
      .rows <- if (.trial) {
        .cuiWithSeed(.cuiStreamSeed(seed, "population", .d), sample.int(.nPop, nSub, replace = TRUE))
      } else {
        seq_len(.nPop)
      }
      .n <- length(.rows)
      .bt <- if (.d == 0 && nBatch >= 2 && .n >= 2 * nBatch) (seq_len(.n) - 1L) %% nBatch + 1L else NULL
      .res <- tryCatch(
        list(ok = TRUE, check = .check,
             out = .cuiRunDraw(.d, .groups, .pd, .rows, if (.trial) .d else 0L,
                                          population, regimens, times, .obsTimes, cmt, .bt, seed)),
        error = function(e) {
          if (.d == 0) stop(e)
          list(ok = FALSE, message = conditionMessage(e))
        })
      # failed draws are not cached, so a rerun retries them
      if (.res$ok) .cuiCacheWrite(.cache, .d, .res)
    }
    if (!.res$ok) {
      .status[.d + 1L] <- "failed"
      .msg[.d + 1L] <- .res$message
      next
    }
    .vals <- c(.vals, lapply(.res$out, `[[`, "value"))
    if (.d == 0) .batch <- do.call(rbind, lapply(.res$out, `[[`, "batch"))
  }
  .nOk <- sum(.status[-1] == "ok")
  if (nDraw > 0 && .nOk / nDraw < minSuccess) {
    stop(sprintf("only %d of %d draws simulated (minSuccess = %g); first error: %s",
                 .nOk, nDraw, minSuccess, .msg[.status != "ok"][1]), call. = FALSE)
  }
  .desc <- sprintf("Simulated: %d subjects x %d regimens; %s", .nPop, nrow(regimens),
                   switch(if (mode == "trial") "trial" else uncertainty,
                          trial = sprintf("%d simulated trials of %d subjects", nDraw, nSub),
                          cov = sprintf("%d parameter draws from the fit covariance", nDraw),
                          bootstrap = sprintf("%d bootstrap parameter draws%s", nDraw, .bootNote),
                          none = "no parameter uncertainty"))
  .ret <- .cuiNew(values = do.call(rbind, .vals), endpoints = lapply(seq_along(endpoints), function(i) {
    .e <- endpoints[[i]]
    .e$summary <- .groups[[.g$epGroup[i]]]$summaries[[.e$name]]
    .e
  }), regimens = regimens, weights = weights, combine = combine, eligible = eligible,
  level = level, tieTol = tieTol,
  draws = data.frame(draw = 0:nDraw, status = .status,
                     source = if (mode == "trial") "trial" else uncertainty,
                     message = .msg, stringsAsFactors = FALSE),
  settings = list(source = "simulation", mode = mode, uncertainty = uncertainty,
                  nDraw = nDraw, nPopulation = .nPop, nSub = nSub, seed = seed,
                  times = times, exposures = .g$exposures, cacheDir = .cache,
                  description = .desc),
  batch = .batch)
  if (!is.null(.ret$mcse) && any(.ret$mcse$se > mcTol)) {
    warning(sprintf(paste0("Monte Carlo standard error of the difference between the two best ",
                           "regimens (%s) exceeds mcTol = %g; use a larger population"),
                    paste(signif(.ret$mcse$se, 2), collapse = ", "), mcTol), call. = FALSE)
  }
  .ret
}
