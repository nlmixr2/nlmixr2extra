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
#' @param fit The endpoints' fit (an `nlmixr2` fit or an `rxUi`).
#' @param eps Endpoints sharing this fit.
#' @param group Group number (keys its random streams).
#' @return List with the simulation model, resolved summaries and parameters.
#' @noRd
.cuiModelGroup <- function(fit, eps, group) {
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
  .sums <- lapply(eps, function(.ep) {
    .s <- .cuiResolveSummary(.ep$summary, .cols, .ep$name)
    if (!is.null(.s$window)) {
      .state <- paste0("cuiInt", length(.lines) + 1L)
      while (.state %in% .cols) .state <- paste0(.state, "x")
      .lines[[length(.lines) + 1L]] <<- bquote(d / dt(.(as.name(.state))) <- .(.s$expr))
      .s$state <- .state
    }
    .s
  })
  names(.sums) <- vapply(eps, `[[`, character(1), "name")
  .sim <- .ui
  for (.l in .lines) {
    .sim <- suppressMessages(eval(bquote(rxode2::model(.(.sim), .(.l), append = TRUE))))
  }
  list(group = group, fit = fit, isFit = .isFit, ui = .sim, theta = .theta, struct = .struct,
       lower = setNames(.ini$lower, .ini$name), upper = setNames(.ini$upper, .ini$name),
       omega = .omega, covs = .ui$allCovs, summaries = .sums)
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
.cuiEvents <- function(regimens, nSub, obsTimes, cmt) {
  .maxT <- max(obsTimes)
  .ev <- lapply(seq_len(nrow(regimens)), function(.r) {
    .ids <- (.r - 1L) * nSub + seq_len(nSub)
    .reg <- regimens[.r, ]
    .obs <- data.frame(id = rep(.ids, each = length(obsTimes)),
                       time = rep(obsTimes, times = nSub),
                       amt = 0, evid = 0L, ii = 0, addl = 0L)
    if (.reg$dose == 0) return(.obs)
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
#' @return Data frame of simulated rows with `.reg` and `.sub` columns.
#' @noRd
.cuiSolveDraw <- function(mg, theta, eta, cov, regimens, obsTimes, cmt) {
  .nSub <- nrow(eta)
  .nReg <- nrow(regimens)
  .p <- as.data.frame(matrix(rep(theta, each = .nSub * .nReg), nrow = .nSub * .nReg,
                             dimnames = list(NULL, names(theta))))
  if (ncol(eta) > 0) {
    .p <- cbind(.p, eta[rep(seq_len(.nSub), times = .nReg), , drop = FALSE])
  }
  if (length(mg$covs) > 0) {
    .p <- cbind(.p, cov[rep(seq_len(.nSub), times = .nReg), mg$covs, drop = FALSE])
  }
  rownames(.p) <- NULL
  .ev <- .cuiEvents(regimens, .nSub, obsTimes, cmt)
  .s <- suppressMessages(rxode2::rxSolve(mg$ui, params = .p, events = .ev, nStud = 1,
                                         thetaMat = NULL, simVariability = FALSE, sigma = NA,
                                         returnType = "data.frame"))
  .id <- as.integer(.s$id)
  .s$.reg <- (.id - 1L) %/% .nSub + 1L
  .s$.sub <- (.id - 1L) %% .nSub + 1L
  .s
}

#' Evaluate one summary at one time for every regimen (and batch)
#'
#' @return List with `value` (per regimen) and `batch` (regimen x batch).
#' @noRd
.cuiEvalSummary <- function(s, sim, time, nReg, nSub, batch, ep) {
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
  } else {
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
      .v <- as.numeric(.v)
    } else if (!is.numeric(.v) || anyNA(.v)) {
      stop(sprintf("endpoint '%s' must give numeric values without NA", ep), call. = FALSE)
    } else if (s$type == "prob" && any(.v < 0 | .v > 1)) {
      stop(sprintf("endpoint '%s': cuiProb() values must be in [0, 1]", ep), call. = FALSE)
    }
  }
  .m <- matrix(.v, nrow = nSub, ncol = nReg)
  .ret <- list(value = colMeans(.m))
  if (!is.null(batch)) {
    .ret$batch <- rowsum(.m, batch, reorder = TRUE) / as.vector(table(batch))
  }
  .ret
}

#' Simulated clinical utility index from nlmixr2 models
#'
#' Simulates each endpoint's model for every candidate regimen in one
#' target population, summarizes the endpoints per regimen, and combines
#' them into a clinical utility index (CUI) as in [cuiFromValues()].
#'
#' The point estimate (draw 0) uses the fits' estimates.  With
#' `uncertainty = "cov"`, draws `1..nDraw` sample the structural thetas
#' from each fit's covariance (residual-error and fixed parameters, and
#' omega, stay at their estimates), giving intervals for the expected CUI
#' and selection frequencies.  Every regimen and every draw sees the same
#' virtual subjects (covariates and standardized random effects), so
#' regimen comparisons are not blurred by simulation noise.  With
#' `mode = "trial"` (as in Winzenborg et al. 2021), parameters stay at
#' their estimates and each draw is a new trial of `nSub` subjects
#' resampled from the population, so the intervals describe trial-to-trial
#' variation.
#'
#' Simulations use model predictions (no residual error).  The Monte Carlo
#' precision of the comparison between the two best regimens is estimated
#' from `nBatch` batches of subjects; a warning suggests a larger
#' population when its standard error exceeds `mcTol`.
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
#' @param uncertainty `"none"` or `"cov"` (parameter uncertainty from the
#'   fits' covariance).
#' @param nDraw Number of uncertainty draws (`"cov"`) or simulated trials
#'   (`mode = "trial"`).
#' @param mode `"expectation"` (expected endpoints in the population) or
#'   `"trial"` (finite simulated trials).
#' @param nSub Subjects per simulated trial in `mode = "trial"`; defaults
#'   to the population size.
#' @param cmt Dosing compartment, or `NULL` for the model default.
#' @param nBatch Number of subject batches for the Monte Carlo standard
#'   error.
#' @param mcTol Warn when the Monte Carlo standard error of the CUI
#'   difference between the two best regimens exceeds this.
#' @param seed Seed, or `NULL` for the current random stream.
#' @param minSuccess Minimum fraction of draws that must simulate.
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
                uncertainty = c("none", "cov"), nDraw = 0L, mode = c("expectation", "trial"),
                nSub = NULL, cmt = NULL, nBatch = 10L, mcTol = 0.01, seed = NULL,
                minSuccess = 0.9, level = 0.95, tieTol = sqrt(.Machine$double.eps)) {
  uncertainty <- match.arg(uncertainty)
  mode <- match.arg(mode)
  combine <- match.arg(combine)
  checkmate::assertCount(nDraw)
  checkmate::assertCount(nBatch)
  checkmate::assertNumber(mcTol, lower = 0)
  checkmate::assertNumber(minSuccess, lower = 0, upper = 1)
  if (inherits(endpoints, "cuiEndpoint")) endpoints <- list(endpoints)
  if (!all(vapply(endpoints, function(e) {
    inherits(e, "cuiEndpoint") && inherits(e$summary, "cuiSummary") &&
      (inherits(e$fit, "nlmixr2FitCore") || inherits(e$fit, "rxUi"))
  }, logical(1)))) {
    stop("each endpoint needs cuiEndpoint(name, fit, summary) with an nlmixr2 fit or rxUi ",
         "and a cuiProb()/cuiMean()/cuiAuc()/cuiThreshold() summary", call. = FALSE)
  }
  .epNames <- vapply(endpoints, `[[`, character(1), "name")
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
    stop("nDraw > 0 needs uncertainty = \"cov\" or mode = \"trial\"", call. = FALSE)
  } else if (uncertainty == "cov" && nDraw < 1) {
    stop("uncertainty = \"cov\" needs nDraw >= 1", call. = FALSE)
  }
  # assessments
  if (is.null(times)) {
    .dur <- unique(regimens$duration)
    .needsTime <- vapply(endpoints, function(e) {
      is.null(e$summary$time) && is.null(e$summary$window)
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
  # model groups: endpoints sharing a fit are simulated together
  .key <- lapply(endpoints, function(e) if (inherits(e$fit, "nlmixr2FitCore")) e$fit$env else e$fit)
  .grp <- integer(length(endpoints))
  for (.i in seq_along(endpoints)) {
    .grp[.i] <- .i
    for (.j in seq_len(.i - 1L)) {
      if (identical(.key[[.i]], .key[[.j]])) {
        .grp[.i] <- .grp[.j]
        break
      }
    }
  }
  .groups <- lapply(unique(.grp), function(g) {
    .cuiModelGroup(endpoints[[which(.grp == g)[1]]]$fit, endpoints[.grp == g], g)
  })
  .miss <- setdiff(unlist(lapply(.groups, `[[`, "covs")), names(population))
  if (length(.miss) > 0) {
    stop(sprintf("model covariate(s) not in 'population': %s", paste(.miss, collapse = ", ")),
         call. = FALSE)
  }
  .covChol <- lapply(.groups, function(mg) {
    if (uncertainty != "cov") return(NULL)
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
  # observation times: assessments, fixed summary times and window ends
  .obsTimes <- times[!is.na(times)]
  for (.mg in .groups) for (.s in .mg$summaries) {
    .obsTimes <- c(.obsTimes, .s$time, .s$window)
  }
  .obsTimes <- sort(unique(unname(.obsTimes)))
  # draws
  .nPop <- nrow(population)
  .etaZ <- function(mg, draw, n) {
    .nEta <- if (is.null(mg$omega)) 0L else nrow(mg$omega)
    if (.nEta == 0L) return(matrix(numeric(0), nrow = n, ncol = 0))
    .z <- .cuiWithSeed(.cuiStreamSeed(seed, "eta", mg$group, draw),
                       matrix(stats::rnorm(n * .nEta), nrow = n))
    .e <- .z %*% chol(mg$omega)
    colnames(.e) <- rownames(mg$omega)
    .e
  }
  .vals <- list()
  .batch <- NULL
  .status <- rep("ok", nDraw + 1L)
  .msg <- rep("", nDraw + 1L)
  for (.d in 0:nDraw) {
    .trial <- mode == "trial" && .d > 0
    .rows <- if (.trial) {
      .cuiWithSeed(.cuiStreamSeed(seed, "population", .d), sample.int(.nPop, nSub, replace = TRUE))
    } else {
      seq_len(.nPop)
    }
    .n <- length(.rows)
    .bt <- if (.d == 0 && nBatch >= 2 && .n >= 2 * nBatch) (seq_len(.n) - 1L) %% nBatch + 1L else NULL
    .res <- tryCatch({
      .out <- list()
      for (.gi in seq_along(.groups)) {
        .mg <- .groups[[.gi]]
        .theta <- if (uncertainty == "cov" && .d > 0) {
          .cuiDrawTheta(.mg, .covChol[[.gi]], seed, .d)
        } else {
          .mg$theta
        }
        .eta <- .etaZ(.mg, if (.trial) .d else 0L, .n)
        .sim <- .cuiSolveDraw(.mg, .theta, .eta, population[.rows, , drop = FALSE],
                              regimens, .obsTimes, cmt)
        for (.ep in names(.mg$summaries)) {
          .s <- .mg$summaries[[.ep]]
          for (.a in names(times)) {
            .t <- if (!is.null(.s$time)) .s$time else times[[.a]]
            .r <- .cuiEvalSummary(.s, .sim, .t, nrow(regimens), .n, .bt, .ep)
            .out[[length(.out) + 1L]] <- list(
              value = data.frame(assessment = .a, regimen = regimens$label, draw = .d,
                                 endpoint = .ep, value = .r$value, stringsAsFactors = FALSE),
              batch = if (!is.null(.r$batch)) {
                data.frame(assessment = .a, regimen = rep(regimens$label, each = nrow(.r$batch)),
                           endpoint = .ep, batch = rep(seq_len(nrow(.r$batch)), nrow(regimens)),
                           value = as.vector(.r$batch), stringsAsFactors = FALSE)
              })
          }
        }
      }
      .out
    }, error = function(e) {
      if (.d == 0) stop(e)
      conditionMessage(e)
    })
    if (is.character(.res)) {
      .status[.d + 1L] <- "failed"
      .msg[.d + 1L] <- .res
      next
    }
    .vals <- c(.vals, lapply(.res, `[[`, "value"))
    if (.d == 0) .batch <- do.call(rbind, lapply(.res, `[[`, "batch"))
  }
  .nOk <- sum(.status[-1] == "ok")
  if (nDraw > 0 && .nOk / nDraw < minSuccess) {
    stop(sprintf("only %d of %d draws simulated (minSuccess = %g); first error: %s",
                 .nOk, nDraw, minSuccess, .msg[.status != "ok"][1]), call. = FALSE)
  }
  .desc <- sprintf("Simulated: %d subjects x %d regimens; %s", .nPop, nrow(regimens),
                   if (mode == "trial") {
                     sprintf("%d simulated trials of %d subjects", nDraw, nSub)
                   } else if (uncertainty == "cov") {
                     sprintf("%d parameter draws from the fit covariance", nDraw)
                   } else {
                     "no parameter uncertainty"
                   })
  .ret <- .cuiNew(values = do.call(rbind, .vals), endpoints = lapply(seq_along(endpoints), function(i) {
    .e <- endpoints[[i]]
    .e$summary <- .groups[[match(.grp[i], unique(.grp))]]$summaries[[.e$name]]
    .e
  }), regimens = regimens, weights = weights, combine = combine, eligible = eligible,
  level = level, tieTol = tieTol,
  draws = data.frame(draw = 0:nDraw, status = .status,
                     source = if (mode == "trial") "trial" else uncertainty,
                     message = .msg, stringsAsFactors = FALSE),
  settings = list(source = "simulation", mode = mode, uncertainty = uncertainty,
                  nDraw = nDraw, nPopulation = .nPop, nSub = nSub, seed = seed,
                  times = times, description = .desc),
  batch = .batch)
  if (!is.null(.ret$mcse) && any(.ret$mcse$se > mcTol)) {
    warning(sprintf(paste0("Monte Carlo standard error of the difference between the two best ",
                           "regimens (%s) exceeds mcTol = %g; use a larger population"),
                    paste(signif(.ret$mcse$se, 2), collapse = ", "), mcTol), call. = FALSE)
  }
  .ret
}
