# Multistart estimation ----
#
# Refit a model from many perturbed starting points so that a suspected local
# optimum can be recognised.  The pieces here are deliberately small and
# testable: `.msPerturbIni()` builds the candidate starting points,
# `.msRunOne()` fits one of them, and `.multistartRun()` orchestrates
# screening, fitting, caching and summarisation.

#' Control options for `multistart()`
#'
#' @param n Number of candidate starting points to generate.  The first
#'   candidate is always the unperturbed starting point, so `n = 1` reproduces
#'   the original fit.
#' @param nFit Number of candidates to fully estimate after screening.  `NULL`
#'   (the default) estimates every candidate.  Ignored when `screen = "none"`.
#' @param sampling How the starting points are drawn around the initial
#'   estimates: `"uniform"` (the default) draws uniformly within `spread`, `"lhs"`
#'   uses a Latin hypercube over the same interval so the range is covered more
#'   evenly, and `"normal"` draws normally with a standard deviation of `spread`.
#' @param spread Fractional spread of the perturbation on the estimation scale.
#'   A parameter with initial estimate `est` is perturbed within
#'   `est +/- spread*max(abs(est), 1)`.
#' @param which Names of the population parameters to perturb; `NULL` (the
#'   default) perturbs every unfixed population parameter.
#' @param perturbOmega Should the between-subject variability estimates be
#'   perturbed as well?
#' @param omegaFold Fold-range for the `omega` perturbation.  A variance `v` is
#'   drawn within `v/omegaFold` and `v*omegaFold`.
#' @param around When starting from a fit, perturb around the fit's `"final"`
#'   estimates (the default, which asks "is this a local optimum?") or around
#'   the `"initial"` estimates the fit started from (which asks "how sensitive
#'   was this fit to where I started?").
#' @param screen Cheap pre-selection of candidates.  `"posthoc"` (the default)
#'   evaluates the objective function at each candidate with an empirical Bayes
#'   step only and fully estimates the best `nFit`; `"none"` estimates every
#'   candidate.
#' @param refitBest Re-run the best start with the full estimation control, so
#'   the returned fit has the covariance step and tables the exploratory runs
#'   skip.
#' @param excludeBoundary Should fits with a parameter at a boundary be excluded
#'   when picking the best start?  Matches [getMinAICFit()].
#' @param keepFits Keep each start's fit in the result.  Setting this to `FALSE`
#'   keeps only the summary and the best fit, which is much smaller.
#' @param seed Integer seed.  The perturbations are drawn from this seed and
#'   each start is estimated with its own derived seed.
#' @param cores Number of starts to estimate at once.  Each estimation already
#'   uses every available thread internally, so the default of `1` is usually
#'   the fastest choice; see the "Parallel estimation" section of
#'   [multistart()].
#' @param cacheDir Directory used to cache the individual starts so that an
#'   interrupted run can be resumed.  `NULL` (the default) derives a name from
#'   the model; `NA` disables caching.
#' @param restart Discard any cached results and start over.
#'
#' @returns A validated list of control options for [multistart()]
#' @family Multistart
#' @seealso [nlmixr2extra::multistart()]
#' @export
multistartControl <- function(n = 10L,
                              nFit = NULL,
                              sampling = c("uniform", "lhs", "normal"),
                              spread = 0.2,
                              which = NULL,
                              perturbOmega = TRUE,
                              omegaFold = 2,
                              around = c("final", "initial"),
                              screen = c("posthoc", "none"),
                              refitBest = TRUE,
                              excludeBoundary = TRUE,
                              keepFits = TRUE,
                              seed = 1234L,
                              cores = 1L,
                              cacheDir = NULL,
                              restart = FALSE) {
  ret <-
    list(
      n = checkmate::assert_integerish(n, lower = 1, len = 1, any.missing = FALSE, null.ok = FALSE, coerce = TRUE),
      nFit = checkmate::assert_integerish(nFit, lower = 1, len = 1, any.missing = FALSE, null.ok = TRUE, coerce = TRUE),
      sampling = match.arg(sampling),
      spread = checkmate::assert_number(spread, lower = 1e-8, upper = 10, finite = TRUE, na.ok = FALSE, null.ok = FALSE),
      which = checkmate::assert_character(which, any.missing = FALSE, min.len = 1, null.ok = TRUE),
      perturbOmega = checkmate::assert_logical(perturbOmega, len = 1, any.missing = FALSE, null.ok = FALSE),
      omegaFold = checkmate::assert_number(omegaFold, lower = 1.000001, upper = 1000, finite = TRUE, na.ok = FALSE, null.ok = FALSE),
      around = match.arg(around),
      screen = match.arg(screen),
      refitBest = checkmate::assert_logical(refitBest, len = 1, any.missing = FALSE, null.ok = FALSE),
      excludeBoundary = checkmate::assert_logical(excludeBoundary, len = 1, any.missing = FALSE, null.ok = FALSE),
      keepFits = checkmate::assert_logical(keepFits, len = 1, any.missing = FALSE, null.ok = FALSE),
      seed = checkmate::assert_integerish(seed, len = 1, any.missing = FALSE, null.ok = FALSE, coerce = TRUE),
      cores = checkmate::assert_integerish(cores, lower = 1, len = 1, any.missing = FALSE, null.ok = FALSE, coerce = TRUE),
      cacheDir = .msAssertCacheDir(cacheDir),
      restart = checkmate::assert_logical(restart, len = 1, any.missing = FALSE, null.ok = FALSE)
    )
  if (is.null(ret$nFit)) {
    ret$nFit <- ret$n
  } else if (ret$nFit > ret$n) {
    warning("'nFit' is larger than 'n'; using nFit = ", ret$n,
            call. = FALSE)
    ret$nFit <- ret$n
  }
  class(ret) <- "multistartControl"
  ret
}

#' @export
rxUiDeparse.multistartControl <- function(object, var) {
  .default <- multistartControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}

# `NA` disables caching, `NULL` asks for a derived name, anything else must be a
# single directory name.  `assert_string` would reject both sentinels.
.msAssertCacheDir <- function(cacheDir) {
  if (is.null(cacheDir)) {
    return(NULL)
  }
  if (length(cacheDir) == 1L && is.na(cacheDir)) {
    return(NA_character_)
  }
  checkmate::assert_string(cacheDir, min.chars = 1)
}

# Multistart ----

#' Estimate a model from many starting points
#'
#' Refits a model from several perturbed sets of initial estimates and collects
#' the results, so that a fit which settled in a local optimum can be
#' recognised.  Use [plot()] on the result for the objective-function waterfall
#' and the parameter-stability plots.
#'
#' @details
#'
#' # Starting points
#'
#' The first candidate is always the unperturbed starting point, so the original
#' fit is always represented in the comparison.  Every other candidate perturbs
#' the unfixed population parameters on the scale the model is estimated on
#' (which for a mu-referenced parameter is usually the log scale), and clips the
#' result to the parameter's declared bounds.  Between-subject variances are
#' perturbed multiplicatively so they stay positive, and the resulting matrix is
#' made positive-definite with [lotri::lotriNearPD()].
#'
#' # Screening
#'
#' Fully estimating every candidate is wasteful when many of them start far from
#' anywhere sensible.  With the default `screen = "posthoc"` each candidate is
#' first evaluated with an empirical-Bayes step only, which costs a small
#' fraction of a full estimation, and only the best `nFit` candidates are then
#' fully estimated.
#'
#' # Parallel estimation
#'
#' Each estimation already runs across every available thread, so estimating
#' several starts at once oversubscribes the machine unless the thread budget is
#' divided.  `multistartControl(cores=)` therefore restricts each worker to a
#' single thread.  Whether that is faster than the serial default depends
#' entirely on the model; a model whose subjects parallelise well is usually
#' better off left serial.  Parallel estimation uses [parallel::mclapply()] and
#' is not available on Windows.
#'
#' # Resuming
#'
#' Each start is cached to `cacheDir` as it completes, so an interrupted run
#' resumes where it left off.  Increasing `n` on a later call re-uses the starts
#' already estimated and only estimates the new ones.  Pass `restart = TRUE` to
#' discard the cache, or `cacheDir = NA` to never write one.
#'
#' Changing anything that alters what a starting point *is* (`sampling`,
#' `spread`, `which`, `perturbOmega`, `omegaFold` or `seed`) gives the run its
#' own cache, so a cached estimation is never re-used for a start it did not
#' come from.  A Latin hypercube is a design over all `n` candidates at once, so
#' growing an `"lhs"` run moves its earlier starting points; those starts are
#' detected and re-estimated rather than reported against the wrong starting
#' point.
#'
#' @param object A nlmixr2 fit, a nlmixr2 model function, or a `rxUi` model
#'   object
#' @param data The data to estimate with; taken from `object` when it is a fit
#' @param est The estimation method; taken from `object` when it is a fit
#' @param estControl The control for `est`; taken from `object` when it is a
#'   fit, and otherwise the method's default
#' @param control A list passed to [multistartControl()]
#' @param ... ignored
#'
#' @returns An object of class `nlmixr2Multistart`, a list with elements
#'   `starts` (one row per candidate starting point), `summary` (one row per
#'   estimated start), `fits`, `best`, and `bestIndex`
#'
#' @family Multistart
#' @author Matthew Fidler
#' @examples
#' \dontrun{
#' # Every start is a full estimation, so this takes a few minutes.
#' fit <- nlmixr2extra::theoFitOde
#'
#' ms <- multistart(fit, control = list(n = 8, spread = 0.3))
#' ms
#'
#' # objective function values, best to worst
#' plot(ms)
#'
#' # how stable each parameter is across the best starts
#' plot(ms, "parameters")
#'
#' # the best fit found, ready to use like any other fit
#' ms$best
#' }
#' @export
multistart <- function(object, ...) {
  UseMethod("multistart")
}

#' @rdname multistart
#' @export
multistart.nlmixr2FitCore <- function(object, ..., data = NULL, est = NULL,
                                      estControl = NULL, control = list()) {
  control <- do.call(multistartControl, control)
  if (control$around == "final") {
    ui <- object$finalUiEnv
  } else {
    ui <- rxode2::rxUiDecompress(object$iniUi)
  }
  if (is.null(data)) data <- nlme::getData(object)
  if (is.null(est)) est <- getFitMethod(object)
  if (is.null(estControl)) estControl <- object$control
  .multistartRun(ui = ui, data = data, est = est, estControl = estControl,
                 control = control, origFit = object)
}

#' @rdname multistart
#' @export
multistart.rxUi <- function(object, data, ..., est = "focei",
                            estControl = NULL, control = list()) {
  control <- do.call(multistartControl, control)
  checkmate::assert_data_frame(data, min.rows = 1)
  if (is.null(estControl)) estControl <- .msDefaultControl(est)
  .multistartRun(ui = rxode2::rxUiDecompress(object), data = data, est = est,
                 estControl = estControl, control = control, origFit = NULL)
}

#' @rdname multistart
#' @export
multistart.function <- function(object, data, ...) {
  multistart(rxode2::as.rxUi(object), data = data, ...)
}

#' @rdname multistart
#' @export
multistart.default <- function(object, ...) {
  stop("'multistart()' needs a nlmixr2 fit or a nlmixr2 model",
       call. = FALSE)
}

# The default control for an estimation method, e.g. foceiControl() for
# "focei".  Falls back to an empty list for a method with no control function,
# which nlmixr2() then fills in itself.
.msDefaultControl <- function(est) {
  .fn <- paste0(est, "Control")
  if (exists(.fn, envir = asNamespace("nlmixr2est"), mode = "function")) {
    do.call(get(.fn, envir = asNamespace("nlmixr2est"), mode = "function"), list())
  } else {
    list()
  }
}

# Perturbation ----

# One Latin hypercube sample on (0, 1): one draw from each of `n` equal strata,
# in random order.
.msLhsUnit <- function(n) {
  sample(((seq_len(n) - 1L) + stats::runif(n)) / n)
}

# `n` draws centred on 0 with a scale of 1, by sampling method.  The caller
# multiplies these by the spread it wants, so one mapping handles all three
# methods.  The uniform and Latin hypercube draws stay inside `[-1, 1]`; the
# normal draws do not, and rely on the caller's clipping to stay sensible.
.msSpreadDraws <- function(n, sampling) {
  if (sampling == "lhs") {
    2 * .msLhsUnit(n) - 1
  } else if (sampling == "normal") {
    stats::rnorm(n)
  } else {
    stats::runif(n, -1, 1)
  }
}

# Clip to the declared bounds, keeping a small margin so an estimate never
# starts exactly on a boundary.  Matches the margin FOCEi uses when it resamples
# after a failure.  `lower` and `upper` are single values; `x` may be a vector.
.msClip <- function(x, lower, upper) {
  margin <- .Machine$double.eps^(1 / 7)
  if (is.finite(lower)) x <- pmax(x, lower + margin)
  if (is.finite(upper)) x <- pmin(x, upper - margin)
  x
}

#' Build the perturbed starting points for a multistart run
#'
#' @param iniDf The `iniDf` of the model to perturb
#' @param control A `multistartControl()` object
#' @param omegaSameMap The model's `omegaSameMap`; when non-`NULL` the model has
#'   repeated `same()` blocks and `omega` is left alone
#' @returns A list of `n` `iniDf` data frames.  The first is `iniDf` unchanged.
#' @noRd
.msPerturbIni <- function(iniDf, control, omegaSameMap = NULL) {
  n <- control$n
  ret <- vector(mode = "list", length = n)
  ret[[1]] <- iniDf
  if (n == 1L) return(ret)
  nPerturb <- n - 1L

  isTheta <- !is.na(iniDf$ntheta)
  isEtaDiag <- !is.na(iniDf$neta1) & iniDf$neta1 == iniDf$neta2
  # `fix` is never perturbed; neither is anything the user left out of `which`
  eligible <- !iniDf$fix
  if (!is.null(control$which)) {
    unknown <- setdiff(control$which, iniDf$name)
    if (length(unknown) > 0) {
      stop("'which' names parameters that are not in the model: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    eligible <- eligible & iniDf$name %in% control$which
  }

  perturbOmega <- control$perturbOmega
  if (perturbOmega && !is.null(omegaSameMap)) {
    warning("the model uses 'same()' variability blocks; leaving 'omega' at its initial estimates",
            call. = FALSE)
    perturbOmega <- FALSE
  }

  thetaRows <- which(isTheta & eligible)
  etaRows <- if (perturbOmega) which(isEtaDiag & eligible) else integer(0)

  changed <- c(thetaRows, etaRows)
  if (length(changed) == 0L) {
    warning("no parameters left to perturb; every start is the same",
            call. = FALSE)
    for (j in seq_len(nPerturb)) ret[[j + 1L]] <- iniDf
    return(ret)
  }

  # Draw candidate-major, so candidate `j` occupies a fixed stretch of the
  # stream: asking for more starts later appends draws instead of shifting the
  # ones already used, which is what lets a cached run be resumed with a larger
  # `n`.  A Latin hypercube is a joint design over all the candidates, so it
  # cannot have that property and is drawn per parameter.
  nEl <- length(changed)
  if (control$sampling == "lhs") {
    z <- vapply(changed, function(i) .msSpreadDraws(nPerturb, "lhs"),
                numeric(nPerturb))
    dim(z) <- c(nPerturb, nEl)
  } else {
    z <- matrix(.msSpreadDraws(nPerturb * nEl, control$sampling),
                nrow = nPerturb, ncol = nEl, byrow = TRUE)
  }

  logFold <- log(control$omegaFold)
  for (j in seq_len(nPerturb)) {
    cur <- iniDf
    for (k in seq_len(nEl)) {
      i <- changed[k]
      est <- iniDf$est[i]
      if (i %in% etaRows) {
        # multiplicative, so a variance stays positive whatever is drawn
        val <- est * exp(logFold * z[j, k])
      } else {
        # the floor keeps a parameter that starts at 0 from being frozen there
        val <- est + control$spread * max(abs(est), 1) * z[j, k]
      }
      cur$est[i] <- .msClip(val, iniDf$lower[i], iniDf$upper[i])
    }
    ret[[j + 1L]] <- .msFixOmega(cur)
  }
  ret
}

# Scaling the variances on their own leaves the covariances implying
# correlations that are too large, so rebuild each variability block and make it
# positive definite again.  Blocks without off-diagonal terms cannot go wrong
# and are left alone.
.msFixOmega <- function(iniDf) {
  isEta <- !is.na(iniDf$neta1)
  isOff <- isEta & iniDf$neta1 != iniDf$neta2
  if (!any(isOff)) return(iniDf)
  # each `condition` (id, occasion, ...) is its own block
  for (cond in unique(iniDf$condition[isOff])) {
    rows <- which(isEta & iniDf$condition %in% cond)
    idx <- sort(unique(c(iniDf$neta1[rows], iniDf$neta2[rows])))
    m <- matrix(0, nrow = length(idx), ncol = length(idx))
    for (i in rows) {
      r <- match(iniDf$neta1[i], idx)
      cc <- match(iniDf$neta2[i], idx)
      m[r, cc] <- m[cc, r] <- iniDf$est[i]
    }
    ev <- try(min(eigen(m, symmetric = TRUE, only.values = TRUE)$values), silent = TRUE)
    if (inherits(ev, "try-error") || !is.finite(ev) || ev > 1e-8) next
    fixed <- try(lotri::lotriNearPD(m), silent = TRUE)
    if (inherits(fixed, "try-error")) next
    for (i in rows) {
      iniDf$est[i] <- fixed[match(iniDf$neta1[i], idx), match(iniDf$neta2[i], idx)]
    }
  }
  iniDf
}

# Estimation helpers ----

# `setQuietFastControl()` forces `covMethod <- 0L`, which is right for FOCEi but
# wrong for a method whose covMethod is a character choice (saem's is one of
# "", "linFim", "fim").  Keep everything else it does.
.msQuietControl <- function(ctl) {
  wasChar <- !is.null(ctl$covMethod) && is.character(ctl$covMethod)
  orig <- ctl$covMethod
  ctl <- setQuietFastControl(ctl)
  if (wasChar) ctl$covMethod <- ""
  else if (is.null(orig)) ctl$covMethod <- NULL
  ctl
}

# The objective function of a fit, taking the row of `objDf` that matches the
# fit's current objective-function type.  saem leaves that row NA until the
# objective is computed, which reading `fit$objf` triggers.
.msOfv <- function(fit) {
  objDf <- fit$objDf
  if (is.null(objDf) || nrow(objDf) == 0L) return(NA_real_)
  w <- which(tolower(rownames(objDf)) == tolower(fit$ofvType))
  if (length(w) != 1L) w <- 1L
  ret <- objDf[w, "OBJF"]
  if (is.na(ret)) {
    # lazily computed for saem; costs a quadrature evaluation
    ret <- try(suppressWarnings(fit$objf), silent = TRUE)
    if (inherits(ret, "try-error") || length(ret) != 1L) return(NA_real_)
    ret <- as.numeric(ret)
  }
  ret
}

# A row of the summary table for one estimated start.  `fit` may be a
# try-error, in which case everything but the failure message is NA.
.msSummaryRow <- function(fit, index, elapsed) {
  if (inherits(fit, "try-error") || !inherits(fit, "nlmixr2FitCore")) {
    msg <- if (inherits(fit, "try-error")) trimws(conditionMessage(attr(fit, "condition"))) else "not estimated"
    return(data.frame(start = index, OBJF = NA_real_, AIC = NA_real_, BIC = NA_real_,
                      converged = NA, boundary = NA, covMethod = NA_character_,
                      message = msg, elapsed = elapsed, stringsAsFactors = FALSE))
  }
  objDf <- fit$objDf
  w <- which(tolower(rownames(objDf)) == tolower(fit$ofvType))
  if (length(w) != 1L) w <- 1L
  conv <- fit$convergence
  msg <- fit$message
  data.frame(
    start = index,
    OBJF = .msOfv(fit),
    AIC = if ("AIC" %in% names(objDf)) objDf[w, "AIC"] else NA_real_,
    BIC = if ("BIC" %in% names(objDf)) objDf[w, "BIC"] else NA_real_,
    converged = if (is.null(conv) || length(conv) != 1L) NA else isTRUE(conv == 0),
    boundary = isBoundaryFit(fit),
    covMethod = if (is.null(fit$covMethod)) NA_character_ else paste(fit$covMethod, collapse = "; "),
    message = if (is.null(msg) || length(msg) != 1L) "" else msg,
    elapsed = elapsed,
    stringsAsFactors = FALSE
  )
}

# The estimates of one fit as a one-row data.frame: fixed effects plus the
# between-subject variances.  Only the "Estimate" column of `parFixedDf` is
# used; the other columns come and go with the model and the covariance step.
.msEstimateRow <- function(fit) {
  if (!inherits(fit, "nlmixr2FitCore")) return(NULL)
  pf <- fit$parFixedDf
  ret <- list()
  if (!is.null(pf) && "Estimate" %in% names(pf)) {
    ret <- as.list(stats::setNames(pf[["Estimate"]], rownames(pf)))
  }
  om <- fit$omega
  if (!is.null(om) && nrow(om) > 0L) {
    ret <- c(ret, stats::setNames(as.list(diag(om)), colnames(om)))
  }
  if (length(ret) == 0L) return(NULL)
  as.data.frame(ret, stringsAsFactors = FALSE, check.names = FALSE)
}

# The starting estimates of one candidate as a one-row data.frame.
.msStartRow <- function(iniDf) {
  keep <- !is.na(iniDf$ntheta) | (!is.na(iniDf$neta1) & iniDf$neta1 == iniDf$neta2)
  as.data.frame(as.list(stats::setNames(iniDf$est[keep], iniDf$name[keep])),
                stringsAsFactors = FALSE, check.names = FALSE)
}

# cbind that tolerates a NULL half; `cbind(data.frame(), NULL)` errors rather
# than returning the non-NULL side.
.msCbind <- function(a, b) {
  if (is.null(b)) a else if (is.null(a)) b else cbind(a, b)
}

# rbind a list of one-row data frames whose columns may differ, padding the
# missing ones with NA.  do.call(rbind, ...) errors on mismatched names.
.msRbindFill <- function(lst) {
  lst <- lst[!vapply(lst, is.null, logical(1))]
  if (length(lst) == 0L) return(NULL)
  nms <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(x) {
    miss <- setdiff(nms, names(x))
    # logical NA rather than NA_real_, so padding a character column still
    # rbinds to a character column
    for (m in miss) x[[m]] <- NA
    x[, nms, drop = FALSE]
  }))
}

# Estimate one start.  Never throws: a failure comes back as a try-error so the
# rest of the run continues.
.msRunOne <- function(ui, iniDf, data, est, estControl, seed, quiet = TRUE) {
  cur <- rxode2::rxUiDecompress(rxode2::.copyUi(ui))
  suppressMessages(rxode2::ini(cur) <- iniDf)
  ctl <- estControl
  # each start gets its own stream; the estimators run inside rxWithSeed(), so
  # this is the only thing that varies their Monte-Carlo draws
  if ("seed" %in% names(ctl)) ctl$seed <- seed
  # do not inherit the previous fit's ETAs, which would correlate the starts
  if ("etaMat" %in% names(ctl)) ctl$etaMat <- NA
  fn <- function() {
    nlmixr2est::nlmixr2(cur, data, est = est, control = ctl)
  }
  if (quiet) {
    try(suppressWarnings(suppressMessages(fn())), silent = TRUE)
  } else {
    try(suppressWarnings(fn()), silent = TRUE)
  }
}

# Cache ----

.msCacheFile <- function(dir, what, index) {
  file.path(dir, sprintf("%s_%04d.rds", what, index))
}

.msCacheSetup <- function(control, ui, est, origFit) {
  if (length(control$cacheDir) == 1L && is.na(control$cacheDir)) return(NULL)
  dir <- control$cacheDir
  if (is.null(dir)) {
    # Everything that changes what a start *is* goes into the key, so changing
    # the spread or the sampling method cannot silently re-use stale fits.  `n`,
    # `nFit` and the run-time options are deliberately left out: varying those
    # is exactly what resuming a cached run means.
    md5 <- digest::digest(list(ui$iniDf, ui$lstExpr, est,
                               control[c("sampling", "spread", "which",
                                         "perturbOmega", "omegaFold", "seed")]))
    dir <- paste0("nlmixr2MultistartCache_", md5)
  }
  if (dir.exists(dir) && control$restart) {
    unlink(dir, recursive = TRUE, force = TRUE)
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  dir
}

.msCacheRead <- function(dir, what, index) {
  if (is.null(dir)) return(NULL)
  f <- .msCacheFile(dir, what, index)
  if (!file.exists(f)) return(NULL)
  ret <- try(readRDS(f), silent = TRUE)
  if (inherits(ret, "try-error")) NULL else ret
}

.msCacheWrite <- function(dir, what, index, value) {
  if (is.null(dir)) return(invisible(NULL))
  saveRDS(value, file = .msCacheFile(dir, what, index))
  invisible(NULL)
}

# Driver ----

# Apply `fn` over `seq_len(n)`, in parallel when asked for it and possible.
.msLapply <- function(n, cores, fn) {
  if (cores <= 1L) return(lapply(seq_len(n), fn))
  if (.Platform$OS.type == "windows") {
    warning("parallel multistart estimation is not available on Windows; estimating serially",
            call. = FALSE)
    return(lapply(seq_len(n), fn))
  }
  # mclapply forks, so the thread count each worker pins for itself (see the
  # setRxThreads() calls in the loops below) stays in that worker's process and
  # the parent's is left alone.
  parallel::mclapply(seq_len(n), fn, mc.cores = min(cores, n),
                     mc.preschedule = FALSE)
}

.multistartRun <- function(ui, data, est, estControl, control, origFit) {
  checkmate::assert_data_frame(data, min.rows = 1)
  iniDf <- ui$iniDf
  cacheDir <- .msCacheSetup(control, ui, est, origFit)
  quietControl <- .msQuietControl(estControl)

  # everything random happens inside one seeded block, so a run is reproducible
  # and does not disturb the caller's RNG state
  starts <- rxode2::rxWithSeed(control$seed, {
    .msPerturbIni(iniDf, control, omegaSameMap = ui$omegaSameMap)
  })
  n <- length(starts)
  seeds <- control$seed + seq_len(n)
  # A cached start that no longer matches the one we just built (which "lhs" can
  # do, since its design changes with `n`) makes the cached estimation for that
  # index meaningless, so drop it rather than report a fit against the wrong
  # starting point.
  for (i in seq_len(n)) {
    cached <- .msCacheRead(cacheDir, "start", i)
    if (!is.null(cached) && !isTRUE(all.equal(cached, starts[[i]]))) {
      unlink(c(.msCacheFile(cacheDir, "screen", i), .msCacheFile(cacheDir, "fit", i)))
    }
    .msCacheWrite(cacheDir, "start", i, starts[[i]])
  }

  # Screening ----
  screenOfv <- rep(NA_real_, n)
  if (control$screen == "posthoc" && control$nFit < n) {
    cli::cli_h1("Multistart screening ({n} candidates)")
    scr <- .msLapply(n, control$cores, function(i) {
      cached <- .msCacheRead(cacheDir, "screen", i)
      if (!is.null(cached)) return(cached)
      if (control$cores <= 1L) cli::cli_alert_info("screening start {i}/{n}")
      else rxode2::setRxThreads(1L)
      f <- .msRunOne(ui, starts[[i]], data, "posthoc", quietControl, seeds[i])
      v <- if (inherits(f, "nlmixr2FitCore")) .msOfv(f) else NA_real_
      if (!is.finite(v)) v <- Inf
      .msCacheWrite(cacheDir, "screen", i, v)
      v
    })
    screenOfv <- vapply(scr, function(x) as.numeric(x)[1], numeric(1))
    toFit <- order(screenOfv)[seq_len(control$nFit)]
  } else {
    toFit <- seq_len(min(control$nFit, n))
  }
  toFit <- sort(toFit)

  # Estimation ----
  cli::cli_h1("Multistart estimation ({length(toFit)} start{?s})")
  res <- .msLapply(length(toFit), control$cores, function(j) {
    i <- toFit[j]
    cached <- .msCacheRead(cacheDir, "fit", i)
    if (!is.null(cached)) {
      if (control$cores <= 1L) cli::cli_alert_success("start {i} loaded from cache")
      return(cached)
    }
    if (control$cores <= 1L) cli::cli_alert_info("estimating start {i} ({j}/{length(toFit)})")
    else rxode2::setRxThreads(1L)
    t0 <- proc.time()[["elapsed"]]
    f <- .msRunOne(ui, starts[[i]], data, est, quietControl, seeds[i])
    out <- list(fit = f, elapsed = proc.time()[["elapsed"]] - t0)
    .msCacheWrite(cacheDir, "fit", i, out)
    out
  })

  fits <- lapply(res, `[[`, "fit")
  names(fits) <- as.character(toFit)

  summary <- .msRbindFill(lapply(seq_along(toFit), function(j) {
    .msCbind(.msSummaryRow(res[[j]]$fit, toFit[j], res[[j]]$elapsed),
             .msEstimateRow(res[[j]]$fit))
  }))
  row.names(summary) <- NULL

  ok <- !is.na(summary$OBJF) & is.finite(summary$OBJF)
  if (!any(ok)) {
    stop("every multistart estimation failed; first message: ",
         summary$message[1], call. = FALSE)
  }
  eligible <- ok
  if (control$excludeBoundary && any(ok & !summary$boundary %in% TRUE)) {
    eligible <- ok & !summary$boundary %in% TRUE
  }
  bestRow <- which(eligible)[which.min(summary$OBJF[eligible])]
  bestIndex <- summary$start[bestRow]
  summary$dOBJF <- summary$OBJF - summary$OBJF[bestRow]
  summary <- summary[, c("start", "OBJF", "dOBJF",
                         setdiff(names(summary), c("start", "OBJF", "dOBJF")))]
  summary <- summary[order(summary$OBJF, na.last = TRUE), , drop = FALSE]
  row.names(summary) <- NULL

  best <- fits[[as.character(bestIndex)]]
  if (control$refitBest) {
    cli::cli_h1("Re-estimating the best start ({bestIndex}) with the full control")
    refit <- .msRunOne(ui, starts[[bestIndex]], data, est, estControl,
                       seeds[bestIndex], quiet = FALSE)
    if (inherits(refit, "nlmixr2FitCore")) {
      best <- refit
    } else {
      warning("could not re-estimate the best start with the full control; ",
              "returning the exploratory fit", call. = FALSE)
    }
  }

  startsDf <- .msRbindFill(lapply(seq_len(n), function(i) {
    .msCbind(data.frame(start = i, seed = seeds[i], screenOFV = screenOfv[i],
                        fitted = i %in% toFit, stringsAsFactors = FALSE),
             .msStartRow(starts[[i]]))
  }))
  row.names(startsDf) <- NULL

  ret <- list(
    starts = startsDf,
    summary = summary,
    fits = if (control$keepFits) fits else NULL,
    best = best,
    bestIndex = bestIndex,
    origFit = origFit,
    est = est,
    control = control,
    cacheDir = cacheDir
  )
  class(ret) <- "nlmixr2Multistart"
  ret
}

# Methods ----

#' @export
print.nlmixr2Multistart <- function(x, ..., n = 10L) {
  nFit <- nrow(x$summary)
  nFail <- sum(is.na(x$summary$OBJF))
  cli::cli_h1(cli::col_red("Multistart ({nrow(x$starts)} candidate{?s}, {nFit} estimated)"))
  cli::cli_li(cli::col_magenta(cli::style_bold("Best start: "), cli::col_yellow("{x$bestIndex}")))
  if (nFail > 0) {
    cli::cli_li(cli::col_red("{nFail} start{?s} failed to estimate"))
  }
  nBoundary <- sum(x$summary$boundary %in% TRUE)
  if (nBoundary > 0) {
    cli::cli_li(cli::col_yellow("{nBoundary} start{?s} finished with a parameter at a boundary"))
  }
  cli::cli_li(cli::col_magenta(cli::style_bold("Objective functions"),
                               cli::col_yellow(" (x$summary)")))
  cols <- intersect(c("start", "OBJF", "dOBJF", "AIC", "BIC", "converged", "boundary"),
                    names(x$summary))
  print(utils::head(x$summary[, cols, drop = FALSE], n))
  cli::cli_li(cli::col_magenta(cli::style_bold("Best fit"),
                               cli::col_yellow(" (x$best)")))
  cli::cli_h1("end")
  invisible(x)
}

#' @export
as.data.frame.nlmixr2Multistart <- function(x, ...) {
  x$summary
}

#' Plot a multistart result
#'
#' @param x A `nlmixr2Multistart` object from [multistart()]
#' @param type `"waterfall"` plots each start's objective function relative to
#'   the best one, worst to best; `"parameters"` plots how each parameter
#'   estimate varies across the best starts.
#' @param kBest Number of starts to show in the `"parameters"` plot
#' @param dOfvMax Upper limit for the waterfall's objective-function axis, for
#'   zooming in when one start is far worse than the rest
#' @param ... ignored
#' @returns A ggplot2 object
#' @family Multistart
#' @export
plot.nlmixr2Multistart <- function(x, type = c("waterfall", "parameters"),
                                   kBest = 20L, dOfvMax = NULL, ...) {
  type <- match.arg(type)
  if (type == "waterfall") {
    .msPlotWaterfall(x, dOfvMax = dOfvMax)
  } else {
    .msPlotParameters(x, kBest = kBest)
  }
}

# The status of each start, as an ordered factor so the fill scale is stable
# whichever statuses happen to be present.
.msStatus <- function(summary) {
  status <- rep("converged", nrow(summary))
  status[!summary$converged %in% TRUE] <- "not converged"
  status[summary$boundary %in% TRUE] <- "boundary issue"
  factor(status, levels = c("converged", "boundary issue", "not converged"))
}

.msStatusColors <- c("converged" = "#2c7fb8",
                     "boundary issue" = "#d95f02",
                     "not converged" = "#999999")

.msPlotWaterfall <- function(x, dOfvMax = NULL) {
  df <- x$summary[!is.na(x$summary$OBJF), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("no successfully estimated starts to plot", call. = FALSE)
  }
  df <- df[order(df$dOBJF), , drop = FALSE]
  df$rank <- seq_len(nrow(df))
  df$status <- .msStatus(df)
  nFail <- sum(is.na(x$summary$OBJF))

  sub <- paste0("Best objective function: ", signif(min(df$OBJF), 6))
  if (nFail > 0) {
    sub <- paste0(sub, "; ", nFail, " start", if (nFail > 1) "s" else "", " failed")
  }
  .plot <-
    ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$dOBJF,
                                     fill = .data$status)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(name = "Status", values = .msStatusColors) +
    ggplot2::scale_x_continuous(breaks = .msIntBreaks) +
    ggplot2::xlab("Start (best to worst)") +
    ggplot2::ylab("Delta objective function") +
    ggplot2::labs(title = "Multistart objective functions", subtitle = sub) +
    rxode2::rxTheme() +
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")
  if (!is.null(dOfvMax)) {
    .plot <- .plot + ggplot2::coord_cartesian(ylim = c(0, dOfvMax))
  }
  .plot
}

# Integer-only axis breaks; the x axis counts starts.
.msIntBreaks <- function(limits) {
  b <- unique(floor(pretty(limits)))
  b[b >= 1]
}

.msPlotParameters <- function(x, kBest = 20L) {
  df <- x$summary[!is.na(x$summary$OBJF), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("no successfully estimated starts to plot", call. = FALSE)
  }
  df <- df[order(df$dOBJF), , drop = FALSE]
  if (kBest < nrow(df)) df <- df[seq_len(kBest), , drop = FALSE]
  df$rank <- seq_len(nrow(df))
  df$status <- .msStatus(df)

  meta <- c("start", "OBJF", "dOBJF", "AIC", "BIC", "converged", "boundary",
            "covMethod", "message", "elapsed", "rank", "status")
  pars <- setdiff(names(df), meta)
  pars <- pars[vapply(df[pars], is.numeric, logical(1))]
  if (length(pars) == 0L) {
    stop("no parameter estimates were kept for these starts", call. = FALSE)
  }

  long <- do.call(rbind, lapply(pars, function(p) {
    data.frame(rank = df$rank, status = df$status, parameter = p,
               value = df[[p]], stringsAsFactors = FALSE)
  }))
  long$parameter <- factor(long$parameter, levels = pars)
  # the best start's value, as the reference line in each panel
  ref <- data.frame(parameter = factor(pars, levels = pars),
                    value = vapply(pars, function(p) df[[p]][1], numeric(1)),
                    stringsAsFactors = FALSE)

  ggplot2::ggplot(long, ggplot2::aes(x = .data$rank, y = .data$value)) +
    ggplot2::geom_hline(data = ref,
                        ggplot2::aes(yintercept = .data$value),
                        linetype = "dashed", color = "grey40") +
    ggplot2::geom_point(ggplot2::aes(color = .data$status)) +
    ggplot2::scale_color_manual(name = "Status", values = .msStatusColors) +
    ggplot2::scale_x_continuous(breaks = .msIntBreaks) +
    ggplot2::facet_wrap("parameter", scales = "free_y") +
    ggplot2::xlab("Start (best to worst)") +
    ggplot2::ylab("Estimate") +
    ggplot2::labs(title = "Multistart parameter stability",
                  subtitle = paste0("Best ", nrow(df), " start",
                                    if (nrow(df) > 1) "s" else "",
                                    "; dashed line is the best start")) +
    rxode2::rxTheme() +
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")
}
