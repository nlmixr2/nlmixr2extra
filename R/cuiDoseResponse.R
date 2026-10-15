#' Binary endpoint for a dose-response clinical utility index
#'
#' @param name Column of the data holding the 0/1 endpoint, as a bare name
#'   (`tox`) or a string (`"tox"`); also the endpoint name.
#' @param model How the probability of `1` at each dose is estimated:
#'   `"empirical"` uses the observed rate at each dose.
#' @param direction `"higher"` when a `1` is desirable, `"lower"` when it is
#'   not (for example toxicity), giving utility `1 - p`.
#' @param ... Not used.
#' @param .name Column name as a string, for programmatic use instead of
#'   `name`.
#' @return A `cuiBinary` endpoint for [cuiDoseResponse()].
#' @family Clinical utility
#' @export
#' @examples
#' cuiBinary(tox, direction = "lower")
#' cuiBinary("eff")
cuiBinary <- function(name, model = "empirical", direction = c("higher", "lower"), ...,
                      .name = NULL) {
  .name <- .cuiCaptureName(substitute(name), missing(name), .name, "name")
  model <- match.arg(model)
  .ret <- cuiEndpoint(.name = .name, utility = cuiIdentity(), direction = match.arg(direction))
  .ret$model <- model
  class(.ret) <- c("cuiBinary", class(.ret))
  .ret
}

#' @export
print.cuiBinary <- function(x, ...) {
  cat(sprintf("cuiBinary '%s': %s model, %s is better\n", x$name, x$model, x$direction))
  invisible(x)
}

#' Evaluate code with a fixed seed, restoring the caller's random stream
#'
#' @param seed Seed, or `NULL` to use the current stream.
#' @param code Code to evaluate.
#' @return The value of `code`.
#' @noRd
.cuiWithSeed <- function(seed, code) {
  if (is.null(seed)) return(code)
  .env <- globalenv()
  .had <- exists(".Random.seed", envir = .env, inherits = FALSE)
  if (.had) .old <- get(".Random.seed", envir = .env, inherits = FALSE)
  on.exit({
    if (.had) {
      assign(".Random.seed", .old, envir = .env)
    } else if (exists(".Random.seed", envir = .env, inherits = FALSE)) {
      rm(".Random.seed", envir = .env)
    }
  })
  set.seed(seed)
  code
}

#' Estimate the probability of each endpoint at each dose
#'
#' @param y Matrix of 0/1 outcomes, one column per endpoint.
#' @param dose Dose of each row.
#' @param doses Sorted unique doses.
#' @param endpoints List of `cuiBinary` endpoints matching the columns.
#' @return Matrix of probabilities, doses x endpoints.
#' @noRd
.cuiDrFit <- function(y, dose, doses, endpoints) {
  .n <- as.vector(table(factor(dose, levels = doses)))
  .p <- rowsum(y, factor(dose, levels = doses), reorder = TRUE) / .n
  for (.k in seq_along(endpoints)) {
    .model <- endpoints[[.k]]$model
    if (.model != "empirical") {
      stop(sprintf("unknown dose-response model '%s'", .model), call. = FALSE)
    }
  }
  .p
}

#' Fit the original data and the stratified patient bootstrap draws
#'
#' @param fitFun Function `(y, dose)` returning a doses x endpoints matrix.
#' @return List with `values` (long table) and `draws` (manifest).
#' @noRd
.cuiDrRun <- function(y, dose, doses, endpoints, nBoot, fitFun) {
  .epNames <- vapply(endpoints, `[[`, character(1), "name")
  .long <- function(p, draw) {
    data.frame(regimen = rep(.cuiNumLabel(doses), times = ncol(p)),
               endpoint = rep(.epNames, each = nrow(p)),
               draw = draw, value = as.vector(p), stringsAsFactors = FALSE)
  }
  .p0 <- fitFun(y, dose)
  .vals <- list(.long(.p0, 0L))
  .status <- rep("ok", nBoot + 1L)
  .msg <- rep("", nBoot + 1L)
  # resample whole patients within each dose: every endpoint uses the same rows
  .rowsByDose <- split(seq_along(dose), factor(dose, levels = doses))
  for (.b in seq_len(nBoot)) {
    .idx <- unlist(lapply(.rowsByDose, function(.r) .r[sample.int(length(.r), replace = TRUE)]),
                   use.names = FALSE)
    .p <- tryCatch(fitFun(y[.idx, , drop = FALSE], dose[.idx]),
                   error = function(e) conditionMessage(e))
    if (is.character(.p)) {
      .status[.b + 1L] <- "failed"
      .msg[.b + 1L] <- .p
    } else {
      .vals[[length(.vals) + 1L]] <- .long(.p, .b)
    }
  }
  list(values = do.call(rbind, .vals),
       draws = data.frame(draw = 0:nBoot, status = .status, source = "bootstrap",
                          message = .msg, stringsAsFactors = FALSE))
}

#' Clinical utility index from binary endpoints in a multiple-dose trial
#'
#' Estimates the probability of each binary endpoint at each dose level and
#' combines them into a clinical utility index (CUI), following CUI-MET
#' (Zhang et al. 2025).  Uncertainty comes from a bootstrap that resamples
#' whole patients within each dose (so the endpoints of a patient stay
#' together), giving percentile intervals and how often each dose has the
#' highest CUI.
#'
#' @param data Data frame with one row per patient.
#' @param endpoints List of [cuiBinary()] endpoints (columns of `data`
#'   holding 0/1 or `TRUE`/`FALSE`).
#' @param dose Numeric dose column, as a bare name or a string.
#' @param id Optional patient id column (bare name or string); ids must be
#'   unique, one row per patient.
#' @param nBoot Number of bootstrap draws (0 for none).
#' @param seed Seed for the bootstrap, or `NULL` to use the current random
#'   stream; the caller's stream is restored afterwards.
#' @param minSuccess Minimum fraction of bootstrap draws that must succeed.
#' @param .dose,.id Column names as strings, for programmatic use instead
#'   of `dose` and `id`.
#' @inheritParams cuiFromValues
#' @return An `nlmixr2Cui` object; `observed` holds the number of patients,
#'   events and observed rate per dose and endpoint.
#' @references Zhang F, Broglio K, Sweeting M, D'Angelo G (2025).  CUI-MET:
#'   Clinical Utility Index Dose Optimization Approach for Multiple-Dose,
#'   Multiple-Outcome Randomized Trial Designs.  arXiv:2505.03633.
#' @family Clinical utility
#' @export
#' @examples
#' set.seed(42)
#' p <- data.frame(dose = 1:5,
#'                 tox = c(0.038, 0.058, 0.106, 0.232, 0.533),
#'                 eff = c(0.022, 0.119, 0.343, 0.568, 0.681),
#'                 tol = c(0.147, 0.206, 0.280, 0.368, 0.466))
#' d <- p[rep(1:5, each = 30), ]
#' d$id <- seq_len(nrow(d))
#' for (v in c("tox", "eff", "tol")) d[[v]] <- rbinom(nrow(d), 1, d[[v]])
#' res <- cuiDoseResponse(d, dose = dose, id = id,
#'   endpoints = list(cuiBinary(tox, direction = "lower"),
#'                    cuiBinary(eff), cuiBinary(tol)),
#'   weights = c(tox = 0.2, eff = 0.5, tol = 0.3), nBoot = 200, seed = 1)
#' res
cuiDoseResponse <- function(data, endpoints, dose, id, weights = NULL,
                            combine = c("additive", "geometric"), eligible = NULL,
                            nBoot = 1000L, seed = NULL, minSuccess = 0.9, level = 0.95,
                            tieTol = sqrt(.Machine$double.eps), .dose = NULL, .id = NULL) {
  checkmate::assertDataFrame(data, min.rows = 1)
  .doseCol <- .cuiCaptureName(substitute(dose), missing(dose), .dose, "dose")
  .idCol <- NULL
  if (!missing(id) || !is.null(.id)) {
    .idCol <- .cuiCaptureName(substitute(id), missing(id), .id, "id")
  }
  checkmate::assertCount(nBoot)
  checkmate::assertNumber(minSuccess, lower = 0, upper = 1)
  if (inherits(endpoints, "cuiEndpoint")) endpoints <- list(endpoints)
  if (!all(vapply(endpoints, inherits, logical(1), what = "cuiBinary"))) {
    stop("'endpoints' must be a list of cuiBinary() endpoints", call. = FALSE)
  }
  .epNames <- vapply(endpoints, `[[`, character(1), "name")
  .missing <- setdiff(c(.doseCol, .idCol, .epNames), names(data))
  if (length(.missing) > 0) {
    stop(sprintf("column(s) not in 'data': %s; available: %s",
                 paste(.missing, collapse = ", "), paste(names(data), collapse = ", ")),
         call. = FALSE)
  }
  .doseVal <- data[[.doseCol]]
  if (!is.numeric(.doseVal) || anyNA(.doseVal)) {
    stop(sprintf("dose column '%s' must be numeric without missing values", .doseCol),
         call. = FALSE)
  }
  if (!is.null(.idCol) && anyDuplicated(data[[.idCol]]) != 0L) {
    stop(sprintf("id column '%s' must be unique (one row per patient)", .idCol), call. = FALSE)
  }
  .y <- vapply(.epNames, function(.e) {
    .v <- data[[.e]]
    if (is.logical(.v)) .v <- as.numeric(.v)
    if (!is.numeric(.v) || anyNA(.v) || !all(.v %in% c(0, 1))) {
      stop(sprintf("endpoint column '%s' must be 0/1 or TRUE/FALSE without missing values", .e),
           call. = FALSE)
    }
    as.numeric(.v)
  }, numeric(nrow(data)))
  .y <- matrix(.y, nrow = nrow(data), dimnames = list(NULL, .epNames))
  .doses <- sort(unique(.doseVal))
  .fitFun <- function(y, dose) .cuiDrFit(y, dose, .doses, endpoints)
  .run <- .cuiWithSeed(seed, .cuiDrRun(.y, .doseVal, .doses, endpoints, nBoot, .fitFun))
  .nOk <- sum(.run$draws$status[-1] == "ok")
  if (nBoot > 0 && .nOk / nBoot < minSuccess) {
    stop(sprintf("only %d of %d bootstrap draws succeeded (minSuccess = %g); first error: %s",
                 .nOk, nBoot, minSuccess,
                 .run$draws$message[.run$draws$status != "ok"][1]), call. = FALSE)
  }
  .models <- vapply(endpoints, `[[`, character(1), "model")
  .ret <- .cuiNew(values = .run$values, endpoints = endpoints,
                  regimens = cuiRegimens(dose = .doses), weights = weights,
                  combine = match.arg(combine), eligible = eligible, level = level,
                  tieTol = tieTol, draws = .run$draws,
                  settings = list(source = "doseResponse", nBoot = nBoot, seed = seed,
                                  models = setNames(.models, .epNames),
                                  description = sprintf("Dose-response: %s; patient bootstrap within dose",
                                                        paste(unique(.models), collapse = ", "))))
  .n <- as.vector(table(factor(.doseVal, levels = .doses)))
  .events <- rowsum(.y, factor(.doseVal, levels = .doses), reorder = TRUE)
  .ret$observed <- data.frame(regimen = rep(.cuiNumLabel(.doses), times = length(.epNames)),
                              endpoint = rep(.epNames, each = length(.doses)),
                              n = rep(.n, times = length(.epNames)),
                              events = as.vector(.events),
                              rate = as.vector(.events) / rep(.n, times = length(.epNames)),
                              stringsAsFactors = FALSE)
  .ret
}
