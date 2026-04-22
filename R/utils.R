#' Apply FUN over X with optional parallel execution and progress reporting
#'
#' When \pkg{future.apply} is available, uses \code{future_lapply()} under the
#' current \code{future::plan()}.  If \pkg{progressr} is also available,
#' progress is reported.  Otherwise falls back to \code{base::lapply()}.
#'
#' @param X       vector or list to iterate over
#' @param FUN     function applied to each element of \code{X}
#' @param ...     additional arguments passed to \code{FUN}
#' @param .label  optional \code{function(x) -> character} producing a per-item
#'   progress label; \code{x} is each element of \code{X}
#' @return list of results in the same order as \code{X}
#' @noRd
.plap <- function(X, FUN, ..., .label = NULL) {
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    return(lapply(X, FUN, ...))
  }

  if (!requireNamespace("progressr", quietly = TRUE)) {
    return(future.apply::future_lapply(
      X,
      FUN,
      ...,
      future.seed = TRUE,
      future.packages = "nlmixr2extra"
    ))
  }

  progressr::with_progress({
    p <- progressr::progressor(steps = length(X))
    future.apply::future_lapply(
      X,
      function(x, ...) {
        on.exit(p(message = if (!is.null(.label)) .label(x) else ""))
        FUN(x, ...)
      },
      ...,
      future.seed = TRUE,
      future.packages = "nlmixr2extra"
    )
  })
}

.validateWorkers <- function(workers) {
  if (is.null(workers) || identical(workers, "auto")) {
    return(invisible(NULL))
  }
  if (
    !is.numeric(workers) ||
      length(workers) != 1L ||
      is.na(workers) ||
      !is.finite(workers) ||
      workers < 1 ||
      workers != as.integer(workers)
  ) {
    cli::cli_abort(
      "{.arg workers} must be NULL, \"auto\", 1, or a positive integer."
    )
  }
  invisible(NULL)
}

#' Temporarily set a future parallel plan for the duration of an expression
#'
#' @param workers \code{NULL} (leave the current plan unchanged),
#'   \code{1} (force sequential), a positive integer (use that many
#'   \code{multisession} workers), or \code{"auto"} (use
#'   \code{future::availableCores(omit = 1)}).
#' @param expr expression to evaluate; the prior plan is always restored on
#'   exit, even if \code{expr} throws an error.
#' @return value of \code{expr}
#' @noRd
.withWorkerPlan <- function(workers, expr) {
  .validateWorkers(workers)
  if (is.null(workers)) {
    return(force(expr))
  }
  if (!requireNamespace("future", quietly = TRUE)) {
    cli::cli_warn(c(
      "!" = "Package {.pkg future} is not installed.",
      "i" = "Ignoring {.arg workers} and running sequentially."
    ))
    return(force(expr))
  }
  if (identical(workers, "auto")) {
    workers <- max(1L, as.integer(future::availableCores(omit = 1L)))
  } else {
    workers <- as.integer(workers)
  }
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  if (workers == 1L) {
    future::plan("sequential")
  } else {
    future::plan("multisession", workers = workers)
  }
  force(expr)
}

#' Make a control step that is quieter and faster
#'
#' @param ctl the control object
#' @return A faster and quieter control object
#' @noRd
setQuietFastControl <- function(ctl) {
  # make estimation steps quieter
  ctl$print <- 0L
  # make estimation steps faster
  ctl$covMethod <- 0L
  ctl$calcTables <- FALSE
  ctl$compress <- FALSE
  ctl
}
