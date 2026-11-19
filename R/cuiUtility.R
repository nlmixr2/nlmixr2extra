#' Utility functions for a clinical utility index endpoint
#'
#' A utility maps an endpoint summary (a probability, a mean, ...) to a
#' desirability between 0 and 1, where 1 is best.  `direction = "lower"`
#' in [cuiEndpoint()] is applied afterwards as `1 - u`.
#'
#' - `cuiIdentity()`: the value is already a utility; it must lie in
#'   `[0, 1]` (for example a probability).
#' - `cuiLinear(lo, hi)`: `clamp((x - lo) / (hi - lo), 0, 1)`.
#' - `cuiStep(cut, side)`: 1 when `x >= cut` (`side = "above"`) or
#'   `x <= cut` (`side = "below"`), otherwise 0.
#'
#' @param lo,hi Values mapped to utility 0 and 1; `lo < hi`.
#' @param cut Threshold of the step.
#' @param side Which side of `cut` has utility 1.
#' @return A `cuiUtility` object.
#' @family Clinical utility
#' @export
#' @examples
#' cuiLinear(0, 50)
#' cuiStep(-1, "above")
cuiIdentity <- function() {
  .cuiUtilityNew("identity")
}

#' @rdname cuiIdentity
#' @export
cuiLinear <- function(lo, hi) {
  checkmate::assertNumber(lo, finite = TRUE)
  checkmate::assertNumber(hi, finite = TRUE)
  if (lo >= hi) {
    stop("'lo' must be less than 'hi'", call. = FALSE)
  }
  .cuiUtilityNew("linear", lo = lo, hi = hi)
}

#' @rdname cuiIdentity
#' @export
cuiStep <- function(cut, side = c("above", "below")) {
  checkmate::assertNumber(cut, finite = TRUE)
  side <- match.arg(side)
  .cuiUtilityNew("step", cut = cut, side = side)
}

.cuiUtilityNew <- function(type, ...) {
  structure(list(type = type, ...), class = "cuiUtility")
}

#' @export
format.cuiUtility <- function(x, ...) {
  switch(x$type,
    identity = "identity",
    linear = sprintf("linear(%s, %s)", format(x$lo), format(x$hi)),
    step = sprintf("step(%s %s)", ifelse(x$side == "above", ">=", "<="), format(x$cut))
  )
}

#' @export
print.cuiUtility <- function(x, ...) {
  cat("cuiUtility: ", format(x), "\n", sep = "")
  invisible(x)
}

#' Apply a utility and endpoint direction to endpoint values
#'
#' @param utility A `cuiUtility`.
#' @param direction `"higher"` or `"lower"`.
#' @param x Numeric endpoint values.
#' @param endpoint Endpoint name, for error messages.
#' @return Utilities in `[0, 1]`.
#' @noRd
.cuiApplyUtility <- function(utility, direction, x, endpoint) {
  .u <- switch(utility$type,
    identity = {
      if (any(x < 0 | x > 1)) {
        stop(sprintf(
          "endpoint '%s' uses cuiIdentity() but has values outside [0, 1]; use cuiLinear() or cuiStep()",
          endpoint
        ), call. = FALSE)
      }
      x
    },
    linear = pmin(pmax((x - utility$lo) / (utility$hi - utility$lo), 0), 1),
    step = {
      if (utility$side == "above") {
        as.numeric(x >= utility$cut)
      } else {
        as.numeric(x <= utility$cut)
      }
    },
    stop(sprintf("unknown utility type '%s'", utility$type), call. = FALSE)
  )
  if (direction == "lower") .u <- 1 - .u
  .u
}

#' Normalize clinical utility weights
#'
#' @param weights Named nonnegative weights, or `NULL` for equal weights.
#' @param endpoints Endpoint names.
#' @return Normalized weights in `endpoints` order.
#' @noRd
.cuiNormWeights <- function(weights, endpoints) {
  if (is.null(weights)) {
    return(setNames(rep(1 / length(endpoints), length(endpoints)), endpoints))
  }
  checkmate::assertNumeric(weights, lower = 0, finite = TRUE, any.missing = FALSE, min.len = 1)
  if (is.null(names(weights)) || !setequal(names(weights), endpoints) ||
    anyDuplicated(names(weights)) != 0L) {
    stop(sprintf(
      "'weights' must be named with exactly the endpoints: %s",
      paste(endpoints, collapse = ", ")
    ), call. = FALSE)
  }
  .w <- weights[endpoints]
  if (sum(.w) <= 0) {
    stop("'weights' must have a positive total", call. = FALSE)
  }
  .w / sum(.w)
}

#' Combine a utility matrix into a clinical utility index
#'
#' @param u Matrix of utilities, one column per endpoint.
#' @param w Normalized weights matching the columns.
#' @param combine `"additive"` or `"geometric"`.
#' @return CUI per row.
#' @noRd
.cuiCombine <- function(u, w, combine) {
  if (combine == "additive") {
    return(drop(u %*% w))
  }
  # geometric: only positive-weight endpoints take part; a zero utility gives 0
  .pos <- w > 0
  .u <- u[, .pos, drop = FALSE]
  .w <- w[.pos] / sum(w[.pos])
  .zero <- rowSums(.u == 0) > 0
  .ret <- numeric(nrow(.u))
  if (any(!.zero)) {
    .ret[!.zero] <- exp(drop(log(.u[!.zero, , drop = FALSE]) %*% .w))
  }
  .ret
}
