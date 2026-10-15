# Categorical slots in fixed order (validated: adjacent CVD and normal-vision
# separation pass); the CUI itself is drawn in neutral ink
.cuiSeriesColors <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
                      "#6250d6", "#e34948")
.cuiInk <- "#0b0b0b"
.cuiMuted <- "#9b9a95"
# ordinal blue ramp: discrete tiles start at step 250 so the lightest stays visible
.cuiSeqRamp <- c("#86b6ef", "#3987e5", "#1c5cab", "#0d366b")

#' Regimen positions for plotting: dose when every regimen has a distinct
#' dose, otherwise the regimen order
#' @noRd
.cuiRegimenX <- function(regimens) {
  .d <- regimens$dose
  if (!anyNA(.d) && anyDuplicated(.d) == 0L) {
    return(list(numeric = TRUE, x = setNames(.d, regimens$label), lab = "Dose"))
  }
  list(numeric = FALSE, x = setNames(seq_len(nrow(regimens)), regimens$label), lab = "Regimen")
}

.cuiXScale <- function(rx) {
  if (rx$numeric) return(ggplot2::scale_x_continuous(name = rx$lab))
  ggplot2::scale_x_continuous(name = rx$lab, breaks = unname(rx$x), labels = names(rx$x))
}

.cuiFacet <- function(n, var = "assessment") {
  if (n > 1) ggplot2::facet_wrap(var)
}

.cuiPickAssessment <- function(x, assessment) {
  if (is.null(assessment)) return(x$assessments)
  .bad <- setdiff(assessment, x$assessments)
  if (length(.bad) > 0) {
    stop(sprintf("unknown assessment(s): %s", paste(.bad, collapse = ", ")), call. = FALSE)
  }
  assessment
}

#' Plot a clinical utility index
#'
#' - `"cui"`: the CUI of each regimen (point estimate, with its percentile
#'   interval when there are uncertainty draws); the recommended regimen is
#'   ringed and ineligible regimens are grey.
#' - `"endpoints"`: each endpoint's utility (after its utility function and
#'   direction) and the CUI, at the point estimate.
#' - `"selection"`: how often each eligible regimen has the highest CUI
#'   across uncertainty draws.
#' - `"doseResponse"`: for [cuiDoseResponse()] results, the observed rate of
#'   each endpoint per dose with the fitted model and its bootstrap
#'   percentile interval.
#'
#' Regimens are placed by dose when every regimen has a distinct dose,
#' otherwise in their order.  Several assessments are shown as panels.
#'
#' @param x An `nlmixr2Cui` object.
#' @param type Which plot, as described above.
#' @param assessment Assessments to show; `NULL` shows all.
#' @param ... Not used.
#' @return A ggplot object.
#' @family Clinical utility
#' @export
#' @examples
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
#' plot(res)
#' plot(res, "endpoints")
plot.nlmixr2Cui <- function(x, type = c("cui", "endpoints", "selection", "doseResponse"),
                            assessment = NULL, ...) {
  type <- match.arg(type)
  .a <- .cuiPickAssessment(x, assessment)
  switch(type,
    cui = .cuiPlotCui(x, .a),
    endpoints = .cuiPlotEndpoints(x, .a),
    selection = .cuiPlotSelection(x, .a),
    doseResponse = .cuiPlotDoseResponse(x)
  )
}

.cuiNOk <- function(x) sum(x$draws$status == "ok" & x$draws$draw > 0)

.cuiPlotCui <- function(x, assessments) {
  .rx <- .cuiRegimenX(x$regimens)
  .d <- x$summary[x$summary$assessment %in% assessments, , drop = FALSE]
  .d$x <- unname(.rx$x[.d$regimen])
  .d$assessment <- factor(.d$assessment, levels = assessments)
  .best <- vapply(assessments, function(a) {
    paste(.d$regimen[.d$assessment == a & .d$best], collapse = ", ")
  }, character(1))
  .sub <- paste0("Recommended: ", paste(if (length(assessments) > 1) {
    paste0(assessments, " ", .best)
  } else {
    .best
  }, collapse = "; "))
  .hasCi <- .cuiNOk(x) > 0
  if (.hasCi) .sub <- paste0(.sub, "; band: ", 100 * x$level, "% percentile interval")
  .plotData(.cuiCuiFigure(.rx, .hasCi, any(!.d$eligible), .sub, length(assessments)), .d)
}

.cuiCuiFigure <- function(rx, hasCi, anyIneligible, subtitle, nAssess) {
  # force every argument: an unevaluated promise would keep the caller's frame
  invisible(list(rx, hasCi, anyIneligible, subtitle, nAssess))
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x, y = .data$est, group = .data$assessment)) +
    (if (hasCi) {
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                           fill = .cuiSeriesColors[1], alpha = 0.18)
    }) +
    ggplot2::geom_line(color = .cuiSeriesColors[1], linewidth = 0.7) +
    ggplot2::geom_point(data = function(d) d[d$eligible, , drop = FALSE],
                        color = .cuiSeriesColors[1], size = 2.5) +
    (if (anyIneligible) {
      ggplot2::geom_point(data = function(d) d[!d$eligible, , drop = FALSE],
                          color = .cuiMuted, size = 2.5)
    }) +
    ggplot2::geom_point(data = function(d) d[d$best, , drop = FALSE], shape = 21, size = 5,
                        stroke = 0.9, color = .cuiInk, fill = NA) +
    .cuiXScale(rx) +
    ggplot2::ylab("Clinical utility index") +
    ggplot2::labs(title = "Clinical utility index", subtitle = subtitle) +
    (.cuiFacet(nAssess)) +
    rxode2::rxTheme()
}

.cuiPlotEndpoints <- function(x, assessments) {
  .rx <- .cuiRegimenX(x$regimens)
  .epNames <- vapply(x$endpoints, `[[`, character(1), "name")
  .u <- x$utility[x$utility$draw == 0 & x$utility$assessment %in% assessments, , drop = FALSE]
  .c <- x$summary[x$summary$assessment %in% assessments, , drop = FALSE]
  .d <- rbind(data.frame(assessment = .u$assessment, regimen = .u$regimen,
                         series = .u$endpoint, value = .u$utility, stringsAsFactors = FALSE),
              data.frame(assessment = .c$assessment, regimen = .c$regimen,
                         series = "CUI", value = .c$est, stringsAsFactors = FALSE))
  .d$x <- unname(.rx$x[.d$regimen])
  .d$series <- factor(.d$series, levels = c(.epNames, "CUI"))
  .d$assessment <- factor(.d$assessment, levels = assessments)
  # label each line at its last regimen, beside the legend, so identity is
  # never color alone
  .d$last <- .d$x == max(.d$x)
  .d$labY <- .d$value
  for (.a in assessments) {
    .i <- which(.d$last & .d$assessment == .a)
    .d$labY[.i] <- .cuiSpreadLabels(.d$value[.i], 0.045)
  }
  .cols <- c(setNames(rep_len(.cuiSeriesColors, length(.epNames)), .epNames), CUI = .cuiInk)
  .sub <- paste0("Weights: ", paste(sprintf("%s %.2f", .epNames, x$weightsNorm), collapse = ", "))
  .plotData(.cuiEndpointsFigure(.rx, .cols, .sub, length(assessments)), .d)
}

.cuiEndpointsFigure <- function(rx, cols, subtitle, nAssess) {
  invisible(list(rx, cols, subtitle, nAssess))
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x, y = .data$value, color = .data$series,
                                         group = .data$series)) +
    ggplot2::geom_line(ggplot2::aes(linewidth = .data$series == "CUI")) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(data = function(d) d[d$last, , drop = FALSE],
                       ggplot2::aes(y = .data$labY, label = .data$series), hjust = -0.2,
                       size = 3.2, show.legend = FALSE) +
    ggplot2::scale_color_manual(name = NULL, values = cols) +
    ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.1), guide = "none") +
    .cuiXScale(rx) +
    ggplot2::scale_y_continuous(name = "Utility (higher is better)", limits = c(0, 1)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = "Endpoint utilities and clinical utility index", subtitle = subtitle) +
    (.cuiFacet(nAssess)) +
    rxode2::rxTheme() +
    ggplot2::theme(legend.position = "bottom", plot.margin = ggplot2::margin(5, 40, 5, 5))
}

#' Spread label positions so neighbors are at least `gap` apart
#' @noRd
.cuiSpreadLabels <- function(y, gap) {
  .o <- order(y)
  .s <- y[.o]
  for (.i in seq_along(.s)[-1]) .s[.i] <- max(.s[.i], .s[.i - 1] + gap)
  # center the shift so labels stay near their lines
  .s <- .s - (mean(.s) - mean(y[.o]))
  .ret <- y
  .ret[.o] <- .s
  .ret
}

.cuiPlotSelection <- function(x, assessments) {
  .n <- .cuiNOk(x)
  if (.n == 0) stop("selection frequencies need uncertainty draws", call. = FALSE)
  .d <- x$summary[x$summary$assessment %in% assessments & x$summary$eligible, , drop = FALSE]
  .d$regimen <- factor(.d$regimen, levels = x$regimens$label)
  .d$assessment <- factor(.d$assessment, levels = assessments)
  .sub <- sprintf("Share of %d draws in which each regimen has the highest CUI", .n)
  .plotData(.cuiSelectionFigure(.sub, length(assessments),
                                if (.cuiRegimenX(x$regimens)$numeric) "Dose" else "Regimen"), .d)
}

.cuiSelectionFigure <- function(subtitle, nAssess, xlab) {
  invisible(list(subtitle, nAssess, xlab))
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$regimen, y = .data$selFreq)) +
    ggplot2::geom_col(fill = .cuiSeriesColors[1], width = 0.6) +
    ggplot2::geom_text(data = function(d) d[d$selFreq > 0, , drop = FALSE],
                       ggplot2::aes(label = sprintf("%.0f%%", 100 * .data$selFreq)),
                       vjust = -0.4, size = 3.2, color = .cuiInk) +
    ggplot2::scale_y_continuous(name = "Selection frequency", limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::xlab(xlab) +
    ggplot2::labs(title = "How often each regimen is selected", subtitle = subtitle) +
    (.cuiFacet(nAssess)) +
    rxode2::rxTheme()
}

.cuiPlotDoseResponse <- function(x) {
  if (is.null(x$observed)) {
    stop("dose-response plots need a cuiDoseResponse() result", call. = FALSE)
  }
  .doses <- x$regimens$dose
  .epNames <- vapply(x$endpoints, `[[`, character(1), "name")
  .obs <- x$observed
  .obs$dose <- .doses[match(.obs$regimen, x$regimens$label)]
  .obs$kind <- "observed"
  .obs$value <- .obs$rate
  # fitted curves: smooth for parametric models, through the doses otherwise
  .grid <- seq(min(.doses), max(.doses), length.out = 101)
  .fit <- do.call(rbind, lapply(.epNames, function(e) {
    .f <- x$doseResponse[[e]]
    .dd <- if (.f$model %in% c("empirical", "fitted")) .doses else .grid
    data.frame(endpoint = e, dose = .dd, value = .cuiDrPredict(.f, .dd), stringsAsFactors = FALSE)
  }))
  .ci <- NULL
  if (.cuiNOk(x) > 0) {
    .v <- x$values[x$values$draw > 0, , drop = FALSE]
    .alpha <- (1 - x$level) / 2
    .ci <- stats::aggregate(value ~ endpoint + regimen, data = .v,
                            FUN = function(v) stats::quantile(v, c(.alpha, 1 - .alpha)))
    .ci <- data.frame(endpoint = .ci$endpoint, regimen = .ci$regimen,
                      lower = .ci$value[, 1], upper = .ci$value[, 2], stringsAsFactors = FALSE)
    .ci$dose <- .doses[match(.ci$regimen, x$regimens$label)]
  }
  .d <- rbind(data.frame(layer = "observed", endpoint = .obs$endpoint, dose = .obs$dose,
                         value = .obs$value, lower = NA_real_, upper = NA_real_),
              data.frame(layer = "fit", endpoint = .fit$endpoint, dose = .fit$dose,
                         value = .fit$value, lower = NA_real_, upper = NA_real_),
              if (!is.null(.ci)) {
                data.frame(layer = "interval", endpoint = .ci$endpoint, dose = .ci$dose,
                           value = NA_real_, lower = .ci$lower, upper = .ci$upper)
              })
  .d$endpoint <- factor(.d$endpoint, levels = .epNames)
  .sub <- paste0("Points: observed rate; line: fitted model",
                 if (!is.null(.ci)) sprintf("; bars: %g%% bootstrap interval", 100 * x$level) else "")
  .plotData(.cuiDoseResponseFigure(.sub, !is.null(.ci)), .d)
}

.cuiDoseResponseFigure <- function(subtitle, hasCi) {
  invisible(list(subtitle, hasCi))
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$dose)) +
    (if (hasCi) {
      ggplot2::geom_linerange(data = function(d) d[d$layer == "interval", , drop = FALSE],
                              ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                              color = .cuiSeriesColors[1], alpha = 0.5, linewidth = 0.8)
    }) +
    ggplot2::geom_line(data = function(d) d[d$layer == "fit", , drop = FALSE],
                       ggplot2::aes(y = .data$value), color = .cuiSeriesColors[1], linewidth = 0.7) +
    ggplot2::geom_point(data = function(d) d[d$layer == "observed", , drop = FALSE],
                        ggplot2::aes(y = .data$value), color = .cuiInk, size = 2) +
    ggplot2::facet_wrap("endpoint") +
    ggplot2::scale_y_continuous(name = "Probability of the event", limits = c(0, 1)) +
    ggplot2::xlab("Dose") +
    ggplot2::labs(title = "Dose-response of each endpoint", subtitle = subtitle) +
    rxode2::rxTheme()
}

#' Sensitivity of the recommended regimen to the endpoint weights
#'
#' Varies the weight of one endpoint at a time over `weights` (the other
#' weights stay as they are) and records which eligible regimen has the
#' highest point-estimate CUI.  This shows how far a preference has to move
#' before the recommendation changes.
#'
#' @param x An `nlmixr2Cui` object.
#' @param weights Weights to try for each endpoint, on the scale of the
#'   weights given to the CUI (unnormalized; equal weights count as 1).
#' @param assessment Assessment to analyze; by default the first.
#' @return An `nlmixr2CuiSensitivity` data frame with columns `endpoint`,
#'   `weight`, `regimen` (the recommendation) and `cui`, with `print()` and
#'   `plot()` methods.
#' @family Clinical utility
#' @export
#' @examples
#' values <- data.frame(regimen = rep(1:5, 3),
#'                      endpoint = rep(c("tox", "eff", "tol"), each = 5),
#'                      value = c(0.038, 0.058, 0.106, 0.232, 0.533,
#'                                0.022, 0.119, 0.343, 0.568, 0.681,
#'                                0.147, 0.206, 0.280, 0.368, 0.466))
#' res <- cuiFromValues(values,
#'   endpoints = list(cuiEndpoint(tox, direction = "lower"),
#'                    cuiEndpoint(eff), cuiEndpoint(tol)),
#'   weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
#' s <- cuiSensitivity(res)
#' s
#' plot(s)
cuiSensitivity <- function(x, weights = seq(0, 5, by = 0.25), assessment = NULL) {
  if (!inherits(x, "nlmixr2Cui")) stop("'x' must be an nlmixr2Cui object", call. = FALSE)
  checkmate::assertNumeric(weights, lower = 0, finite = TRUE, any.missing = FALSE, min.len = 2)
  .a <- if (is.null(assessment)) x$assessments[1] else .cuiPickAssessment(x, assessment)
  if (length(.a) != 1L) stop("give a single assessment", call. = FALSE)
  .epNames <- vapply(x$endpoints, `[[`, character(1), "name")
  .w0 <- if (is.null(x$weights)) setNames(rep(1, length(.epNames)), .epNames) else x$weights[.epNames]
  .el <- x$eligible$regimen[x$eligible$assessment == .a & x$eligible$eligible]
  .el <- x$regimens$label[x$regimens$label %in% .el]
  .u <- x$utility[x$utility$draw == 0 & x$utility$assessment == .a, , drop = FALSE]
  .m <- matrix(.u$utility[match(paste(rep(.el, each = length(.epNames)), .epNames, sep = "\r"),
                                paste(.u$regimen, .u$endpoint, sep = "\r"))],
               ncol = length(.epNames), byrow = TRUE, dimnames = list(.el, .epNames))
  .ret <- do.call(rbind, lapply(.epNames, function(e) {
    do.call(rbind, lapply(weights, function(w) {
      .w <- .w0
      .w[[e]] <- w
      if (sum(.w) <= 0) {
        return(data.frame(endpoint = e, weight = w, regimen = NA_character_, cui = NA_real_,
                          stringsAsFactors = FALSE))
      }
      .c <- .cuiCombine(.m, .w / sum(.w), x$combine)
      .i <- which(.c >= max(.c) - x$tieTol)[1]
      data.frame(endpoint = e, weight = w, regimen = .el[.i], cui = .c[.i],
                 stringsAsFactors = FALSE)
    }))
  }))
  attr(.ret, "current") <- .w0
  attr(.ret, "assessment") <- .a
  attr(.ret, "regimens") <- x$regimens
  class(.ret) <- c("nlmixr2CuiSensitivity", "data.frame")
  .ret
}

#' @export
print.nlmixr2CuiSensitivity <- function(x, ...) {
  .cur <- attr(x, "current")
  cat(sprintf("Recommended regimen as one endpoint weight varies (assessment %s)\n",
              attr(x, "assessment")))
  for (.e in unique(x$endpoint)) {
    .s <- x[x$endpoint == .e, , drop = FALSE]
    .run <- rle(ifelse(is.na(.s$regimen), "none", .s$regimen))
    .end <- cumsum(.run$lengths)
    .start <- c(1, .end[-length(.end)] + 1)
    .rng <- ifelse(.start == .end, .cuiNumLabel(.s$weight[.start]),
                   paste0(.cuiNumLabel(.s$weight[.start]), "-", .cuiNumLabel(.s$weight[.end])))
    .txt <- sprintf("%s: %s", .rng, .run$values)
    cat(sprintf("  %s (now %s): %s\n", .e, .cuiNumLabel(.cur[[.e]]), paste(.txt, collapse = "; ")))
  }
  invisible(x)
}

#' @rdname cuiSensitivity
#' @param ... Not used.
#' @export
plot.nlmixr2CuiSensitivity <- function(x, ...) {
  .reg <- attr(x, "regimens")
  .cur <- attr(x, "current")
  .d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  .d$endpoint <- factor(.d$endpoint, levels = rev(unique(.d$endpoint)))
  .d$current <- abs(.d$weight - .cur[as.character(.d$endpoint)]) < 1e-12
  .numeric <- .cuiRegimenX(.reg)$numeric
  if (.numeric) .d$dose <- .reg$dose[match(.d$regimen, .reg$label)]
  .step <- min(diff(sort(unique(.d$weight))))
  .plotData(.cuiSensitivityFigure(.numeric, attr(x, "assessment"), .step,
                                  if (.numeric) NULL else sort(unique(stats::na.omit(.d$regimen)))),
            .d)
}

.cuiSensitivityFigure <- function(numeric, assessment, step, levels) {
  invisible(list(numeric, assessment, step, levels))
  .fill <- if (numeric) {
    ggplot2::scale_fill_gradientn(name = "Recommended dose", colours = .cuiSeqRamp)
  } else {
    ggplot2::scale_fill_manual(name = "Recommended regimen",
                               values = setNames(rep_len(.cuiSeriesColors, length(levels)), levels))
  }
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$weight, y = .data$endpoint)) +
    ggplot2::geom_tile(ggplot2::aes(fill = if (numeric) .data$dose else .data$regimen),
                       width = step, height = 0.8, color = "#fcfcfb", linewidth = 0.4) +
    ggplot2::geom_point(data = function(d) d[d$current, , drop = FALSE], shape = 21, size = 3.5,
                        stroke = 1, color = .cuiInk, fill = "#fcfcfb") +
    .fill +
    ggplot2::xlab("Weight of the endpoint (others unchanged)") +
    ggplot2::ylab(NULL) +
    ggplot2::labs(title = "Sensitivity of the recommendation to the weights",
                  subtitle = sprintf("Assessment %s; ring: current weight", assessment)) +
    rxode2::rxTheme() +
    ggplot2::theme(legend.position = "bottom")
}
