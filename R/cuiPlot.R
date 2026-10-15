# Visual grammar shared with the ACoP 2026 CUI workshop app: the CUI is a
# thick navy line, its peak a rust point labeled "max at ...", benefit
# endpoints are solid and safety endpoints dashed or dotted, and reference
# regimens are dashed vertical lines.  Endpoint hues are the validated
# categorical slots in fixed order (the linetypes are a second encoding);
# `colors` in plot() can replace them.
.cuiSeriesColors <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
                      "#6250d6", "#e34948")
.cuiNavy <- "#0B1F3A"
.cuiPeakCol <- "#B4532A"
.cuiMuted <- "#9b9a95"
.cuiSafetyLty <- c("dashed", "dotted", "dotdash", "longdash", "twodash")
# ordinal ramp ending in the CUI navy: discrete tiles start light but visible
.cuiSeqRamp <- c("#86b6ef", "#3987e5", "#1c5cab", "#0B1F3A")

.cuiTheme <- function() {
  ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", color = .cuiNavy),
      plot.subtitle = ggplot2::element_text(color = "#4A4A4A", size = 11),
      legend.position = "bottom"
    )
}

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
#' - `"cui"`: the CUI of each regimen (point estimate line, with its
#'   percentile band when there are uncertainty draws); the peak
#'   (recommended regimen) is marked and ineligible regimens are grey.
#' - `"endpoints"`: each endpoint's utility (after its utility function and
#'   direction) with the CUI, as point-estimate lines with percentile bands
#'   from the uncertainty draws; endpoints with `direction = "lower"` are
#'   dashed or dotted.
#' - `"selection"`: how often each eligible regimen has the highest CUI
#'   across uncertainty draws.
#' - `"doseResponse"`: for [cuiDoseResponse()] results, the observed rate of
#'   each endpoint per dose with the fitted model and its bootstrap
#'   percentile interval.
#'
#' Regimens are placed by dose when every regimen has a distinct dose,
#' otherwise in their order.  Several assessments are shown as panels.  The
#' figures follow the style of the ACoP 2026 CUI workshop app.
#'
#' @param x An `nlmixr2Cui` object.
#' @param type Which plot, as described above.
#' @param assessment Assessments to show; `NULL` shows all.
#' @param reference Regimens marked with dashed vertical lines, as regimen
#'   labels or a named numeric vector of positions (names are the line
#'   labels, for example `c("150 QD" = 150, "200 BID" = 400)`).  By default
#'   the eligible regimens when eligibility is restricted.  `character(0)`
#'   draws none.
#' @param colors Optional named colors replacing the defaults, by endpoint
#'   name and `"CUI"` / `"Peak"`.
#' @param bands Draw percentile bands from the uncertainty draws (when
#'   there are any) around the lines.
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
#' plot(res, "endpoints", reference = c("low" = 2, "high" = 4))
plot.nlmixr2Cui <- function(x, type = c("cui", "endpoints", "selection", "doseResponse"),
                            assessment = NULL, reference = NULL, colors = NULL, bands = TRUE,
                            ...) {
  type <- match.arg(type)
  checkmate::assertFlag(bands)
  .a <- .cuiPickAssessment(x, assessment)
  .cols <- .cuiPlotColors(x, colors)
  switch(type,
    cui = .cuiPlotLines(x, .a, reference, .cols, endpoints = FALSE, bands = bands),
    endpoints = .cuiPlotLines(x, .a, reference, .cols, endpoints = TRUE, bands = bands),
    selection = .cuiPlotSelection(x, .a),
    doseResponse = .cuiPlotDoseResponse(x, .cols, bands)
  )
}

.cuiNOk <- function(x) sum(x$draws$status == "ok" & x$draws$draw > 0)

#' Endpoint, CUI and peak colors, with user overrides
#' @noRd
.cuiPlotColors <- function(x, colors) {
  .epNames <- vapply(x$endpoints, `[[`, character(1), "name")
  .cols <- c(setNames(rep_len(.cuiSeriesColors, length(.epNames)), .epNames),
             CUI = .cuiNavy, Peak = .cuiPeakCol)
  if (!is.null(colors)) {
    checkmate::assertCharacter(colors, names = "unique", any.missing = FALSE)
    .bad <- setdiff(names(colors), names(.cols))
    if (length(.bad) > 0) {
      stop(sprintf("unknown name(s) in 'colors': %s", paste(.bad, collapse = ", ")), call. = FALSE)
    }
    .cols[names(colors)] <- colors
  }
  .cols
}

#' Reference lines: positions and labels on the plotting axis
#' @noRd
.cuiReference <- function(x, reference, rx, assessments) {
  if (is.null(reference)) {
    .el <- x$eligible[x$eligible$assessment %in% assessments, , drop = FALSE]
    if (all(.el$eligible)) return(NULL)
    reference <- unique(.el$regimen[.el$eligible])
  }
  if (length(reference) == 0) return(NULL)
  if (is.numeric(reference)) {
    if (is.null(names(reference))) names(reference) <- .cuiNumLabel(reference)
    return(data.frame(x = unname(reference), label = names(reference), stringsAsFactors = FALSE))
  }
  .bad <- setdiff(reference, names(rx$x))
  if (length(.bad) > 0) {
    stop(sprintf("unknown regimen(s) in 'reference': %s", paste(.bad, collapse = ", ")),
         call. = FALSE)
  }
  data.frame(x = unname(rx$x[reference]), label = reference, stringsAsFactors = FALSE)
}

#' CUI (and optionally endpoint) lines with the peak and reference regimens
#' @noRd
.cuiPlotLines <- function(x, assessments, reference, cols, endpoints, bands = TRUE) {
  .rx <- .cuiRegimenX(x$regimens)
  .epNames <- vapply(x$endpoints, `[[`, character(1), "name")
  .dir <- vapply(x$endpoints, `[[`, character(1), "direction")
  .s <- x$summary[x$summary$assessment %in% assessments, , drop = FALSE]
  .hasCi <- bands && .cuiNOk(x) > 0
  .line <- data.frame(layer = "line", assessment = .s$assessment, regimen = .s$regimen,
                      series = "CUI", value = .s$est,
                      lower = if (.hasCi) .s$lower else NA_real_,
                      upper = if (.hasCi) .s$upper else NA_real_,
                      eligible = .s$eligible, label = NA_character_, stringsAsFactors = FALSE)
  if (endpoints) {
    .u <- x$utility[x$utility$draw == 0 & x$utility$assessment %in% assessments, , drop = FALSE]
    .ep <- data.frame(layer = "line", assessment = .u$assessment, regimen = .u$regimen,
                      series = .u$endpoint, value = .u$utility, lower = NA_real_,
                      upper = NA_real_, eligible = TRUE, label = NA_character_,
                      stringsAsFactors = FALSE)
    if (.hasCi) {
      # percentile band of each endpoint's utility over the uncertainty draws
      .ud <- x$utility[x$utility$draw > 0 & x$utility$assessment %in% assessments, , drop = FALSE]
      .k <- paste(.ud$assessment, .ud$regimen, .ud$endpoint, sep = "\r")
      .alpha <- (1 - x$level) / 2
      .lo <- tapply(.ud$utility, .k, stats::quantile, probs = .alpha, names = FALSE)
      .hi <- tapply(.ud$utility, .k, stats::quantile, probs = 1 - .alpha, names = FALSE)
      .ek <- paste(.ep$assessment, .ep$regimen, .ep$series, sep = "\r")
      .ep$lower <- as.numeric(.lo[.ek])
      .ep$upper <- as.numeric(.hi[.ek])
    }
    .line <- rbind(.ep, .line)
  }
  .units <- x$regimens$units[match(.s$regimen, x$regimens$label)]
  .pk <- .s[.s$best, , drop = FALSE]
  .pkUnits <- .units[.s$best]
  .peak <- data.frame(layer = "peak", assessment = .pk$assessment, regimen = .pk$regimen,
                      series = "Peak", value = .pk$est, lower = NA_real_, upper = NA_real_,
                      eligible = TRUE,
                      label = paste0("max at ", .pk$regimen,
                                     ifelse(is.na(.pkUnits), "", paste0(" ", .pkUnits))),
                      stringsAsFactors = FALSE)
  .d <- rbind(.line, .peak)
  .d$x <- unname(.rx$x[.d$regimen])
  # default reference lines are each panel's own eligible regimens
  .refs <- lapply(assessments, function(a) {
    .r <- .cuiReference(x, reference, .rx, a)
    if (is.null(.r)) return(NULL)
    data.frame(layer = "ref", assessment = a, regimen = NA_character_, series = NA_character_,
               value = NA_real_, lower = NA_real_, upper = NA_real_, eligible = TRUE,
               label = .r$label, x = .r$x, stringsAsFactors = FALSE)
  })
  .ref <- do.call(rbind, .refs)
  if (!is.null(.ref)) .d <- rbind(.d, .ref)
  .lv <- c(if (endpoints) .epNames, "CUI", "Peak")
  .d$series <- factor(.d$series, levels = .lv)
  .d$assessment <- factor(.d$assessment, levels = assessments)
  .lty <- setNames(rep("solid", length(.epNames)), .epNames)
  .safe <- which(.dir == "lower")
  .lty[.safe] <- rep_len(.cuiSafetyLty, length(.safe))
  .lty <- c(if (endpoints) .lty, CUI = "solid", Peak = "blank")
  .lwd <- c(if (endpoints) setNames(rep(1.05, length(.epNames)), .epNames), CUI = 1.4, Peak = 0)
  .sub <- if (endpoints) {
    paste0("Weights: ", paste(sprintf("%s %.2f", .epNames, x$weightsNorm), collapse = ", "))
  } else {
    ""
  }
  if (.hasCi) {
    .band <- sprintf("bands: %g%% percentile interval", 100 * x$level)
    .sub <- if (nzchar(.sub)) paste0(.sub, "; ", .band) else paste0("B", substring(.band, 2))
  }
  .plotData(.cuiLinesFigure(.rx, cols[.lv], .lty, .lwd, .hasCi,
                            any(!.s$eligible), !is.null(.ref), .sub,
                            if (endpoints) "Desirability / CUI" else "CUI",
                            if (endpoints) "Endpoint utilities and clinical utility index"
                            else "Clinical utility index", length(assessments)), .d)
}

.cuiLinesFigure <- function(rx, cols, lty, lwd, hasCi, anyIneligible, hasRef, subtitle, ylab,
                            title, nAssess) {
  # force every argument: an unevaluated promise would keep the caller's frame
  invisible(list(rx, cols, lty, lwd, hasCi, anyIneligible, hasRef, subtitle, ylab, title, nAssess))
  .lineRows <- function(d) d[d$layer == "line", , drop = FALSE]
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$x, y = .data$value)) +
    (if (hasRef) {
      list(
        ggplot2::geom_vline(data = function(d) d[d$layer == "ref", , drop = FALSE],
                            ggplot2::aes(xintercept = .data$x), linetype = "dashed",
                            color = "grey50"),
        ggplot2::geom_text(data = function(d) d[d$layer == "ref", , drop = FALSE],
                           ggplot2::aes(x = .data$x, y = 1, label = .data$label), vjust = -0.4,
                           size = 3, color = "#5C5C5C")
      )
    }) +
    (if (hasCi) {
      list(
        ggplot2::geom_ribbon(data = function(d) d[d$layer == "line" & !is.na(d$lower), , drop = FALSE],
                             ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$series,
                                          group = interaction(.data$series, .data$assessment)),
                             alpha = 0.13, show.legend = FALSE),
        ggplot2::scale_fill_manual(values = cols, guide = "none")
      )
    }) +
    ggplot2::geom_line(data = .lineRows,
                       ggplot2::aes(color = .data$series, linetype = .data$series,
                                    linewidth = .data$series,
                                    group = interaction(.data$series, .data$assessment))) +
    ggplot2::geom_point(data = function(d) d[d$layer == "line" & d$series == "CUI" & d$eligible, ,
                                            drop = FALSE],
                        color = cols[["CUI"]], size = 1.1, show.legend = FALSE) +
    (if (anyIneligible) {
      ggplot2::geom_point(data = function(d) d[d$layer == "line" & d$series == "CUI" & !d$eligible, ,
                                              drop = FALSE],
                          color = .cuiMuted, size = 2, show.legend = FALSE)
    }) +
    ggplot2::geom_point(data = function(d) d[d$layer == "peak", , drop = FALSE],
                        ggplot2::aes(color = .data$series), size = 3.2) +
    ggplot2::geom_text(data = function(d) d[d$layer == "peak", , drop = FALSE],
                       ggplot2::aes(y = pmin(1, .data$value + 0.04), label = .data$label),
                       color = cols[["Peak"]], fontface = "bold", size = 3.5) +
    # the legend guide lives on the scale: `+ guides()` keeps the plot frames
    ggplot2::scale_color_manual(
      values = cols, breaks = names(cols), name = NULL, drop = FALSE,
      guide = ggplot2::guide_legend(
        nrow = 1,
        override.aes = list(linetype = unname(lty),
                            shape = c(rep(NA, length(cols) - 1L), 16),
                            linewidth = c(unname(lwd)[-length(lwd)], 0.8))
      )
    ) +
    ggplot2::scale_linetype_manual(values = lty, breaks = names(cols), guide = "none") +
    ggplot2::scale_linewidth_manual(values = lwd, guide = "none") +
    .cuiXScale(rx) +
    ggplot2::scale_y_continuous(name = ylab, limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off") +
    ggplot2::labs(title = title, subtitle = if (nzchar(subtitle)) subtitle) +
    (.cuiFacet(nAssess)) +
    .cuiTheme()
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
    ggplot2::geom_col(fill = .cuiNavy, width = 0.6) +
    ggplot2::geom_text(data = function(d) d[d$selFreq > 0, , drop = FALSE],
                       ggplot2::aes(label = sprintf("%.0f%%", 100 * .data$selFreq)),
                       vjust = -0.4, size = 3.2, color = .cuiNavy) +
    ggplot2::scale_y_continuous(name = "Selection frequency", limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::xlab(xlab) +
    ggplot2::labs(title = "How often each regimen is selected", subtitle = subtitle) +
    (.cuiFacet(nAssess)) +
    .cuiTheme()
}

.cuiPlotDoseResponse <- function(x, cols, bands = TRUE) {
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
  if (bands && .cuiNOk(x) > 0) {
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
  .plotData(.cuiDoseResponseFigure(.sub, !is.null(.ci), cols[.epNames]), .d)
}

.cuiDoseResponseFigure <- function(subtitle, hasCi, cols) {
  invisible(list(subtitle, hasCi, cols))
  ggplot2::ggplot(mapping = ggplot2::aes(x = .data$dose)) +
    (if (hasCi) {
      ggplot2::geom_linerange(data = function(d) d[d$layer == "interval", , drop = FALSE],
                              ggplot2::aes(ymin = .data$lower, ymax = .data$upper,
                                           color = .data$endpoint), alpha = 0.5,
                              linewidth = 0.8, show.legend = FALSE)
    }) +
    ggplot2::geom_line(data = function(d) d[d$layer == "fit", , drop = FALSE],
                       ggplot2::aes(y = .data$value, color = .data$endpoint), linewidth = 1.05,
                       show.legend = FALSE) +
    ggplot2::geom_point(data = function(d) d[d$layer == "observed", , drop = FALSE],
                        ggplot2::aes(y = .data$value), color = .cuiNavy, size = 2) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::facet_wrap("endpoint") +
    ggplot2::scale_y_continuous(name = "Probability of the event", limits = c(0, 1)) +
    ggplot2::xlab("Dose") +
    ggplot2::labs(title = "Dose-response of each endpoint", subtitle = subtitle) +
    .cuiTheme()
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
  .numeric <- .cuiRegimenX(.reg)$numeric
  if (.numeric) .d$dose <- .reg$dose[match(.d$regimen, .reg$label)]
  .step <- min(diff(sort(unique(.d$weight))))
  # the current weights as their own rows, so the ring shows off the grid too
  .d$current <- FALSE
  .now <- data.frame(endpoint = factor(names(.cur), levels = levels(.d$endpoint)),
                     weight = unname(.cur), regimen = NA_character_, cui = NA_real_,
                     current = TRUE, stringsAsFactors = FALSE)
  if (.numeric) .now$dose <- NA_real_
  .d <- rbind(.d, .now)
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
    ggplot2::geom_tile(data = function(d) d[!d$current, , drop = FALSE],
                       ggplot2::aes(fill = if (numeric) .data$dose else .data$regimen),
                       width = step, height = 0.8, color = "white", linewidth = 0.4) +
    ggplot2::geom_point(data = function(d) d[d$current, , drop = FALSE], shape = 21, size = 3.5,
                        stroke = 1, color = .cuiNavy, fill = "white") +
    .fill +
    ggplot2::xlab("Weight of the endpoint (others unchanged)") +
    ggplot2::ylab(NULL) +
    ggplot2::labs(title = "Sensitivity of the recommendation to the weights",
                  subtitle = sprintf("Assessment %s; ring: current weight", assessment)) +
    .cuiTheme() +
    ggplot2::theme(legend.position = "bottom")
}
