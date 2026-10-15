#' Interactive clinical utility index desk tool
#'
#' A Shiny app around an `nlmixr2Cui` result (from [cui()],
#' [cuiDoseResponse()] or [cuiFromValues()]), modeled on the ACoP 2026 CUI
#' workshop app.  Weights are hierarchical, as in [cuiGroupWeights()]: one
#' slider per endpoint group (for example efficacy and safety) and, under
#' it, one slider per endpoint for its share within the group.  Every change
#' only re-weights the stored endpoint values (`update()`), so the app
#' responds at once; nothing is simulated or refit.
#'
#' The figures are the package's own (`plot()`, with lines and percentile
#' bands), so the app and the reports match.  Tabs show the CUI and
#' endpoint utilities, the selection frequencies (with uncertainty draws),
#' the weight sensitivity ([cuiSensitivity()]) and, for [cuiDoseResponse()]
#' results, the dose-response fits.  The current weights and results can
#' be downloaded as CSV.
#'
#' In `mode = "desk"` recommendations are always restricted to the eligible
#' regimens (for example approved doses and durations, per assessment);
#' `mode = "explorer"` lets the user switch the restriction off.
#'
#' @param x An `nlmixr2Cui` object.
#' @param mode `"explorer"` or `"desk"`, as described above.
#' @param eligible Eligible regimens, replacing those of `x` (see
#'   [cuiFromValues()]).
#' @param presets Optional named list of weight presets; each is a named
#'   vector of endpoint weights, or a list with `group` and `within` as in
#'   [cuiGroupWeights()].
#' @param title App title.
#' @param reference,colors Passed to [plot.nlmixr2Cui()].
#' @param maxWeight Upper end of the weight sliders.
#' @return A `shiny.appobj`; printing it (or [shiny::runApp()]) starts the
#'   app.
#' @family Clinical utility
#' @export
#' @examples
#' \dontrun{
#' values <- data.frame(regimen = rep(c(0, 150, 400), 3),
#'                      endpoint = rep(c("pain", "dysp", "hotFlash"), each = 3),
#'                      value = c(0.2, 0.45, 0.74, 0.3, 0.39, 0.54, 0.07, 0.24, 0.42))
#' res <- cuiFromValues(values, list(cuiEndpoint(pain), cuiEndpoint(dysp),
#'                                   cuiEndpoint(hotFlash, direction = "lower")),
#'                      regimens = cuiRegimens(dose = c(0, 150, 400), units = "mg",
#'                                             label = c("0", "150", "400")))
#' cuiApp(res, presets = list(trial = c(pain = 1, dysp = 0, hotFlash = 1)))
#' }
cuiApp <- function(x, mode = c("explorer", "desk"), eligible = NULL, presets = NULL,
                   title = "Clinical utility index", reference = NULL, colors = NULL,
                   maxWeight = 5) {
  if (!requireNamespace("shiny", quietly = TRUE) || !requireNamespace("bslib", quietly = TRUE)) {
    stop("cuiApp() needs the 'shiny' and 'bslib' packages", call. = FALSE)
  }
  if (!inherits(x, "nlmixr2Cui")) stop("'x' must be an nlmixr2Cui object", call. = FALSE)
  mode <- match.arg(mode)
  checkmate::assertString(title)
  checkmate::assertNumber(maxWeight, lower = 0, finite = TRUE)
  if (!is.null(eligible)) x <- update(x, eligible = eligible)
  .spec <- .cuiAppSpec(x, presets, maxWeight)
  shiny::shinyApp(ui = .cuiAppUi(x, .spec, mode, title),
                  server = .cuiAppServer(x, .spec, mode, reference, colors))
}

#' Endpoints, groups, slider ids and starting weights of the app
#' @noRd
.cuiAppSpec <- function(x, presets, maxWeight) {
  .epNames <- vapply(x$endpoints, `[[`, character(1), "name")
  .grp <- .cuiEndpointGroups(x)
  .groups <- unique(.grp)
  .raw <- if (is.null(x$weights)) setNames(rep(1, length(.epNames)), .epNames) else x$weights[.epNames]
  .presets <- lapply(presets, .cuiAppPreset, x = x, epNames = .epNames, grp = .grp)
  if (!is.null(presets) && (is.null(names(presets)) || any(names(presets) == ""))) {
    stop("'presets' must be a named list", call. = FALSE)
  }
  .start <- .cuiAppPreset(.raw, x, .epNames, .grp)
  if (max(c(unlist(.start), unlist(.presets))) > maxWeight) {
    stop("a starting or preset weight exceeds 'maxWeight'", call. = FALSE)
  }
  # start on the preset that matches the starting weights, if any
  .match <- Filter(function(p) {
    isTRUE(all.equal(cuiGroupWeights(x, p$group, p$within), cuiGroupWeights(x, .start$group, .start$within)))
  }, .presets)
  list(epNames = .epNames, grp = .grp, groups = .groups,
       startPreset = if (length(.match) > 0) names(.match)[1] else "custom",
       groupId = setNames(paste0("cuiGroup", seq_along(.groups)), .groups),
       epId = setNames(paste0("cuiWithin", seq_along(.epNames)), .epNames),
       start = .start, presets = .presets, maxWeight = maxWeight)
}

#' A preset as group and within-group weights
#'
#' Endpoint weights map exactly: a group weighs the sum of its endpoints'
#' weights and the endpoints share it in proportion.
#' @noRd
.cuiAppPreset <- function(p, x, epNames, grp) {
  if (is.list(p)) {
    .w <- cuiGroupWeights(x, group = p$group, within = p$within)
  } else {
    checkmate::assertNumeric(p, lower = 0, finite = TRUE, any.missing = FALSE, names = "unique")
    .bad <- setdiff(names(p), epNames)
    if (length(.bad) > 0) {
      stop(sprintf("unknown endpoint(s) in a preset: %s", paste(.bad, collapse = ", ")), call. = FALSE)
    }
    .w <- setNames(rep(1, length(epNames)), epNames)
    .w[names(p)] <- p
  }
  .g <- vapply(unique(grp), function(k) sum(.w[grp == k]), numeric(1))
  # within-group weights are the endpoint weights themselves (a zero group keeps 1s)
  .within <- .w
  for (.k in names(.g)) if (.g[[.k]] == 0) .within[grp == .k] <- 1
  list(group = .g, within = .within)
}

.cuiAppUi <- function(x, spec, mode, title) {
  .sliders <- lapply(spec$groups, function(g) {
    .eps <- spec$epNames[spec$grp == g]
    shiny::div(
      class = "cui-group",
      shiny::sliderInput(spec$groupId[[g]], paste0(.cuiAppTitle(g), " weight"), 0, spec$maxWeight,
                         unname(spec$start$group[[g]]), step = 0.05),
      shiny::div(
        class = "cui-within",
        lapply(.eps, function(e) {
          shiny::sliderInput(spec$epId[[e]], e, 0, spec$maxWeight,
                             unname(spec$start$within[[e]]), step = 0.05, ticks = FALSE)
        })
      )
    )
  })
  .hasDraws <- .cuiNOk(x) > 0
  .tabs <- list(
    bslib::nav_panel(
      "Clinical utility",
      bslib::card(
        full_screen = TRUE,
        bslib::card_header(shiny::uiOutput("caption")),
        bslib::card_body(
          shiny::plotOutput("mainPlot", height = "55vh",
                            hover = shiny::hoverOpts("mainHover", delay = 150)),
          shiny::uiOutput("hoverInfo")
        )
      ),
      bslib::card(bslib::card_header("Regimens"), shiny::tableOutput("regimenTable"))
    ),
    if (.hasDraws) {
      bslib::nav_panel("Selection",
                       bslib::card(full_screen = TRUE,
                                   shiny::plotOutput("selectionPlot", height = "55vh")))
    },
    bslib::nav_panel(
      "Weight sensitivity",
      bslib::card(full_screen = TRUE, shiny::plotOutput("sensitivityPlot", height = "45vh"),
                  shiny::verbatimTextOutput("sensitivityText"))
    ),
    if (!is.null(x$observed)) {
      bslib::nav_panel("Dose-response",
                       bslib::card(full_screen = TRUE,
                                   shiny::plotOutput("doseResponsePlot", height = "55vh")))
    }
  )
  .tabs <- Filter(Negate(is.null), .tabs)
  bslib::page_sidebar(
    title = title,
    fillable = FALSE,
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#1B365D",
                            "navbar-bg" = "#1B365D"),
    sidebar = bslib::sidebar(
      width = 320,
      shiny::tags$style(shiny::HTML(paste(
        ".cui-within { border-left: 3px solid #d8dee6; padding-left: 10px; margin-left: 4px; }",
        ".cui-group { margin-bottom: 10px; }",
        ".cui-within .irs--shiny .irs-bar { background: #8aa1bd; border-color: #8aa1bd; }",
        sep = "\n"))),
      if (length(spec$presets) > 0) {
        shiny::selectInput("preset", "Weight preset",
                           choices = c(names(spec$presets), "Custom" = "custom"),
                           selected = spec$startPreset)
      },
      shiny::radioButtons("combine", "CUI form",
                          choices = c(Additive = "additive", Geometric = "geometric"),
                          selected = x$combine, inline = TRUE),
      if (length(x$assessments) > 1) {
        shiny::selectInput("assessment", "Assessment", choices = x$assessments)
      },
      .sliders,
      shiny::checkboxInput("showEndpoints", "Show endpoint curves", value = TRUE),
      shiny::checkboxInput("bands", "Show uncertainty bands", value = TRUE),
      if (mode == "explorer") {
        shiny::checkboxInput("eligibleOnly", "Recommend eligible regimens only",
                             value = !all(x$eligible$eligible))
      },
      shiny::downloadButton("downloadWeights", "Weights (CSV)", class = "btn-sm"),
      shiny::downloadButton("downloadResults", "Results (CSV)", class = "btn-sm")
    ),
    do.call(bslib::navset_tab, c(list(id = "tabs"), .tabs))
  )
}

.cuiAppTitle <- function(g) paste0(toupper(substring(g, 1, 1)), substring(g, 2))

.cuiAppServer <- function(x, spec, mode, reference, colors) {
  force(x)
  force(spec)
  force(mode)
  force(reference)
  force(colors)
  function(input, output, session) {
    .slider <- function(id) {
      .v <- input[[id]]
      if (is.null(.v)) NA_real_ else .v
    }
    groupWeights <- shiny::reactive(vapply(spec$groupId, .slider, numeric(1)))
    withinWeights <- shiny::reactive(vapply(spec$epId, .slider, numeric(1)))
    # slider drags re-weight after a short pause
    weightsNow <- shiny::reactive({
      .g <- groupWeights()
      .w <- withinWeights()
      if (anyNA(.g) || anyNA(.w)) return(NULL)
      cuiGroupWeights(x, group = .g, within = .w)
    })
    weights <- shiny::debounce(weightsNow, 250)
    assessment <- shiny::reactive({
      if (is.null(input$assessment)) x$assessments[1] else input$assessment
    })
    eligibleSpec <- shiny::reactive({
      if (mode == "desk" || isTRUE(input$eligibleOnly)) x$eligibleSpec else NULL
    })
    result <- shiny::reactive({
      .w <- weights()
      shiny::req(.w)
      shiny::validate(shiny::need(sum(.w) > 0, "Give at least one weight above zero."))
      update(x, weights = .w, combine = if (is.null(input$combine)) x$combine else input$combine,
             eligible = eligibleSpec())
    })
    mainPlot <- shiny::reactive({
      plot(result(), if (isTRUE(input$showEndpoints)) "endpoints" else "cui",
           assessment = assessment(), reference = reference, colors = colors,
           bands = isTRUE(input$bands))
    })
    output$mainPlot <- shiny::renderPlot(mainPlot())
    output$caption <- shiny::renderUI({
      .r <- result()
      .s <- .r$summary[.r$summary$assessment == assessment(), , drop = FALSE]
      .b <- .s[.s$best, , drop = FALSE]
      .w <- .r$weightsNorm
      shiny::HTML(sprintf(
        "<b>Recommended: %s</b> &nbsp;(CUI %.3f) &nbsp;&middot;&nbsp; weights: %s",
        paste(.b$regimen, collapse = ", "), .b$est[1],
        paste(sprintf("%s %.2f", names(.w), .w), collapse = ", ")))
    })
    output$hoverInfo <- shiny::renderUI({
      .h <- input$mainHover
      if (is.null(.h)) return(shiny::helpText("Hover over a line to read its value."))
      .d <- mainPlot()$data
      .d <- .d[.d$layer == "line" & .d$assessment == assessment(), , drop = FALSE]
      .np <- shiny::nearPoints(.d, .h, xvar = "x", yvar = "value", threshold = 15, maxpoints = 3)
      if (nrow(.np) == 0) return(shiny::helpText("Hover over a line to read its value."))
      shiny::helpText(paste(sprintf("%s at %s: %.3f%s", .np$series, .np$regimen, .np$value,
                                    ifelse(is.na(.np$lower), "",
                                           sprintf(" (%.3f to %.3f)", .np$lower, .np$upper))),
                            collapse = "; "))
    })
    regimenTable <- shiny::reactive({
      .r <- result()
      .s <- .r$summary[.r$summary$assessment == assessment() & .r$summary$eligible, , drop = FALSE]
      .tab <- data.frame(Regimen = .s$regimen, CUI = .s$est, stringsAsFactors = FALSE)
      if (.cuiNOk(.r) > 0) {
        .tab$Lower <- .s$lower
        .tab$Upper <- .s$upper
        .tab$Selected <- sprintf("%.0f%%", 100 * .s$selFreq)
      }
      .v <- .r$values[.r$values$draw == 0 & .r$values$assessment == assessment(), , drop = FALSE]
      for (.e in spec$epNames) {
        .tab[[.e]] <- .v$value[.v$endpoint == .e][match(.s$regimen, .v$regimen[.v$endpoint == .e])]
      }
      .tab$Recommended <- ifelse(.s$best, "yes", "")
      .tab
    })
    output$regimenTable <- shiny::renderTable(regimenTable(), digits = 3)
    output$selectionPlot <- shiny::renderPlot(plot(result(), "selection", assessment = assessment()))
    sensitivity <- shiny::reactive(cuiSensitivity(result(), assessment = assessment()))
    output$sensitivityPlot <- shiny::renderPlot(plot(sensitivity()))
    output$sensitivityText <- shiny::renderPrint(print(sensitivity()))
    output$doseResponsePlot <- shiny::renderPlot(
      plot(result(), "doseResponse", colors = colors, bands = isTRUE(input$bands)))
    # presets set the sliders; moving a slider away from a preset marks it custom
    shiny::observeEvent(input$preset, {
      .p <- spec$presets[[input$preset]]
      if (is.null(.p)) return()
      for (.g in spec$groups) shiny::updateSliderInput(session, spec$groupId[[.g]], value = .p$group[[.g]])
      for (.e in spec$epNames) shiny::updateSliderInput(session, spec$epId[[.e]], value = .p$within[[.e]])
    }, ignoreInit = TRUE)
    shiny::observeEvent(weightsNow(), {
      if (is.null(input$preset) || input$preset == "custom") return()
      .p <- spec$presets[[input$preset]]
      .now <- weightsNow()
      if (is.null(.now)) return()
      if (max(abs(.now - cuiGroupWeights(x, group = .p$group, within = .p$within))) > 0.03) {
        shiny::updateSelectInput(session, "preset", selected = "custom")
      }
    }, ignoreInit = TRUE)
    weightTable <- shiny::reactive({
      .r <- result()
      data.frame(endpoint = spec$epNames, group = spec$grp,
                 groupWeight = unname(groupWeights()[spec$grp]),
                 withinWeight = unname(withinWeights()),
                 weight = unname(weights()), normalized = unname(.r$weightsNorm),
                 combine = .r$combine, assessment = assessment(),
                 eligibleOnly = !is.null(eligibleSpec()),
                 source = if (is.null(.r$settings$description)) .r$settings$source
                          else .r$settings$description,
                 stringsAsFactors = FALSE)
    })
    output$downloadWeights <- shiny::downloadHandler(
      filename = function() "cui-weights.csv",
      content = function(file) utils::write.csv(weightTable(), file, row.names = FALSE)
    )
    output$downloadResults <- shiny::downloadHandler(
      filename = function() "cui-results.csv",
      content = function(file) utils::write.csv(regimenTable(), file, row.names = FALSE)
    )
  }
}

#' Write a deployable CUI desk tool
#'
#' Saves `x` (which keeps no fits or data, only endpoint summaries) and the
#' app settings to `dir` with an `app.R` that calls [cuiApp()], ready for
#' `shiny::runApp(dir)` or deployment (for example to shinyapps.io or
#' Posit Connect).  The deployed app needs 'nlmixr2extra', 'shiny' and
#' 'bslib'.
#'
#' @param x An `nlmixr2Cui` object.
#' @param dir Directory to write (created if needed).
#' @param ... Arguments for [cuiApp()], such as `mode`, `presets` or
#'   `eligible`.
#' @param overwrite Replace an existing `app.R` / `cui.rds`.
#' @return `dir`, invisibly.
#' @family Clinical utility
#' @export
#' @examples
#' values <- data.frame(regimen = rep(1:3, 2), endpoint = rep(c("eff", "ae"), each = 3),
#'                      value = c(0.2, 0.5, 0.7, 0.05, 0.1, 0.3))
#' res <- cuiFromValues(values, list(cuiEndpoint(eff), cuiEndpoint(ae, direction = "lower")))
#' d <- cuiAppExport(res, file.path(tempdir(), "cuiDesk"), mode = "desk")
#' list.files(d)
cuiAppExport <- function(x, dir, ..., overwrite = FALSE) {
  if (!inherits(x, "nlmixr2Cui")) stop("'x' must be an nlmixr2Cui object", call. = FALSE)
  checkmate::assertString(dir, min.chars = 1)
  checkmate::assertFlag(overwrite)
  .args <- list(...)
  if (length(.args) > 0 && (is.null(names(.args)) || any(names(.args) == ""))) {
    stop("arguments for cuiApp() must be named", call. = FALSE)
  }
  .app <- file.path(dir, "app.R")
  .rds <- file.path(dir, "cui.rds")
  if (!overwrite && (file.exists(.app) || file.exists(.rds))) {
    stop(sprintf("'%s' already has an app; use overwrite = TRUE", dir), call. = FALSE)
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(x = x, args = .args), .rds)
  writeLines(c(
    "# Clinical utility index desk tool written by nlmixr2extra::cuiAppExport()",
    "library(nlmixr2extra)",
    ".cui <- readRDS(\"cui.rds\")",
    "do.call(cuiApp, c(list(.cui$x), .cui$args))"
  ), .app)
  invisible(dir)
}
