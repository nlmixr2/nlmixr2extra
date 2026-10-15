# caches default to the working directory; keep them out of the package
withr::local_dir(withr::local_tempdir())

.cur <- loadNamespace("nlmixr2extra")

.appValues <- function() {
  data.frame(regimen = rep(c(0, 150, 400), 4),
             endpoint = rep(c("dysm", "nmpp", "hf", "bmd"), each = 3),
             value = c(0.21, 0.45, 0.74, 0.365, 0.50, 0.56, 0.07, 0.24, 0.42, 0.10, 0.12, 0.40))
}
.appRes <- function(...) {
  cuiFromValues(.appValues(),
                list(cuiEndpoint(dysm), cuiEndpoint(nmpp),
                     cuiEndpoint(hf, direction = "lower"), cuiEndpoint(bmd, direction = "lower")),
                regimens = cuiRegimens(dose = c(0, 150, 400), units = "mg",
                                       label = c("0", "150", "400")), ...)
}

# the app's page as HTML
.uiHtml <- function(x, mode = "explorer", presets = NULL) {
  .spec <- .cur$.cuiAppSpec(x, presets, 5)
  as.character(htmltools::renderTags(.cur$.cuiAppUi(x, .spec, mode, "CUI"))$html)
}

test_that("endpoint groups and hierarchical weights", {
  expect_equal(cuiEndpoint(e)$group, "efficacy")
  expect_equal(cuiEndpoint(e, direction = "lower")$group, "safety")
  expect_equal(cuiEndpoint(e, group = "tolerability")$group, "tolerability")
  expect_equal(cuiBinary(e, group = "pk")$group, "pk")
  .r <- .appRes()
  expect_equal(.cur$.cuiEndpointGroups(.r), c("efficacy", "efficacy", "safety", "safety"))
  .w <- cuiGroupWeights(.r, group = c(efficacy = 2, safety = 1), within = c(dysm = 3, nmpp = 1))
  expect_equal(.w, c(dysm = 1.5, nmpp = 0.5, hf = 0.5, bmd = 0.5))
  # a group's endpoints together carry the group weight
  expect_equal(sum(.w[c("dysm", "nmpp")]), 2)
  expect_equal(cuiGroupWeights(.r), c(dysm = 0.5, nmpp = 0.5, hf = 0.5, bmd = 0.5))
  expect_equal(cuiGroupWeights(.r, within = c(hf = 0, bmd = 0))[c("hf", "bmd")], c(hf = 0, bmd = 0))
  expect_error(cuiGroupWeights(.r, group = c(nope = 1)), "unknown")
  expect_error(cuiGroupWeights(.r, within = c(dysm = -1)))
  expect_equal(update(.r, weights = .w)$weightsNorm, .w / sum(.w))
})

test_that("presets map exactly onto group and within-group sliders", {
  .r <- .appRes()
  .spec <- .cur$.cuiAppSpec(.r, list(dce = c(dysm = 1.48, nmpp = 1.49, hf = 3.62, bmd = 0.49),
                                      hier = list(group = c(efficacy = 3), within = c(hf = 2))), 5)
  .p <- .spec$presets$dce
  expect_equal(cuiGroupWeights(.r, .p$group, .p$within),
               c(dysm = 1.48, nmpp = 1.49, hf = 3.62, bmd = 0.49))
  expect_equal(unname(.spec$presets$hier$group), c(3, 1))
  # the starting weights (all 1) match no preset here; a matching preset is selected
  expect_equal(.spec$startPreset, "custom")
  .r2 <- .appRes(weights = c(dysm = 1.48, nmpp = 1.49, hf = 3.62, bmd = 0.49))
  expect_equal(.cur$.cuiAppSpec(.r2, list(a = c(dysm = 1), dce = c(dysm = 1.48, nmpp = 1.49,
                                                                    hf = 3.62, bmd = 0.49)), 5)$startPreset,
               "dce")
  expect_error(.cur$.cuiAppSpec(.r, list(c(dysm = 1)), 5), "named list")
  expect_error(.cur$.cuiAppSpec(.r, list(a = c(nope = 1)), 5), "unknown endpoint")
  expect_error(.cur$.cuiAppSpec(.r, list(a = c(dysm = 9)), 5), "maxWeight")
})

test_that("the app builds its sliders and tabs", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  .res <- .appRes(weights = c(dysm = 1, nmpp = 1, hf = 2, bmd = 1))
  .pre <- list(trial = c(dysm = 1, nmpp = 1, hf = 1, bmd = 1))
  expect_s3_class(cuiApp(.res, presets = .pre), "shiny.appobj")
  .html <- .uiHtml(.res, presets = .pre)
  for (.id in c("cuiGroup1", "cuiGroup2", "cuiWithin1", "cuiWithin4", "preset", "combine",
                "mainPlot", "regimenTable", "sensitivityPlot", "downloadWeights")) {
    expect_match(.html, .id, fixed = TRUE)
  }
  expect_match(.html, "Efficacy weight", fixed = TRUE)
  expect_match(.html, "Safety weight", fixed = TRUE)
  # no draws, no dose-response data: no selection or dose-response tabs
  expect_false(grepl("selectionPlot", .html, fixed = TRUE))
  expect_false(grepl("doseResponsePlot", .html, fixed = TRUE))
  expect_error(cuiApp(1), "nlmixr2Cui")
})

test_that("the server re-weights, restricts eligibility and downloads", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  .r <- .appRes(eligible = c("150", "400"))
  shiny::testServer(cuiApp(.r), {
    session$setInputs(cuiGroup1 = 2, cuiGroup2 = 1, cuiWithin1 = 3, cuiWithin2 = 1,
                      cuiWithin3 = 1, cuiWithin4 = 1, combine = "additive",
                      showEndpoints = TRUE, bands = TRUE, eligibleOnly = FALSE)
    session$elapse(500)
    .w <- cuiGroupWeights(.r, c(efficacy = 2, safety = 1), c(dysm = 3, nmpp = 1, hf = 1, bmd = 1))
    expect_equal(result()$weightsNorm, .w / sum(.w))
    # explorer with the restriction switched off: every regimen is eligible
    expect_true(all(result()$eligible$eligible))
    expect_equal(nrow(regimenTable()), 3)
    session$setInputs(eligibleOnly = TRUE)
    expect_equal(regimenTable()$Regimen, c("150", "400"))
    expect_equal(mainPlot()$data$label[mainPlot()$data$layer == "peak"],
                 paste0("max at ", result()$summary$regimen[result()$summary$best], " mg"))
    session$setInputs(combine = "geometric")
    expect_equal(result()$combine, "geometric")
    .f <- output$downloadWeights
    .wt <- utils::read.csv(.f)
    expect_equal(.wt$weight, unname(.w))
    expect_equal(.wt$group, c("efficacy", "efficacy", "safety", "safety"))
    expect_equal(nrow(utils::read.csv(output$downloadResults)), 2)
    expect_match(output$sensitivityText, "now")
    # all weights zero is a message, not a crash
    session$setInputs(cuiGroup1 = 0, cuiGroup2 = 0)
    session$elapse(500)
    expect_error(result(), "at least one weight")
  })
  # desk mode always restricts
  shiny::testServer(cuiApp(.r, mode = "desk"), {
    session$setInputs(cuiGroup1 = 1, cuiGroup2 = 1, cuiWithin1 = 1, cuiWithin2 = 1,
                      cuiWithin3 = 1, cuiWithin4 = 1, combine = "additive")
    session$elapse(500)
    expect_equal(regimenTable()$Regimen, c("150", "400"))
  })
})

test_that("draws and dose-response data add their tabs", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  .d <- withr::with_seed(1, {
    .d <- data.frame(dose = rep(c(0, 50, 100), each = 20))
    .d$eff <- stats::rbinom(60, 1, 0.2 + 0.004 * .d$dose)
    .d$tox <- stats::rbinom(60, 1, 0.05 + 0.002 * .d$dose)
    .d
  })
  .r <- cuiDoseResponse(.d, list(cuiBinary(eff), cuiBinary(tox, direction = "lower")),
                        dose = dose, nBoot = 20, seed = 1)
  .html <- .uiHtml(.r)
  expect_match(.html, "selectionPlot", fixed = TRUE)
  expect_match(.html, "doseResponsePlot", fixed = TRUE)
})

test_that("cuiAppExport() writes a runnable app", {
  .r <- .appRes()
  .dir <- cuiAppExport(.r, "deskApp", mode = "desk", title = "Desk")
  expect_setequal(list.files(.dir), c("app.R", "cui.rds"))
  .saved <- readRDS(file.path(.dir, "cui.rds"))
  expect_equal(.saved$args, list(mode = "desk", title = "Desk"))
  expect_equal(.saved$x$summary, .r$summary)
  expect_silent(parse(file.path(.dir, "app.R")))
  expect_error(cuiAppExport(.r, "deskApp"), "overwrite")
  expect_invisible(cuiAppExport(.r, "deskApp", overwrite = TRUE))
  expect_error(cuiAppExport(.r, "other", "desk"), "named")
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  .app <- withr::with_dir(.dir, eval(parse(file.path("app.R"))))
  expect_s3_class(.app, "shiny.appobj")
})
