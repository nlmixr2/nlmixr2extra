.cur <- loadNamespace("nlmixr2extra")

.cuiPlotValues <- function() {
  data.frame(regimen = rep(1:5, 3),
             endpoint = rep(c("tox", "eff", "tol"), each = 5),
             value = c(0.038, 0.058, 0.106, 0.232, 0.533,
                       0.022, 0.119, 0.343, 0.568, 0.681,
                       0.147, 0.206, 0.280, 0.368, 0.466))
}
.cuiPlotEps <- function() {
  list(cuiEndpoint(tox, direction = "lower"), cuiEndpoint(eff), cuiEndpoint(tol))
}

.cuiPlotDr <- function() {
  .d <- withr::with_seed(42, {
    .p <- data.frame(dose = c(10, 25, 50, 100, 200), tox = c(0.04, 0.06, 0.11, 0.23, 0.53),
                     eff = c(0.02, 0.12, 0.34, 0.57, 0.68))
    .d <- .p[rep(1:5, each = 30), ]
    for (.v in c("tox", "eff")) .d[[.v]] <- stats::rbinom(nrow(.d), 1, .d[[.v]])
    .d
  })
  cuiDoseResponse(.d, list(cuiBinary(tox, model = "exponential", direction = "lower"),
                           cuiBinary(eff, model = "logitQuadratic")),
                  dose = dose, nBoot = 30, seed = 1)
}

test_that("every plot type builds and keeps no caller frame or extra data", {
  .r <- .cuiPlotDr()
  for (.t in c("cui", "endpoints", "selection", "doseResponse")) {
    .fig <- plot(.r, .t)
    expect_s3_class(.fig, "ggplot")
    expect_identical(.figureHeldData(.fig), character(0))
    invisible(ggplot2::ggplot_build(.fig))
    expect_identical(.figureHeldData(.fig), character(0))
  }
  .s <- cuiSensitivity(.r)
  expect_identical(.figureHeldData(plot(.s)), character(0))
  .fromCaller <- local({
    .callerData <- data.frame(x = seq_len(1000))
    plot(.r, "endpoints")
  })
  expect_identical(.figureHeldData(.fromCaller), character(0))
})

test_that("plot data reflect the CUI, eligibility and assessments", {
  .r <- cuiFromValues(.cuiPlotValues(), .cuiPlotEps(), weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
  .p <- plot(.r)
  expect_equal(.p$data$est, .r$summary$est)
  expect_equal(.p$data$x, 1:5)
  expect_match(.p$labels$subtitle, "Recommended: 5")
  .e <- plot(update(.r, eligible = c("1", "2", "3", "4")))
  expect_match(.e$labels$subtitle, "Recommended: 4")
  expect_equal(nrow(ggplot2::layer_data(.e, 2)), 4)
  .ep <- plot(.r, "endpoints")
  expect_equal(levels(.ep$data$series), c("tox", "eff", "tol", "CUI"))
  expect_equal(.ep$data$value[.ep$data$series == "tox"], 1 - .cuiPlotValues()$value[1:5])
  # non-numeric regimens are placed in order
  .lab <- cuiFromValues(transform(.cuiPlotValues(), regimen = letters[regimen]), .cuiPlotEps())
  expect_equal(plot(.lab)$data$x, 1:5)
  expect_error(plot(.r, "selection"), "uncertainty draws")
  expect_error(plot(.r, "doseResponse"), "cuiDoseResponse")
  .two <- cuiFromValues(rbind(cbind(.cuiPlotValues(), assessment = "m3"),
                              cbind(.cuiPlotValues(), assessment = "m6")), .cuiPlotEps())
  expect_equal(levels(plot(.two)$data$assessment), c("m3", "m6"))
  expect_equal(levels(plot(.two, assessment = "m6")$data$assessment), "m6")
  expect_error(plot(.two, assessment = "m9"), "unknown assessment")
})

test_that("end labels are spread apart", {
  .y <- .cur$.cuiSpreadLabels(c(0.62, 0.63, 0.30, 0.80), 0.045)
  expect_true(all(diff(sort(.y)) >= 0.045 - 1e-12))
  expect_equal(.y[3], 0.30, tolerance = 0.05)
})

test_that("sensitivity agrees with the recommendation at the current weights", {
  .r <- cuiFromValues(.cuiPlotValues(), .cuiPlotEps(), weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
  .s <- cuiSensitivity(.r, weights = c(0, 0.2, 0.3, 0.5, 2, 5))
  for (.e in c("tox", "eff", "tol")) {
    .now <- .r$weights[[.e]]
    expect_equal(.s$regimen[.s$endpoint == .e & .s$weight == .now], "5")
  }
  # matches update() with the swept weight
  .row <- .s[.s$endpoint == "tox" & .s$weight == 2, ]
  .u <- update(.r, weights = c(tox = 2, eff = 0.5, tol = 0.3))
  expect_equal(.row$regimen, .u$summary$regimen[.u$summary$best])
  expect_equal(.row$cui, max(.u$summary$est))
  expect_output(print(.s), "tox \\(now 0.2\\)")
  # all weights zero for an endpoint that is the only weighted one
  .one <- update(.r, weights = c(tox = 0, eff = 1, tol = 0))
  .z <- cuiSensitivity(.one, weights = c(0, 1))
  expect_true(is.na(.z$regimen[.z$endpoint == "eff" & .z$weight == 0]))
  expect_error(cuiSensitivity(.r, weights = 1), "length")
  expect_error(cuiSensitivity(1), "nlmixr2Cui")
  # labels instead of doses use a categorical fill
  .lab <- cuiFromValues(transform(.cuiPlotValues(), regimen = letters[regimen]), .cuiPlotEps())
  expect_s3_class(plot(cuiSensitivity(.lab)), "ggplot")
})
