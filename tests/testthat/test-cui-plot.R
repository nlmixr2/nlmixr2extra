# caches default to the working directory; keep them out of the package
withr::local_dir(withr::local_tempdir())

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
  .line <- .p$data[.p$data$layer == "line", ]
  expect_equal(.line$value, .r$summary$est)
  expect_equal(.line$x, 1:5)
  expect_equal(.p$data$label[.p$data$layer == "peak"], "max at 5")
  expect_false(any(.p$data$layer == "ref"))
  # restricted eligibility: the peak moves and eligible regimens get reference lines
  .e <- plot(update(.r, eligible = c("1", "2", "3", "4")))
  expect_equal(.e$data$label[.e$data$layer == "peak"], "max at 4")
  expect_equal(.e$data$x[.e$data$layer == "ref"], 1:4)
  .ref <- plot(.r, reference = c("150 QD" = 2, "200 BID" = 4))
  expect_equal(.ref$data$label[.ref$data$layer == "ref"], c("150 QD", "200 BID"))
  expect_false(any(plot(update(.r, eligible = c("1", "2")), reference = character(0))$data$layer == "ref"))
  expect_error(plot(.r, reference = "nope"), "unknown regimen")
  .ep <- plot(.r, "endpoints")
  expect_equal(levels(.ep$data$series), c("tox", "eff", "tol", "CUI", "Peak"))
  expect_equal(.ep$data$value[.ep$data$series %in% "tox"], 1 - .cuiPlotValues()$value[1:5])
  # safety endpoints dashed, benefit solid, CUI navy; colors can be replaced
  .sc <- .ep$scales$get_scales("linetype")
  expect_equal(unname(.sc$palette(5)), c("dashed", "solid", "solid", "solid", "blank"))
  .cc <- plot(.r, "endpoints", colors = c(eff = "#1B9E77", CUI = "#000000"))
  .pal <- .cc$scales$get_scales("colour")$palette(5)
  expect_equal(unname(.pal[c(2, 4)]), c("#1B9E77", "#000000"))
  expect_error(plot(.r, colors = c(nope = "red")), "unknown name")
  # units label the peak
  .u <- cuiFromValues(.cuiPlotValues(), .cuiPlotEps(),
                      regimens = cuiRegimens(dose = 1:5, units = "mg", label = as.character(1:5)))
  expect_equal(plot(.u)$data$label[plot(.u)$data$layer == "peak"], "max at 4 mg")
  # non-numeric regimens are placed in order
  .lab <- cuiFromValues(transform(.cuiPlotValues(), regimen = letters[regimen]), .cuiPlotEps())
  .lp <- plot(.lab)
  expect_equal(.lp$data$x[.lp$data$layer == "line"], 1:5)
  expect_error(plot(.r, "selection"), "uncertainty draws")
  expect_error(plot(.r, "doseResponse"), "cuiDoseResponse")
  .two <- cuiFromValues(rbind(cbind(.cuiPlotValues(), assessment = "m3"),
                              cbind(.cuiPlotValues(), assessment = "m6")), .cuiPlotEps())
  expect_equal(levels(plot(.two)$data$assessment), c("m3", "m6"))
  expect_equal(sum(plot(.two)$data$layer == "peak"), 2)
  expect_equal(levels(plot(.two, assessment = "m6")$data$assessment), "m6")
  expect_error(plot(.two, assessment = "m9"), "unknown assessment")
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

test_that("lines carry percentile bands for the CUI and every endpoint", {
  .r <- .cuiPlotDr()
  .ep <- plot(.r, "endpoints")
  .line <- .ep$data[.ep$data$layer == "line", ]
  expect_false(anyNA(.line$lower))
  expect_true(all(.line$lower <= .line$value + 1e-12 | .line$series != "CUI"))
  # endpoint bands are the utility percentiles over the draws
  .u <- .r$utility[.r$utility$draw > 0 & .r$utility$endpoint == "eff" & .r$utility$regimen == "100", ]
  expect_equal(.line$upper[.line$series == "eff" & .line$regimen == "100"],
               unname(stats::quantile(.u$utility, 0.975)))
  expect_true(any(vapply(.ep$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1))))
  .nb <- plot(.r, "endpoints", bands = FALSE)
  expect_false(any(vapply(.nb$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1))))
  expect_match(.ep$labels$subtitle, "bands: 95% percentile interval")
  expect_false(any(vapply(plot(.r, "doseResponse", bands = FALSE)$layers,
                          function(l) inherits(l$geom, "GeomLinerange"), logical(1))))
  expect_identical(.figureHeldData(.ep), character(0))
})

test_that("default reference lines follow each panel and the current weight is ringed", {
  .two <- cuiFromValues(rbind(cbind(.cuiPlotValues(), assessment = "A"),
                              cbind(.cuiPlotValues(), assessment = "B")), .cuiPlotEps(),
                        eligible = data.frame(assessment = c("A", "A", "B", "B", "B"),
                                              regimen = c("1", "2", "3", "4", "5")))
  .ref <- plot(.two)$data
  .ref <- .ref[.ref$layer == "ref", ]
  expect_equal(.ref$x[.ref$assessment == "A"], 1:2)
  expect_equal(.ref$x[.ref$assessment == "B"], 3:5)
  .r <- cuiFromValues(.cuiPlotValues(), .cuiPlotEps(), weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
  .p <- plot(cuiSensitivity(.r))
  .cur <- .p$data[.p$data$current, ]
  expect_equal(sort(.cur$weight), c(0.2, 0.3, 0.5))
  expect_identical(.figureHeldData(.p), character(0))
})
