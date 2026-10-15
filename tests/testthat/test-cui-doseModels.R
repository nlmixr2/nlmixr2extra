.cur <- loadNamespace("nlmixr2extra")

.doses <- c(0, 10, 25, 50, 100)
.n <- rep(30, 5)
.fitP <- function(x, ...) .cur$.cuiDrModel(cuiBinary(e, ...), x, .n, .doses)$p

test_that("logit models without constraints match glm()", {
  .x <- c(2, 5, 9, 14, 20)
  .d <- data.frame(dose = .doses, x = .x, n = .n)
  .lin <- stats::glm(cbind(x, n - x) ~ dose, binomial, data = .d)
  expect_equal(.fitP(.x, model = "logitLinear"), unname(stats::fitted(.lin)), tolerance = 1e-5)
  .quad <- stats::glm(cbind(x, n - x) ~ dose + I(dose^2), binomial, data = .d)
  expect_equal(.fitP(.x, model = "logitQuadratic"), unname(stats::fitted(.quad)),
               tolerance = 1e-5)
})

test_that("monotone logit fits allow flat fits and follow the requested direction", {
  .down <- c(20, 16, 12, 8, 4)
  .pooled <- rep(sum(.down) / sum(.n), 5)
  expect_equal(.fitP(.down, model = "logitLinear", monotone = "increasing"), .pooled,
               tolerance = 1e-5)
  expect_equal(.fitP(.down, model = "logitLinear", monotone = "decreasing"),
               .fitP(.down, model = "logitLinear"), tolerance = 1e-5)
  # a peaked response: the constrained quadratic is nondecreasing over the range
  .peak <- c(5, 15, 22, 15, 8)
  .q <- .fitP(.peak, model = "logitQuadratic", monotone = "increasing")
  expect_true(all(diff(.q) >= -1e-8))
  .qd <- .fitP(.peak, model = "logitQuadratic", monotone = "decreasing")
  expect_true(all(diff(.qd) <= 1e-8))
  # the CUI-MET formulations: exponentiated slope and derivative penalty
  .pl <- .fitP(.down, model = "logitLinear", monotone = "increasing", penalty = TRUE)
  expect_true(all(diff(.pl) >= 0))
  expect_equal(.pl, .pooled, tolerance = 1e-3)
  .pq <- .fitP(.peak, model = "logitQuadratic", monotone = "increasing", penalty = TRUE)
  expect_true(all(diff(.pq) >= -1e-4))
})

test_that("Emax and exponential fits match DoseFinding", {
  skip_if_not_installed("DoseFinding")
  .xe <- c(3, 8, 14, 17, 18)
  expect_equal(.fitP(.xe, model = "emax"), .fitP(.xe, model = "emax", backend = "DoseFinding"),
               tolerance = 1e-3)
  .xx <- c(2, 3, 5, 9, 18)
  expect_equal(.fitP(.xx, model = "exponential"),
               .fitP(.xx, model = "exponential", backend = "DoseFinding"), tolerance = 1e-3)
  # the fit is the curve through the per-dose log-odds
  .f <- .cur$.cuiDrModel(cuiBinary(e, model = "emax"), .xe, .n, .doses)$fit
  expect_equal(.f$model, "emax")
  expect_length(.f$coef, 3)
  expect_equal(.cur$.cuiDrPredict(.f, .doses), .fitP(.xe, model = "emax"))
  # a dose scale does not change the fit (multiplicative scaling)
  .f10 <- .cur$.cuiDrModel(cuiBinary(e, model = "emax"), .xe, .n, .doses * 10)
  expect_equal(.f10$p, .fitP(.xe, model = "emax"), tolerance = 1e-6)
})

test_that("monotone Emax and exponential fits constrain the sign of the effect", {
  .down <- c(20, 16, 12, 8, 4)
  for (.m in c("emax", "exponential")) {
    .p <- .fitP(.down, model = .m, monotone = "increasing")
    expect_equal(max(.p) - min(.p), 0, tolerance = 1e-8)
    .pd <- .fitP(.down, model = .m, monotone = "decreasing")
    expect_true(all(diff(.pd) <= 1e-8))
  }
})

test_that("dose groups with all or no events", {
  .x <- c(0, 6, 12, 16, 18)
  .g <- .fitP(.x, model = "emax")
  expect_true(all(is.finite(.g)))
  # the boundary group gets almost no weight: the other doses fit as without it
  .without <- .cur$.cuiDrModel(cuiBinary(e, model = "emax"), .x[-1], .n[-1], .doses[-1])$fit
  .rescaled <- .without
  expect_equal(.g[-1], .cur$.cuiDrPredict(.rescaled, .doses[-1]), tolerance = 1e-3)
  expect_error(.fitP(.x, model = "emax", boundary = "error"), "all or no events")
  .c <- .fitP(.x, model = "emax", boundary = "correction")
  expect_gt(abs(.c[1] - .g[1]), 1e-4)
  expect_true(all(is.finite(.fitP(c(0, 0, 0, 30, 30), model = "logitLinear"))))
})

test_that("options are validated and printed", {
  expect_error(cuiBinary(e, model = "emax", monotone = "increasing", backend = "DoseFinding"),
               "monotone")
  expect_error(cuiBinary(e, model = "nope"))
  expect_error(cuiBinary(e, penalty = NA))
  expect_output(print(cuiBinary(tox, model = "emax", monotone = "increasing")),
                "emax model \\(increasing\\)")
})

# CUI-MET (Zhang et al. 2025) Table 1 models with Table 2 probabilities
.cuimetSim <- function(p, n, seed) {
  withr::with_seed(seed, {
    .d <- data.frame(dose = rep(1:5, each = n))
    for (.v in names(p)) .d[[.v]] <- stats::rbinom(nrow(.d), 1, rep(p[[.v]], each = n))
    .d
  })
}

test_that("CUI-MET example 3: models, decreasing tolerability and recommendations", {
  .p <- list(tox = c(0.105, 0.152, 0.233, 0.374, 0.583),
             eff = c(0.223, 0.387, 0.470, 0.442, 0.312),
             tol = c(0.705, 0.607, 0.500, 0.393, 0.295))
  .d <- .cuimetSim(.p, 3000, 7)
  .eps <- list(cuiBinary(tox, model = "exponential", direction = "lower"),
               cuiBinary(eff, model = "logitQuadratic"),
               cuiBinary(tol, model = "logitLinear"))
  .r <- cuiDoseResponse(.d, .eps, dose = dose, nBoot = 0,
                        weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
  expect_equal(.r$summary$est, c(0.502, 0.545, 0.538, 0.464, 0.328), tolerance = 0.02)
  expect_equal(.r$summary$regimen[.r$summary$best], "2")
  .w2 <- update(.r, weights = c(tox = 0.4, eff = 0.2, tol = 0.4))
  expect_equal(.w2$summary$regimen[.w2$summary$best], "1")
  # tolerability falls with dose: a decreasing constraint changes nothing
  .tolDec <- cuiDoseResponse(.d, cuiBinary(tol, model = "logitLinear", monotone = "decreasing"),
                             dose = dose, nBoot = 0)
  .tolFree <- cuiDoseResponse(.d, cuiBinary(tol, model = "logitLinear"), dose = dose, nBoot = 0)
  expect_equal(.tolDec$summary$est, .tolFree$summary$est, tolerance = 1e-5)
  expect_equal(names(.r$doseResponse), c("tox", "eff", "tol"))
  expect_equal(.r$doseResponse$eff$model, "logitQuadratic")
})

test_that("CUI-MET example 1 recommendations with modeled endpoints and a bootstrap", {
  .p <- list(tox = c(0.038, 0.058, 0.106, 0.232, 0.533),
             eff = c(0.022, 0.119, 0.343, 0.568, 0.681),
             tol = c(0.147, 0.206, 0.280, 0.368, 0.466))
  .eps <- list(cuiBinary(tox, model = "exponential", direction = "lower"),
               cuiBinary(eff, model = "logitQuadratic"),
               cuiBinary(tol, model = "logitLinear"))
  .big <- cuiDoseResponse(.cuimetSim(.p, 3000, 3), .eps, dose = dose, nBoot = 0,
                          weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
  expect_equal(.big$summary$regimen[.big$summary$best], "5")
  expect_equal(update(.big, weights = c(tox = 0.3, eff = 0.2, tol = 0.5))$summary$regimen[
    update(.big, weights = c(tox = 0.3, eff = 0.2, tol = 0.5))$summary$best], "4")
  # 30 per dose as in the paper: zero-event groups are common in bootstraps
  .small <- cuiDoseResponse(.cuimetSim(.p, 30, 4), .eps, dose = dose, nBoot = 50, seed = 1,
                            weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
  expect_equal(sum(.small$draws$status == "ok"), 51)
  expect_equal(sum(.small$summary$selFreq), 1)
})

test_that("logit fits match glm() for doses far from 0", {
  .n100 <- rep(100, 3)
  for (.dd in list(c(80, 90, 100), c(200, 250, 300))) {
    .x <- c(10, 40, 15)
    .d <- data.frame(dose = .dd, x = .x, n = .n100)
    .quad <- stats::glm(cbind(x, n - x) ~ dose + I(dose^2), binomial, data = .d)
    .f <- .cur$.cuiDrModel(cuiBinary(e, model = "logitQuadratic"), .x, .n100, .dd)
    expect_equal(.f$p, unname(stats::fitted(.quad)), tolerance = 1e-5)
    expect_equal(.cur$.cuiDrPredict(.f$fit, .dd), .f$p)
  }
  .dd <- c(98, 99, 100)
  .x <- c(10, 30, 60)
  .lin <- stats::glm(cbind(x, n - x) ~ dose, binomial,
                     data = data.frame(dose = .dd, x = .x, n = .n100))
  expect_equal(.cur$.cuiDrModel(cuiBinary(e, model = "logitLinear"), .x, .n100, .dd)$p,
               unname(stats::fitted(.lin)), tolerance = 1e-5)
})

test_that("too few dose levels and negative doses are errors", {
  .one <- data.frame(dose = 5, eff = c(0, 1, 1))
  expect_equal(cuiDoseResponse(.one, cuiBinary(eff), dose = dose, nBoot = 0)$summary$est, 2 / 3)
  expect_error(cuiDoseResponse(.one, cuiBinary(eff, model = "logitLinear"), dose = dose,
                               nBoot = 0), "more dose levels")
  .two <- data.frame(dose = rep(c(1, 2), each = 3), eff = c(0, 1, 1, 1, 1, 0))
  expect_error(cuiDoseResponse(.two, cuiBinary(eff, model = "logitQuadratic", monotone = "increasing"),
                               dose = dose, nBoot = 0), "logitQuadratic")
  .neg <- data.frame(dose = rep(c(-1, 0, 1), each = 3), eff = rep(c(0, 1, 1), 3))
  expect_error(cuiDoseResponse(.neg, cuiBinary(eff, model = "emax"), dose = dose, nBoot = 0),
               "nonnegative")
})
