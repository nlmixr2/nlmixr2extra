.cur <- loadNamespace("nlmixr2extra")

# a dose-response trial with a patient covariate
.trial <- function(n = 40, seed = 1) {
  withr::with_seed(seed, {
    .d <- data.frame(id = seq_len(5 * n), dose = rep(c(0, 10, 25, 50, 100), each = n),
                     age = round(stats::runif(5 * n, 20, 70)),
                     site = factor(sample(c("a", "b", "c"), 5 * n, replace = TRUE)))
    .d$tox <- stats::rbinom(nrow(.d), 1, stats::plogis(-3 + 0.03 * .d$dose))
    .d$eff <- stats::rbinom(nrow(.d), 1, stats::plogis(-2 + 0.04 * .d$dose + 0.02 * (.d$age - 45)))
    .d
  })
}

test_that("glm and lm adapters predict at the estimates and at parameter draws", {
  .d <- .trial()
  .fit <- stats::glm(eff ~ dose + age, binomial, data = .d)
  .nd <- data.frame(dose = c(0, 50), age = c(30, 60))
  expect_equal(cuiModelPredict(.fit, .nd), unname(stats::predict(.fit, .nd, type = "response")))
  # a draw equal to the estimates reproduces predict()
  expect_equal(cuiModelPredict(.fit, .nd, draw = stats::coef(.fit)), cuiModelPredict(.fit, .nd))
  .dr <- withr::with_seed(1, cuiModelDraws(.fit, 2000))
  expect_length(.dr, 2000)
  expect_equal(colMeans(do.call(rbind, .dr)), stats::coef(.fit), tolerance = 0.05)
  expect_equal(cuiModelVars(.fit), c("dose", "age"))
  # factors use the fit's levels and contrasts
  .ff <- stats::glm(eff ~ dose + site, binomial, data = .d)
  .ndf <- data.frame(dose = 10, site = factor("c", levels = c("a", "b", "c")))
  expect_equal(cuiModelPredict(.ff, .ndf, draw = stats::coef(.ff)), cuiModelPredict(.ff, .ndf))
  .lm <- stats::lm(age ~ dose, data = .d)
  expect_equal(cuiModelPredict(.lm, data.frame(dose = 5), draw = stats::coef(.lm)),
               unname(stats::predict(.lm, data.frame(dose = 5))))
  expect_equal(nrow(cuiModelRefit(.fit, .d[1:100, ])$data), 100)
  expect_error(cuiModelDraws(structure(list(), class = "unknownModel"), 2), "unknownModel")
})

test_that("lme4 and glmmTMB adapters give population-level predictions", {
  skip_if_not_installed("lme4")
  .d <- .trial(60)
  .g <- suppressMessages(lme4::glmer(eff ~ dose + (1 | site), binomial, data = .d))
  .nd <- data.frame(dose = c(0, 50))
  expect_equal(cuiModelPredict(.g, .nd), cuiModelPredict(.g, .nd, draw = lme4::fixef(.g)),
               tolerance = 1e-8)
  expect_equal(cuiModelVars(.g), "dose")
  expect_length(withr::with_seed(1, cuiModelDraws(.g, 3)), 3)
  skip_if_not_installed("glmmTMB")
  .t <- glmmTMB::glmmTMB(eff ~ dose + (1 | site), binomial, data = .d)
  expect_equal(cuiModelPredict(.t, .nd), cuiModelPredict(.t, .nd, draw = glmmTMB::fixef(.t)$cond),
               tolerance = 1e-6)
  expect_equal(cuiModelVars(.t), "dose")
  expect_length(withr::with_seed(1, cuiModelDraws(.t, 3)), 3)
})

test_that("a fitted glm in cuiDoseResponse() equals the built-in logit model", {
  .d <- .trial()
  .fit <- stats::glm(tox ~ dose, binomial, data = .d)
  .a <- cuiDoseResponse(.d, cuiBinary(tox, model = .fit, direction = "lower"), dose = dose,
                        nBoot = 20, seed = 1)
  .b <- cuiDoseResponse(.d, cuiBinary(tox, model = "logitLinear", direction = "lower"),
                        dose = dose, nBoot = 20, seed = 1)
  # same model, same resamples: bootstrap refits of the glm agree too
  expect_equal(.a$cui$cui, .b$cui$cui, tolerance = 1e-4)
  expect_output(print(cuiBinary(tox, model = .fit)), "fitted glm")
  expect_error(cuiBinary(tox, model = .fit, monotone = "increasing"), "built-in")
})

test_that("a fitted model with covariates is averaged over the patients", {
  .d <- .trial()
  .fit <- stats::glm(eff ~ dose + age, binomial, data = .d)
  .r <- cuiDoseResponse(.d, cuiBinary(eff, model = .fit), dose = dose, nBoot = 0)
  .want <- vapply(c(0, 10, 25, 50, 100), function(dd) {
    mean(stats::predict(.fit, transform(.d, dose = dd), type = "response"))
  }, numeric(1))
  expect_equal(.r$summary$est, .want)
  .noDose <- stats::glm(eff ~ age, binomial, data = .d)
  expect_error(cuiDoseResponse(.d, cuiBinary(eff, model = .noDose), dose = dose, nBoot = 0),
               "does not use the dose column")
  .gauss <- stats::glm(age ~ dose, gaussian, data = .d)
  expect_error(cuiDoseResponse(.d, cuiBinary(eff, model = .gauss), dose = dose, nBoot = 0),
               "probabilities")
})

test_that("other model classes work once they have methods", {
  # a toy model class: probability = plogis(a + b * dose)
  .toy <- structure(list(a = -1, b = 0.02), class = "cuiToyModel")
  registerS3method("cuiModelPredict", "cuiToyModel", function(fit, newdata, draw = NULL, ...) {
    .ab <- if (is.null(draw)) c(fit$a, fit$b) else draw
    stats::plogis(.ab[1] + .ab[2] * newdata$dose)
  }, envir = asNamespace("nlmixr2extra"))
  registerS3method("cuiModelVars", "cuiToyModel", function(fit, ...) "dose",
                   envir = asNamespace("nlmixr2extra"))
  registerS3method("cuiModelDraws", "cuiToyModel", function(fit, n, ...) {
    lapply(seq_len(n), function(i) c(fit$a, fit$b) + stats::rnorm(2, sd = 0.01))
  }, envir = asNamespace("nlmixr2extra"))
  registerS3method("cuiModelRefit", "cuiToyModel", function(fit, data, ...) fit,
                   envir = asNamespace("nlmixr2extra"))
  .d <- .trial(10)
  .r <- cuiDoseResponse(.d, cuiBinary(eff, model = .toy), dose = dose, nBoot = 3, seed = 1)
  expect_equal(.r$summary$est, stats::plogis(-1 + 0.02 * c(0, 10, 25, 50, 100)))
  .c <- cui(cuiEndpoint(eff, .toy, cuiResponse(dose = dose)), cuiRegimens(dose = c(0, 50)),
            population = 10, uncertainty = "cov", nDraw = 5, seed = 1, mcTol = 1)
  expect_equal(.c$summary$est, stats::plogis(-1 + 0.02 * c(0, 50)))
  expect_equal(sum(.c$draws$status == "ok"), 6)
})

test_that("fitted models in cui(): population averaging, draws and chained exposure", {
  skip_on_cran()
  .d <- .trial()
  .fit <- stats::glm(eff ~ dose + age, binomial, data = .d)
  .pop <- cuiPopulation(.d, n = 50, covs = "age", id = "id", seed = 2)
  .reg <- cuiRegimens(dose = c(0, 25, 100))
  .r <- cui(cuiEndpoint(eff, .fit, cuiResponse(dose = dose)), .reg, .pop, seed = 1, mcTol = 1)
  .want <- vapply(c(0, 25, 100), function(dd) {
    mean(stats::predict(.fit, data.frame(dose = dd, age = .pop$age), type = "response"))
  }, numeric(1))
  expect_equal(.r$summary$est, .want)
  expect_output(print(.r), "response\\(dose = dose\\)")
  .u <- cui(cuiEndpoint(eff, .fit, cuiResponse(dose = "dose")), .reg, .pop,
            uncertainty = "cov", nDraw = 20, seed = 1, mcTol = 1)
  expect_equal(.u$summary$est, .want)
  expect_true(all(.u$summary$lower < .u$summary$upper))
  expect_error(cui(cuiEndpoint(eff, .fit, cuiProb(peff, time = 1)), .reg, .pop), "cuiResponse")
  expect_error(cui(cuiEndpoint(eff, .fit, cuiResponse(dose = dose)), .reg, 10),
               "not in 'population': age")
  .nod <- stats::glm(eff ~ age, binomial, data = .d)
  expect_error(cui(cuiEndpoint(eff, .nod, cuiResponse(dose = dose)), .reg, .pop),
               "does not use the dose column")
  # a fitted exposure-response model of a chained Cavg
  .pkModel <- function() {
    ini({
      tka <- log(1.5)
      tcl <- log(3)
      tv <- log(30)
      eta.cl ~ 0.2
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  .pk <- suppressMessages(rxode2::rxode2(.pkModel))
  .er <- data.frame(CAVG = seq(0, 20, length.out = 50))
  .er$y <- 2 * .er$CAVG
  .lmEr <- stats::lm(y ~ CAVG, data = .er)
  expect_equal(format(cuiResponse()), "response()")
  .ch <- cui(list(cuiEndpoint(lin, .lmEr, cuiResponse(), utility = cuiLinear(0, 50)),
                  cuiEndpoint(cavg, .pk, cuiMean(cp, window = c(0, 24)), utility = cuiLinear(0, 50))),
             cuiRegimens(dose = c(200, 400), ii = 12), population = 100,
             exposure = cuiExposure(.pk, CAVG = cuiMean(cp, window = c(0, 24))), seed = 1, mcTol = 1)
  .v <- .ch$values
  expect_equal(.v$value[.v$endpoint == "lin"], 2 * .v$value[.v$endpoint == "cavg"],
               tolerance = 1e-8)
  expect_error(cui(cuiEndpoint(eff, .fit, cuiResponse(dose = dose)), .reg, .pop,
                   exposure = cuiExposure(.fit, X = cuiMean(cp, time = 1))), "nlmixr2 fit or an rxUi")
})

test_that("brms fits use posterior draws", {
  skip_on_cran()
  skip_if_not_installed("brms")
  .d <- .trial(30)
  .b <- suppressMessages(suppressWarnings(
    brms::brm(eff ~ dose, data = .d, family = brms::bernoulli(), chains = 1, iter = 600,
              refresh = 0, seed = 1, silent = 2)))
  .nd <- data.frame(dose = c(0, 100))
  .mean <- cuiModelPredict(.b, .nd)
  expect_equal(.mean, colMeans(brms::posterior_epred(.b, newdata = .nd)))
  .ids <- withr::with_seed(1, cuiModelDraws(.b, 4))
  expect_length(.ids, 4)
  expect_equal(cuiModelPredict(.b, .nd, draw = .ids[[1]]),
               as.numeric(brms::posterior_epred(.b, newdata = .nd, draw_ids = .ids[[1]])))
  expect_equal(cuiModelVars(.b), "dose")
  .c <- cui(cuiEndpoint(eff, .b, cuiResponse(dose = dose)), cuiRegimens(dose = c(0, 50, 100)),
            population = 20, uncertainty = "cov", nDraw = 10, seed = 1, mcTol = 1)
  expect_equal(.c$summary$est, colMeans(brms::posterior_epred(.b, newdata = data.frame(dose = c(0, 50, 100)))))
  .r <- cuiDoseResponse(.d, cuiBinary(eff, model = .b), dose = dose, nBoot = 0)
  expect_equal(.r$summary$est, colMeans(brms::posterior_epred(.b, newdata = data.frame(dose = c(0, 10, 25, 50, 100)))))
  # random-effect terms are not needed for population-level predictions
  expect_equal(all.vars(.cur$.cuiNoBars(y ~ dose + (1 | site) + age)), c("y", "dose", "age"))
})

test_that("mixed-model draws keep transformed terms and factor levels of the fit", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("glmmTMB")
  .d <- .trial(60)
  .d$grp <- factor(rep(1:8, length.out = nrow(.d)))
  .nd <- data.frame(dose = c(0, 50, 100), age = c(30, 45, 60),
                    site = factor(c("b", "c", "c"), levels = c("b", "c")))
  .g <- suppressMessages(lme4::glmer(eff ~ poly(dose, 2) + scale(age) + site + (1 | grp),
                                     binomial, data = .d))
  expect_equal(cuiModelPredict(.g, .nd, draw = lme4::fixef(.g)), cuiModelPredict(.g, .nd),
               tolerance = 1e-6)
  .t <- glmmTMB::glmmTMB(eff ~ poly(dose, 2) + scale(age) + site + (1 | grp), binomial,
                         data = .d)
  expect_equal(cuiModelPredict(.t, .nd, draw = glmmTMB::fixef(.t)$cond),
               cuiModelPredict(.t, .nd), tolerance = 1e-6)
  # a single level of a factor covariate
  .one <- transform(.nd, site = factor("c"))
  expect_equal(cuiModelPredict(.g, .one, draw = lme4::fixef(.g)), cuiModelPredict(.g, .one),
               tolerance = 1e-6)
  expect_equal(cuiModelPredict(.t, .one, draw = glmmTMB::fixef(.t)$cond),
               cuiModelPredict(.t, .one), tolerance = 1e-6)
})

test_that("refits work for models built inside functions, offsets are handled", {
  .d <- .trial()
  .make <- function(dd) {
    fam <- stats::binomial
    form <- tox ~ dose
    stats::glm(form, family = fam, data = dd)
  }
  .r <- cuiDoseResponse(.d, cuiBinary(tox, model = .make(.d)), dose = dose, nBoot = 5, seed = 1)
  expect_equal(sum(.r$draws$status == "ok"), 6)
  .d$t <- 1 + .d$age / 100
  .off <- stats::glm(tox ~ dose, poisson, data = .d, offset = log(t))
  expect_true("t" %in% cuiModelVars(.off))
  expect_error(cuiModelPredict(.off, data.frame(dose = 1, t = 1), draw = stats::coef(.off)),
               "offset")
})

test_that("a fit shared by endpoints with and without a dose column needs every input", {
  .d <- .trial()
  .fit <- stats::glm(eff ~ dose, binomial, data = .d)
  expect_error(cui(list(cuiEndpoint(a, .fit, cuiResponse(dose = dose)),
                        cuiEndpoint(b, .fit, cuiResponse())),
                   cuiRegimens(dose = c(0, 50)), 10, seed = 1),
               "not in 'population': dose")
})
