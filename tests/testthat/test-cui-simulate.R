skip_on_cran()

# caches default to the working directory; keep them out of the package
withr::local_dir(withr::local_tempdir())

.cur <- loadNamespace("nlmixr2extra")

# one-compartment oral model with weight on clearance and two probability outputs
.cuiModel <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl) * (WT / 70)^0.75
    v <- exp(tv + eta.v)
    cp <- linCmt()
    peff <- expit(-2 + 0.8 * cp)
    pae <- expit(-6 + 0.4 * cp)
    cp ~ add(add.sd)
  })
}

# same structure without random effects or covariates, for analytic checks
.cuiFixedModel <- function() {
  ini({
    tka <- log(1.5)
    tcl <- log(3)
    tv <- log(30)
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

.ui <- suppressMessages(rxode2::rxode2(.cuiModel))
.uiFixed <- suppressMessages(rxode2::rxode2(.cuiFixedModel))
.pop <- cuiPopulation(nlmixr2data::theo_sd, n = 100, covs = "WT", seed = 1)

.conc <- function(t, dose, ka, k, v) dose * ka / (v * (ka - k)) * (exp(-k * t) - exp(-ka * t))
.auc <- function(t, dose, ka, k, v) {
  dose * ka / (v * (ka - k)) * ((1 - exp(-k * t)) / k - (1 - exp(-ka * t)) / ka)
}

test_that("supplied thetas, etas and covariates reproduce the analytic solution", {
  .mg <- .cur$.cuiModelGroup(.ui, list(cuiEndpoint(eff, .ui, cuiProb(peff, time = 4))), 1L)
  .eta <- matrix(c(0, 0.5, -0.3, 0, 0.2, -0.2, 0, 0.1, 0), nrow = 3,
                 dimnames = list(NULL, c("eta.ka", "eta.cl", "eta.v")))
  # a discrete covariate would make rxode2 regroup subjects if passed as iCov
  .cov <- data.frame(WT = c(70, 50, 70))
  .reg <- cuiRegimens(dose = c(100, 300))
  .sim <- .cur$.cuiSolveDraw(.mg, .mg$theta, .eta, .cov, .reg, c(2, 4), NULL)
  .th <- .mg$theta
  for (.r in 1:2) for (.i in 1:3) {
    .ka <- exp(.th[["tka"]] + .eta[.i, 1])
    .cl <- exp(.th[["tcl"]] + .eta[.i, 2]) * (.cov$WT[.i] / 70)^0.75
    .v <- exp(.th[["tv"]] + .eta[.i, 3])
    .row <- .sim[.sim$.reg == .r & .sim$.sub == .i & .sim$time == 4, ]
    expect_equal(.row$cp, .conc(4, .reg$dose[.r], .ka, .cl / .v, .v), tolerance = 1e-6,
                 ignore_attr = TRUE)
  }
})

test_that("window averages and integrals match the analytic multiple-dose solution", {
  .ka <- 1.5
  .k <- 3 / 30
  .doseTimes <- c(0, 12, 24, 36)
  .aucTo <- function(T) sum(.auc(pmax(T - .doseTimes, 0), 100, .ka, .k, 30))
  .r <- cui(list(cuiEndpoint(cavg, .uiFixed, cuiMean(cp, window = c(0, 48)), utility = cuiLinear(0, 100)),
                 cuiEndpoint(auc, .uiFixed, cuiAuc(cp, window = c(12, 36)), utility = cuiLinear(0, 1e4)),
                 cuiEndpoint(c0, .uiFixed, cuiMean(cp, time = 0), utility = cuiLinear(0, 100))),
            regimens = cuiRegimens(dose = 100, ii = 12, duration = 48), population = 3)
  .v <- .r$values
  expect_equal(.v$value[.v$endpoint == "cavg"], .aucTo(48) / 48, tolerance = 1e-4)
  expect_equal(.v$value[.v$endpoint == "auc"], .aucTo(36) - .aucTo(12), tolerance = 1e-4)
  # an observation at time 0 is matched (absolute tolerance)
  expect_equal(.v$value[.v$endpoint == "c0"], 0)
})

test_that("summaries take names, strings, expressions and constants", {
  .reg <- cuiRegimens(dose = c(100, 400))
  .run <- function(...) {
    cui(list(...), regimens = .reg, population = .pop, seed = 1, mcTol = 1)$values
  }
  .a <- .run(cuiEndpoint(e, .ui, cuiProb(peff, time = 4)))
  .b <- .run(cuiEndpoint(e, .ui, cuiProb("peff", time = 4)))
  .col <- "peff"
  .c <- .run(cuiEndpoint(e, .ui, cuiProb(.x = .col, time = 4)))
  .wrap <- function(v) cuiEndpoint(e, .ui, cuiProb(.x = v, time = 4))
  .d <- .run(.wrap("peff"))
  expect_equal(.a, .b)
  expect_equal(.a, .c)
  expect_equal(.a, .d)
  # expressions, strings parsed as expressions, function calls and constants
  cut <- 5
  .t1 <- .run(cuiEndpoint(e, .ui, cuiThreshold(cp > cut, time = 2)))
  .t2 <- .run(cuiEndpoint(e, .ui, cuiThreshold("cp > 5", time = 2)))
  .t3 <- .run(cuiEndpoint(e, .ui, cuiThreshold(abs(cp - 0) > .env$cut, time = 2)))
  expect_equal(.t1, .t2)
  expect_equal(.t1, .t3)
  # .env$ is frozen when the summary is made; plain constants when cui() runs
  .frozen <- cuiThreshold(cp > .env$cut, time = 2)
  .late <- cuiThreshold(cp > cut, time = 2)
  cut <- 1000
  expect_equal(.run(cuiEndpoint(e, .ui, .frozen)), .t1)
  expect_true(all(.run(cuiEndpoint(e, .ui, .late))$value == 0))
  # a variable named like a model column: the column wins, with a warning
  cl <- 1
  expect_warning(.run(cuiEndpoint(e, .ui, cuiThreshold(cl > 1, time = 2))), "model column")
  expect_no_warning(.run(cuiEndpoint(e, .ui, cuiThreshold(.data$cl > .env$cl, time = 2))))
  expect_error(.run(cuiEndpoint(e, .ui, cuiProb(nope, time = 2))), "not a model column")
  expect_error(.run(cuiEndpoint(e, .ui, cuiThreshold(cp > nope2, time = 2))), "not a model column or a constant")
  expect_error(.run(cuiEndpoint(e, .ui, cuiThreshold(cp, time = 2))), "TRUE/FALSE")
  expect_error(.run(cuiEndpoint(e, .ui, cuiProb(cp, time = 2))), "\\[0, 1\\]")
  expect_error(cuiProb(peff, 4), "by name")
  expect_error(cuiProb(peff, .x = "peff"), "only one")
  expect_error(cuiProb(), "required")
  expect_error(cuiMean(cp, time = 1, window = c(0, 1)), "not both")
  expect_error(cuiAuc(cp), "window")
  expect_error(cuiAuc(cp, window = c(2, 1)), "from < to")
  expect_equal(format(cuiMean(cp, window = c(0, 24))), "mean(cp) over [0, 24]")
})

test_that("every regimen and run sees the same virtual subjects", {
  .eps <- list(cuiEndpoint(e, .ui, cuiProb(peff, time = 4)))
  .a <- cui(.eps, regimens = cuiRegimens(dose = c(100, 100, 400), label = c("x", "y", "z")),
            population = .pop, seed = 3, mcTol = 1)
  .s <- .a$summary
  expect_equal(.s$est[.s$regimen == "x"], .s$est[.s$regimen == "y"])
  # regimens do not change each other's subjects
  .b <- cui(.eps, regimens = cuiRegimens(dose = c(400, 100), label = c("z", "x")),
            population = .pop, seed = 3, mcTol = 1)
  expect_equal(.b$summary$est[.b$summary$regimen == "x"], .s$est[.s$regimen == "x"])
  # identical regimens give a zero Monte Carlo difference
  .c <- cui(.eps, regimens = cuiRegimens(dose = c(100, 100), label = c("x", "y")),
            population = .pop, seed = 3)
  expect_equal(.c$mcse$se, 0)
})

test_that("endpoints sharing a fit are simulated once per draw", {
  .ui2 <- suppressMessages(rxode2::rxode2(.cuiFixedModel))
  .calls <- 0
  .real <- .cur$.cuiSolveDraw
  local_mocked_bindings(.cuiSolveDraw = function(...) {
    .calls <<- .calls + 1
    .real(...)
  })
  .r <- cui(list(cuiEndpoint(a, .ui, cuiProb(peff, time = 4)),
                 cuiEndpoint(b, .ui, cuiProb(pae, time = 2), direction = "lower"),
                 cuiEndpoint(c, .ui2, cuiMean(cp, time = 2), utility = cuiLinear(0, 10))),
            regimens = cuiRegimens(dose = c(100, 400)), population = .pop, seed = 1, mcTol = 1)
  expect_equal(.calls, 2)
})

test_that("trial mode resamples subjects per trial around the expected point estimate", {
  .eps <- list(cuiEndpoint(eff, .ui, cuiProb(peff, time = 4)),
               cuiEndpoint(ae, .ui, cuiProb(pae, time = 2), direction = "lower"))
  .reg <- cuiRegimens(dose = c(100, 200, 400, 800))
  .t <- cui(.eps, regimens = .reg, population = .pop, mode = "trial", nDraw = 10, nSub = 30,
            seed = 1, mcTol = 1)
  .e <- cui(.eps, regimens = .reg, population = .pop, seed = 1, mcTol = 1)
  expect_equal(.t$summary$est, .e$summary$est)
  expect_equal(sum(.t$summary$selFreq), 1)
  expect_equal(sort(unique(.t$cui$draw)), 0:10)
  expect_output(print(.t), "10 simulated trials of 30 subjects")
  .t2 <- cui(.eps, regimens = .reg, population = .pop, mode = "trial", nDraw = 10, nSub = 30,
             seed = 1, mcTol = 1)
  expect_equal(.t$summary, .t2$summary)
})

test_that("failed draws are recorded and the Monte Carlo warning fires", {
  .eps <- list(cuiEndpoint(eff, .ui, cuiProb(peff, time = 4)))
  .reg <- cuiRegimens(dose = c(100, 400))
  .calls <- 0
  .real <- .cur$.cuiSolveDraw
  local_mocked_bindings(.cuiSolveDraw = function(...) {
    .calls <<- .calls + 1
    if (.calls == 3) stop("solver failed")
    .real(...)
  })
  .t <- cui(.eps, regimens = .reg, population = .pop, mode = "trial", nDraw = 10, nSub = 30,
            seed = 1, mcTol = 1, cacheDir = NA)
  expect_equal(.t$draws$status[3], "failed")
  expect_equal(.t$draws$message[3], "solver failed")
  expect_equal(sum(.t$summary$selFreq), 1)
  .calls <- 0
  expect_error(cui(.eps, regimens = .reg, population = .pop, mode = "trial", nDraw = 4,
                   nSub = 30, seed = 1, mcTol = 1, cacheDir = NA), "only 3 of 4 draws")
  .calls <- -100
  expect_warning(cui(.eps, regimens = .reg, population = .pop, seed = 1, mcTol = 1e-12),
                 "Monte Carlo standard error")
})

test_that("inputs are validated", {
  .ep <- cuiEndpoint(eff, .ui, cuiProb(peff, time = 4))
  .reg <- cuiRegimens(dose = c(100, 400))
  expect_error(cui(list(cuiEndpoint(eff)), .reg, .pop), "needs cuiEndpoint")
  expect_error(cui(.ep, .reg, 10), "covariate")
  expect_error(cui(.ep, .reg, .pop, nDraw = 3), "needs uncertainty")
  expect_error(cui(.ep, .reg, .pop, uncertainty = "cov"), "nDraw >= 1")
  expect_error(cui(.ep, .reg, .pop, uncertainty = "cov", nDraw = 2), "covariance")
  expect_error(cui(.ep, .reg, .pop, mode = "trial"), "nDraw >= 1")
  expect_error(cui(cuiEndpoint(eff, .ui, cuiProb(peff)), .reg, .pop), "times")
  expect_error(cui(.ep, .reg, .pop, times = c(a = 1, a = 2)), "unique")
  expect_error(cui(.ep, data.frame(dose = 1), .pop), "cuiRegimens")
  expect_error(cuiPopulation(nlmixr2data::theo_sd, covs = "nope"), "not in 'data'")
  expect_error(cuiPopulation(nlmixr2data::theo_sd, id = "subj"), "id column")
  expect_equal(nrow(cuiPopulation(nlmixr2data::theo_sd)), 12)
})

test_that("assessment times and the result object", {
  .reg <- cuiRegimens(dose = c(100, 400), duration = 24)
  .r <- cui(list(cuiEndpoint(eff, .ui, cuiProb(peff)),
                 cuiEndpoint(ae, .ui, cuiProb(pae, time = 2), direction = "lower")),
            regimens = .reg, population = .pop, times = c(early = 2, late = 12), seed = 1,
            mcTol = 1)
  expect_equal(.r$assessments, c("early", "late"))
  .v <- .r$values
  # a summary with its own time is the same at every assessment
  expect_equal(.v$value[.v$endpoint == "ae" & .v$assessment == "early"],
               .v$value[.v$endpoint == "ae" & .v$assessment == "late"])
  expect_false(isTRUE(all.equal(.v$value[.v$endpoint == "eff" & .v$assessment == "early"],
                                .v$value[.v$endpoint == "eff" & .v$assessment == "late"])))
  # without times, the common duration is the assessment
  .d <- cui(cuiEndpoint(eff, .ui, cuiProb(peff)), regimens = .reg, population = .pop,
            seed = 1, mcTol = 1)
  .d24 <- cui(cuiEndpoint(eff, .ui, cuiProb(peff)), regimens = .reg, population = .pop,
              times = c(all = 24), seed = 1, mcTol = 1)
  expect_equal(.d$values, .d24$values)
  expect_output(print(.r), "prob\\(peff\\)")
  # no environments or fits are kept
  .envs <- function(x) {
    if (is.environment(x)) return(TRUE)
    if (is.list(x)) return(any(vapply(x, .envs, logical(1))))
    FALSE
  }
  expect_false(.envs(unclass(.r)))
})

test_that("parameter uncertainty from a fit's covariance", {
  skip_if_not_installed("nlmixr2data")
  .fitModel <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      peff <- expit(-2 + 0.8 * cp)
      cp ~ add(add.sd)
    })
  }
  # tableControl(keep =) must not leak into the simulation
  .fit <- suppressMessages(nlmixr2est::nlmixr2(.fitModel, nlmixr2data::theo_sd, "focei",
                                               control = nlmixr2est::foceiControl(print = 0),
                                               table = nlmixr2est::tableControl(keep = "WT")))
  .eps <- list(cuiEndpoint(eff, .fit, cuiProb(peff, time = 4)),
               cuiEndpoint(c2, .fit, cuiMean(cp, time = 2), utility = cuiLinear(0, 15),
                           direction = "lower"))
  .reg <- cuiRegimens(dose = c(100, 300, 600))
  expect_no_warning(.u <- cui(.eps, .reg, population = 200, uncertainty = "cov", nDraw = 20,
                              seed = 1, mcTol = 1))
  .n <- cui(.eps, .reg, population = 200, seed = 1, mcTol = 1)
  expect_equal(.u$summary$est, .n$summary$est)
  expect_equal(sum(.u$draws$status == "ok"), 21)
  expect_true(all(.u$summary$lower <= .u$summary$upper))
  expect_equal(sum(.u$summary$selFreq), 1)
  # draws differ, but only through the structural thetas
  .mg <- .cur$.cuiModelGroup(.fit, .eps, 1L)
  expect_equal(.mg$struct, c("tka", "tcl", "tv"))
  .ch <- chol(.fit$cov[.mg$struct, .mg$struct])
  .t1 <- .cur$.cuiDrawTheta(.mg, .ch, 1, 1)
  expect_equal(.t1[["add.sd"]], .mg$theta[["add.sd"]])
  expect_false(isTRUE(all.equal(.t1[["tcl"]], .mg$theta[["tcl"]])))
  expect_equal(.cur$.cuiDrawTheta(.mg, .ch, 1, 1), .t1)
})

test_that("seed = NULL follows the current random stream", {
  .ep <- cuiEndpoint(c1, .ui, cuiMean(cp, time = 4), utility = cuiLinear(0, 20))
  .reg <- cuiRegimens(dose = c(100, 200))
  set.seed(1)
  .a <- cui(.ep, .reg, .pop, mcTol = 1)
  set.seed(2)
  .b <- cui(.ep, .reg, .pop, mcTol = 1)
  set.seed(1)
  .c <- cui(.ep, .reg, .pop, mcTol = 1)
  expect_false(isTRUE(all.equal(.a$values, .b$values)))
  expect_equal(.a$values, .c$values)
})

test_that("one subject and one regimen, and zero-variance etas", {
  .r <- cui(cuiEndpoint(c1, .uiFixed, cuiMean(cp, time = 4), utility = cuiLinear(0, 20)),
            regimens = cuiRegimens(dose = 100), population = 1)
  expect_equal(nrow(.r$summary), 1)
  .zero <- function() {
    ini({
      tcl <- log(3)
      tv <- log(30)
      eta.cl ~ 0.1
      eta.v ~ fix(0)
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  .uz <- suppressMessages(rxode2::rxode2(.zero))
  expect_equal(.cur$.cuiOmegaFactor(.uz$omega)[2, ], c(0, 0))
  .z <- cui(cuiEndpoint(v, .uz, cuiMean(v, time = 1), utility = cuiLinear(0, 100)),
            regimens = cuiRegimens(dose = 100), population = 20, seed = 1, mcTol = 1)
  expect_equal(.z$summary$est, 30 / 100)
})

test_that("cui() caches by default with a seed, not without", {
  .ep <- cuiEndpoint(eff, .ui, cuiProb(peff, time = 4))
  .reg <- cuiRegimens(dose = c(100, 400))
  .r <- cui(.ep, .reg, .pop, seed = 9, mcTol = 1)
  expect_true(dir.exists(.r$settings$cacheDir))
  expect_true(startsWith(basename(.r$settings$cacheDir), "nlmixr2CuiCache_"))
  .n <- length(list.files(pattern = "^nlmixr2CuiCache_"))
  expect_null(cui(.ep, .reg, .pop, mcTol = 1)$settings$cacheDir)
  expect_null(cui(.ep, .reg, .pop, seed = 9, mcTol = 1, cacheDir = NA)$settings$cacheDir)
  expect_equal(length(list.files(pattern = "^nlmixr2CuiCache_")), .n)
})
