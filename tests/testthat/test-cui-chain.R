skip_on_cran()

.cur <- loadNamespace("nlmixr2extra")

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

# exposure-response model driven by a chained Cavg, with no compartments
.erModel <- function() {
  ini({
    te0 <- -2
    tslope <- 0.5
    add.sd <- 0.1
  })
  model({
    lin <- 2 * CAVG
    peff <- expit(te0 + tslope * CAVG)
    lin ~ add(add.sd)
  })
}

.pk <- suppressMessages(rxode2::rxode2(.pkModel))
.er <- suppressMessages(rxode2::rxode2(.erModel))

test_that("sampling() records the original ids of each resample", {
  .d <- nlmixr2data::theo_sd
  .s <- withr::with_seed(1, .cur$sampling(.d))
  .ids <- attr(.s, "bootIds")
  expect_length(.ids, 12)
  .wtNew <- .s$WT[!duplicated(.s$ID)]
  .wtOrig <- .d$WT[match(.ids, .d$ID)]
  expect_equal(.wtNew, .wtOrig)
  .d$grp <- ifelse(.d$ID <= 6, "a", "b")
  .st <- withr::with_seed(2, .cur$sampling(.d, performStrat = TRUE, stratVar = "grp"))
  .ids2 <- attr(.st, "bootIds")
  expect_equal(.st$WT[!duplicated(.st$ID)], .d$WT[match(.ids2, .d$ID)])
})

test_that("a chained exposure drives exposure-response models per subject", {
  .calls <- 0
  .real <- .cur$.cuiSolveDraw
  local_mocked_bindings(.cuiSolveDraw = function(...) {
    .calls <<- .calls + 1
    .real(...)
  })
  .r <- cui(list(cuiEndpoint(lin, .er, cuiMean(lin, time = 1), utility = cuiLinear(0, 50)),
                 cuiEndpoint(eff, .er, cuiProb(peff, time = 1)),
                 cuiEndpoint(cavg, .pk, cuiMean(cp, window = c(0, 24)), utility = cuiLinear(0, 50))),
            regimens = cuiRegimens(dose = c(200, 400), ii = 12),
            population = 200, exposure = cuiExposure(.pk, CAVG = cuiMean(cp, window = c(0, 24))),
            seed = 1, mcTol = 1)
  # the PK model is simulated once and feeds both its own endpoint and the exposure
  expect_equal(.calls, 2)
  .v <- .r$values
  .get <- function(ep, reg) .v$value[.v$endpoint == ep & .v$regimen == reg]
  expect_equal(.get("lin", "200 q12"), 2 * .get("cavg", "200 q12"))
  # linear PK: exposure doubles with the dose
  expect_equal(.get("lin", "400 q12"), 2 * .get("lin", "200 q12"))
  # predictions are made per subject, not at the mean exposure
  .atMean <- plogis(-2 + 0.5 * .get("cavg", "200 q12"))
  expect_gt(abs(.get("eff", "200 q12") - .atMean), 1e-3)
  expect_equal(.r$settings$exposures, "CAVG")
  expect_output(print(cuiExposure(.pk, CAVG = cuiMean(cp, window = c(0, 24)))), "CAVG = mean")
})

test_that("chained exposures are validated", {
  .ep <- cuiEndpoint(eff, .er, cuiProb(peff, time = 1))
  .x <- cuiExposure(.pk, CAVG = cuiMean(cp, window = c(0, 24)))
  .reg <- cuiRegimens(dose = 100)
  expect_error(cui(.ep, .reg, 10), "CAVG")
  .pop <- structure(data.frame(id = 1:10, CAVG = 1), class = c("cuiPopulation", "data.frame"))
  expect_error(cui(.ep, .reg, .pop, exposure = .x), "also in 'population'")
  expect_error(cuiExposure(.pk, CAVG = cuiMean(cp)), "explicit 'time' or 'window'")
  expect_error(cuiExposure(.pk, cuiMean(cp, time = 1)), "named")
  expect_error(cuiExposure(.pk, CAVG = 1), "summary")
  expect_error(cuiExposure(1, CAVG = cuiMean(cp, time = 1)), "nlmixr2 fit or an rxUi")
  expect_error(cui(.ep, .reg, 10, exposure = list(.x, .x)), "unique")
  expect_error(cui(.ep, .reg, 10, exposure = "x"), "cuiExposure")
  # an exposure model that needs a chained exposure itself
  .x2 <- cuiExposure(.er, LIN = cuiMean(lin, time = 1))
  expect_error(cui(.ep, .reg, 10, exposure = list(.x, .x2)), "cannot itself use")
})

test_that("bootstrap ensembles are paired by resample or used independently", {
  .item <- function(r) list(ok = TRUE, resample = r)
  .e1 <- list(items = lapply(c("a", "b", "c", "d"), .item))
  .e2 <- list(items = lapply(c("c", "a", "x", "b"), .item))
  .p <- .cur$.cuiBootPair(list(.e1, .e2))
  expect_true(.p$paired)
  expect_equal(.p$draws, matrix(c(1L, 2L, 3L, 2L, 4L, 1L), ncol = 2))
  .e3 <- list(items = lapply(c("p", "q", "r"), .item))
  .q <- .cur$.cuiBootPair(list(.e1, .e3))
  expect_false(.q$paired)
  expect_equal(.q$draws, matrix(rep(1:3, 2), ncol = 2))
  .single <- .cur$.cuiBootPair(list(.e1))
  expect_equal(.single$draws, matrix(1:4, ncol = 1))
})

test_that("bootstrap refits give the parameter draws", {
  withr::local_dir(withr::local_tempdir())
  .mod <- function() {
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
  .fit <- suppressMessages(nlmixr2est::nlmixr2(.mod, nlmixr2data::theo_sd, "focei",
                                               control = nlmixr2est::foceiControl(print = 0)))
  .boot <- suppressMessages(suppressWarnings(
    withr::with_seed(1, bootstrapFit(.fit, nboot = 3, returnType = "fitList"))))
  for (.i in 1:3) {
    .info <- .boot[[.i]]$bootstrapResample
    expect_equal(.info$parentMd5, .fit$bootstrapMd5)
    expect_length(.info$ids, 12)
  }
  expect_setequal(vapply(.boot, function(f) f$bootstrapResample$index, numeric(1)), 1:3)
  .eps <- list(cuiEndpoint(eff, .fit, cuiProb(peff, time = 4)))
  .reg <- cuiRegimens(dose = c(100, 300))
  .b <- cui(.eps, .reg, 100, uncertainty = cuiBootstrap(.boot), seed = 1, mcTol = 1)
  expect_equal(sum(.b$draws$status == "ok"), 4)
  expect_output(print(.b), "3 bootstrap parameter draws")
  # draw b simulates with bootstrap fit b (same virtual subjects)
  .ord <- order(vapply(.boot, function(f) f$bootstrapResample$index, numeric(1)))
  .b2 <- cui(list(cuiEndpoint(eff, .boot[[.ord[2]]], cuiProb(peff, time = 4))), .reg, 100,
             seed = 1, mcTol = 1)
  expect_equal(.b$values$value[.b$values$draw == 2], .b2$values$value)
  # point estimate unchanged by the uncertainty source
  .n <- cui(.eps, .reg, 100, seed = 1, mcTol = 1)
  expect_equal(.b$summary$est, .n$summary$est)
  expect_equal(nrow(cui(.eps, .reg, 100, uncertainty = cuiBootstrap(.boot), nDraw = 2,
                        seed = 1, mcTol = 1)$draws), 3)
  expect_error(cui(.eps, .reg, 100, uncertainty = cuiBootstrap(.boot), nDraw = 9, seed = 1),
               "only 3 bootstrap draws")
  # failed refits, the wrong return type, and ensembles of other fits
  expect_error(cuiBootstrap(c(.boot, list(NA, NA))), "only 3 of 5")
  expect_equal(cuiBootstrap(c(.boot, list(NA)), minSuccess = 0.5)$ensembles[[1]]$nFailed, 1)
  expect_error(cuiBootstrap(.fit), "returnType")
  expect_error(cui(.eps, .reg, 100, uncertainty = cuiBootstrap(.boot, .boot), seed = 1),
               "same fit")
  .other <- suppressMessages(rxode2::rxode2(.mod))
  expect_error(cui(list(cuiEndpoint(eff, .other, cuiProb(peff, time = 4))), .reg, 100,
                   uncertainty = cuiBootstrap(.boot), seed = 1), "does not belong")
  expect_output(print(cuiBootstrap(.boot)), "3 of 3 refits")
})

test_that("draws are cached, resumed and restarted", {
  withr::local_dir(withr::local_tempdir())
  .eps <- list(cuiEndpoint(c1, .pk, cuiMean(cp, time = 2), utility = cuiLinear(0, 20)))
  .reg <- cuiRegimens(dose = c(100, 200))
  .run <- function(nDraw, ...) {
    cui(.eps, .reg, 50, mode = "trial", nDraw = nDraw, nSub = 20, seed = 1, mcTol = 1, ...)
  }
  .calls <- 0
  .real <- .cur$.cuiSolveDraw
  local_mocked_bindings(.cuiSolveDraw = function(...) {
    .calls <<- .calls + 1
    .real(...)
  })
  .a <- .run(3, cacheDir = "cc")
  expect_equal(.calls, 4)
  expect_true(file.exists(file.path("cc", "draw_00003.rds")))
  .calls <- 0
  .b <- .run(3, cacheDir = "cc")
  expect_equal(.calls, 0)
  expect_equal(.a$summary, .b$summary)
  # a larger nDraw simulates only the new draws, and matches an uncached run
  .calls <- 0
  .c <- .run(5, cacheDir = "cc")
  expect_equal(.calls, 2)
  .fresh <- .run(5)
  expect_equal(.c$summary, .fresh$summary)
  expect_equal(.c$cui, .fresh$cui)
  .calls <- 0
  .run(3, cacheDir = "cc", restart = TRUE)
  expect_equal(.calls, 4)
  # a cache belongs to one simulation
  expect_error(cui(.eps, .reg, 50, mode = "trial", nDraw = 2, nSub = 20, seed = 2,
                   cacheDir = "cc"), "different simulation")
  # NULL derives the directory name; NA (the default) writes nothing
  .run(1, cacheDir = NULL)
  expect_length(list.files(pattern = "^nlmixr2CuiCache_"), 1)
  expect_equal(sort(list.files()), sort(c("cc", list.files(pattern = "^nlmixr2CuiCache_"))))
})
