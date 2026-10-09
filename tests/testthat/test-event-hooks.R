skip_if_not(exists("rxEventEmit", envir = asNamespace("rxode2"), inherits = FALSE),
            "rxode2 has no event bus")

.evRec <- new.env()
.evListen <- function(env = parent.frame()) {
  .evRec$ev <- list()
  rxode2::rxEventListen("nlmixr2extra-test", function(event, ...) {
    .evRec$ev[[length(.evRec$ev) + 1L]] <- list(event = event, p = list(...))
  })
  withr::defer(rxode2::rxEventUnlisten("nlmixr2extra-test"), envir = env)
}
.evNames <- function() vapply(.evRec$ev, function(e) e$event, character(1))
.evFake <- function() structure(list(env = new.env()), class = c("nlmixr2FitCore", "list"))

test_that(".extraEventExit picks the event from what the driver returned", {
  .evListen()
  fit <- .evFake()
  other <- .evFake()
  .ex <- function(res, ...) {
    .extraEventEnter()
    .extraEventExit(res, fit, quote(drv(fit)), "k", "drv", ...)
  }
  .ex(fit)
  .ex(other)
  .ex(data.frame(a = 1))
  .ex(1, update = TRUE)
  .ex(NULL)
  expect_identical(.evNames(), c("fitUpdate", "fitComplete", "fitResult", "fitUpdate"))
  expect_identical(.evRec$ev[[1]]$p$what, "k")
  expect_identical(.evRec$ev[[2]]$p$source, "k")
  expect_identical(.evRec$ev[[2]]$p$object, fit)
  expect_identical(.evRec$ev[[3]]$p$kind, "k")
  ## a non-fit input never emits, and the depth is always restored
  .extraEventEnter()
  .extraEventExit(data.frame(a = 1), list(), NULL, "k", "drv")
  expect_length(.evRec$ev, 4L)
  expect_identical(rxode2::rxEventDepth(), 0L)
})

test_that("summaries never carry fits", {
  .e <- new.env()
  assign("objDf", data.frame(OBJF = 12.5), envir = .e)
  fit <- structure(list(env = .e), class = c("nlmixr2FitCore", "list"))
  s <- .extraEventSummary(structure(list(best = fit, tab = data.frame(x = 1), fn = sum,
                                         nest = list(fit)), class = "myResult"))
  expect_identical(s$best$OBJF, 12.5)
  expect_null(s$fn)
  expect_identical(s$nest[[1]]$OBJF, 12.5)
  expect_identical(attr(s, "nlmixr2extraClass"), "myResult")
  expect_false(any(vapply(unlist(s, recursive = TRUE), function(x) inherits(x, "nlmixr2FitCore"), TRUE)))
})

.evModel <- function() {
  ini({
    tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5)
    eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}
.evFit <- function() rxode2::rxEventScope(suppressMessages(suppressWarnings(
  nlmixr2(.evModel, nlmixr2data::theo_sd, est = "focei", control = list(print = 0))
)))

test_that("bootstrapFit: one fitUpdate, no fitComplete from the bootstrap fits", {
  skip_on_cran()
  withr::local_dir(withr::local_tempdir())
  fit <- .evFit()
  .evListen()
  suppressMessages(suppressWarnings(bootstrapFit(fit, nboot = 2, restart = TRUE)))
  expect_identical(.evNames(), "fitUpdate")
  expect_identical(.evRec$ev[[1]]$p$what, "bootstrap")
  expect_identical(rxode2::rxEventDepth(), 0L)
})

test_that("preconditionFit: one fitUpdate (covariance)", {
  skip_on_cran()
  fit <- .evFit()
  .evListen()
  suppressMessages(suppressWarnings(preconditionFit(fit)))
  expect_identical(.evNames(), "fitUpdate")
  expect_identical(.evRec$ev[[1]]$p$what, "covariance")
})

test_that("multistart and profile: one fitResult each, with no fit inside", {
  skip_on_cran()
  fit <- .evFit()
  .evListen()
  suppressMessages(suppressWarnings(multistart(fit, control = list(
    n = 2, spread = 0.3, screen = "none", cacheDir = NA, refitBest = FALSE
  ))))
  suppressMessages(profile(fit, which = data.frame(tka = log(c(1.4, 1.6))), method = "fixed"))
  expect_identical(.evNames(), c("fitResult", "fitResult"))
  expect_identical(vapply(.evRec$ev, function(e) e$p$kind, ""), c("multistart", "profile"))
  .has <- function(x) inherits(x, "nlmixr2FitCore") ||
    (is.list(x) && !is.data.frame(x) && any(vapply(x, .has, TRUE)))
  expect_false(.has(.evRec$ev[[1]]$p$result))
  expect_true(object.size(.evRec$ev[[1]]$p$result) < 1e6)
})

test_that("linearize: one fitComplete linked to the input fit", {
  skip_on_cran()
  fit <- .evFit()
  .evListen()
  lin <- suppressMessages(suppressWarnings(linearize(fit)))
  expect_identical(.evNames(), "fitComplete")
  expect_identical(.evRec$ev[[1]]$p$source, "linearize")
  expect_identical(.evRec$ev[[1]]$p$object, fit)
})

test_that("a driver error restores the depth and emits nothing", {
  .evListen()
  expect_error(bootstrapFit(.evFake(), nboot = 2, stratVar = "nope"))
  expect_length(.evRec$ev, 0L)
  expect_identical(rxode2::rxEventDepth(), 0L)
})

test_that("the event names the user's fit even when the driver reassigns it", {
  skip_on_cran()
  ## linearize() refits a posthoc fit with focei into its own `fit`
  fit <- rxode2::rxEventScope(suppressMessages(suppressWarnings(
    nlmixr2(.evModel, nlmixr2data::theo_sd, est = "posthoc", control = list(print = 0))
  )))
  .evListen()
  suppressMessages(suppressWarnings(linearize(fit)))
  expect_identical(.evNames(), "fitComplete")
  expect_identical(.evRec$ev[[1]]$p$object, fit)
})
