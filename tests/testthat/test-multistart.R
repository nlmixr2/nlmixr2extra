skip_on_cran()

# A small, fast model used by every test that needs a real fit.  linCmt() keeps
# each start to about a second.
.msTestModel <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    linCmt() ~ add(add.sd)
  })
}

.msTestUi <- function() {
  rxode2::rxUiDecompress(rxode2::as.rxUi(.msTestModel()))
}

.msTestFit <- function() {
  suppressMessages(suppressWarnings(
    nlmixr2est::nlmixr2(.msTestModel(), nlmixr2data::theo_sd, est = "focei",
                        control = nlmixr2est::foceiControl(print = 0))
  ))
}

# The estimates of every start, as a matrix, for comparing perturbations.
.msEstMatrix <- function(starts) {
  do.call(rbind, lapply(starts, function(x) stats::setNames(x$est, x$name)))
}

test_that("multistartControl() validates and defaults", {
  ctl <- multistartControl()
  expect_s3_class(ctl, "multistartControl")
  expect_equal(ctl$n, 10L)
  # nFit defaults to n
  expect_equal(ctl$nFit, ctl$n)
  expect_equal(ctl$sampling, "uniform")
  expect_equal(ctl$screen, "posthoc")
  expect_null(ctl$cacheDir)
  # NA is the "no cache" sentinel and must survive validation
  expect_true(is.na(multistartControl(cacheDir = NA)$cacheDir))
  expect_equal(multistartControl(cacheDir = "abc")$cacheDir, "abc")

  expect_error(multistartControl(n = 0))
  expect_error(multistartControl(spread = -1))
  expect_error(multistartControl(sampling = "nope"))
  expect_error(multistartControl(omegaFold = 1))
  expect_error(multistartControl(cores = 0))
  expect_error(multistartControl(cacheDir = ""))

  # nFit is clamped to n rather than silently over-fitting
  expect_warning(ctl <- multistartControl(n = 3, nFit = 10), "larger than")
  expect_equal(ctl$nFit, 3L)

  dep <- rxode2::rxUiDeparse(multistartControl(n = 3), "x")
  expect_match(deparse1(dep), "^x <- multistartControl\\(")
})

test_that(".msPerturbIni() respects fix, bounds and 'which'", {
  iniDf <- .msTestUi()$iniDf
  ctl <- multistartControl(n = 6, spread = 0.3)
  starts <- rxode2::rxWithSeed(1234, .msPerturbIni(iniDf, ctl))

  expect_length(starts, 6L)
  # the first candidate is always the unperturbed starting point
  expect_equal(starts[[1]], iniDf)
  m <- .msEstMatrix(starts)
  expect_false(any(is.na(m)))
  # every other candidate really did move
  expect_true(all(apply(m[-1, , drop = FALSE], 1, function(r) any(r != m[1, ]))))

  # add.sd is bounded below at 0, and a wide spread must not cross it
  wide <- rxode2::rxWithSeed(1, .msPerturbIni(iniDf, multistartControl(n = 200, spread = 5)))
  addSd <- vapply(wide, function(x) x$est[x$name == "add.sd"], numeric(1))
  expect_true(all(addSd > 0))
  # variances are perturbed multiplicatively, so they stay positive too
  etaKa <- vapply(wide, function(x) x$est[x$name == "eta.ka"], numeric(1))
  expect_true(all(etaKa > 0))

  # fixed parameters are never touched
  fixed <- iniDf
  fixed$fix[fixed$name == "tcl"] <- TRUE
  starts <- rxode2::rxWithSeed(1, .msPerturbIni(fixed, multistartControl(n = 5, spread = 1)))
  expect_equal(unique(vapply(starts, function(x) x$est[x$name == "tcl"], numeric(1))),
               fixed$est[fixed$name == "tcl"])

  # 'which' restricts the perturbation to the named parameters
  starts <- rxode2::rxWithSeed(1, .msPerturbIni(iniDf, multistartControl(n = 4, which = c("tka", "tv"))))
  m <- .msEstMatrix(starts)
  untouched <- setdiff(colnames(m), c("tka", "tv"))
  for (p in untouched) {
    expect_equal(unname(unique(m[, p])), unname(m[1, p]), info = p)
  }
  expect_true(all(m[-1, "tka"] != m[1, "tka"]))

  expect_error(.msPerturbIni(iniDf, multistartControl(which = "nosuchparam")), "not in the model")
})

test_that(".msPerturbIni() is reproducible and, for uniform and normal, independent of n", {
  iniDf <- .msTestUi()$iniDf
  g <- function(n, sampling = "uniform", seed = 1234) {
    rxode2::rxWithSeed(seed, .msEstMatrix(.msPerturbIni(iniDf, multistartControl(n = n, sampling = sampling))))
  }
  expect_equal(g(4), g(4))
  expect_false(isTRUE(all.equal(g(4), g(4, seed = 99))))

  # asking for more starts later must not move the ones already estimated, or a
  # cached run could not be resumed
  for (samp in c("uniform", "normal")) {
    expect_equal(g(2, samp), g(6, samp)[1:2, , drop = FALSE], info = samp)
  }
  # a Latin hypercube is a joint design over all n candidates, so it does move;
  # .multistartRun() detects that and re-estimates those starts
  expect_false(isTRUE(all.equal(g(2, "lhs"), g(6, "lhs")[1:2, , drop = FALSE])))
})

test_that(".msPerturbIni() covers each Latin hypercube stratum exactly once", {
  iniDf <- .msTestUi()$iniDf
  n <- 21L
  starts <- rxode2::rxWithSeed(3, .msPerturbIni(iniDf, multistartControl(n = n, sampling = "lhs", spread = 0.5)))
  tka <- vapply(starts[-1], function(x) x$est[x$name == "tka"], numeric(1))
  est <- iniDf$est[iniDf$name == "tka"]
  half <- 0.5 * max(abs(est), 1)
  # map back onto (0, 1) and check one draw landed in each of the n-1 strata
  u <- (tka - (est - half)) / (2 * half)
  expect_equal(sort(floor(u * length(tka))), seq_along(tka) - 1L)
})

test_that(".msPerturbIni() keeps a correlated omega positive definite", {
  corMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka + eta.cl ~ c(0.6, 0.05, 0.3)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  iniDf <- suppressMessages(rxode2::rxUiDecompress(rxode2::as.rxUi(corMod)))$iniDf
  starts <- rxode2::rxWithSeed(1, .msPerturbIni(iniDf, multistartControl(n = 30, omegaFold = 8, spread = 0.5)))
  minEigen <- vapply(starts, function(x) {
    e <- x[!is.na(x$neta1), ]
    m <- matrix(0, 2, 2)
    for (i in seq_len(nrow(e))) {
      m[e$neta1[i], e$neta2[i]] <- m[e$neta2[i], e$neta1[i]] <- e$est[i]
    }
    min(eigen(m, symmetric = TRUE, only.values = TRUE)$values)
  }, numeric(1))
  expect_true(all(minEigen > 0))
})

test_that(".msPerturbIni() leaves a same() block alone", {
  sameMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.iov1 ~ 0.2
      eta.iov2 ~ same()
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka + eta.iov1*(OCC == 1) + eta.iov2*(OCC == 2))
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  ui <- suppressWarnings(suppressMessages(rxode2::rxUiDecompress(rxode2::as.rxUi(sameMod))))
  expect_false(is.null(ui$omegaSameMap))
  expect_warning(
    starts <- rxode2::rxWithSeed(1, .msPerturbIni(ui$iniDf, multistartControl(n = 3, spread = 0.3),
                                                 omegaSameMap = ui$omegaSameMap)),
    "same"
  )
  m <- .msEstMatrix(starts)
  for (p in c("eta.ka", "eta.iov1", "eta.iov2")) {
    expect_equal(unname(unique(m[, p])), unname(m[1, p]), info = p)
  }
  # the population parameters were still perturbed
  expect_true(all(m[-1, "tka"] != m[1, "tka"]))
})

test_that(".msSummaryRow() reports a failed start instead of aborting", {
  failed <- try(stop("boom"), silent = TRUE)
  row <- .msSummaryRow(failed, 3L, 1.5)
  expect_equal(row$start, 3L)
  expect_true(is.na(row$OBJF))
  expect_match(row$message, "boom")
  expect_equal(row$elapsed, 1.5)
  expect_null(.msEstimateRow(failed))
})

test_that(".msQuietControl() leaves a character covMethod alone", {
  focei <- .msQuietControl(nlmixr2est::foceiControl())
  expect_equal(focei$covMethod, 0L)
  expect_equal(focei$print, 0L)
  expect_false(focei$calcTables)

  saem <- .msQuietControl(nlmixr2est::saemControl())
  # setQuietFastControl() would put 0L here, which saem cannot use
  expect_true(is.character(saem$covMethod))
  expect_equal(saem$print, 0L)
  expect_false(saem$calcTables)
})

test_that("multistart() estimates every start and picks the best", {
  fit <- .msTestFit()
  ms <- multistart(fit, control = list(n = 3, spread = 0.3, screen = "none",
                                       cacheDir = NA, refitBest = FALSE))
  expect_s3_class(ms, "nlmixr2Multistart")
  expect_equal(nrow(ms$starts), 3L)
  expect_true(all(ms$starts$fitted))
  expect_equal(nrow(ms$summary), 3L)
  expect_true(all(c("start", "OBJF", "dOBJF", "AIC", "BIC", "converged",
                    "boundary", "elapsed") %in% names(ms$summary)))
  # the final estimates are carried alongside the objective functions
  expect_true(all(c("tka", "tcl", "tv", "add.sd") %in% names(ms$summary)))
  # sorted best first, with dOBJF measured from the best
  expect_false(is.unsorted(ms$summary$OBJF))
  expect_equal(min(ms$summary$dOBJF), 0)
  expect_equal(ms$bestIndex, ms$summary$start[1])
  expect_s3_class(ms$best, "nlmixr2FitCore")
  expect_length(ms$fits, 3L)
  expect_equal(as.data.frame(ms), ms$summary)
  expect_output(print(ms), "start")
  expect_invisible(print(ms))

  # a well behaved model finds the same optimum from every start
  expect_lt(max(ms$summary$dOBJF), 1)
})

test_that("multistart() screening fits only the best candidates", {
  fit <- .msTestFit()
  ms <- multistart(fit, control = list(n = 4, nFit = 2, spread = 0.4,
                                       screen = "posthoc", cacheDir = NA,
                                       refitBest = FALSE))
  expect_equal(nrow(ms$starts), 4L)
  # every candidate is screened ...
  expect_true(all(is.finite(ms$starts$screenOFV)))
  # ... but only nFit of them are estimated
  expect_equal(sum(ms$starts$fitted), 2L)
  expect_equal(nrow(ms$summary), 2L)
  # and the ones estimated are the ones that screened best
  expect_setequal(ms$starts$start[ms$starts$fitted],
                  ms$starts$start[order(ms$starts$screenOFV)][1:2])
})

test_that("multistart() works from a model and data", {
  ms <- multistart(.msTestModel(), nlmixr2data::theo_sd,
                   control = list(n = 2, screen = "none", cacheDir = NA,
                                  refitBest = FALSE))
  expect_s3_class(ms, "nlmixr2Multistart")
  expect_equal(nrow(ms$summary), 2L)
  expect_null(ms$origFit)
})

test_that("multistart() rejects things that are not models", {
  expect_error(multistart(1), "nlmixr2 fit or a nlmixr2 model")
})

test_that("multistart() resumes from its cache", {
  withr::with_tempdir({
    fit <- .msTestFit()
    args <- list(spread = 0.3, screen = "none", cacheDir = "msCache",
                 refitBest = FALSE)
    m1 <- multistart(fit, control = c(list(n = 2), args))
    expect_equal(sort(list.files("msCache", pattern = "^fit_")),
                 c("fit_0001.rds", "fit_0002.rds"))

    m2 <- multistart(fit, control = c(list(n = 4), args))
    # only the two new starts were added
    expect_equal(sort(list.files("msCache", pattern = "^fit_")),
                 c("fit_0001.rds", "fit_0002.rds", "fit_0003.rds", "fit_0004.rds"))
    # and the cached starts came back unchanged
    keep <- setdiff(names(m1$starts), "fitted")
    expect_equal(m1$starts[, keep], m2$starts[1:2, keep])
    expect_equal(sort(m1$summary$OBJF),
                 sort(m2$summary$OBJF[m2$summary$start <= 2]))

    # restart = TRUE throws the cache away
    m3 <- multistart(fit, control = c(list(n = 2, restart = TRUE), args))
    expect_equal(sort(list.files("msCache", pattern = "^fit_")),
                 c("fit_0001.rds", "fit_0002.rds"))
    expect_equal(sort(m1$summary$OBJF), sort(m3$summary$OBJF))
  })
})

test_that("multistart() drops a cached start that no longer matches", {
  withr::with_tempdir({
    fit <- .msTestFit()
    args <- list(sampling = "lhs", screen = "none", cacheDir = "lhsCache",
                 refitBest = FALSE)
    m1 <- multistart(fit, control = c(list(n = 2), args))
    m2 <- multistart(fit, control = c(list(n = 4), args))
    # the Latin hypercube design changed, so start 2 is a different start
    keep <- setdiff(names(m1$starts), "fitted")
    expect_false(isTRUE(all.equal(m1$starts[2, keep], m2$starts[2, keep],
                                  check.attributes = FALSE)))
    # and it was re-estimated rather than reported against the wrong start
    expect_equal(nrow(m2$summary), 4L)
  })
})

test_that("plot() draws the waterfall and parameter-stability plots", {
  fit <- .msTestFit()
  ms <- multistart(fit, control = list(n = 3, spread = 0.3, screen = "none",
                                       cacheDir = NA, refitBest = FALSE))
  expect_s3_class(plot(ms), "ggplot")
  expect_s3_class(plot(ms, "waterfall"), "ggplot")
  expect_s3_class(plot(ms, "waterfall", dOfvMax = 1), "ggplot")
  expect_s3_class(plot(ms, "parameters"), "ggplot")
  expect_s3_class(plot(ms, "parameters", kBest = 2), "ggplot")
  expect_error(plot(ms, "nope"))

  # nothing to plot when every start failed
  empty <- ms
  empty$summary$OBJF <- NA_real_
  expect_error(plot(empty), "no successfully estimated starts")
  expect_error(plot(empty, "parameters"), "no successfully estimated starts")
})
