skip_on_cran()

# covarSearchAuto() covariate-selection fixes (Issue 103) ----

# Two bugs were fixed here:
#
# 1. When a covariate passed the inclusion threshold, forwardSearch()/
#    backwardSearch() tried to recover the best fit with
#    `covSearchRes[[which.min(resTable$pchisqr)]][[1]]`. The elements of
#    covSearchRes are model objects (rxUi/environments), not lists, so `[[1]]`
#    raised "wrong arguments for subsetting an environment" and the search
#    crashed. The best model is now re-fit to recover its fit object.
#
# 2. The forward inclusion test had an inverted sign: with
#    `dObjf <- fit$objf - x$objf` an improving candidate has dObjf > 0, but the
#    code computed a p-value only when dObjf < 0. As a result an improving
#    covariate always received pchisqr == 1 and was never selected, while a
#    worsening candidate was the one that got a p-value (and triggered bug 1).

# Simulate a data set with a strong, well-scaled covariate effect on clearance
# so that the covariate is unambiguously selected -- this exercises the re-fit
# line that used to crash and the corrected forward inclusion test.
simCovModel <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.0
    tv <- 3.45
    b.cl <- 0.5
    eta.ka ~ 0.3
    eta.cl ~ 0.02
    eta.v ~ 0.05
    add.sd <- 0.2
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + b.cl * SCWT + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

makeCovData <- function() {
  set.seed(3)
  obsT <- seq(0.5, 60, length.out = 12)
  ev <- do.call(rbind, lapply(seq_len(50), function(id) {
    z <- round(stats::rnorm(1), 3)
    e <- as.data.frame(rxode2::et(amt = 320) |> rxode2::et(obsT))
    e$ID <- id
    e$SCWT <- z
    e
  }))
  s <- rxode2::rxSolve(simCovModel(), ev, returnType = "data.frame")
  dose <- do.call(rbind, lapply(seq_len(50), function(id) {
    data.frame(ID = id, time = 0, DV = NA, amt = 320, evid = 1,
               SCWT = ev$SCWT[ev$ID == id][1])
  }))
  obs <- data.frame(ID = s$id, time = s$time, DV = s$sim, amt = NA, evid = 0,
                    SCWT = s$SCWT)
  d <- rbind(dose, obs)
  d[order(d$ID, d$time, -d$evid), ]
}

baseCovModel <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.0
    tv <- 3.45
    eta.ka ~ 0.3
    eta.cl ~ 0.1
    eta.v ~ 0.1
    add.sd <- 0.3
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

test_that("covarSearchAuto completes and selects a real covariate (Issue 103)", {
  d <- makeCovData()
  fit <- suppressWarnings(
    nlmixr2(baseCovModel(), d, est = "focei",
            control = nlmixr2est::foceiControl(print = 0))
  )

  res <- suppressWarnings(
    covarSearchAuto(fit, varsVec = c("cl", "v"), covarsVec = "SCWT",
                    searchType = "forward", restart = TRUE)
  )

  # returns the documented structure ...
  expect_true(is.list(res))
  expect_true(all(c("summaryTable", "resFwd") %in% names(res)))

  # ... the re-fit line recovered a genuine fit object (bug 1: no more
  # "wrong arguments for subsetting an environment")
  expect_s3_class(res$resFwd[[1]], "nlmixr2FitData")

  # ... and the strong SCWT-on-cl effect was selected (bug 2: improving
  # covariates are now included instead of always getting pchisqr == 1)
  included <- res$summaryTable[res$summaryTable$included == "yes", ]
  expect_true(nrow(included) >= 1)
  expect_true(any(unlist(included$covar) == "SCWT" & unlist(included$var) == "cl"))
})

# A covariate coefficient is added at exactly 0, which carries no scale; FOCEi
# then nudges it to foceiControl(zeroTheta) = 0.001 and steps by that, which
# pins a unit-scale covariate (SCWT above) at zero.  .scmCovIni() starts each
# new coefficient at 0.1/max(|X|) instead, so the step follows the covariate's
# own units -- it must not overshoot a covariate measured in the tens.
.cur <- loadNamespace("nlmixr2extra")

test_that(".scmCovIni() scales new covariate coefficients by the data", {
  d <- data.frame(ID = 1:4, SCWT = c(-2, -1, 1, 2), WT = c(50, 70, 90, 100))
  ui <- nlmixr2extra::buildupatedUI(baseCovModel(), varsVec = c("cl", "v"),
                                    covarsVec = c("SCWT", "WT"),
                                    indep = FALSE, add = TRUE)
  ini <- .cur$.scmCovIni(ui, d)$iniDf
  est <- setNames(ini$est, ini$name)
  expect_equal(est[["cov_SCWT_cl"]], 0.1 / 2)
  expect_equal(est[["cov_WT_v"]], 0.1 / 100)
  # the starting term stays within 0.1 for every subject
  expect_lte(max(abs(est[["cov_WT_v"]] * d$WT)), 0.1)
  # structural parameters are untouched
  expect_equal(est[["tcl"]], 1.0)
})

test_that(".scmCovIni() leaves estimated or unusable coefficients alone", {
  ui <- nlmixr2extra::buildupatedUI(baseCovModel(), varsVec = "cl",
                                    covarsVec = "SCWT", indep = FALSE, add = TRUE)
  # already estimated: kept
  ini <- .cur$.scmCovIni(rxode2::ini(ui, cov_SCWT_cl = 0.3),
                         data.frame(SCWT = c(-2, 2)))$iniDf
  expect_equal(ini$est[ini$name == "cov_SCWT_cl"], 0.3)
  # covariate missing, constant zero, or non-numeric: left at 0
  for (d in list(data.frame(AGE = 1:3),
                 data.frame(SCWT = c(0, 0, 0)),
                 data.frame(SCWT = c("a", "b", "c")))) {
    ini <- .cur$.scmCovIni(ui, d)$iniDf
    expect_equal(ini$est[ini$name == "cov_SCWT_cl"], 0)
  }
})

test_that("covarSearchAuto selects a covariate measured in the tens", {
  # true effect 0.03 per kg on log(CL); a coefficient started at 0.1 (rather
  # than 0.1/max(WT)) would put the starting CL off by exp(0.1*110)
  wtModel <- function() {
    ini({
      tka <- 0.45
      tcl <- -1.0
      tv <- 3.45
      b.cl <- 0.03
      eta.ka ~ 0.3
      eta.cl ~ 0.02
      eta.v ~ 0.05
      add.sd <- 0.15
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + b.cl * WT + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  set.seed(11)
  obsT <- seq(0.5, 60, length.out = 12)
  ev <- do.call(rbind, lapply(seq_len(60), function(id) {
    e <- as.data.frame(rxode2::et(amt = 320) |> rxode2::et(obsT))
    e$ID <- id
    e$WT <- round(stats::runif(1, 50, 110))
    e
  }))
  s <- rxode2::rxSolve(wtModel(), ev, returnType = "data.frame")
  dose <- do.call(rbind, lapply(seq_len(60), function(id) {
    data.frame(ID = id, time = 0, DV = NA, amt = 320, evid = 1,
               WT = ev$WT[ev$ID == id][1])
  }))
  obs <- data.frame(ID = s$id, time = s$time, DV = s$sim, amt = NA, evid = 0,
                    WT = s$WT)
  d <- rbind(dose, obs)
  d <- d[order(d$ID, d$time, -d$evid), ]

  fit <- suppressWarnings(
    nlmixr2(baseCovModel(), d, est = "focei",
            control = nlmixr2est::foceiControl(print = 0))
  )
  res <- suppressWarnings(
    covarSearchAuto(fit, varsVec = "cl", covarsVec = "WT",
                    searchType = "forward", restart = TRUE)
  )
  included <- res$summaryTable[res$summaryTable$included == "yes", ]
  expect_true(any(unlist(included$covar) == "WT" & unlist(included$var) == "cl"))
  expect_equal(unname(unlist(included$covarEffect)[1]), 0.03, tolerance = 0.2)
})
