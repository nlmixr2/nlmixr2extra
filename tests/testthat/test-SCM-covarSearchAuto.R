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
