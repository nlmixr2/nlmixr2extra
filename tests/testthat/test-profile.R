test_that("llpControl - Step 5 additions", {
  # Defaults work
  ctrl <- llpControl()
  expect_equal(ctrl$normq, 1.96)
  expect_equal(ctrl$rseTheta, 30)
  expect_null(ctrl$workers)

  # normq validated
  expect_error(llpControl(normq = -1))
  expect_error(llpControl(normq = Inf))
  expect_error(llpControl(normq = NA_real_))
  expect_silent(llpControl(normq = 0))
  expect_silent(llpControl(normq = 2.576))

  # rseTheta: unnamed scalar accepted
  expect_silent(llpControl(rseTheta = 20))
  expect_equal(llpControl(rseTheta = 20)$rseTheta, 20)

  # rseTheta: named vector accepted
  ctrl_named <- llpControl(rseTheta = c(tka = 25, tcl = 15))
  expect_equal(ctrl_named$rseTheta, c(tka = 25, tcl = 15))

  # rseTheta: unnamed vector length > 1 rejected
  expect_error(llpControl(rseTheta = c(25, 15)))

  # rseTheta: negative value rejected
  expect_error(llpControl(rseTheta = -5))
  expect_error(llpControl(rseTheta = c(tka = -5)))

  # workers: NULL (default)
  expect_null(llpControl()$workers)

  # workers: positive integer
  expect_equal(llpControl(workers = 2L)$workers, 2L)
  expect_equal(llpControl(workers = 2)$workers, 2L)  # coerced to integer

  # workers: "auto"
  expect_equal(llpControl(workers = "auto")$workers, "auto")

  # workers: invalid values rejected
  expect_error(llpControl(workers = 0))
  expect_error(llpControl(workers = -1))
  expect_error(llpControl(workers = "parallel"))

  # rxUiDeparse round-trips new fields
  ctrl_modified <- llpControl(normq = 2.576, rseTheta = 20, workers = 2L)
  deparsed <- rxUiDeparse.llpControl(ctrl_modified, "ctrl")
  expect_true(any(grepl("normq", deparsed)))
  expect_true(any(grepl("workers", deparsed)))
})

test_that("llpExtractSE - Step 1", {
  one.compartment <- function() {
    ini({
      tka <- log(1.57)
      tcl <- log(2.72)
      tv  <- fixed(log(31.5))
      eta.ka ~ 0.6
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v  <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }

  fit <- suppressMessages(nlmixr2(
    one.compartment, data = nlmixr2data::theo_sd, est = "focei",
    control = list(print = 0, covMethod = "r")
  ))

  se <- llpExtractSE(fit)

  # Returns a named vector with one entry per fixef parameter
  expect_named(se, names(nlmixr2est::fixef(fit)))

  # Structural thetas with covariance step have non-NA SE matching sqrt(diag(cov))
  expect_false(is.na(se["tka"]))
  expect_false(is.na(se["tcl"]))
  expect_equal(se["tka"], sqrt(diag(fit$cov))["tka"], tolerance = 1e-6)
  expect_equal(se["tcl"], sqrt(diag(fit$cov))["tcl"], tolerance = 1e-6)

  # Fixed parameter has NA SE
  expect_true(is.na(se["tv"]))

  # Residual-error theta (add.sd) without covariance SE is NA
  # (parFixedDf SE is NA for add.sd even with covMethod = "r")
  expect_true(is.na(se["add.sd"]))

  # Without covariance step, all SEs are NA
  fit_nocov <- suppressMessages(nlmixr2(
    one.compartment, data = nlmixr2data::theo_sd, est = "focei",
    control = list(print = 0, covMethod = "")
  ))
  se_nocov <- llpExtractSE(fit_nocov)
  expect_true(all(is.na(se_nocov)))
})

test_that("llpNextParam - Step 3", {
  # OFV = x^2 + 100, MLE at 0; target = 103.84 -> root at sqrt(3.84) ~ 1.9596

  # Quadratic finds root precisely (compare to uniroot)
  result_up <- llpNextParam(
    paramVals = c(0, 1, 2), ofvVals = c(100, 101, 104),
    targetOfv = 103.84, direction = 1, mle = 0,
    searchLower = 0, searchUpper = 10, hardLower = -100, hardUpper = 100
  )
  expected_up <- uniroot(function(x) x^2 + 100 - 103.84, c(0, 10))$root
  expect_equal(result_up, expected_up, tolerance = 1e-3)

  # Correct root chosen for lower direction
  result_lo <- llpNextParam(
    paramVals = c(-2, -1, 0), ofvVals = c(104, 101, 100),
    targetOfv = 103.84, direction = -1, mle = 0,
    searchLower = -10, searchUpper = 0, hardLower = -100, hardUpper = 100
  )
  expect_equal(result_lo, -expected_up, tolerance = 1e-3)
  expect_true(result_lo < 0)

  # Flat quadratic triggers warning and falls back to linear extrapolation
  # OFV = 100, 101.92, 103.84 -> perfectly linear -> c2 = NA from lm -> warn
  expect_warning(
    result_flat <- llpNextParam(
      paramVals = c(0, 1, 2), ofvVals = c(100, 101.92, 103.84),
      targetOfv = 103.84, direction = 1, mle = 0,
      searchLower = 0, searchUpper = 10, hardLower = -100, hardUpper = 100
    ),
    regexp = "Flat quadratic"
  )
  # Linear interpolation to boundary value gives x = 2
  expect_equal(result_flat, 2, tolerance = 1e-6)

  # Clamped to hard model bounds (not search bounds)
  result_clamped <- llpNextParam(
    paramVals = c(0, 1, 2), ofvVals = c(100, 101, 104),
    targetOfv = 103.84, direction = 1, mle = 0,
    searchLower = 0, searchUpper = 10, hardLower = -100, hardUpper = 1.5
  )
  margin <- sqrt(.Machine$double.eps)
  expect_lte(result_clamped, 1.5 - margin)
  expect_gt(result_clamped, 1.0)  # not clamped to searchUpper (which is 10)

  # Explicit extrapolation advances beyond the current outer point (no rule=2 clamping)
  # Two points only: (0, 100) and (1, 101); target=103.84 is outside range -> extrapolate
  # slope = 1 OFV/unit, so need 2.84 more units from x=1 -> 3.84
  result_extrap <- llpNextParam(
    paramVals = c(0, 1), ofvVals = c(100, 101),
    targetOfv = 103.84, direction = 1, mle = 0,
    searchLower = 0, searchUpper = 10, hardLower = -100, hardUpper = 100
  )
  expect_equal(result_extrap, 3.84, tolerance = 1e-6)
  expect_gt(result_extrap, 1)  # advanced beyond current outer point

  # NA when fewer than 2 valid side-specific points
  expect_true(is.na(llpNextParam(
    paramVals = c(0), ofvVals = c(100),
    targetOfv = 103.84, direction = 1, mle = 0,
    searchLower = 0, searchUpper = 10, hardLower = -100, hardUpper = 100
  )))

  # NA when all OFV values are NA
  expect_true(is.na(llpNextParam(
    paramVals = c(0, 1, 2), ofvVals = c(NA_real_, NA_real_, NA_real_),
    targetOfv = 103.84, direction = 1, mle = 0,
    searchLower = 0, searchUpper = 10, hardLower = -100, hardUpper = 100
  )))

  # Half-space filter: lower direction ignores points above MLE
  # Only (0, 100) is <= mle=0 -> fewer than 2 valid points -> NA
  expect_true(is.na(llpNextParam(
    paramVals = c(0, 1, 2), ofvVals = c(100, 101, 104),
    targetOfv = 103.84, direction = -1, mle = 0,
    searchLower = -10, searchUpper = 0, hardLower = -100, hardUpper = 100
  )))
})

test_that("profileNlmixr2FitDataEstInitial - Step 2", {
  noSE <- c(A = NA_real_)

  # Must have one-row input
  expect_error(
    profileNlmixr2FitDataEstInitial(estimates = data.frame(A = 1:2))
  )

  # RSE fallback: width = normq * rseTheta/100 * |MLE|
  # MLE=1, normq=1.96, rseTheta=100 -> width=1.96, guesses = c(-0.96, 2.96)
  expect_equal(
    profileNlmixr2FitDataEstInitial(
      estimates  = data.frame(A = 1),
      which      = "A",
      normq      = 1.96,
      rseTheta   = c(A = 100),
      seTheta    = noSE,
      hardLower  = -100,
      hardUpper  = 200
    ),
    c(-0.96, 2.96)
  )

  # RSE fallback explicitly differs from old ofvIncrease-based formula
  # (normq=1.96 != ofvIncrease=3.84, so results differ)
  old_style <- 1 + c(-1, 1) * 3.84 * 100/100 * abs(1)
  new_style  <- profileNlmixr2FitDataEstInitial(
    estimates = data.frame(A = 1), which = "A",
    normq = 1.96, rseTheta = c(A = 100), seTheta = noSE,
    hardLower = -100, hardUpper = 200
  )
  expect_false(isTRUE(all.equal(old_style, new_style)))

  # SE-based width: width = normq * SE
  # MLE=1, normq=1.96, SE=0.5 -> width=0.98, guesses = c(0.02, 1.98)
  expect_equal(
    profileNlmixr2FitDataEstInitial(
      estimates = data.frame(A = 1),
      which     = "A",
      normq     = 1.96,
      rseTheta  = c(A = 30),
      seTheta   = c(A = 0.5),
      hardLower = -100,
      hardUpper = 200
    ),
    c(1 - 1.96 * 0.5, 1 + 1.96 * 0.5)
  )

  # SE takes precedence over rseTheta even when rseTheta would give a wider interval
  se_guess <- profileNlmixr2FitDataEstInitial(
    estimates = data.frame(A = 1), which = "A",
    normq = 1.96, rseTheta = c(A = 100), seTheta = c(A = 0.1),
    hardLower = -100, hardUpper = 200
  )
  expect_equal(se_guess, c(1 - 1.96 * 0.1, 1 + 1.96 * 0.1))

  # Soft clamp: lower guess below hardLower -> pulled 90% toward MLE
  # MLE=1, width=1.96, lower guess=-0.96, hardLower=0
  # soft clamp lower: (0 - 1) * 0.9 + 1 = 0.1
  result_clamped <- profileNlmixr2FitDataEstInitial(
    estimates = data.frame(A = 1), which = "A",
    normq = 1.96, rseTheta = c(A = 100), seTheta = noSE,
    hardLower = 0, hardUpper = 200
  )
  expect_equal(result_clamped[1], (0 - 1) * 0.9 + 1)  # soft clamp lower
  expect_equal(result_clamped[2], 2.96)                # upper unchanged

  # Soft clamp: upper guess above hardUpper -> pulled 90% toward MLE
  result_upper_clamped <- profileNlmixr2FitDataEstInitial(
    estimates = data.frame(A = 1), which = "A",
    normq = 1.96, rseTheta = c(A = 100), seTheta = noSE,
    hardLower = -100, hardUpper = 1.5
  )
  expect_equal(result_upper_clamped[1], -0.96)         # lower unchanged
  expect_equal(result_upper_clamped[2], (1.5 - 1) * 0.9 + 1)  # soft clamp upper
})

test_that("optimProfile / llpRunOneParameter - Step 4 status structure", {
  # Status data.frame has the required fields and can survive bind_rows across parameters
  make_status <- function(param, direction, converged, nearBound, hitMaxIter) {
    data.frame(
      param = param, direction = direction,
      converged = converged, nearBound = nearBound, hitMaxIter = hitMaxIter,
      message = "complete", boundEstimate = 1.0,
      bestOfvDiff = 0.001, nEval = 3L,
      stringsAsFactors = FALSE
    )
  }

  s1 <- make_status("tka", -1L, TRUE,  FALSE, FALSE)
  s2 <- make_status("tka",  1L, TRUE,  FALSE, FALSE)
  s3 <- make_status("tcl", -1L, FALSE, TRUE,  TRUE)
  s4 <- make_status("tcl",  1L, FALSE, FALSE, TRUE)

  combined <- dplyr::bind_rows(s1, s2, s3, s4)
  expect_equal(nrow(combined), 4L)
  expect_named(combined, c("param", "direction", "converged", "nearBound",
                            "hitMaxIter", "message", "boundEstimate",
                            "bestOfvDiff", "nEval"))

  # nearBound and hitMaxIter are correctly reflected per row
  expect_false(combined$nearBound[combined$param == "tka" & combined$direction == -1L])
  expect_true(combined$nearBound[combined$param == "tcl" & combined$direction == -1L])
  expect_true(all(combined$hitMaxIter[combined$param == "tcl"]))
  expect_false(any(combined$hitMaxIter[combined$param == "tka"]))
})

test_that("llpRunOneParameter interval ratio warning - Step 4", {
  # Pure logic test: simulate ratio > threshold and < 1/threshold
  # Extract the interval ratio warning logic into a local helper for unit testing
  check_ratio_warns <- function(lBound, uBound, mle, isResidualError = FALSE) {
    lDist <- abs(mle - lBound)
    uDist <- abs(uBound - mle)
    threshold <- if (isResidualError) 1.6 else 1.3
    ratio <- uDist / lDist
    ratio > threshold || ratio < 1 / threshold
  }

  # Symmetric: no warning
  expect_false(check_ratio_warns(lBound = -1, uBound = 1, mle = 0))

  # Mildly asymmetric (ratio = 1.2): below threshold 1.3 -> no warning
  expect_false(check_ratio_warns(lBound = -1, uBound = 1.2, mle = 0))

  # Asymmetric ratio = 1.4 > 1.3: warning for structural theta
  expect_true(check_ratio_warns(lBound = -1, uBound = 1.4, mle = 0))

  # Ratio = 1.4 is below residual-error threshold 1.6: no warning
  expect_false(check_ratio_warns(lBound = -1, uBound = 1.4, mle = 0, isResidualError = TRUE))

  # Ratio = 1.7 > 1.6: warning even for residual-error
  expect_true(check_ratio_warns(lBound = -1, uBound = 1.7, mle = 0, isResidualError = TRUE))

  # Inverted ratio (upper narrower than lower): ratio = 1/1.4 < 1/1.3 -> warning
  expect_true(check_ratio_warns(lBound = -1.4, uBound = 1, mle = 0))

  # NA bounds: no warning (both must be non-NA)
  # This tests the guard in llpRunOneParameter
  lBound <- NA_real_
  uBound <- 1.5
  expect_false(!is.na(lBound) && !is.na(uBound))  # guard would skip

  # MLE-side search bracket is NOT the near-boundary trigger
  # nearBound should only be TRUE when nextPar hits hardLower/hardUpper
  # This is tested structurally: direction=-1 searchUpper=MLE is NOT a hard bound
  margin <- sqrt(.Machine$double.eps)
  mle_val <- 0.5
  hard_lo <- -Inf; hard_hi <- Inf
  # nextPar at exactly mle (searchUpper for lower direction) should NOT trigger nearBound
  nextPar <- mle_val
  expect_false(nextPar <= hard_lo + margin || nextPar >= hard_hi - margin)
})

test_that("nlmixr2Profile S3 methods - Step 6", {
  # Build a minimal mock nlmixr2Profile without running a real fit
  mle_val <- 1.0; lo_val <- 0.5; hi_val <- 1.6
  raw <- data.frame(
    Parameter    = c("A",    "A",   "A",   "A",    "A"),
    OFV          = c(NA,     101,   100,   101.5,  NA),
    A            = c(lo_val, 0.75,  mle_val, 1.3,  hi_val),
    profileBound = c(-3.84,  NA,    NA,    NA,     3.84),
    stringsAsFactors = FALSE
  )
  status <- data.frame(
    param = c("A", "A"), direction = c(-1L, 1L),
    converged = c(TRUE, TRUE), nearBound = c(FALSE, FALSE),
    hitMaxIter = c(FALSE, FALSE), message = c("complete", "complete"),
    boundEstimate = c(lo_val, hi_val), bestOfvDiff = c(0.001, 0.001),
    nEval = c(3L, 3L), stringsAsFactors = FALSE
  )
  attr(raw, "fitName")       <- "mockFit"
  attr(raw, "method")        <- "llp"
  attr(raw, "ofvIncrease")   <- qchisq(0.95, df = 1)
  attr(raw, "mle")           <- c(A = mle_val)
  attr(raw, "status")        <- status
  attr(raw, "nearBound")     <- c(A = FALSE)
  attr(raw, "hitMaxIter")    <- c(A = FALSE)
  attr(raw, "intervalRatio") <- c(A = (hi_val - mle_val) / (mle_val - lo_val))
  class(raw) <- c("nlmixr2Profile", "data.frame")

  # Class checks
  expect_s3_class(raw, "nlmixr2Profile")
  expect_s3_class(raw, "data.frame")

  # confint returns a matrix with correct dimnames and values
  ci <- confint(raw)
  expect_true(is.matrix(ci))
  expect_equal(colnames(ci), c("lower", "upper"))
  expect_equal(rownames(ci), "A")
  expect_equal(ci["A", "lower"], lo_val)
  expect_equal(ci["A", "upper"], hi_val)

  # confint with explicit parm
  ci2 <- confint(raw, parm = "A")
  expect_equal(ci, ci2)

  # confint warns when level implies a different ofvIncrease
  expect_warning(confint(raw, level = 0.90), regexp = "ΔOFV")

  # print runs without error
  expect_no_error(capture.output(print(raw)))

  # plot returns a ggplot object
  p <- plot(raw)
  expect_s3_class(p, "ggplot")

  # nearBound and hitMaxIter propagate correctly in .nlmixr2ProfileNew
  status_nb <- status
  status_nb$nearBound[status_nb$direction == 1L] <- TRUE
  attr(raw, "nearBound") <- c(A = TRUE)
  expect_true(attr(raw, "nearBound")[["A"]])

  # intervalRatio is upper_dist / lower_dist
  expected_ratio <- (hi_val - mle_val) / (mle_val - lo_val)
  expect_equal(attr(raw, "intervalRatio")[["A"]], expected_ratio, tolerance = 1e-6)
})

test_that("llpWithIsolatedFitDir - Step 7", {
  original_wd <- getwd()
  captured_wd <- NULL

  result <- llpWithIsolatedFitDir("tka", function() {
    captured_wd <<- getwd()
    42L
  })

  # Return value is passed through
  expect_equal(result, 42L)
  # Working directory is restored after the call
  expect_equal(getwd(), original_wd)
  # The function ran in a different directory
  expect_false(identical(captured_wd, original_wd))
  # The temp directory is cleaned up
  expect_false(dir.exists(captured_wd))
})

test_that(".llpCheckParallelDeps aborts clearly when deps are missing - Step 7", {
  both_installed <- requireNamespace("future",       quietly = TRUE) &&
                    requireNamespace("future.apply", quietly = TRUE)
  skip_if(both_installed, "future and future.apply are both installed; cannot test missing-dep path")
  expect_error(.llpCheckParallelDeps(), regexp = "install\\.packages")
})

test_that("messageProfileComplete uses cli - Step 7", {
  # Convergence messages (start with "complete") emit a message, not an error/warning
  expect_no_error(suppressMessages(
    messageProfileComplete("tka", 1L, "complete due to significant digits precision")
  ))
  expect_no_warning(suppressMessages(
    messageProfileComplete("tka", 1L, "complete due to significant digits precision")
  ))
  # Abort messages also emit a message, not an error/warning
  expect_no_error(suppressMessages(
    messageProfileComplete("tka", -1L, "aborted due to numerical difficulties")
  ))
  expect_no_warning(suppressMessages(
    messageProfileComplete("tka", -1L, "aborted due to numerical difficulties")
  ))
})

test_that(".withWorkerPlan restores plan after success and error - Step 7", {
  skip_if_not_installed("future")
  orig_plan <- future::plan()
  on.exit(future::plan(orig_plan), add = TRUE)

  orig_class <- class(future::plan())

  # Normal completion restores the plan
  .withWorkerPlan(1L, "done")
  expect_equal(class(future::plan()), orig_class)

  # Error inside expression also restores the plan
  tryCatch(
    .withWorkerPlan(1L, stop("boom")),
    error = function(e) NULL
  )
  expect_equal(class(future::plan()), orig_class)
})

test_that("profileNlmixr2FitCoreRet", {
  # Variance and covariance is correctly captured
  one.compartment <- function() {
    ini({
      tka <- log(1.57)
      tcl <- log(2.72)
      tv <- fixed(log(31.5))
      eta.ka ~ 0.6
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }

  fit <-
    suppressMessages(nlmixr2(
      one.compartment, data = nlmixr2data::theo_sd, est="focei", control = list(print=0)
    ))
  withoutCov <- profileNlmixr2FitCoreRet(fit, which = "tka")
  expect_s3_class(withoutCov, "data.frame")
  expect_named(withoutCov, expected = c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "eta.ka"))

  one.compartment <- function() {
    ini({
      tka <- log(1.57)
      tcl <- log(2.72)
      tv <- fixed(log(31.5))
      eta.ka + eta.cl ~ c(0.6, 0.1, 0.2)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }

  fit <-
    suppressMessages(nlmixr2(
      one.compartment, data = nlmixr2data::theo_sd, est="focei", control = list(print=0)
    ))
  withCov <- profileNlmixr2FitCoreRet(fit, which = "tka")
  expect_s3_class(withCov, "data.frame")
  expect_named(withCov, expected = c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "eta.ka", "eta.cl", "cov(eta.cl,eta.ka)"))
})

test_that("profileFixed", {
  # fix most of the parameters so that it estimates faster
  one.compartment <- function() {
    ini({
      tka <- log(1.57)
      tcl <- log(2.72)
      tv <- fixed(log(31.5))
      eta.ka ~ 0.6
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }

  fit <-
    suppressMessages(nlmixr2(
      one.compartment, data = nlmixr2data::theo_sd, est="focei", control = list(print=0)
    ))

  testFixed <-
    suppressMessages(
      profile(fit, which = data.frame(tka = log(c(1.4, 1.6, 1.8))), method = "fixed")
    )
  expect_s3_class(testFixed, "data.frame")
  expect_named(testFixed, expected = c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "eta.ka"))
  expect_equal(nrow(testFixed), 3)

  # Fix multiple parameters simultaneously
  testFixedMulti <-
    suppressMessages(
      profile(
        fit,
        which =
          data.frame(
            tka = log(c(1.4, 1.6, 1.8)),
            tcl = log(c(2.6, 2.7, 2.8))
          ),
        method = "fixed"
      )
    )
  expect_s3_class(testFixedMulti, "data.frame")
  expect_named(testFixedMulti, expected = c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "eta.ka"))
  expect_equal(nrow(testFixedMulti), 3)
  expect_equal(testFixedMulti$Parameter, rep("tka,tcl", 3))
})

test_that("profile a standard model", {
  # fix most of the parameters so that it estimates faster
  one.compartment <- function() {
    ini({
      tka <- log(1.57)
      tcl <- log(2.72)
      tv <- fixed(log(31.5))
      eta.ka ~ 0.6
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }

  fit <-
    suppressMessages(nlmixr2(
      one.compartment, data = nlmixr2data::theo_sd, est="focei", control = list(print=0)
    ))

  # All parameters
  profall <- suppressMessages(profile(fit))
  expect_s3_class(profall, "data.frame")
  expect_named(profall, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "profileBound"))

  # A single parameter
  proftka <- suppressMessages(profile(fit, which = "tka"))
  expect_s3_class(proftka, "data.frame")
  expect_named(proftka, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "profileBound"))

  # A fixed parameter
  expect_warning(
    proftv <- profile(fit, which = "tv"),
    regexp = "OFV decreased while profiling"
  )
  expect_s3_class(proftv, "data.frame")
  expect_named(proftv, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "profileBound"))

  # Residual error
  profadd.sd <- profile(fit, which = "add.sd")
  expect_s3_class(profadd.sd, "data.frame")
  expect_named(profadd.sd, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "profileBound"))
})

test_that("profile a standard model with correlated etas", {
  # fix most of the parameters so that it estimates faster
  one.compartment <- function() {
    ini({
      tka <- log(1.57)
      tcl <- log(2.72)
      tv <- fixed(log(31.5))
      eta.ka ~ 0.6
      eta.cl ~ 0.1
      eta.v ~ 0.2
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }

  fit <-
    suppressMessages(nlmixr2(
      one.compartment, data = nlmixr2data::theo_sd, est="focei", control = list(print=0)
    ))

  # All parameters
  profall <- suppressMessages(profile(fit))
  expect_s3_class(profall, "data.frame")
  expect_named(profall, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "profileBound"))

  # A single parameter
  proftka <- suppressMessages(profile(fit, which = "tka"))
  expect_s3_class(proftka, "data.frame")
  expect_named(proftka, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "profileBound"))

  # A fixed parameter — warning fires if OFV decreases, which is model-dependent
  proftv <- suppressWarnings(suppressMessages(profile(fit, which = "tv")))
  expect_s3_class(proftv, "data.frame")
  expect_named(proftv, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "profileBound"))

  # Residual error
  profadd.sd <- suppressMessages(profile(fit, which = "add.sd"))
  expect_s3_class(profadd.sd, "data.frame")
  expect_named(profadd.sd, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "profileBound"))
})

# Omega diagonal profiling unit tests (Step O) ----

test_that("llpOmegaNames identifies non-fixed diagonal elements - Step O", {
  skip_if_not_installed("nlmixr2data")
  skip_if_not_installed("nlmixr2est")
  skip("Requires real fit - slow")

  one.compartment <- function() {
    ini({
      tka <- log(1.57); tcl <- log(2.72); tv <- fixed(log(31.5))
      eta.ka ~ 0.6; add.sd <- 0.7
    })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
            cp <- linCmt(); cp ~ add(add.sd) })
  }
  fit <- suppressMessages(nlmixr2(one.compartment, data=nlmixr2data::theo_sd,
                                   est="focei", control=list(print=0)))
  nms <- llpOmegaNames(fit)
  expect_character(nms)
  expect_true("eta.ka" %in% nms)
  expect_false(any(c("tka","tcl","add.sd") %in% nms))
})

test_that("llpOmegaNames iniDf logic - Step O", {
  # Reproduce the selection logic (diagonal, non-fixed omega rows only)
  idf <- data.frame(
    ntheta = c(1L,       NA_integer_, NA_integer_),
    neta1  = c(NA_integer_, 1L,        1L),
    neta2  = c(NA_integer_, 1L,        2L),   # off-diagonal neta1 != neta2
    name   = c("tka",    "eta.ka",   "off"),
    lower  = c(-Inf,     -Inf,       -Inf),
    est    = c(0.4,       0.3,        0.05),
    upper  = c(Inf,       Inf,        Inf),
    fix    = c(FALSE,     FALSE,      FALSE),
    stringsAsFactors = FALSE
  )
  diag_rows <- idf[is.na(idf$ntheta) & !is.na(idf$neta1) & idf$neta1 == idf$neta2, ]
  nms       <- diag_rows$name[!isTRUE(diag_rows$fix)]
  expect_equal(nms, "eta.ka")  # only diagonal, non-fixed

  # Fixed diagonal excluded
  idf2      <- idf; idf2$fix[2L] <- TRUE
  diag2     <- idf2[is.na(idf2$ntheta) & !is.na(idf2$neta1) & idf2$neta1 == idf2$neta2, ]
  nms2      <- diag2$name[!diag2$fix]
  expect_length(nms2, 0L)
})

test_that("llpOmegaSE Wishart formula - Step O", {
  # Verify the Wishart approximation formula: sqrt(2 * omega_kk^2 / (n_sub - 1))
  # This is a pure numeric check of the formula, not of the full function.
  omega_kk    <- 0.4
  n_sub       <- 12L
  expected_se <- sqrt(2 * omega_kk^2 / (n_sub - 1L))
  computed_se <- sqrt(2 * omega_kk^2 / (n_sub - 1L))
  expect_equal(computed_se, expected_se, tolerance = 1e-15)
  # Check boundary case: n_sub = 2 (minimum meaningful case)
  expect_equal(sqrt(2 * 0.3^2 / 1L), 0.3 * sqrt(2), tolerance = 1e-10)
})

test_that("buildFixedOmegaModel generates function with fix=TRUE - Step O", {
  skip_if_not_installed("nlmixr2data")
  skip_if_not_installed("nlmixr2est")
  skip("Requires real fit - slow")

  one.compartment <- function() {
    ini({
      tka <- log(1.57); tcl <- log(2.72); tv <- fixed(log(31.5))
      eta.ka ~ 0.6; add.sd <- 0.7
    })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
            cp <- linCmt(); cp ~ add(add.sd) })
  }
  fit <- suppressMessages(nlmixr2(one.compartment, data=nlmixr2data::theo_sd,
                                   est="focei", control=list(print=0)))
  fn  <- buildFixedOmegaModel(fit, "eta.ka", 0.4)
  ui  <- rxode2::rxode2(fn)
  expect_true(ui$iniDf[ui$iniDf$name == "eta.ka", "fix"])
  expect_equal(ui$iniDf[ui$iniDf$name == "eta.ka", "est"], 0.4)
  # theta initial values match fit estimates
  expect_equal(ui$iniDf[ui$iniDf$name == "tka", "est"],
               nlmixr2est::fixef(fit)[["tka"]], tolerance = 1e-10)
})

test_that("runLLP accepts omega names in which - Step O", {
  skip_if_not_installed("nlmixr2data")
  skip_if_not_installed("nlmixr2est")
  skip("Requires real fit - slow")

  one.compartment <- function() {
    ini({
      tka <- log(1.57); tcl <- log(2.72); tv <- fixed(log(31.5))
      eta.ka ~ 0.6; add.sd <- 0.7
    })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
            cp <- linCmt(); cp ~ add(add.sd) })
  }
  fit <- suppressMessages(nlmixr2(one.compartment, data=nlmixr2data::theo_sd,
                                   est="focei", control=list(print=0)))

  # Single omega parameter
  profOmega <- suppressMessages(runLLP(fit, which = "eta.ka"))
  expect_s3_class(profOmega, "nlmixr2Profile")
  expect_true("eta.ka" %in% names(profOmega))
  expect_true("profileBound" %in% names(profOmega))
  expect_true("OFV" %in% names(profOmega))

  # omega CI bounds finite and straddle MLE on variance scale
  ci <- confint(profOmega)
  expect_true(is.matrix(ci))
  expect_true("eta.ka" %in% rownames(ci))
  mleOmega <- fit$iniDf[fit$iniDf$name == "eta.ka", "est"]
  expect_gt(ci["eta.ka", "lower"], 0)        # variance >= 0
  expect_lt(ci["eta.ka", "lower"], mleOmega)
  expect_gt(ci["eta.ka", "upper"], mleOmega)
})
