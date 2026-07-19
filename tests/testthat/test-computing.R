skip_on_cran()

test_that("Function to return pop mean, pop std of a given covariate (.popMeanStd)", {
  funstring1 <- .popMeanStd(nlmixr2data::theo_sd,"WT")[[1]]
  funstring2 <- .popMeanStd(nlmixr2data::theo_sd,"WT")[[2]]

  expect_equal(funstring1, 69.58333, tolerance=0.1)
  expect_equal(funstring2, 9.098519, tolerance=0.1)
})

test_that("Function to return pop mean, pop std of a given covariate", {
  expect_error(.popMeanStd(theo_sd,"BMI"))
})

test_that("Function to return normalization column of a covariates (.normalizeDf)", {
  funstring1 <- mean(.normalizeDf(nlmixr2data::theo_sd,"WT")$normalized_WT)
  funstring2 <- sd(.normalizeDf(nlmixr2data::theo_sd,"WT")$normalized_WT)

  expect_equal(funstring1,0,tolerance=0.2)
  expect_equal(funstring2,1,tolerance=0.2)
})

test_that("Function to return normalization column of a covariates", {
  expect_error(.normalizeDf(nlmixr2data::theo_sd,"BMI"))
})

# foldgen() loop fix (Issue 3) ----

test_that("foldgen assigns folds to ALL strata classes, not just the first", {
  # With the bug (1:seq_along(numInClass)), R's ':' operator used only the
  # first element of seq_along()'s vector result, so the loop ran exactly
  # once. Only the first stratum class ever received fold assignments; all
  # others remained as 0L. After the fix (seq_along(numInClass)), every
  # class gets complete fold coverage.
  set.seed(42)
  n <- 40
  df <- data.frame(
    ID    = seq_len(n),
    TIME  = 0,
    DV    = stats::rnorm(n),
    AMT   = 0,
    EVID  = 0,
    CMT   = 1,
    STRAT = rep(c("A", "B"), each = n / 2)
  )
  result <- foldgen(df, nfold = 5, stratVar = "STRAT")
  # All rows must have a valid fold number
  expect_true(all(result$fold >= 1L & result$fold <= 5L))
  # Both strata classes must have rows assigned to all 5 folds
  folds_A <- result$fold[result$STRAT == "A"]
  folds_B <- result$fold[result$STRAT == "B"]
  expect_equal(sort(unique(folds_A)), 1:5)
  expect_equal(sort(unique(folds_B)), 1:5)
})

# extractVars() loop fix (Issue 5) ----

# optimUnisampling() bug fix ----

test_that("optimUnisampling respects N argument (was hardcoded to 1000)", {
  set.seed(42)
  result50  <- optimUnisampling(xvec = c(6.7, 140), N = 50,  medValue = 71.7)
  set.seed(42)
  result200 <- optimUnisampling(xvec = c(6.7, 140), N = 200, medValue = 71.7)

  expect_length(result50, 50)
  expect_length(result200, 200)
})

test_that("optimUnisampling respects floorT=FALSE (was always TRUE in recursive path)", {
  set.seed(42)
  result <- optimUnisampling(xvec = c(6.7, 140), N = 200, medValue = 71.7, floorT = FALSE)

  expect_length(result, 200)
  # With floorT=FALSE the values should not all be integers
  expect_false(all(result == floor(result)))
})

# extractVars() loop fix (Issue 5) ----

test_that("extractVars checks all message strings, not just the last", {
  # With the bug (for (i in length(res))), only res[[length(res)]] was ever
  # inspected. A non-empty message in an earlier element was silently ignored
  # and the function incorrectly returned "". After the fix (seq_along(res)),
  # every element is checked.
  fitlist <- list(
    list(message = ""),
    list(message = "converged"),
    list(message = "")
  )
  result <- nlmixr2extra:::extractVars(fitlist, id = "message")
  expect_false(identical(result, ""))
  expect_true("converged" %in% result)
})
