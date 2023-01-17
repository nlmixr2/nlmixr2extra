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
