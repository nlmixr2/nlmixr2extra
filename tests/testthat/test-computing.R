test_that("Function to return pop mean, pop std of a given covariate", {
  funstring1 <- .popMeanStd(nlmixr2data::theo_sd,"WT")[[1]]
  funstring2 <- .popMeanStd(nlmixr2data::theo_sd,"WT")[[2]]
  
  expect_equal(funstring1,69.58333,tolerance=0.1)
  expect_equal(funstring2,9.098519,tolerance=0.1)
})

test_that("Function to return pop mean, pop std of a given covariate", {
  expect_error(.popMeanStd(theo_sd,"BMI"))
})

# ==== .normalizeDf

test_that("Function to return normalization column of a covariates", {
  funstring1 <- mean(.normalizeDf(nlmixr2data::theo_sd,"WT")$normalized_WT)
  funstring2 <- sd(.normalizeDf(nlmixr2data::theo_sd,"WT")$normalized_WT)
  
  expect_equal(funstring1,0,tolerance=0.2)
  expect_equal(funstring2,1,tolerance=0.2)
})

test_that("Function to return normalization column of a covariates", {
  expect_error(.normalizeDf(nlmixr2data::theo_sd,"BMI"))
})

test_that("single column covarsVec works (#15)", {
  expect_equal(
    normalizedData(data = data.frame(ID=1:3, A=factor(1:3)), covarsVec = "A"),
    data.frame(ID=1:3, A=factor(1:3))
  )
})
