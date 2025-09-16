test_that("isBoundaryFit, getMinAICFit, listModelsTested work", {
  # fit models for testing
  d_noec50 <-
    data.frame(
      conc = c(rep(0, 10), rep(1:20, each = 10)),
      DV = c(rnorm(n = 10, mean = 1, sd = 1e-5), rnorm(n = 200, mean = 5, sd = 1e-5)),
      TIME = 0
    )

  modEmax <- function() {
    ini({
      e0 = 1
      emax = 5
      ec50 = c(0, 1.1)
      addSd = 0.5
    })
    model({
      effect <- e0 + emax*conc/(ec50 + conc)
      effect ~ add(addSd)
    })
  }

  modStep <- function() {
    ini({
      e0 = 1
      emax = 5
      addSd = 1e-5
    })
    model({
      effect <- e0 + emax*(conc > 0)
      effect ~ add(addSd)
    })
  }

  modLinear <- function() {
    ini({
      e0 = 1
      slope = 5
      addSd = 1
    })
    model({
      effect <- e0 + slope*conc
      effect ~ add(addSd)
    })
  }

  fitEmaxBoundaryIssue <-
    suppressMessages(
      nlmixr2est::nlmixr2(modEmax, data = d_noec50, est = "focei", control = list(print = 0))
    )
  fitStep <-
    suppressMessages(
      nlmixr2est::nlmixr2(modStep, data = d_noec50, est = "focei", control = list(print = 0))
    )
  fitLinear <-
    suppressMessages(
      nlmixr2est::nlmixr2(modLinear, data = d_noec50, est = "focei", control = list(print = 0))
    )
  fitError <- try(stop("Model did not converge"), silent = TRUE)
  # Boundary conditions are accurately caught
  expect_true(isBoundaryFit(fitEmaxBoundaryIssue))
  expect_false(isBoundaryFit(fitStep))
  # The correct model is returned
  expect_message(
    minAICfit <- getMinAICFit(fitEmaxBoundaryIssue, fitStep, fitLinear, fitError),
    regexp = "Removing model with a parameter at the boundary"
  )
  expect_equal(minAICfit, fitLinear)

  # getMinAICFit gives NULL when expected
  expect_warning(
    getMinAICFit(fitEmaxBoundaryIssue, fitEmaxBoundaryIssue, fitError),
    regexp = "No model fits in input after excluding models"
  )
  expect_warning(
    getMinAICFit(list()),
    regexp = "No model fits in input list"
  )

  # listModelsTested ----
  tabTested <- listModelsTested(list(Emax = fitEmaxBoundaryIssue, Step = fitStep, Linear = fitLinear, "Model error" = fitError), caption = "Listing of models tested.")
  expect_equal(
    tabTested$Description,
    c("Emax", "Step", "Linear", "Model error")
  )
  # system differences may cause different AIC values over time, so cannot test
  # those directly
  expect_type(tabTested$AIC, "double")
  expect_equal(tabTested$AIC[4], NA_real_)
  expect_equal(tabTested$dAIC[c(1, 4)], c("-", "-"))
  expect_equal(tabTested$Exclude, c("parameter at boundary", "", "", ""))
  expect_equal(
    attr(tabTested, "caption"),
    "Listing of models tested. Abbreviations: AIC = Akaike's Information Criterion; dAIC = change from minimum AIC"
  )
})
