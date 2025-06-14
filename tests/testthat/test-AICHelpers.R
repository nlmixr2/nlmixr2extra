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

  fitEmaxBoundaryIssue <- nlmixr2est::nlmixr2(modEmax, data = d_noec50, est = "focei", control = list(print = 0))
  fitStep <- nlmixr2est::nlmixr2(modStep, data = d_noec50, est = "focei", control = list(print = 0))
  fitLinear <- nlmixr2est::nlmixr2(modLinear, data = d_noec50, est = "focei", control = list(print = 0))
  # Boundary conditions are accurately caught
  expect_true(isBoundaryFit(fitEmaxBoundaryIssue))
  expect_false(isBoundaryFit(fitStep))
  # The correct model is returned
  expect_message(
    minAICfit <- getMinAICFit(fitEmaxBoundaryIssue, fitStep, fitLinear),
    regexp = "Removing model with a parameter at the boundary"
  )
  expect_equal(minAICfit, fitLinear)

  # getMinAICFit gives NULL when expected
  expect_warning(
    getMinAICFit(fitEmaxBoundaryIssue, fitEmaxBoundaryIssue),
    regexp = "No model fits in input after excluding models"
  )
  expect_warning(
    getMinAICFit(list()),
    regexp = "No model fits in input list"
  )

  # listModelsTested ----
  expect_equal(
    listModelsTested(list(Emax = fitEmaxBoundaryIssue, Step = fitStep, Linear = fitLinear), caption = "Listing of models tested."),
    structure(
      data.frame(
        Description = c("Emax", "Step", "Linear"),
        AIC = c(-1832.41345872027, 545.60562057833, 521.880978234241),
        dAIC = c("-", "23.72", "0"),
        Exclude = c("parameter at boundary", "", "")
      ),
      caption = "Listing of models tested. Abbreviations: AIC = Akaike's Information Criterion; dAIC = change from minimum AIC"
    ),
    tolerance = 0.1
  )
})
