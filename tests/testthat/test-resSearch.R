test_that("Test residual search on linearized models", {
  skip_on_cran()
  one.cmpt.adderr <- function() {
    ini({
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            tka <- log(1.56) #  Ka
            eta.cl ~ 0.3
            eta.v ~ 0.1
            eta.ka ~ 0.6
            add.sd <- 0.7
        })
    model({
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd)
        })
  }

  fit <- nlmixr(one.cmpt.adderr, nlmixr2data::theo_sd, est = "focei")
  suppressWarnings({
    linFit <- linearize(fit)
    resRes <- resSearch(linFit)
  })

  expect_named(resRes, c("summary", "originalFit"))
  expect_s3_class(resRes$originalFit, "nlmixr2Linearize")
  # the base fit plus one row per residual model tried
  expect_equal(resRes$summary$search, c("base fit", "prop", "combined2", "combined1"))
  expect_named(resRes$summary, c("OBJF", "AIC", "BIC", "search"))
  expect_true(all(is.finite(resRes$summary$OBJF)))
  # theo_sd is fit best by the additive (base) model; proportional-only is
  # worse (about 117 vs 120)
  expect_lt(
    resRes$summary$OBJF[resRes$summary$search == "base fit"],
    resRes$summary$OBJF[resRes$summary$search == "prop"]
  )
})
