skip_on_cran()

test_that("Add covariates and lasso string to ui ", {
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  ui1 <- nlmixr(one.cmt)
  ui <- ui1
  varsVec <- c("ka","cl","v")
  covarsVec <- c("WT","BMI")

  funstring1 <- intersect((.lassoUicovariate(ui,varsVec,covarsVec))$iniDf$name,
                          c("cov_WT_ka","cov_WT_cl","cov_WT_v","cov_BMI_ka","cov_BMI_cl","cov_BMI_v"))
  funstring2 <- c("cov_WT_ka","cov_WT_cl","cov_WT_v","cov_BMI_ka","cov_BMI_cl","cov_BMI_v")
  funstring3 <- .lassoUicovariate(ui,varsVec,covarsVec)$funTxt
  funstring4 <- "tvalue <- 0.05\nabssum <- sum(abs(cov_BMI_ka) + abs(cov_WT_ka) + abs(cov_WT_cl) + abs(cov_BMI_cl) + abs(cov_BMI_v) + abs(cov_WT_v))\nratio <- abssum/tvalue\nfactor <- exp(1 - ratio)\nka <- exp(tka + eta.ka + cov_WT_ka * factor * WT + cov_BMI_ka * factor * BMI)\ncl <- exp(tcl + eta.cl + cov_BMI_cl * factor * BMI + cov_WT_cl * factor * WT)\nv <- exp(tv + eta.v + cov_WT_v * factor * WT + cov_BMI_v * factor * BMI)\nlinCmt() ~ add(add.sd)"
  expect_equal(funstring1, funstring2)
  expect_equal(funstring3, funstring4)
})

test_that("Add covariates and lasso string to ui ", {
  two.compartment <- function() {
    ini({
      tcl <- log(53.4)
      tv1 <- log(73.6)
      tv2 <- log(320)
      tQ <- log(191)
      eta.cl ~ 0.43^2
      eta.v1 ~ 0.48^2
      eta.v2 ~ 0.49^2
      eta.Q ~ 0.36^2
      prop.sd <- 0.44^2
    })
    model({
      cl <- exp(tcl + eta.cl)
      v1 <- exp(tv1 + eta.v1)
      v2 <- exp(tv2 + eta.v2)
      Q <- exp(tQ+eta.Q)
      linCmt() ~ prop(prop.sd)
    })
  }

  ui2 <- nlmixr(two.compartment)
  ui <- ui2
  varsVec <- "ka"
  covarsVec <- c("WT","BMI")
  expect_error(.lassoUicovariate(ui,varsVec,covarsVec))

})



