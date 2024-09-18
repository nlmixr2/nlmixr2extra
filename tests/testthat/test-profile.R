test_that("profileNlmixr2FitDataEstInitial", {
  # Must have one-row input
  expect_error(
    profileNlmixr2FitDataEstInitial(estimates = data.frame(A = 1:2))
  )

  expect_equal(
    profileNlmixr2FitDataEstInitial(
      estimates = data.frame(A = 1),
      which = "A",
      ofvIncrease = 1.92,
      rseTheta = c(A=100),
      lower = -100, upper = 200
    ),
    c(-0.92, 2.92)
  )
  # Bounds are respected
  expect_equal(
    profileNlmixr2FitDataEstInitial(
      estimates = data.frame(A = 1),
      which = "A",
      ofvIncrease = 1.92,
      rseTheta = c(A=100),
      lower = 0, upper = 200
    ),
    c(sqrt(.Machine$double.eps), 2.92)
  )
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
  expect_named(profall, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "eta.ka", "profileBound"))

  # A single parameter
  proftka <- suppressMessages(profile(fit, which = "tka"))
  expect_s3_class(proftka, "data.frame")
  expect_named(proftka, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "eta.ka", "profileBound"))

  # A fixed parameter
  expect_warning(
    proftv <- profile(fit, which = "tv"),
    regexp = "OFV decreased while profiling"
  )
  expect_s3_class(proftv, "data.frame")
  expect_named(proftv, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "eta.ka", "profileBound"))

  # Residual error
  profadd.sd <- profile(fit, which = "add.sd")
  expect_s3_class(profadd.sd, "data.frame")
  expect_named(profadd.sd, c("Parameter", "OFV", "tka", "tcl", "tv", "add.sd", "eta.ka", "profileBound"))
})
