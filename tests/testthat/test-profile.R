test_that("profileNlmixr2FitDatEstBound", {
  estimatesSetup <- data.frame(A = 0:4)

  # Outside range, bound not present
  expect_equal(
    profileNlmixr2FitDatEstBound(estimates = estimatesSetup, which="A", newEst = 4.5, bound = 4.25, direction = 1, tol = 0.001, atol = 0.005),
    4.25
  )
  # Inside range, bound not present
  expect_equal(
    profileNlmixr2FitDatEstBound(estimates = estimatesSetup, which="A", newEst = 4.5, bound = 5, direction = 1, tol = 0.001, atol = 0.005),
    4.5
  )
  # Inside range, bound already present
  expect_equal(
    profileNlmixr2FitDatEstBound(estimates = estimatesSetup, which="A", newEst = 0.5, bound = 1, direction = 1, tol = 0.001, atol = 0.005),
    0.5
  )
  # Outside range, bound already present
  expect_equal(
    profileNlmixr2FitDatEstBound(estimates = estimatesSetup, which="A", newEst = 1.5, bound = 1, direction = 1, tol = 0.001, atol = 0.005),
    0.9995
  )
  # Non-zero bound is already present, coming from below
  expect_equal(
    profileNlmixr2FitDatEstBound(estimates = estimatesSetup, which="A", newEst = 1, bound = 1, direction = 1, tol = 0.001, atol = 0.005),
    0.9995
  )
  # Non-zero bound is already present, coming from above
  expect_equal(
    profileNlmixr2FitDatEstBound(estimates = estimatesSetup, which="A", newEst = 1, bound = 1, direction = -1, tol = 0.001, atol = 0.005),
    1.0005
  )
  # Zero bound is already present, coming from below
  expect_equal(
    profileNlmixr2FitDatEstBound(estimates = estimatesSetup, which="A", newEst = 1, bound = 0, direction = 1, tol = 0.001, atol = 0.005),
    -0.0025
  )
  # Zero bound is already present, coming from above
  expect_equal(
    profileNlmixr2FitDatEstBound(estimates = estimatesSetup, which="A", newEst = 0, bound = 0, direction = -1, tol = 0.001, atol = 0.005),
    0.0025
  )
})

test_that("profileNlmixr2FitDataEstInitial", {
  # Must have one-row input
  expect_error(
    profileNlmixr2FitDataEstInitial(estimates = data.frame(A = 1:2))
  )

  expect_equal(
    profileNlmixr2FitDataEstInitial(estimates = data.frame(A = 1), which = "A", normQuantile = 1.96, rse_theta = c(A=100)),
    c(-0.96, 2.96)
  )
})

test_that("profileNlmixr2MultiParam", {
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
      one.compartment, data = theo_sd, est="focei", control = list(print=0)
    ))

  testMultiParam <-
    suppressMessages(
      profileNlmixr2MultiParam(fit, which = data.frame(tka = log(c(1.4, 1.6, 1.8))))
    )
  expect_s3_class(testMultiParam, "data.frame")
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
  profall <- profile(fit)
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
})
