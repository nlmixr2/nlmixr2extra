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

test_that("profileNlmixr2FitDataNewEst", {
  # Go the range in the down direction
  estimateSetup <- data.frame(A = 1:2, OFV = c(1, 5))
  expect_equal(
    profileNlmixr2FitDataNewEst(
      estimates = estimateSetup, which = "A", direction = -1, bound = Inf, ofvIncrease = 3, method = "linapprox"
    ),
    0
  )
  # Go the range in the up direction
  estimateSetup <- data.frame(A = 1:2, OFV = c(1, 5))
  expect_equal(
    profileNlmixr2FitDataNewEst(
      estimates = estimateSetup, which = "A", direction = 1, bound = Inf, ofvIncrease = 3, method = "linapprox"
    ),
    2
  )
})
