# Tests for thread-safety and memory-safety fixes

# preCondInv zero-eigenvalue guard (Issue 1) ----

test_that("preCondInv errors on singular matrix (zero eigenvalue)", {
  # A singular matrix has a zero eigenvalue. Before the fix, preCondInv()
  # would silently return Inf values (1/0) in the preconditioner. After
  # the fix it errors with an informative message.
  singular_mat <- matrix(c(1, 1, 1, 1), 2, 2)
  # Confirm the matrix is indeed singular
  expect_equal(det(singular_mat), 0)
  expect_error(
    preCondInv(singular_mat),
    regexp = "singular or near-singular"
  )
})

test_that("preCondInv errors on near-singular matrix (eigenvalue < 1e-10)", {
  # A near-singular matrix has an eigenvalue smaller than the tolerance
  # (1e-10). Inverting it would produce numerically meaningless values.
  eps <- 1e-11
  near_singular <- matrix(c(1, 1, 1, 1 + eps), 2, 2)
  expect_error(
    preCondInv(near_singular),
    regexp = "singular or near-singular"
  )
})

test_that("preCondInv still works correctly on a well-conditioned matrix", {
  # Confirm that the zero-eigenvalue guard does not break the normal path.
  well_conditioned <- diag(c(2.0, 3.0, 5.0))
  result <- preCondInv(well_conditioned)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3L, 3L))
  expect_true(all(is.finite(result)))
})

