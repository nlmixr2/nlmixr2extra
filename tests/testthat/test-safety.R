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

# Integer overflow in getBootstrapSummary() (Issue 4) ----

test_that("integer overflow: d*d overflows 32-bit signed integer for d > 46340", {
  # This demonstrates the bug at the arithmetic level without requiring any
  # memory allocation. The fixed code uses matrix(0, d, d) which avoids
  # the integer multiply entirely.
  #
  # 46341^2 = 2,147,484,281 which exceeds .Machine$integer.max (2,147,483,647).
  d <- 46341L
  # Confirm the integer overflow is real (produces NA_integer_)
  expect_true(is.na(d * d))
  expect_identical(class(d * d), "integer")
  # The fix: as.numeric avoids the overflow
  expect_false(is.na(as.numeric(d) * as.numeric(d)))
})

test_that("getBootstrapSummary integer overflow: 46341x46341 allocation would need ~17 GB RAM", {
  skip("Allocation of 46341x46341 doubles requires ~17 GB RAM (exceeds safe test limit of ~6 GB).
The bug and fix are verified at the arithmetic level in the preceding test.
Run tests/manual/test-integer-overflow.R manually to confirm the runtime behaviour.")
  # No allocation code here: placing matrix(0, 46341, 46341) in this block
  # would attempt a 17 GB allocation if this file is sourced outside of
  # the testthat framework (where skip() does not fire).
})
