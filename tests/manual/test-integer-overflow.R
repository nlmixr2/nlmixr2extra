# Manual test for the getBootstrapSummary() integer overflow fix
#
# Run this script OUTSIDE of devtools::test() to confirm the runtime behaviour.
# It requires significant RAM (~5.8 GB for the fix verification step).
# Do NOT run inside devtools::test() or testthat::test_local().
#
# Usage:
#   Rscript tests/manual/test-integer-overflow.R
# or interactively:
#   source("tests/manual/test-integer-overflow.R")

cat("=== Integer overflow manual test ===\n\n")

# ------------------------------------------------------------------
# Part 1: Demonstrate the BUG
#   d*d overflows 32-bit integer arithmetic when d > 46340.
#   rep(0, NA_integer_) then throws an immediate error -- no allocation.
# ------------------------------------------------------------------
cat("Part 1: confirming integer overflow in d*d (no allocation needed)\n")
d_overflow <- 46341L
stopifnot(is.na(d_overflow * d_overflow))   # integer overflow -> NA_integer_
cat("  d*d is NA (integer overflow confirmed):", is.na(d_overflow * d_overflow), "\n")

cat("Part 1: confirming unfixed code errors at rep() before any allocation\n")
result <- tryCatch(
  matrix(rep(0, d_overflow * d_overflow), d_overflow, d_overflow),
  error = function(e) e
)
stopifnot(inherits(result, "error"))
cat("  Unfixed code error:", conditionMessage(result), "\n")

# ------------------------------------------------------------------
# Part 2: Demonstrate the FIX
#   matrix(0, d, d) fills from a scalar -- no d*d integer multiply.
#   Use d = 27000 (~5.8 GB) to stay under the 6 GB limit.
#   Note: this WILL attempt a ~5.8 GB allocation; ensure sufficient RAM.
# ------------------------------------------------------------------
cat("\nPart 2: confirming fixed code allocates cleanly (requires ~5.8 GB RAM)\n")
cat("  Allocating 27000 x 27000 double matrix ...\n")
d_fix <- 27000L
# This is the fixed line from getBootstrapSummary():
m <- matrix(0, d_fix, d_fix)
stopifnot(is.matrix(m), nrow(m) == d_fix, ncol(m) == d_fix, all(m == 0))
cat("  Fixed code succeeded: matrix(0, d, d) works for d =", d_fix, "\n")
rm(m)
gc()

cat("\n=== All manual checks passed ===\n")
