skip_on_cran()

# backwardSearch() undefined 'fit' variable fix (Issue 7) ----

test_that("backwardSearch errors with correct message when fitorig is not a nlmixr2 fit", {
  # Before the fix, the function body referenced the undefined variable 'fit'
  # instead of the parameter 'fitorig'. R would look up 'fit' in the parent
  # frame, either finding a stale value from a previous call or throwing the
  # cryptic error "object 'fit' not found". After the fix, passing a non-fit
  # object throws the informative validation message below.
  not_a_fit <- list(a = 1, b = 2)
  expect_error(
    nlmixr2extra:::backwardSearch(
      varsVec   = "cl",
      covarsVec = "WT",
      fitorig   = not_a_fit,
      outputDir = tempdir()
    ),
    regexp = "fitorig.*needs to be a nlmixr2 fit"
  )
})
