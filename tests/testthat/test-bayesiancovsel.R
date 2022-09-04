test_that(".horseshoePrior", {
  v1 <- .horseshoePrior(5)
  expect_equal(
    v1[[1]]$prior,
    c(
      "horseshoe(df = 1, scale_global = tau0, df_global = 1)",
      "normal(0, 10)"
    )
  )
})
