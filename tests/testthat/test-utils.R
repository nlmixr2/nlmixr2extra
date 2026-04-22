skip_on_cran()

.cur <- loadNamespace("nlmixr2extra")

# =============================================================================
# setQuietFastControl
# =============================================================================

test_that("setQuietFastControl: sets print=0", {
  ctl <- list(
    print = 100,
    covMethod = "r,s",
    calcTables = TRUE,
    compress = TRUE
  )
  out <- .cur$setQuietFastControl(ctl)
  expect_equal(out$print, 0L)
})

test_that("setQuietFastControl: sets covMethod=0", {
  ctl <- list(
    print = 100,
    covMethod = "r,s",
    calcTables = TRUE,
    compress = TRUE
  )
  out <- .cur$setQuietFastControl(ctl)
  expect_equal(out$covMethod, 0L)
})

test_that("setQuietFastControl: sets calcTables=FALSE", {
  ctl <- list(
    print = 100,
    covMethod = "r,s",
    calcTables = TRUE,
    compress = TRUE
  )
  out <- .cur$setQuietFastControl(ctl)
  expect_false(out$calcTables)
})

test_that("setQuietFastControl: sets compress=FALSE", {
  ctl <- list(
    print = 100,
    covMethod = "r,s",
    calcTables = TRUE,
    compress = TRUE
  )
  out <- .cur$setQuietFastControl(ctl)
  expect_false(out$compress)
})

test_that("setQuietFastControl: preserves fields it does not touch", {
  ctl <- list(
    print = 100,
    covMethod = "r,s",
    calcTables = TRUE,
    compress = TRUE,
    eval.max = 999,
    grad.eps = 0.001
  )
  out <- .cur$setQuietFastControl(ctl)
  expect_equal(out$eval.max, 999)
  expect_equal(out$grad.eps, 0.001)
})

# =============================================================================
# .plap — parallel lapply wrapper
# =============================================================================

test_that(".plap: applies FUN to every element and preserves order", {
  res <- .cur$.plap(1:5, function(x) x * 10)
  expect_equal(unlist(res), c(10, 20, 30, 40, 50))
})

test_that(".plap: returns a list", {
  res <- .cur$.plap(1:3, function(x) x)
  expect_type(res, "list")
})

test_that(".plap: passes extra ... arguments to FUN", {
  res <- .cur$.plap(1:3, function(x, y) x + y, y = 100)
  expect_equal(unlist(res), c(101, 102, 103))
})

test_that(".plap: empty input returns empty list", {
  res <- .cur$.plap(integer(0), function(x) x)
  expect_equal(res, list())
})

test_that(".plap: .label function is accepted without error", {
  expect_no_error(
    .cur$.plap(1:4, function(x) x * 2, .label = function(x) paste("item", x))
  )
})

test_that(".plap: result length equals input length", {
  n <- 7
  res <- .cur$.plap(seq_len(n), function(x) x)
  expect_length(res, n)
})

test_that(".plap: FUN can return complex objects", {
  res <- .cur$.plap(1:3, function(x) list(val = x, sq = x^2))
  expect_equal(res[[2]]$val, 2)
  expect_equal(res[[3]]$sq, 9)
})

test_that(".plap: future.apply path preserves order", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  plan_before <- future::plan()
  on.exit(future::plan(plan_before), add = TRUE)
  future::plan("sequential")

  res <- .cur$.plap(1:5, function(x) x * 10)

  expect_equal(unlist(res), c(10, 20, 30, 40, 50))
})

# =============================================================================
# .withWorkerPlan
# =============================================================================

test_that(".withWorkerPlan: NULL workers evaluates expr and returns its value", {
  expect_equal(.cur$.withWorkerPlan(NULL, 1 + 1), 2)
})

test_that(".withWorkerPlan: NULL workers works without future installed", {
  # NULL path short-circuits before any requireNamespace("future") check
  result <- .cur$.withWorkerPlan(NULL, "hello")
  expect_equal(result, "hello")
})

test_that(".withWorkerPlan: NULL workers passes side effects through", {
  x <- 0L
  .cur$.withWorkerPlan(NULL, {
    x <- 99L
  })
  expect_equal(x, 99L)
})

test_that(".withWorkerPlan: workers=1 evaluates expr and returns value", {
  skip_if_not_installed("future")
  expect_equal(.cur$.withWorkerPlan(1L, sum(1:10)), 55L)
})

test_that(".withWorkerPlan: workers=1 restores original plan on clean exit", {
  skip_if_not_installed("future")
  plan_before <- class(future::plan())
  on.exit(future::plan("sequential"), add = TRUE) # safety net
  .cur$.withWorkerPlan(1L, NULL)
  expect_equal(class(future::plan()), plan_before)
})

test_that(".withWorkerPlan: plan restored even when expr throws an error", {
  skip_if_not_installed("future")
  plan_before <- class(future::plan())
  on.exit(future::plan("sequential"), add = TRUE)
  try(.cur$.withWorkerPlan(1L, stop("intentional error")), silent = TRUE)
  expect_equal(class(future::plan()), plan_before)
})

test_that(".withWorkerPlan: error from expr is propagated to caller", {
  skip_if_not_installed("future")
  expect_error(.cur$.withWorkerPlan(1L, stop("boom")), "boom")
})

test_that(".withWorkerPlan: workers='auto' evaluates expr without error", {
  skip_if_not_installed("future")
  expect_no_error(.cur$.withWorkerPlan("auto", TRUE))
})

test_that(".withWorkerPlan: workers='auto' returns expr value", {
  skip_if_not_installed("future")
  result <- .cur$.withWorkerPlan("auto", 42L)
  expect_equal(result, 42L)
})

test_that(".withWorkerPlan: rejects invalid workers values", {
  bad_workers <- list(0, -1, NA_real_, NaN, Inf, c(1, 2), 1.5, "bad")

  for (workers in bad_workers) {
    expect_error(.cur$.withWorkerPlan(workers, TRUE), "workers")
  }
})

test_that(".withWorkerPlan: workers=2 restores original plan", {
  skip_if_not_installed("future")
  plan_before <- class(future::plan())
  on.exit(future::plan("sequential"), add = TRUE)

  .cur$.withWorkerPlan(2L, NULL)

  expect_equal(class(future::plan()), plan_before)
})
