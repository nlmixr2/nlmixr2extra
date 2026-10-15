.cur <- loadNamespace("nlmixr2extra")

# patients from CUI-MET (Zhang et al. 2025) Table 2 example 1 probabilities
.cuimetData <- function(n, seed = 1) {
  withr::with_seed(seed, {
    .p <- data.frame(dose = 1:5,
                     tox = c(0.038, 0.058, 0.106, 0.232, 0.533),
                     eff = c(0.022, 0.119, 0.343, 0.568, 0.681),
                     tol = c(0.147, 0.206, 0.280, 0.368, 0.466))
    .d <- .p[rep(1:5, each = n), ]
    .d$id <- seq_len(nrow(.d))
    for (.v in c("tox", "eff", "tol")) .d[[.v]] <- stats::rbinom(nrow(.d), 1, .d[[.v]])
    rownames(.d) <- NULL
    .d
  })
}

.cuimetEps <- function() {
  list(cuiBinary(tox, direction = "lower"), cuiBinary(eff), cuiBinary(tol))
}

test_that("the point estimate is the observed rate and matches cuiFromValues()", {
  .d <- .cuimetData(30)
  .r <- cuiDoseResponse(.d, .cuimetEps(), dose = dose, id = id, nBoot = 0,
                        weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
  .obs <- as.data.frame(.r, type = "observed")
  expect_equal(.obs$rate[.obs$endpoint == "eff"],
               as.vector(tapply(.d$eff, .d$dose, mean)))
  expect_equal(.obs$n, rep(30, 15))
  .v <- data.frame(regimen = .obs$regimen, endpoint = .obs$endpoint, value = .obs$rate)
  .f <- cuiFromValues(.v, .cuimetEps(), weights = c(tox = 0.2, eff = 0.5, tol = 0.3),
                      regimens = cuiRegimens(dose = 1:5))
  expect_equal(.r$summary$est, .f$summary$est)
  expect_true(all(is.na(.r$summary$selFreq)))
  expect_output(print(.r), "Dose-response: empirical")
})

test_that("large trials recover the CUI-MET Table 2 recommendations", {
  .d <- .cuimetData(4000)
  .r <- cuiDoseResponse(.d, .cuimetEps(), dose = dose, nBoot = 20, seed = 3,
                        weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
  expect_equal(.r$summary$est, c(0.248, 0.310, 0.434, 0.548, 0.574), tolerance = 0.02)
  expect_equal(.r$summary$regimen[.r$summary$best], "5")
  .w2 <- update(.r, weights = c(tox = 0.3, eff = 0.2, tol = 0.5))
  expect_equal(.w2$summary$regimen[.w2$summary$best], "4")
  expect_equal(.w2$summary$regimen[.w2$summary$mostSelected], "4")
})

test_that("the bootstrap resamples whole patients within each dose", {
  .d <- .cuimetData(25)
  # eff is a copy of tox: if patients stay together the rates agree in every draw
  .d$eff <- .d$tox
  .r <- cuiDoseResponse(.d, .cuimetEps(), dose = "dose", nBoot = 50, seed = 2)
  .u <- as.data.frame(.r, type = "utility")
  .tox <- .u[.u$endpoint == "tox", ]
  .eff <- .u[.u$endpoint == "eff", ]
  expect_equal(.tox$value, .eff$value)
  # stratified: each dose keeps 25 patients, so rates are multiples of 1/25
  expect_equal(.u$value * 25, round(.u$value * 25))
  expect_equal(length(unique(.u$draw)), 51)
  expect_equal(sum(.r$summary$selFreq), 1)
  # unequal dose groups keep their own sizes
  .un <- .d[!(.d$dose == 1 & .d$id %% 5 != 0), ]
  .ru <- cuiDoseResponse(.un, .cuimetEps(), dose = dose, nBoot = 20, seed = 2)
  .uu <- as.data.frame(.ru, type = "utility")
  .d1 <- .uu$value[.uu$regimen == "1"]
  expect_equal(.d1 * 5, round(.d1 * 5))
})

test_that("seeds are reproducible and restore the caller's stream", {
  .d <- .cuimetData(20)
  .a <- cuiDoseResponse(.d, .cuimetEps(), dose = dose, nBoot = 30, seed = 10)
  .b <- cuiDoseResponse(.d, .cuimetEps(), dose = dose, nBoot = 30, seed = 10)
  .c <- cuiDoseResponse(.d, .cuimetEps(), dose = dose, nBoot = 30, seed = 11)
  expect_equal(.a$summary, .b$summary)
  expect_false(isTRUE(all.equal(.a$cui, .c$cui)))
  set.seed(99)
  .before <- stats::runif(1)
  set.seed(99)
  .x <- cuiDoseResponse(.d, .cuimetEps(), dose = dose, nBoot = 5, seed = 1)
  expect_equal(stats::runif(1), .before)
})

test_that("columns are given as names or strings", {
  .d <- .cuimetData(10)
  names(.d)[names(.d) == "dose"] <- "DOSE"
  .a <- cuiDoseResponse(.d, list(cuiBinary("tox", direction = "lower"), cuiBinary(eff)),
                        dose = DOSE, id = "id", nBoot = 0)
  .col <- "DOSE"
  .b <- cuiDoseResponse(.d, list(cuiBinary(.name = "tox", direction = "lower"),
                                 cuiBinary(.name = "eff")),
                        .dose = .col, .id = "id", nBoot = 0)
  expect_equal(.a$summary, .b$summary)
  expect_error(cuiDoseResponse(.d, cuiBinary(eff), nBoot = 0), "'dose' is required")
  expect_output(print(cuiBinary(tox, direction = "lower")), "empirical")
})

test_that("inputs are validated", {
  .d <- .cuimetData(10)
  expect_error(cuiDoseResponse(.d, cuiBinary(nope), dose = dose, nBoot = 0), "not in 'data'")
  expect_error(cuiDoseResponse(.d, list(cuiEndpoint(eff)), dose = dose, nBoot = 0), "cuiBinary")
  .bad <- .d
  .bad$eff[1] <- 2
  expect_error(cuiDoseResponse(.bad, cuiBinary(eff), dose = dose, nBoot = 0), "0/1")
  .bad$eff[1] <- NA
  expect_error(cuiDoseResponse(.bad, cuiBinary(eff), dose = dose, nBoot = 0), "0/1")
  .dup <- .d
  .dup$id[2] <- .dup$id[1]
  expect_error(cuiDoseResponse(.dup, cuiBinary(eff), dose = dose, id = id, nBoot = 0), "unique")
  .chr <- .d
  .chr$dose <- as.character(.chr$dose)
  expect_error(cuiDoseResponse(.chr, cuiBinary(eff), dose = dose, nBoot = 0), "numeric")
  .lgl <- .d
  .lgl$eff <- .lgl$eff == 1
  expect_equal(cuiDoseResponse(.lgl, cuiBinary(eff), dose = dose, nBoot = 0)$summary,
               cuiDoseResponse(.d, cuiBinary(eff), dose = dose, nBoot = 0)$summary)
})

test_that("failed draws are recorded and minSuccess is enforced", {
  .d <- .cuimetData(10)
  .y <- as.matrix(.d[, c("tox", "eff")])
  .eps <- list(cuiBinary(tox), cuiBinary(eff))
  .calls <- 0
  .flaky <- function(y, dose) {
    .calls <<- .calls + 1
    if (.calls %in% c(3, 5)) stop("no convergence")
    .cur$.cuiDrFit(y, dose, 1:5, .eps)
  }
  .run <- .cur$.cuiDrRun(.y, .d$dose, 1:5, .eps, nBoot = 10, fitFun = .flaky)
  expect_equal(.run$draws$status, c("ok", "ok", "failed", "ok", "failed", rep("ok", 6)))
  expect_equal(.run$draws$message[3], "no convergence")
  expect_false(any(.run$values$draw %in% c(2, 4)))
  .r <- .cur$.cuiNew(.run$values, .eps, cuiRegimens(dose = 1:5), NULL, "additive", NULL,
                     0.95, 1e-8, .run$draws, list())
  expect_equal(sum(.r$summary$selFreq), 1)
  expect_output(print(.r), "8 ok, 2 failed")
})

test_that("too few successful draws is an error", {
  .d <- .cuimetData(10)
  .real <- .cur$.cuiDrFit
  .calls <- 0
  local_mocked_bindings(.cuiDrFit = function(...) {
    .calls <<- .calls + 1
    if (.calls > 1 && .calls %% 2 == 0) stop("boom")
    .real(...)
  })
  expect_error(cuiDoseResponse(.d, cuiBinary(eff), dose = dose, nBoot = 10, seed = 1),
               "only 5 of 10 bootstrap draws succeeded.*boom")
  .calls <- 0
  expect_s3_class(cuiDoseResponse(.d, cuiBinary(eff), dose = dose, nBoot = 10, seed = 1,
                                  minSuccess = 0.5), "nlmixr2Cui")
})

test_that("unsorted, non-integer doses keep rates aligned", {
  .d <- data.frame(dose = c(10, 0.5, 2.5, 10, 0.5, 2.5, 0.5),
                   eff = c(1, 0, 1, 1, 1, 0, 0))
  .r <- cuiDoseResponse(.d, cuiBinary(eff), dose = dose, nBoot = 5, seed = 1)
  expect_equal(.r$summary$regimen, c("0.5", "2.5", "10"))
  expect_equal(.r$summary$est, c(1 / 3, 1 / 2, 1))
  .obs <- as.data.frame(.r, type = "observed")
  expect_equal(.obs$n, c(3, 2, 2))
  .near <- data.frame(dose = c(0.3, 0.1 + 0.2), eff = c(1, 0))
  expect_error(cuiDoseResponse(.near, cuiBinary(eff), dose = dose, nBoot = 0), "15 digits")
})

test_that("the seed helper removes a new .Random.seed and restores after errors", {
  .env <- globalenv()
  .had <- exists(".Random.seed", envir = .env, inherits = FALSE)
  .old <- if (.had) get(".Random.seed", envir = .env)
  withr::defer(if (.had) assign(".Random.seed", .old, envir = .env))
  if (.had) rm(".Random.seed", envir = .env)
  .cur$.cuiWithSeed(1, stats::runif(1))
  expect_false(exists(".Random.seed", envir = .env, inherits = FALSE))
  set.seed(5)
  .want <- stats::runif(1)
  set.seed(5)
  expect_error(.cur$.cuiWithSeed(1, stop("inside")), "inside")
  expect_equal(stats::runif(1), .want)
})
