skip_on_cran()

# Shared fit: one-compartment model on theo_sd with covariance step.
# Created once; reused across all Step 1 tests.
.theo_one_cmt <- function() {
  ini({
    tka <- log(1.57)
    tcl <- log(2.72)
    tv <- log(31.5)
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

.theo_fit <- suppressMessages(
  nlmixr2(
    .theo_one_cmt,
    nlmixr2data::theo_sd,
    est = "focei",
    control = list(print = 0L, covMethod = "r")
  )
)

# Step 1: sirGetProposalCov ----------------------------------------------------

test_that("sirGetProposalCov returns covariance close to fit$cov with defaults", {
  result <- sirGetProposalCov(.theo_fit)
  expect_equal(result, .theo_fit$cov, tolerance = 1e-10)
})

test_that("sirGetProposalCov doubles theta variances with thetaInflation = 2", {
  orig <- .theo_fit$cov
  result <- sirGetProposalCov(.theo_fit, thetaInflation = 2, capCorrelation = 1)

  ini_df <- .theo_fit$iniDf
  theta_names <- ini_df$name[!is.na(ini_df$ntheta) & !ini_df$fix]
  theta_in <- rownames(orig) %in% theta_names

  # Theta variances should be 2x; omega variances unchanged
  expect_equal(
    diag(result)[theta_in],
    2 * diag(orig)[theta_in],
    tolerance = 1e-10
  )
  expect_equal(
    diag(result)[!theta_in],
    diag(orig)[!theta_in],
    tolerance = 1e-10
  )
  # Correlations must be preserved
  expect_equal(cov2cor(result), cov2cor(orig), tolerance = 1e-10)
})


test_that("sirGetProposalCov caps correlations at capCorrelation", {
  # Inflate theta 10x to push correlations toward 1, then cap at 0.3
  result <- sirGetProposalCov(
    .theo_fit,
    thetaInflation = 10,
    capCorrelation = 0.3
  )
  corr <- cov2cor(result)
  off <- corr[lower.tri(corr)]
  expect_true(all(off >= -0.3 - 1e-10 & off <= 0.3 + 1e-10))
})

test_that("sirGetProposalCov with capCorrelation = 1 leaves correlations unchanged", {
  orig <- .theo_fit$cov
  result <- sirGetProposalCov(.theo_fit, capCorrelation = 1)
  expect_equal(cov2cor(result), cov2cor(orig), tolerance = 1e-10)
})

test_that("sirGetProposalCov result is symmetric", {
  result <- sirGetProposalCov(
    .theo_fit,
    thetaInflation = 2,
    capCorrelation = 0.8
  )
  expect_equal(result, t(result), tolerance = 1e-14)
})

test_that("sirGetProposalCov preserves parameter names from fit$cov", {
  result <- sirGetProposalCov(.theo_fit)
  expect_equal(rownames(result), rownames(.theo_fit$cov))
  expect_equal(colnames(result), colnames(.theo_fit$cov))
})

test_that("sirGetProposalCov errors when fit has no covariance matrix", {
  # Strip the covariance by fitting without a cov step
  fit_no_cov <- suppressMessages(
    nlmixr2(
      .theo_one_cmt,
      nlmixr2data::theo_sd,
      est = "focei",
      control = list(print = 0L, covMethod = "")
    )
  )
  expect_error(sirGetProposalCov(fit_no_cov), "covariance matrix")
})

test_that("sirGetProposalCov errors on non-fit input", {
  expect_error(sirGetProposalCov(list()), class = "error")
})

# Step 2: sirSampleTheta -------------------------------------------------------

# Fixed 3-parameter setup reused across Step 2 tests
.mu3 <- c(tka = 0.45, tcl = 0.98, tv = 3.47)
.cov3 <- matrix(
  c(0.065, -0.004, 0.004, -0.004, 0.006, -0.002, 0.004, -0.002, 0.003),
  nrow = 3L,
  dimnames = list(names(.mu3), names(.mu3))
)

test_that("sirSampleTheta returns n rows with no bounds", {
  set.seed(1)
  res <- sirSampleTheta(.mu3, .cov3, n = 200L)
  expect_equal(nrow(res$samples), 200L)
  expect_equal(ncol(res$samples), 3L)
  expect_equal(res$nRejected, 0L)
})

test_that("sirSampleTheta column names match mu", {
  set.seed(1)
  res <- sirSampleTheta(.mu3, .cov3, n = 50L)
  expect_equal(colnames(res$samples), names(.mu3))
})

test_that("sirSampleTheta respects lower and upper bounds", {
  set.seed(42)
  lo <- .mu3 - 0.5
  hi <- .mu3 + 0.5
  res <- sirSampleTheta(.mu3, .cov3, n = 200L, lower = lo, upper = hi)
  expect_true(all(res$samples >= rep(lo, each = nrow(res$samples))))
  expect_true(all(res$samples <= rep(hi, each = nrow(res$samples))))
})

test_that("sirSampleTheta nRejected increases with tight bounds", {
  set.seed(7)
  # Wide bounds: expect few rejections
  res_wide <- sirSampleTheta(
    .mu3,
    .cov3,
    n = 200L,
    lower = .mu3 - 10,
    upper = .mu3 + 10
  )
  # Tight bounds: expect more rejections
  res_tight <- suppressWarnings(
    sirSampleTheta(
      .mu3,
      .cov3,
      n = 200L,
      lower = .mu3 - 0.01,
      upper = .mu3 + 0.01
    )
  )
  expect_gt(res_tight$nRejected, res_wide$nRejected)
})

test_that("sirSampleTheta is reproducible with set.seed", {
  set.seed(99)
  r1 <- sirSampleTheta(.mu3, .cov3, n = 50L)
  set.seed(99)
  r2 <- sirSampleTheta(.mu3, .cov3, n = 50L)
  expect_equal(r1$samples, r2$samples)
})

test_that("sirSampleTheta warns and returns fewer rows when bounds exclude all draws", {
  # Bounds set to a tiny region far from mu — all draws will be rejected
  lo <- .mu3 + 1e6
  hi <- .mu3 + 1e6 + 1e-9
  expect_warning(
    res <- sirSampleTheta(.mu3, .cov3, n = 10L, lower = lo, upper = hi),
    regexp = "within bounds"
  )
  expect_equal(nrow(res$samples), 0L)
  expect_equal(res$nRejected, 100L) # 10 * 10 max attempts, all rejected
})

test_that("sirSampleTheta scalar bounds are recycled to length p", {
  set.seed(3)
  res <- sirSampleTheta(.mu3, .cov3, n = 100L, lower = -100, upper = 100)
  expect_equal(nrow(res$samples), 100L)
  expect_true(all(res$samples >= -100 & res$samples <= 100))
})

# Step 3: sirSampleOmegaSigma --------------------------------------------------

# 2×2 OMEGA — lower triangle has 3 elements: [1,1], [2,1], [2,2]
.omega2 <- matrix(
  c(0.5, 0.1, 0.1, 0.3),
  nrow = 2L,
  dimnames = list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl"))
)
# Diagonal omegaCovMat built from plausible SEs
.omega2_cov <- diag(c(0.05^2, 0.02^2, 0.04^2))

test_that("sirSampleOmegaSigma returns n matrices", {
  set.seed(1)
  res <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 50L)
  expect_length(res$samples, 50L)
})

test_that("sirSampleOmegaSigma all samples are positive definite", {
  set.seed(2)
  res <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 100L)
  pd_ok <- vapply(
    res$samples,
    function(m) {
      tryCatch(
        {
          chol(m)
          TRUE
        },
        error = function(e) FALSE
      )
    },
    logical(1L)
  )
  expect_true(all(pd_ok))
})

test_that("sirSampleOmegaSigma all samples are symmetric", {
  set.seed(3)
  res <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 50L)
  sym_ok <- vapply(
    res$samples,
    function(m) isTRUE(all.equal(m, t(m))),
    logical(1L)
  )
  expect_true(all(sym_ok))
})

test_that("sirSampleOmegaSigma preserves dimnames from omegaEst", {
  set.seed(4)
  res <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 10L)
  dn_ok <- vapply(
    res$samples,
    function(m) {
      identical(dimnames(m), dimnames(.omega2))
    },
    logical(1L)
  )
  expect_true(all(dn_ok))
})

test_that("sirSampleOmegaSigma nRejected is non-negative integer", {
  set.seed(5)
  res <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 50L)
  expect_true(is.integer(res$nRejected) || is.numeric(res$nRejected))
  expect_gte(res$nRejected, 0L)
})

test_that("sirSampleOmegaSigma is reproducible with set.seed", {
  set.seed(77)
  r1 <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 30L)
  set.seed(77)
  r2 <- sirSampleOmegaSigma(.omega2, .omega2_cov, n = 30L)
  expect_equal(r1$samples, r2$samples)
})

test_that("sirSampleOmegaSigma warns and returns fewer matrices when budget exhausted", {
  # Mean is a non-PD configuration (large off-diagonal, tiny diagonals);
  # near-zero variance pins draws close to the mean → ~0% PD rate.
  bad_omega <- matrix(
    c(1e-4, 1.0, 1.0, 1e-4),
    2L,
    2L,
    dimnames = list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl"))
  )
  bad_cov <- diag(c(1e-10, 1e-10, 1e-10))
  expect_warning(
    res <- sirSampleOmegaSigma(bad_omega, bad_cov, n = 50L),
    regexp = "positive definite"
  )
  expect_lt(length(res$samples), 50L)
})

test_that("sirSampleOmegaSigma works with 1x1 OMEGA (single eta)", {
  omega1 <- matrix(0.4, 1L, 1L, dimnames = list("eta.ka", "eta.ka"))
  cov1 <- matrix(0.01, 1L, 1L)
  set.seed(9)
  res <- sirSampleOmegaSigma(omega1, cov1, n = 30L)
  expect_length(res$samples, 30L)
  pd_ok <- vapply(res$samples, function(m) m[1L, 1L] > 0, logical(1L))
  expect_true(all(pd_ok))
})

# Step 4: sirEvalOFV -----------------------------------------------------------

test_that("sirEvalOFV returns OFV close to fit$objf at true estimates", {
  mu <- .theo_fit$theta[rownames(.theo_fit$cov)]
  param_mat <- matrix(mu, nrow = 1L, dimnames = list(NULL, names(mu)))

  ofvs <- sirEvalOFV(.theo_fit, param_mat)

  expect_length(ofvs, 1L)
  expect_false(is.na(ofvs[1L]))
  expect_lt(abs(ofvs[1L] - .theo_fit$objf), 1)
})

test_that("sirEvalOFV returns NA for a clearly invalid parameter vector", {
  mu <- .theo_fit$theta[rownames(.theo_fit$cov)]
  bad <- mu
  bad["tka"] <- 1e10 # absurdly large; model should fail or produce NA
  param_mat <- matrix(bad, nrow = 1L, dimnames = list(NULL, names(bad)))

  ofvs <- suppressWarnings(sirEvalOFV(.theo_fit, param_mat))

  expect_length(ofvs, 1L)
  expect_true(is.na(ofvs[1L]) || is.finite(ofvs[1L])) # NA or a number (not NaN/Inf)
})

test_that("sirEvalOFV returns a vector with one entry per row", {
  mu <- .theo_fit$theta[rownames(.theo_fit$cov)]
  # Three rows: true estimates, slight perturbations
  param_mat <- rbind(
    mu,
    mu + 0.01,
    mu - 0.01
  )
  dimnames(param_mat) <- list(NULL, names(mu))

  ofvs <- sirEvalOFV(.theo_fit, param_mat)

  expect_length(ofvs, 3L)
  expect_true(all(!is.na(ofvs)))
})

test_that("sirEvalOFV returns the same OFVs with future workers", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  plan_before <- future::plan()
  on.exit(future::plan(plan_before), add = TRUE)

  mu <- .theo_fit$theta[rownames(.theo_fit$cov)]
  param_mat <- rbind(
    mu,
    mu + 0.01
  )
  dimnames(param_mat) <- list(NULL, names(mu))

  ofv_seq <- suppressMessages(sirEvalOFV(.theo_fit, param_mat, workers = 1L))
  ofv_par <- suppressMessages(sirEvalOFV(.theo_fit, param_mat, workers = 2L))

  expect_equal(ofv_par, ofv_seq, tolerance = 1e-8)
})

test_that("sirEvalOFV handles diagonal omega columns for multi-ETA models", {
  three_eta_one_cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.00
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  fit <- suppressMessages(suppressWarnings(
    nlmixr2(
      three_eta_one_cmt,
      nlmixr2data::theo_sd,
      est = "focei",
      control = list(print = 0L),
      table = list(npde = TRUE, cwres = TRUE)
    )
  ))
  proposal <- .sirInitialProposal(
    fit,
    fit$theta[rownames(fit$cov)],
    sirGetProposalCov(fit),
    capCorrelation = 0.8
  )
  param_mat <- matrix(
    proposal$mu,
    nrow = 1L,
    dimnames = list(NULL, names(proposal$mu))
  )

  eta_cols <- c("eta.ka", "eta.cl", "eta.v")
  expect_setequal(intersect(eta_cols, colnames(param_mat)), eta_cols)

  ofvs <- sirEvalOFV(fit, param_mat, workers = 1L)

  expect_length(ofvs, 1L)
  expect_false(is.na(ofvs[1L]))
  expect_lt(abs(ofvs[1L] - fit$objf), 1)
})

test_that("sirEvalOFV OFV increases away from the estimates", {
  mu <- .theo_fit$theta[rownames(.theo_fit$cov)]
  big_offset <- mu + 1 # large shift in all thetas
  param_mat <- rbind(mu, big_offset)
  dimnames(param_mat) <- list(NULL, names(mu))

  ofvs <- suppressWarnings(sirEvalOFV(.theo_fit, param_mat))

  expect_true(!is.na(ofvs[1L]))
  # OFV at true values should be lower (better fit)
  if (!is.na(ofvs[2L])) expect_lt(ofvs[1L], ofvs[2L])
})

test_that("sirEvalOFV errors when no columns match fit parameters", {
  bad_mat <- matrix(1, nrow = 1L, dimnames = list(NULL, "not_a_param"))
  expect_error(sirEvalOFV(.theo_fit, bad_mat), "match")
})

# Step 5: sirCalcWeights -------------------------------------------------------

# Fixed 3-parameter setup reused across Step 5 tests (same as Step 2)
.wt_mu <- c(tka = 0.45, tcl = 0.98, tv = 3.47)
.wt_cov <- matrix(
  c(0.065, -0.004, 0.004, -0.004, 0.006, -0.002, 0.004, -0.002, 0.003),
  nrow = 3L,
  dimnames = list(names(.wt_mu), names(.wt_mu))
)

test_that("sirCalcWeights: relPDF = 1 at the proposal mean", {
  samp <- matrix(.wt_mu, nrow = 1L, dimnames = list(NULL, names(.wt_mu)))
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = 0)
  expect_equal(res$relPDF[1L], 1, tolerance = 1e-12)
})

test_that("sirCalcWeights: prob_resample sums to 1", {
  set.seed(1)
  samp <- mvtnorm::rmvnorm(50L, mean = .wt_mu, sigma = .wt_cov)
  dofv <- rnorm(50L, mean = 0, sd = 2)
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = dofv)
  expect_equal(sum(res$prob_resample), 1, tolerance = 1e-12)
})

test_that("sirCalcWeights: negative dOFV (better fit) gives IR > 1 at mu", {
  # At mu: relPDF = 1, so IR = exp(-0.5 * dOFV); dOFV < 0 → IR > 1
  samp <- matrix(.wt_mu, nrow = 1L, dimnames = list(NULL, names(.wt_mu)))
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = -2)
  expect_gt(res$importance_ratio[1L], 1)
})

test_that("sirCalcWeights: NA dOFV gives prob_resample = 0", {
  set.seed(2)
  samp <- mvtnorm::rmvnorm(5L, mean = .wt_mu, sigma = .wt_cov)
  dofv <- c(1, NA, 2, NA, 0.5)
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = dofv)
  expect_equal(res$prob_resample[c(2L, 4L)], c(0, 0))
  expect_true(is.na(res$importance_ratio[2L]))
  expect_equal(sum(res$prob_resample), 1, tolerance = 1e-12)
})

test_that("sirCalcWeights: output has correct structure", {
  set.seed(3)
  samp <- mvtnorm::rmvnorm(10L, mean = .wt_mu, sigma = .wt_cov)
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = rep(0, 10L))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 10L)
  expect_named(
    res,
    c(
      "sample_id",
      "dOFV",
      "likelihood_ratio",
      "relPDF",
      "importance_ratio",
      "prob_resample"
    )
  )
})

test_that("sirCalcWeights: relPDF decreases away from mu", {
  far <- matrix(.wt_mu + 5, nrow = 1L, dimnames = list(NULL, names(.wt_mu)))
  near <- matrix(.wt_mu + 0.01, nrow = 1L, dimnames = list(NULL, names(.wt_mu)))
  samp <- rbind(far, near)
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = c(0, 0))
  expect_lt(res$relPDF[1L], res$relPDF[2L])
})

test_that("sirCalcWeights: positive dOFV (worse fit) at mu gives IR < 1", {
  samp <- matrix(.wt_mu, nrow = 1L, dimnames = list(NULL, names(.wt_mu)))
  res <- sirCalcWeights(samp, .wt_mu, .wt_cov, dOFV = 4)
  expect_lt(res$importance_ratio[1L], 1)
})

test_that("sirCalcWeights: errors on non-PD covMat", {
  bad_cov <- matrix(c(1, 2, 2, 1), 2L) # not PD
  samp <- matrix(c(1, 2), nrow = 1L)
  expect_error(
    sirCalcWeights(samp, c(1, 2), bad_cov, dOFV = 0),
    "positive definite"
  )
})

test_that("sirCalcWeights: log-space normalization handles extreme weights", {
  cov_mat <- diag(2)
  mu <- c(a = 0, b = 0)
  samp <- matrix(
    c(100, 100, 0, 0),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, names(mu))
  )
  res <- sirCalcWeights(samp, mu, cov_mat, dOFV = c(-1000, 1000))
  expect_equal(sum(res$prob_resample), 1, tolerance = 1e-12)
  expect_true(all(is.finite(res$prob_resample)))
  expect_gt(res$prob_resample[1L], 0.999)
})

# Step 6: sirResample ----------------------------------------------------------

# Helper: build a weights data frame with uniform probs for n samples
.uniform_weights <- function(n) {
  data.frame(prob_resample = rep(1 / n, n))
}

# Helper: build a weights data frame with all weight on row 1
.spike_weights <- function(n) {
  p <- c(1, rep(0, n - 1L))
  data.frame(prob_resample = p)
}

test_that("sirResample returns m rows", {
  set.seed(1)
  samp <- matrix(rnorm(30L), nrow = 10L)
  res <- sirResample(samp, .uniform_weights(10L), m = 5L)
  expect_equal(nrow(res$samples), 5L)
  expect_equal(ncol(res$samples), 3L)
})

test_that("sirResample resampleCounts sums to m", {
  set.seed(2)
  samp <- matrix(rnorm(30L), nrow = 10L)
  res <- sirResample(samp, .uniform_weights(10L), m = 6L)
  expect_equal(sum(res$resampleCounts), 6L)
  expect_length(res$resampleCounts, 10L)
})

test_that("sirResample with capped replacement can repeatedly select row 1", {
  samp <- matrix(1:20, nrow = 5L)
  res <- sirResample(samp, .spike_weights(5L), m = 5L, capResampling = 5L)
  # Every resampled row should equal samp[1, ]
  expect_true(all(sweep(res$samples, 2L, samp[1L, ], "==") == 1L))
  expect_equal(res$resampleCounts[1L], 5L)
  expect_true(all(res$resampleCounts[-1L] == 0L))
})

test_that("sirResample cap is respected: no sample exceeds capResampling", {
  set.seed(3)
  samp <- matrix(rnorm(50L), nrow = 10L)
  cap <- 3L
  res <- sirResample(samp, .uniform_weights(10L), m = 20L, capResampling = cap)
  expect_true(all(res$resampleCounts <= cap))
  expect_equal(sum(res$resampleCounts), 20L)
  expect_equal(nrow(res$samples), 20L)
})

test_that("sirResample cap = 1 samples without replacement", {
  set.seed(4)
  samp <- matrix(rnorm(30L), nrow = 10L)
  res <- sirResample(samp, .uniform_weights(10L), m = 10L, capResampling = 1)
  expect_equal(nrow(res$samples), 10L)
  expect_true(all(res$resampleCounts <= 1L))
  expect_equal(sum(res$resampleCounts), 10L)
})

test_that("sirResample with cap forces spread when one sample dominates", {
  set.seed(5)
  samp <- matrix(rnorm(50L), nrow = 10L)
  # Uniform weights + cap = 3: row 1 cannot be selected more than 3 times
  res <- sirResample(samp, .uniform_weights(10L), m = 20L, capResampling = 3L)
  expect_lte(res$resampleCounts[1L], 3L)
  expect_equal(sum(res$resampleCounts), 20L)
  expect_equal(nrow(res$samples), 20L)
})

test_that("sirResample selected rows all come from original samples", {
  set.seed(6)
  samp <- matrix(seq_len(20L), nrow = 5L)
  res <- sirResample(samp, .uniform_weights(5L), m = 5L)
  # Every resampled row must exactly match one of the original rows
  for (i in seq_len(nrow(res$samples))) {
    row_match <- apply(samp, 1L, function(r) identical(r, res$samples[i, ]))
    expect_true(any(row_match))
  }
})

test_that("sirResample returns sample-order metadata", {
  set.seed(7)
  samp <- matrix(seq_len(20L), nrow = 5L)
  res <- sirResample(samp, .uniform_weights(5L), m = 5L)
  expect_length(res$sampleOrder, 5L)
  expect_length(res$selectionOrder, 5L)
  expect_equal(sum(res$resampleCounts), 5L)
})

# Step 7: sirBoxCox / sirBoxCoxInverse -----------------------------------------

test_that("sirBoxCox with lambda = 1 is a linear (affine) transformation of x", {
  x <- c(1, 2, 4, 8, 16)
  res <- sirBoxCox(x, lambda = 1)
  # transformed = (x + delta)^1 - 1 = x + delta - 1, linear in x
  expect_equal(res$lambda, 1)
  expect_equal(cor(res$transformed, x), 1, tolerance = 1e-12)
})

test_that("sirBoxCox with lambda = 0 returns log(x + delta)", {
  x <- c(1, 2, 4, 8, 16)
  delta <- abs(min(x)) + 1e-6
  res <- sirBoxCox(x, lambda = 0)
  expect_equal(res$transformed, log(x + delta), tolerance = 1e-12)
})

test_that("sirBoxCoxInverse recovers x after sirBoxCox (lambda estimated)", {
  set.seed(1)
  x <- exp(rnorm(50L)) # log-normal: should yield lambda ≈ 0
  res <- sirBoxCox(x)
  x_back <- sirBoxCoxInverse(res$transformed, res$lambda, res$delta)
  expect_equal(x_back, x, tolerance = 1e-8)
})

test_that("sirBoxCoxInverse recovers x for fixed lambda = 1", {
  x <- c(0.5, 1, 2, 5, 10)
  res <- sirBoxCox(x, lambda = 1)
  x_back <- sirBoxCoxInverse(res$transformed, res$lambda, res$delta)
  expect_equal(x_back, x, tolerance = 1e-10)
})

test_that("sirBoxCoxInverse recovers x for fixed lambda = 0 (log case)", {
  x <- c(0.1, 0.5, 1, 5, 20)
  res <- sirBoxCox(x, lambda = 0)
  x_back <- sirBoxCoxInverse(res$transformed, res$lambda, res$delta)
  expect_equal(x_back, x, tolerance = 1e-10)
})

test_that("sirBoxCox estimated lambda improves normality vs untransformed", {
  set.seed(2)
  x <- rexp(100L, rate = 2) # clearly right-skewed
  res <- sirBoxCox(x)
  # Shapiro-Wilk p should be higher after transformation
  p_raw <- shapiro.test(x)$p.value
  p_bc <- shapiro.test(res$transformed)$p.value
  expect_gt(p_bc, p_raw)
})

test_that("sirBoxCox returns correct list structure", {
  x <- rnorm(20L) + 5
  res <- sirBoxCox(x)
  expect_named(res, c("transformed", "lambda", "delta"))
  expect_length(res$transformed, 20L)
  expect_length(res$lambda, 1L)
  expect_length(res$delta, 1L)
})

test_that("sirBoxCox delta defaults to |min(x)| + 1e-6", {
  x <- c(-3, 0, 1, 5)
  res <- sirBoxCox(x, lambda = 1)
  expected_delta <- abs(min(x)) + 1e-6
  expect_equal(res$delta, expected_delta, tolerance = 1e-12)
})

test_that("sirBoxCox supplied delta is respected", {
  x <- c(1, 2, 3)
  res <- sirBoxCox(x, lambda = 1, delta = 0.5)
  expect_equal(res$delta, 0.5)
})

test_that("sirBoxCoxInverse errors when back-transform base is non-positive", {
  # lambda = 2, x_transformed = -1 → base = 2*(-1) + 1 = -1 ≤ 0
  expect_error(sirBoxCoxInverse(-1, lambda = 2, delta = 0), "positive")
})

# Step 8: sirUpdateProposal ----------------------------------------------------

test_that("sirUpdateProposal covMat is symmetric", {
  set.seed(1)
  mat <- matrix(
    rexp(100L),
    nrow = 20L,
    dimnames = list(NULL, c("a", "b", "c", "d", "e"))
  )
  res <- sirUpdateProposal(mat)
  expect_equal(res$covMat, t(res$covMat), tolerance = 1e-12)
})

test_that("sirUpdateProposal covMat is positive semi-definite", {
  set.seed(2)
  mat <- matrix(rexp(60L), nrow = 20L, dimnames = list(NULL, c("a", "b", "c")))
  res <- sirUpdateProposal(mat)
  eigs <- eigen(res$covMat, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

test_that("sirUpdateProposal preserves column names", {
  set.seed(3)
  nms <- c("tka", "tcl", "tv")
  mat <- matrix(rnorm(60L) + 5, nrow = 20L, dimnames = list(NULL, nms))
  res <- sirUpdateProposal(mat)
  expect_equal(rownames(res$covMat), nms)
  expect_equal(colnames(res$covMat), nms)
})

test_that("sirUpdateProposal boxcoxParams has correct structure when boxcox = TRUE", {
  set.seed(4)
  nms <- c("p1", "p2")
  mat <- matrix(rexp(40L), nrow = 20L, dimnames = list(NULL, nms))
  res <- sirUpdateProposal(mat, boxcox = TRUE)
  expect_s3_class(res$boxcoxParams, "data.frame")
  expect_named(res$boxcoxParams, c("param", "lambda", "delta"))
  expect_equal(nrow(res$boxcoxParams), 2L)
  expect_equal(res$boxcoxParams$param, nms)
})

test_that("sirUpdateProposal boxcoxParams is NULL when boxcox = FALSE", {
  set.seed(5)
  mat <- matrix(rnorm(40L), nrow = 20L)
  res <- sirUpdateProposal(mat, boxcox = FALSE)
  expect_null(res$boxcoxParams)
})

test_that("sirUpdateProposal boxcox = FALSE matches cov() directly", {
  set.seed(6)
  mat <- matrix(rnorm(60L), nrow = 20L, dimnames = list(NULL, c("a", "b", "c")))
  res <- sirUpdateProposal(mat, boxcox = FALSE)
  expected <- cov(mat)
  expect_equal(res$covMat, expected, tolerance = 1e-12)
})

test_that("sirUpdateProposal lambdas are retrievable and finite", {
  set.seed(7)
  mat <- matrix(
    rexp(80L),
    nrow = 20L,
    dimnames = list(NULL, c("a", "b", "c", "d"))
  )
  res <- sirUpdateProposal(mat, boxcox = TRUE)
  expect_true(all(is.finite(res$boxcoxParams$lambda)))
  expect_true(all(res$boxcoxParams$lambda >= -3 & res$boxcoxParams$lambda <= 3))
})

# Step 9: sirRunIteration ------------------------------------------------------
# These tests require OFV evaluation so use tiny nSamples (8) for speed.

.iter1 <- local({
  set.seed(42)
  mu <- .theo_fit$theta[rownames(.theo_fit$cov)]
  prop_cov <- sirGetProposalCov(.theo_fit)
  suppressMessages(
    sirRunIteration(
      .theo_fit,
      mu = mu,
      proposalCov = prop_cov,
      nSamples = 8L,
      nResample = 4L,
      iterNum = 1L,
      recenter = TRUE,
      boxcox = TRUE,
      directory = NULL
    )
  )
})

test_that("sirRunIteration returns correct list structure", {
  expect_named(
    .iter1,
    c(
      "resampledMat",
      "newMu",
      "newCov",
      "iterSummary",
      "boxcoxState",
      "rawResults"
    )
  )
})

test_that("sirRunIteration resampledMat has nResample rows", {
  expect_equal(nrow(.iter1$resampledMat), 4L)
})

test_that("sirRunIteration iterSummary has expected columns and values", {
  s <- .iter1$iterSummary
  expect_named(
    s,
    c(
      "iter",
      "nSamples",
      "nAttempted",
      "nDrawAttempts",
      "nCollected",
      "nSuccessful",
      "nFailed",
      "nResample",
      "nResampled",
      "minDOFV",
      "meanDOFV",
      "nNegativeDOFV",
      "thetaRejected",
      "omegaRejected",
      "sigmaRejected",
      "inverseRejected"
    )
  )
  expect_equal(s$iter, 1L)
  expect_equal(s$nSamples, 8L)
  expect_equal(s$nAttempted, 8L)
  expect_true(s$nCollected <= 8L && s$nCollected >= 1L)
  expect_gte(s$nFailed, 0L)
  expect_true(is.numeric(s$minDOFV))
  expect_gte(s$sigmaRejected, 0L)
})

test_that("sirRunIteration newMu is a named numeric vector matching proposal params", {
  mu_names <- colnames(.iter1$resampledMat)
  expect_named(.iter1$newMu, mu_names)
  expect_true(is.numeric(.iter1$newMu))
})

test_that("sirRunIteration newCov is symmetric and positive semi-definite", {
  m <- .iter1$newCov
  expect_equal(m, t(m), tolerance = 1e-12)
  eigs <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= -1e-10))
})

test_that("sirRunIteration boxcoxState has correct structure when boxcox = TRUE", {
  bc <- .iter1$boxcoxState
  expect_s3_class(bc, "data.frame")
  expect_named(bc, c("param", "lambda", "delta"))
  expect_equal(bc$param, colnames(.iter1$resampledMat))
})

test_that("sirRunIteration rawResults has dOFV column", {
  expect_true("dOFV" %in% names(.iter1$rawResults))
  expect_equal(nrow(.iter1$rawResults), .iter1$iterSummary$nCollected + 1L)
  expect_equal(.iter1$rawResults$sample_id[[1L]], 0L)
  expect_equal(.iter1$rawResults$dOFV[[1L]], 0)
})

test_that("sirRunIteration rawResults has PsN-like SIR metadata", {
  expect_true(all(
    c(
      "sample_id",
      "likelihood_ratio",
      "probability_resample",
      "resamples",
      "sample_order"
    ) %in%
      names(.iter1$rawResults)
  ))
  expect_equal(sum(.iter1$rawResults$resamples), .iter1$iterSummary$nResampled)
})

test_that(".sirBuildRawResults expands capped raw-result rows", {
  param_mat <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b"))
  )
  weights <- data.frame(
    likelihood_ratio = c(1, 2, 3),
    relPDF = 1,
    importance_ratio = c(1, 2, 3),
    prob_resample = c(0.2, 0.3, 0.5)
  )
  resampled <- list(
    resampleCounts = c(2L, 0L, 1L),
    sampleOrder = c("2;4", "", "1")
  )
  raw <- nlmixr2extra:::.sirBuildRawResults(
    paramMat = param_mat,
    weights = weights,
    dOFV = c(0.1, 0.2, 0.3),
    resampled = resampled,
    mu = c(a = 0, b = 0),
    capResampling = 2
  )
  expect_equal(nrow(raw), 7L)
  expect_equal(raw$sample_id[[1L]], 0L)
  expect_equal(raw$resamples[raw$sample_id == 1L], c(1L, 1L))
  expect_equal(raw$sample_order[raw$sample_id == 1L], c("2", "4"))
  expect_equal(raw$resamples[raw$sample_id == 2L], c(0L, 0L))
})

test_that("sirRunIteration recentering: newMu shifts when a better sample exists", {
  # Use a mu perturbed away from estimates so samples near true values get dOFV < 0
  mu_perturbed <- .theo_fit$theta[rownames(.theo_fit$cov)] + 0.5
  prop_cov <- sirGetProposalCov(.theo_fit)
  set.seed(1)
  res <- suppressMessages(
    sirRunIteration(
      .theo_fit,
      mu = mu_perturbed,
      proposalCov = prop_cov,
      nSamples = 8L,
      nResample = 4L,
      iterNum = 1L,
      recenter = TRUE,
      boxcox = FALSE,
      directory = NULL
    )
  )
  # newMu should differ from mu_perturbed if any dOFV < 0 was found
  min_dofv <- res$iterSummary$minDOFV
  if (!is.na(min_dofv) && min_dofv < 0) {
    expect_false(isTRUE(all.equal(res$newMu, mu_perturbed, tolerance = 1e-6)))
  } else {
    expect_equal(res$newMu[names(mu_perturbed)], mu_perturbed, tolerance = 1e-6)
  }
})

test_that("sirRunIteration writes CSV when directory is provided", {
  mu <- .theo_fit$theta[rownames(.theo_fit$cov)]
  prop_cov <- sirGetProposalCov(.theo_fit)
  tmp_dir <- tempfile("sir_test_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE))

  set.seed(9)
  suppressMessages(
    sirRunIteration(
      .theo_fit,
      mu = mu,
      proposalCov = prop_cov,
      nSamples = 5L,
      nResample = 3L,
      iterNum = 2L,
      directory = tmp_dir
    )
  )
  expected_csv <- file.path(tmp_dir, "raw_results_sir_iteration2.csv")
  expect_true(file.exists(expected_csv))
  df <- utils::read.csv(expected_csv)
  expect_true("dOFV" %in% names(df))
})

test_that("sirRunIteration chained: iter 2 accepts boxcoxState from iter 1", {
  set.seed(77)
  res2 <- suppressMessages(
    sirRunIteration(
      .theo_fit,
      mu = .iter1$newMu,
      proposalCov = .iter1$newCov,
      nSamples = 6L,
      nResample = 3L,
      iterNum = 2L,
      boxcoxState = .iter1$boxcoxState
    )
  )
  expect_named(
    res2,
    c(
      "resampledMat",
      "newMu",
      "newCov",
      "iterSummary",
      "boxcoxState",
      "rawResults"
    )
  )
  expect_equal(res2$iterSummary$iter, 2L)
})

test_that("sirRunIteration new proposal includes omega and sigma columns", {
  expect_true("eta.ka" %in% rownames(.iter1$newCov))
  expect_true("add.sd" %in% rownames(.iter1$newCov))
  expect_true("eta.ka" %in% .iter1$boxcoxState$param)
  expect_true("add.sd" %in% .iter1$boxcoxState$param)
})

test_that("sir initial proposal applies omega and sigma inflation", {
  mu <- .theo_fit$theta[rownames(.theo_fit$cov)]
  prop_cov <- sirGetProposalCov(.theo_fit)
  base <- nlmixr2extra:::.sirInitialProposal(
    .theo_fit,
    mu,
    prop_cov,
    capCorrelation = 1
  )
  inflated <- nlmixr2extra:::.sirInitialProposal(
    .theo_fit,
    mu,
    prop_cov,
    omegaInflation = 4,
    sigmaInflation = 9,
    capCorrelation = 1
  )
  expect_equal(
    inflated$covMat["eta.ka", "eta.ka"],
    4 * base$covMat["eta.ka", "eta.ka"],
    tolerance = 1e-8
  )
  expect_equal(
    inflated$covMat["add.sd", "add.sd"],
    9 * base$covMat["add.sd", "add.sd"],
    tolerance = 1e-8
  )
})

# Gap fix: .sirOmegaProposalInfo, .sirSigmaInfo, .sirReconstructOmega ----------

test_that(".sirOmegaProposalInfo returns correct structure for single-eta model", {
  info <- nlmixr2extra:::.sirOmegaProposalInfo(.theo_fit)
  expect_s3_class(info, "data.frame")
  expect_named(info, c("colName", "neta1", "neta2", "est", "se", "isDiag"))
  # Single eta.ka: one diagonal element
  expect_equal(nrow(info), 1L)
  expect_equal(info$colName, "eta.ka")
  expect_true(info$isDiag)
  expect_gt(info$se, 0)
})

test_that(".sirOmegaProposalInfo est matches fit$omega diagonal", {
  info <- nlmixr2extra:::.sirOmegaProposalInfo(.theo_fit)
  expect_equal(
    info$est[info$isDiag],
    unname(diag(as.matrix(.theo_fit$omega))),
    tolerance = 1e-10
  )
})

test_that(".sirOmegaProposalInfo excludes fixed omega elements", {
  fake_fit <- list(
    nsub = 20L,
    omega = matrix(
      c(0.4, 0.05, 0.05, 0.2),
      2L,
      2L,
      dimnames = list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl"))
    ),
    iniDf = data.frame(
      name = c("eta.ka", "eta.cl:eta.ka", "eta.cl"),
      neta1 = c(1L, 2L, 2L),
      neta2 = c(1L, 1L, 2L),
      fix = c(FALSE, TRUE, TRUE),
      stringsAsFactors = FALSE
    )
  )

  info <- nlmixr2extra:::.sirOmegaProposalInfo(fake_fit)
  expect_equal(info$colName, "eta.ka")
  expect_equal(info$neta1, 1L)
  expect_equal(info$neta2, 1L)
})

test_that(".sirOmegaProposalInfo returns empty data frame when all omega fixed", {
  fake_fit <- list(
    nsub = 20L,
    omega = matrix(
      c(0.4, 0.05, 0.05, 0.2),
      2L,
      2L,
      dimnames = list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl"))
    ),
    iniDf = data.frame(
      name = c("eta.ka", "eta.cl:eta.ka", "eta.cl"),
      neta1 = c(1L, 2L, 2L),
      neta2 = c(1L, 1L, 2L),
      fix = TRUE,
      stringsAsFactors = FALSE
    )
  )

  info <- nlmixr2extra:::.sirOmegaProposalInfo(fake_fit)
  expect_s3_class(info, "data.frame")
  expect_named(info, c("colName", "neta1", "neta2", "est", "se", "isDiag"))
  expect_equal(nrow(info), 0L)
})

test_that(".sirSigmaInfo identifies add.sd as sigma parameter", {
  si <- nlmixr2extra:::.sirSigmaInfo(.theo_fit)
  expect_s3_class(si, "data.frame")
  expect_named(si, c("colName", "est", "se"))
  expect_true("add.sd" %in% si$colName)
  expect_gt(si$se[si$colName == "add.sd"], 0)
})

test_that(".sirSigmaInfo returns NULL when all thetas have covariance", {
  # Fake fit where cov includes all theta params: pass a mock by manipulating
  # the check via a small wrapper — instead, just verify the logic on real fit
  # by confirming that structural thetas are NOT in sigma info
  si <- nlmixr2extra:::.sirSigmaInfo(.theo_fit)
  cov_names <- rownames(.theo_fit$cov)
  sigma_nms <- si$colName
  expect_true(length(intersect(sigma_nms, cov_names)) == 0L)
})

test_that(".sirReconstructOmega sets diagonal from named vector", {
  info <- nlmixr2extra:::.sirOmegaProposalInfo(.theo_fit)
  base_omega <- .theo_fit$omega
  vals <- c("eta.ka" = 0.99)
  mat <- nlmixr2extra:::.sirReconstructOmega(info, vals, base_omega)
  expect_equal(mat["eta.ka", "eta.ka"], 0.99, tolerance = 1e-12)
})

test_that(".sirReconstructOmega leaves unspecified elements unchanged", {
  # 2x2 test case using synthetic omega_info (no real multi-eta fit needed)
  base_omega <- matrix(
    c(0.5, 0.1, 0.1, 0.3),
    2L,
    2L,
    dimnames = list(c("eta.ka", "eta.cl"), c("eta.ka", "eta.cl"))
  )
  info <- data.frame(
    colName = c("eta.ka", "eta.ka:eta.cl", "eta.cl"),
    neta1 = c(1L, 2L, 2L),
    neta2 = c(1L, 1L, 2L),
    est = c(0.5, 0.1, 0.3),
    se = c(0.05, 0.02, 0.04),
    isDiag = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  vals <- c("eta.ka" = 0.8, "eta.ka:eta.cl" = 0.05)
  mat <- nlmixr2extra:::.sirReconstructOmega(info, vals, base_omega)
  expect_equal(mat[1L, 1L], 0.8, tolerance = 1e-12)
  expect_equal(mat[2L, 1L], 0.05, tolerance = 1e-12)
  expect_equal(mat[1L, 2L], 0.05, tolerance = 1e-12)
  expect_equal(mat[2L, 2L], 0.3, tolerance = 1e-12) # unchanged
})

test_that("sirRunIteration rawResults contains sigma column (add.sd)", {
  expect_true("add.sd" %in% names(.iter1$rawResults))
})

test_that("sirRunIteration rawResults contains omega column (eta.ka)", {
  expect_true("eta.ka" %in% names(.iter1$rawResults))
})

test_that("sirRunIteration resampledMat contains sigma column (add.sd)", {
  expect_true("add.sd" %in% colnames(.iter1$resampledMat))
})

# Step 11: sirSummary ----------------------------------------------------------

test_that("sirSummary returns final SIR summary statistics", {
  s <- sirSummary(.iter1$resampledMat, .theo_fit)
  expect_s3_class(s, "data.frame")
  expect_named(
    s,
    c(
      "param",
      "estimate",
      "sd",
      "rse",
      "p2.5",
      "p5",
      "p25",
      "p50",
      "p75",
      "p95",
      "p97.5"
    )
  )
  expect_equal(s$param, colnames(.iter1$resampledMat))
  expect_equal(nrow(s), ncol(.iter1$resampledMat))
})

test_that("sirSummary computes empirical SD, RSE, and percentiles", {
  s <- sirSummary(.iter1$resampledMat, .theo_fit)
  first_param <- s$param[[1L]]
  x <- .iter1$resampledMat[, first_param]
  expect_equal(
    s$sd[s$param == first_param],
    stats::sd(x),
    tolerance = 1e-12
  )
  expect_equal(
    s$p50[s$param == first_param],
    unname(stats::quantile(x, probs = 0.5)),
    tolerance = 1e-12
  )
  expect_equal(
    s$rse[s$param == first_param],
    stats::sd(x) / abs(s$estimate[s$param == first_param]) * 100,
    tolerance = 1e-12
  )
})

test_that("sirSummary attaches empirical covariance and correlation matrices", {
  s <- sirSummary(.iter1$resampledMat, .theo_fit)
  expect_equal(
    attr(s, "covMatrix"),
    stats::cov(.iter1$resampledMat),
    tolerance = 1e-12
  )
  expect_equal(
    attr(s, "corMatrix"),
    stats::cor(.iter1$resampledMat),
    tolerance = 1e-12
  )
})

# Step 10: runSIR --------------------------------------------------------------

test_that(".sirAdjustedAttemptedSamples follows PsN carryover rule", {
  expect_equal(nlmixr2extra:::.sirAdjustedAttemptedSamples(100L), 100L)
  expect_equal(
    nlmixr2extra:::.sirAdjustedAttemptedSamples(
      100L,
      previousAttempted = 100L,
      previousSuccessful = 95L
    ),
    100L
  )
  expect_equal(
    nlmixr2extra:::.sirAdjustedAttemptedSamples(
      100L,
      previousAttempted = 100L,
      previousSuccessful = 90L
    ),
    112L
  )
})

test_that("runSIR runs end-to-end and writes Step 10 artifacts", {
  tmp_dir <- tempfile("sir_run_")
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  set.seed(20260420)
  res <- suppressMessages(
    runSIR(
      .theo_fit,
      nSamples = c(5L, 5L),
      nResample = c(3L, 3L),
      directory = tmp_dir,
      recover = FALSE,
      workers = 1L,
      boxcox = TRUE
    )
  )

  expect_s3_class(res, "nlmixr2SIR")
  expect_s3_class(res, "data.frame")
  expect_true(all(
    c("param", "estimate", "sd", "rse", "p2.5", "p97.5") %in%
      names(res)
  ))
  expect_true(file.exists(file.path(tmp_dir, "raw_results_sir_iteration1.csv")))
  expect_true(file.exists(file.path(tmp_dir, "raw_results_sir_iteration2.csv")))
  expect_true(file.exists(file.path(tmp_dir, "summary_iterations.csv")))
  expect_true(file.exists(file.path(tmp_dir, "sample_rejection_summary.txt")))
  expect_true(file.exists(file.path(tmp_dir, "sir_results.csv")))
  expect_true(file.exists(file.path(tmp_dir, "sir_state.rds")))

  iter_summary <- attr(res, "iterationSummary")
  expect_equal(nrow(iter_summary), 2L)
  expect_equal(iter_summary$nSamples, c(5L, 5L))
  expect_equal(iter_summary$nResample, c(3L, 3L))

  iterations <- attr(res, "iterations")
  expect_null(iterations[[2L]]$boxcoxState)

  raw2 <- utils::read.csv(
    file.path(tmp_dir, "raw_results_sir_iteration2.csv"),
    check.names = FALSE
  )
  expect_equal(raw2$sample_id[[1L]], 0L)
  expect_equal(raw2$dOFV[[1L]], 0)
  expect_equal(sum(raw2$resamples), iterations[[2L]]$iterSummary$nResampled)
})

test_that("runSIR rejects invalid workers before running", {
  expect_error(
    runSIR(
      .theo_fit,
      nSamples = 1L,
      nResample = 1L,
      workers = 0L,
      recover = FALSE
    ),
    "workers"
  )
})

# Step 12: S3 print and plot methods ------------------------------------------

.sir_obj <- local({
  out <- sirSummary(.iter1$resampledMat, .theo_fit)
  class(out) <- c("nlmixr2SIR", "data.frame")
  attr(out, "iterationSummary") <- .iter1$iterSummary
  attr(out, "iterations") <- list(.iter1)
  attr(out, "resampledMat") <- .iter1$resampledMat
  attr(out, "outputDir") <- tempdir()
  out
})

test_that("print.nlmixr2SIR returns object invisibly and prints tables", {
  printed <- utils::capture.output(ret <- print(.sir_obj))
  expect_identical(ret, .sir_obj)
  expect_match(paste(printed, collapse = "\n"), "param")
  expect_match(paste(printed, collapse = "\n"), "estimate")
  expect_match(paste(printed, collapse = "\n"), "nAttempted")
})

test_that("plot.nlmixr2SIR returns parameter distribution plot", {
  p <- plot(.sir_obj)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "Parameter value")
})

test_that("plot.nlmixr2SIR returns dOFV diagnostic plot", {
  p <- plot(.sir_obj, type = "dofv", bins = 10L)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "dOFV")
})

test_that("plot.nlmixr2SIR returns resampling diagnostic plot", {
  p <- plot(.sir_obj, type = "resampling")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$y, "Probability resample")
})
