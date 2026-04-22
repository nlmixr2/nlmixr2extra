# Sampling Importance Resampling (SIR) ----------------------------------------
#
# Reference: Dosne et al. (2013), PAGE 22, Abstract 2907.
# Algorithm mirrors PsN's sir tool (psn_ref/sir, sir.pm, sir_userguide.tex).

# Helpers ----------------------------------------------------------------------

# Apply bounds elementwise across rows of a matrix; returns logical vector.
.sirInBounds <- function(mat, lower, upper) {
  apply(mat, 1L, function(x) all(x >= lower & x <= upper))
}

.sirCapCovCorrelation <- function(covMat, capCorrelation = 0.8) {
  checkmate::assertMatrix(covMat, mode = "numeric")
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)

  cov_mat <- (covMat + t(covMat)) / 2
  if (capCorrelation >= 1 || ncol(cov_mat) < 2L) {
    return(cov_mat)
  }

  param_names <- colnames(cov_mat)
  sd_vals <- sqrt(pmax(diag(cov_mat), 0))
  corr_mat <- suppressWarnings(cov2cor(cov_mat))
  corr_mat[!is.finite(corr_mat)] <- 0
  diag(corr_mat) <- 1

  lo <- lower.tri(corr_mat)
  corr_mat[lo] <- pmax(-capCorrelation, pmin(capCorrelation, corr_mat[lo]))
  corr_mat[upper.tri(corr_mat)] <- t(corr_mat)[upper.tri(corr_mat)]

  capped <- outer(sd_vals, sd_vals) * corr_mat
  dimnames(capped) <- list(param_names, param_names)
  capped
}

.sirEnsurePosDef <- function(covMat, minEigen = sqrt(.Machine$double.eps)) {
  checkmate::assertMatrix(covMat, mode = "numeric")
  dim_names <- dimnames(covMat)
  cov_mat <- (covMat + t(covMat)) / 2
  eig <- eigen(cov_mat, symmetric = TRUE)
  if (all(is.finite(eig$values)) && min(eig$values) >= minEigen) {
    dimnames(cov_mat) <- dim_names
    return(cov_mat)
  }

  values <- pmax(eig$values, minEigen)
  out <- eig$vectors %*% diag(values, nrow = length(values)) %*% t(eig$vectors)
  out <- (out + t(out)) / 2
  dimnames(out) <- dim_names
  out
}

# Step 1 -----------------------------------------------------------------------

#' Extract and inflate the proposal covariance matrix for SIR
#'
#' Pulls the parameter covariance matrix from a fitted nlmixr2 model and
#' optionally inflates variances and caps correlations to widen the proposal
#' distribution for the first SIR iteration.
#'
#' Inflation is applied per parameter type: `thetaInflation` scales the
#' variance of every THETA parameter, `omegaInflation` scales every OMEGA
#' element.  Off-diagonal elements are recomputed so that correlations are
#' preserved after inflation.  `capCorrelation` then hard-clips all absolute
#' pairwise correlations.
#'
#' @param fit An nlmixr2 fit with a successful covariance step.
#' @param thetaInflation Non-negative scalar. Variance multiplier for THETA
#'   parameters. Default `1` (no inflation).
#' @param omegaInflation Non-negative scalar. Variance multiplier for OMEGA
#'   elements. Default `1`.
#' @param sigmaInflation Not used in nlmixr2 (residual error parameters are
#'   THETAs). Retained for API compatibility with PsN. Default `1`.
#' @param capCorrelation Numeric in \[0, 1\]. All absolute pairwise
#'   correlations are capped at this value after inflation. Set to `1` to
#'   disable. Default `0.8`.
#' @return A named, symmetric covariance matrix with the same row/column order
#'   as `fit$cov`.
#' @noRd
sirGetProposalCov <- function(
  fit,
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8
) {
  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertNumber(thetaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(omegaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(sigmaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)

  cov_mat <- fit$cov
  if (is.null(cov_mat) || nrow(cov_mat) == 0L) {
    cli::cli_abort(c(
      "No covariance matrix found in {.arg fit}.",
      "i" = "Re-run the model with a successful covariance step.",
      "i" = "Example: {.code foceiControl(covMethod = \"r\")}"
    ))
  }

  param_names <- rownames(cov_mat)
  ini_df <- fit$iniDf
  theta_names <- ini_df$name[!is.na(ini_df$ntheta) & !ini_df$fix]

  # One inflation factor per row/col of the covariance matrix
  is_theta <- param_names %in% theta_names
  inflation <- ifelse(is_theta, thetaInflation, omegaInflation)

  # Preserve correlations: new SD_i = old SD_i * sqrt(inflation_i)
  # => new cov_ij = old corr_ij * new SD_i * new SD_j
  sd_orig <- sqrt(diag(cov_mat))
  sd_new <- sd_orig * sqrt(inflation)
  corr_mat <- cov2cor(cov_mat)

  # Reconstruct covariance from (capped) correlations and inflated SDs
  new_cov <- outer(sd_new, sd_new) * corr_mat
  dimnames(new_cov) <- list(param_names, param_names)

  .sirCapCovCorrelation(new_cov, capCorrelation = capCorrelation)
}

# Step 2 -----------------------------------------------------------------------

#' Sample THETA vectors from a (truncated) multivariate normal
#'
#' Draws `n` parameter vectors from a multivariate normal distribution and
#' discards any that violate the supplied bounds.  Sampling is repeated in
#' batches until `n` valid vectors are collected or the maximum attempt
#' budget is exhausted.
#'
#' @param mu Named numeric vector of THETA means (length p).
#' @param covMat p×p covariance matrix.
#' @param n Positive integer. Number of valid samples to collect.
#' @param lower Numeric vector (length p, or scalar recycled). Lower bounds.
#'   Default `-Inf` (no lower truncation).
#' @param upper Numeric vector (length p, or scalar recycled). Upper bounds.
#'   Default `Inf` (no upper truncation).
#' @return A named list:
#'   \describe{
#'     \item{`samples`}{Numeric matrix, at most `n` rows × p columns, with
#'       column names from `mu`.  Fewer than `n` rows if the budget was
#'       exhausted before `n` valid draws were found.}
#'     \item{`nRejected`}{Integer. Total number of draws that violated at
#'       least one bound.}
#'   }
#' @noRd
sirSampleTheta <- function(
  mu,
  covMat,
  n,
  lower = rep(-Inf, length(mu)),
  upper = rep(Inf, length(mu))
) {
  p <- length(mu)
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, min.len = 1L)
  checkmate::assertMatrix(covMat, mode = "numeric", nrows = p, ncols = p)
  checkmate::assertCount(n, positive = TRUE)
  if (length(lower) == 1L) {
    lower <- rep(lower, p)
  }
  if (length(upper) == 1L) {
    upper <- rep(upper, p)
  }
  checkmate::assertNumeric(lower, len = p, any.missing = FALSE)
  checkmate::assertNumeric(upper, len = p, any.missing = FALSE)

  n <- as.integer(n)
  max_attempts <- 10L * n
  param_names <- names(mu)

  collected <- matrix(NA_real_, nrow = n, ncol = p)
  n_filled <- 0L
  n_rejected <- 0L
  n_attempted <- 0L

  while (n_filled < n && n_attempted < max_attempts) {
    batch_n <- min(n - n_filled, max_attempts - n_attempted)
    draws <- mvtnorm::rmvnorm(batch_n, mean = mu, sigma = covMat)
    n_attempted <- n_attempted + batch_n

    ok <- .sirInBounds(draws, lower, upper)
    good <- draws[ok, , drop = FALSE]
    n_rejected <- n_rejected + sum(!ok)

    n_take <- min(nrow(good), n - n_filled)
    if (n_take > 0L) {
      idx <- seq(n_filled + 1L, n_filled + n_take)
      collected[idx, ] <- good[seq_len(n_take), ]
      n_filled <- n_filled + n_take
    }
  }

  if (n_filled < n) {
    cli::cli_warn(c(
      "Only {n_filled} of {n} requested samples were within bounds \\
       after {max_attempts} draw attempts.",
      "i" = "Consider widening the parameter bounds or reducing {.arg n}."
    ))
    collected <- collected[seq_len(n_filled), , drop = FALSE]
  }

  colnames(collected) <- param_names
  list(samples = collected, nRejected = n_rejected)
}

# Step 3 -----------------------------------------------------------------------

#' Sample OMEGA matrices from a multivariate normal, retaining only PD draws
#'
#' Vectorizes the lower triangle (including diagonal) of `omegaEst`, draws from
#' a multivariate normal with covariance `omegaCovMat`, reconstructs symmetric
#' matrices, and discards any that are not positive definite.
#'
#' @param omegaEst Square symmetric matrix. Current OMEGA estimate.
#' @param omegaCovMat Square numeric matrix of size `p × p`, where
#'   `p = length(lower.tri(omegaEst, diag = TRUE))`. Covariance of the
#'   lower-triangle elements (e.g., constructed from SEs in `fit$parFixedDf`).
#' @param n Positive integer. Number of valid (positive-definite) samples to
#'   collect.
#' @return A named list:
#'   \describe{
#'     \item{`samples`}{List of at most `n` positive-definite matrices with the
#'       same `dimnames` as `omegaEst`. Fewer than `n` entries if the budget was
#'       exhausted.}
#'     \item{`nRejected`}{Integer. Draws that failed the Cholesky check.}
#'   }
#' @noRd
sirSampleOmegaSigma <- function(omegaEst, omegaCovMat, n) {
  checkmate::assertMatrix(omegaEst, mode = "numeric")
  if (nrow(omegaEst) != ncol(omegaEst)) {
    cli::cli_abort("{.arg omegaEst} must be a square matrix.")
  }
  dim_omega <- nrow(omegaEst)
  lt_idx <- which(lower.tri(omegaEst, diag = TRUE))
  mu <- omegaEst[lt_idx]
  p <- length(mu)
  checkmate::assertMatrix(omegaCovMat, mode = "numeric", nrows = p, ncols = p)
  checkmate::assertCount(n, positive = TRUE)

  n <- as.integer(n)
  max_attempts <- 10L * n

  collected <- vector("list", n)
  n_filled <- 0L
  n_rejected <- 0L
  n_attempted <- 0L

  while (n_filled < n && n_attempted < max_attempts) {
    batch_n <- min(n - n_filled, max_attempts - n_attempted)
    draws <- mvtnorm::rmvnorm(batch_n, mean = mu, sigma = omegaCovMat)
    n_attempted <- n_attempted + batch_n

    for (i in seq_len(nrow(draws))) {
      mat <- matrix(0, dim_omega, dim_omega)
      mat[lt_idx] <- draws[i, ]
      mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
      dimnames(mat) <- dimnames(omegaEst)

      is_pd <- tryCatch(
        {
          chol(mat)
          TRUE
        },
        error = function(e) FALSE
      )
      if (is_pd) {
        n_filled <- n_filled + 1L
        collected[[n_filled]] <- mat
        if (n_filled == n) break
      } else {
        n_rejected <- n_rejected + 1L
      }
    }
  }

  if (n_filled < n) {
    cli::cli_warn(c(
      "Only {n_filled} of {n} OMEGA samples were positive definite \\
       after {max_attempts} draw attempts.",
      "i" = "Consider revising {.arg omegaCovMat} or reducing {.arg n}."
    ))
    collected <- collected[seq_len(n_filled)]
  }

  list(samples = collected, nRejected = n_rejected)
}

# Step 4 -----------------------------------------------------------------------

#' Evaluate OFV for a matrix of sampled parameter vectors
#'
#' For each row of `paramSamples`, seeds the model's ini block with the
#' supplied values (theta and/or lower-triangle omega elements) and calls
#' `nlmixr2(est = "focei")` with `maxOuterIterations = 0`, making it the
#' nlmixr2 equivalent of NONMEM `MAXEVAL=0`.  THETA column names must match
#' parameter names in `fit$iniDf`; omega columns use the SIR lower-triangle
#' proposal names.
#'
#' Parameters NOT present as columns in `paramSamples` retain their estimated
#' values from `fit`.
#'
#' @param fit An nlmixr2 fit object (carries the data internally).
#' @param paramSamples Numeric matrix, one row per sample.  Column names may
#'   include THETA names from `fit$iniDf` and SIR omega lower-triangle names.
#' @param workers Passed to `.withWorkerPlan()`: `NULL` (keep current plan),
#'   `1` (force sequential), a positive integer, or `"auto"`.
#' @return Named numeric vector of length `nrow(paramSamples)`.  Entries are
#'   `NA_real_` for rows that produced an error during evaluation.
#' @noRd
sirEvalOFV <- function(fit, paramSamples, workers = NULL) {
  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertMatrix(
    paramSamples,
    mode = "numeric",
    min.rows = 1L,
    min.cols = 1L
  )

  ini_df <- fit$iniDf
  theta_names <- ini_df$name[!is.na(ini_df$ntheta)]
  omega_info <- .sirOmegaProposalInfo(fit)

  col_names <- colnames(paramSamples)
  theta_cols <- intersect(col_names, theta_names)
  omega_cols <- intersect(col_names, omega_info$colName)
  has_omega <- length(omega_cols) > 0L

  if (length(theta_cols) == 0L && !has_omega) {
    cli::cli_abort(
      "No column names in {.arg paramSamples} match any parameter in {.arg fit}."
    )
  }

  base_omega <- fit$omega
  base_theta <- fit$theta

  eval_one <- function(i) {
    row <- paramSamples[i, ]

    theta_vals <- base_theta
    if (length(theta_cols) > 0L) {
      theta_vals[theta_cols] <- row[theta_cols]
    }

    ini_args <- as.list(theta_vals)

    if (has_omega) {
      omega_mat <- .sirReconstructOmega(omega_info, row, base_omega)
      eta_names <- rownames(omega_mat)
      lt_vals <- unlist(
        lapply(seq_len(nrow(omega_mat)), function(r) {
          omega_mat[r, seq_len(r)]
        }),
        use.names = FALSE
      )
      lhs <- paste(eta_names, collapse = " + ")
      lt_txt <- format(lt_vals, scientific = TRUE, digits = 17, trim = TRUE)
      rhs <- if (length(lt_vals) == 1L) {
        lt_txt
      } else {
        paste0("c(", paste(lt_txt, collapse = ", "), ")")
      }
      omega_expr <- str2lang(paste(lhs, "~", rhs))
      ini_args <- c(ini_args, list(omega_expr))
    }

    tryCatch(
      {
        model_new <- suppressMessages(
          do.call(rxode2::ini, c(list(x = fit), ini_args))
        )
        f <- suppressMessages(
          nlmixr2est::nlmixr2(
            model_new,
            est = "focei",
            control = nlmixr2est::foceiControl(
              calcTables = FALSE,
              covMethod = "",
              compress = FALSE,
              maxOuterIterations = 0L,
              print = 0L
            )
          )
        )
        f$objf
      },
      error = function(e) NA_real_
    )
  }

  results <- .withWorkerPlan(workers, {
    .plap(
      seq_len(nrow(paramSamples)),
      eval_one,
      .label = function(i) sprintf("sample %d", i)
    )
  })

  unlist(results, use.names = FALSE)
}

# Step 5 -----------------------------------------------------------------------

#' Compute relative PDFs, importance ratios, and resampling probabilities
#'
#' Given a matrix of sampled parameter vectors, their delta-OFV values relative
#' to the original fit, and the proposal distribution (mean + covariance),
#' computes the importance ratio for each sample and normalises them into
#' resampling probabilities.
#'
#' The relative PDF is the multivariate normal density evaluated at each sample
#' *relative* to the density at `mu`, so it equals 1 at the proposal mean and
#' decreases for samples further away.  The likelihood ratio converts the OFV
#' difference back to a probability scale.  Samples with `NA` `dOFV` (failed
#' evaluations) receive `prob_resample = 0` and are excluded from resampling.
#'
#' @param samples Numeric matrix, one row per sample (same parameter space as
#'   `mu` / `covMat`).
#' @param mu Named numeric vector. Proposal mean (length = `ncol(samples)`).
#' @param covMat Positive-definite covariance matrix of the proposal (same
#'   dimension as `length(mu)`).
#' @param dOFV Numeric vector of length `nrow(samples)`.
#'   `dOFV[i] = OFV_sample[i] - OFV_original`.  May contain `NA`.
#' @return Data frame with one row per sample and columns:
#'   \describe{
#'     \item{`sample_id`}{Integer row index.}
#'     \item{`dOFV`}{As supplied.}
#'     \item{`likelihood_ratio`}{`exp(-0.5 * dOFV)`.}
#'     \item{`relPDF`}{Relative proposal density (1 at `mu`).}
#'     \item{`importance_ratio`}{`likRatio / relPDF`; `NA` when `dOFV` is `NA`.}
#'     \item{`prob_resample`}{Normalised resampling probability; 0 for `NA` rows.}
#'   }
#' @noRd
sirCalcWeights <- function(samples, mu, covMat, dOFV) {
  n <- nrow(samples)
  p <- ncol(samples)
  checkmate::assertMatrix(samples, mode = "numeric", min.rows = 1L)
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, len = p)
  checkmate::assertMatrix(covMat, mode = "numeric", nrows = p, ncols = p)
  checkmate::assertNumeric(dOFV, len = n)

  L <- tryCatch(
    chol(covMat),
    error = function(e) {
      cli::cli_abort("{.arg covMat} is not positive definite.")
    }
  )

  # log(relPDF_i) = -0.5 * ||L^{-T}(x_i - mu)||^2; equals 0 when x_i = mu.
  log_rel_pdf <- vapply(
    seq_len(n),
    function(i) {
      z <- backsolve(L, samples[i, ] - mu)
      -0.5 * sum(z^2)
    },
    numeric(1L)
  )

  log_lik_ratio <- -0.5 * dOFV
  log_ir <- log_lik_ratio - log_rel_pdf

  valid <- is.finite(log_ir)
  if (!any(valid)) {
    cli::cli_abort(c(
      "No finite SIR importance weights could be computed.",
      "i" = "All OFV evaluations may have failed, or the proposal density is numerically degenerate."
    ))
  }

  max_log_ir <- max(log_ir[valid])
  scaled <- pmax(exp(log_ir[valid] - max_log_ir), .Machine$double.xmin)
  scaled_sum <- sum(scaled)
  if (!is.finite(scaled_sum) || scaled_sum <= 0) {
    cli::cli_abort(c(
      "SIR importance weights could not be normalized.",
      "i" = "Check failed OFV evaluations, proposal covariance, and resampling diagnostics."
    ))
  }

  prob <- numeric(n)
  prob[valid] <- scaled / scaled_sum

  rel_pdf <- exp(log_rel_pdf)
  lik_ratio <- exp(log_lik_ratio)
  ir <- exp(log_ir)
  ir[!valid] <- NA_real_

  data.frame(
    sample_id = seq_len(n),
    dOFV = dOFV,
    likelihood_ratio = lik_ratio,
    relPDF = rel_pdf,
    importance_ratio = ir,
    prob_resample = prob
  )
}

# Step 6 -----------------------------------------------------------------------

#' Weighted resample of parameter vectors with optional appearance cap
#'
#' Draws `m` rows from `samples` using `weights$prob_resample` as
#' probabilities.  `capResampling = 1` follows PsN's default and samples
#' without replacement.  Values greater than one allow limited replacement so
#' no original sample appears more than `capResampling` times.
#'
#' @param samples Numeric matrix. Original sampled parameter vectors (one row
#'   per sample).
#' @param weights Data frame returned by `sirCalcWeights()`.  Must contain a
#'   `prob_resample` column of length `nrow(samples)`.
#' @param m Positive integer. Number of resampled vectors to return.
#' @param capResampling Positive number. Default `1` samples without
#'   replacement; values `> 1` allow limited replacement.
#' @return A named list:
#'   \describe{
#'     \item{`samples`}{Numeric matrix with `m` rows selected from the
#'       originals.}
#'     \item{`resampleCounts`}{Integer vector of length `nrow(samples)`.
#'       `resampleCounts[i]` is how many times original row `i` was selected.}
#'     \item{`sampleOrder`}{Character vector of length `nrow(samples)`.
#'       Semicolon-separated draw positions for each original row.}
#'     \item{`selectionOrder`}{Integer vector of selected row indices, in draw
#'       order.}
#'   }
#' @noRd
sirResample <- function(samples, weights, m, capResampling = 1) {
  n <- nrow(samples)
  checkmate::assertMatrix(samples, mode = "numeric", min.rows = 1L)
  checkmate::assertDataFrame(weights)
  checkmate::assertNumeric(
    weights$prob_resample,
    len = n,
    lower = 0,
    any.missing = FALSE
  )
  checkmate::assertCount(m, positive = TRUE)
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)

  m <- as.integer(m)
  prob <- weights$prob_resample
  valid <- is.finite(prob) & prob > 0
  cap <- as.integer(floor(capResampling))

  if (!any(valid)) {
    cli::cli_abort(c(
      "No samples have non-zero SIR resampling probability.",
      "i" = "Check OFV failures and importance weight diagnostics."
    ))
  }
  if (m > sum(valid) * cap) {
    cli::cli_abort(c(
      "Cannot draw {m} SIR resamples with {.arg capResampling} = {cap}.",
      "i" = "Only {sum(valid)} samples have non-zero probability."
    ))
  }

  if (cap <= 1L) {
    idx <- sample.int(n, m, replace = FALSE, prob = prob)
  } else {
    expanded_idx <- rep(seq_len(n), each = cap)
    expanded_prob <- rep(prob, each = cap)
    keep <- expanded_prob > 0 & is.finite(expanded_prob)
    idx <- sample(
      expanded_idx[keep],
      m,
      replace = FALSE,
      prob = expanded_prob[keep]
    )
  }

  resampled <- samples[idx, , drop = FALSE]
  rownames(resampled) <- NULL
  resample_counts <- tabulate(idx, nbins = n)
  sample_order <- vapply(
    seq_len(n),
    function(i) {
      pos <- which(idx == i)
      if (length(pos) == 0L) "" else paste(pos, collapse = ";")
    },
    character(1L)
  )

  list(
    samples = resampled,
    resampleCounts = resample_counts,
    sampleOrder = sample_order,
    selectionOrder = idx
  )
}

.sirBuildRawResults <- function(
  paramMat,
  weights,
  dOFV,
  resampled,
  mu,
  capResampling = 1
) {
  checkmate::assertMatrix(paramMat, mode = "numeric", min.rows = 1L)
  checkmate::assertDataFrame(weights, nrows = nrow(paramMat))
  checkmate::assertNumeric(dOFV, len = nrow(paramMat))
  checkmate::assertList(resampled)
  checkmate::assertNumeric(mu, len = ncol(paramMat))
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)

  param_names <- colnames(paramMat)
  cap <- as.integer(floor(capResampling))
  raw_df <- data.frame(
    sample_id = seq_len(nrow(paramMat)),
    as.data.frame(paramMat, check.names = FALSE),
    check.names = FALSE
  )
  raw_df$dOFV <- dOFV
  raw_df$deltaofv <- dOFV
  raw_df$likelihood_ratio <- weights$likelihood_ratio
  raw_df$relPDF <- weights$relPDF
  raw_df$importance_ratio <- weights$importance_ratio
  raw_df$IR <- weights$importance_ratio
  raw_df$probability_resample <- weights$prob_resample
  raw_df$prob <- weights$prob_resample

  expanded <- lapply(seq_len(nrow(raw_df)), function(i) {
    out <- raw_df[rep(i, cap), , drop = FALSE]
    orders <- resampled$sampleOrder[[i]]
    orders <- if (identical(orders, "")) {
      character(0L)
    } else {
      strsplit(orders, ";", fixed = TRUE)[[1L]]
    }
    n_selected <- min(length(orders), cap)
    out$resamples <- 0L
    out$sample_order <- ""
    if (n_selected > 0L) {
      idx <- seq_len(n_selected)
      out$resamples[idx] <- 1L
      out$sample_order[idx] <- orders[idx]
    }
    out
  })
  expanded <- do.call(rbind, expanded)
  rownames(expanded) <- NULL

  mu_mat <- matrix(
    mu[param_names],
    nrow = 1L,
    dimnames = list(NULL, param_names)
  )
  mu_df <- data.frame(
    sample_id = 0L,
    as.data.frame(mu_mat, check.names = FALSE),
    check.names = FALSE
  )
  mu_df$dOFV <- 0
  mu_df$deltaofv <- 0
  mu_df$likelihood_ratio <- 1
  mu_df$relPDF <- 1
  mu_df$importance_ratio <- NA_real_
  mu_df$IR <- NA_real_
  mu_df$probability_resample <- 0
  mu_df$prob <- 0
  mu_df$resamples <- 0L
  mu_df$sample_order <- ""

  out <- rbind(mu_df, expanded)
  rownames(out) <- NULL
  out
}

# Step 7 -----------------------------------------------------------------------

#' Box-Cox transform a vector, optionally estimating lambda
#'
#' Applies the Box-Cox power transformation after shifting `x` by `delta` to
#' ensure strict positivity.  When `lambda` is `NULL` it is estimated by
#' maximising the Pearson correlation between the sorted transformed values and
#' their expected normal order statistics (normal scores), using
#' `stats::optimize()` over the interval \[−3, 3\].
#'
#' Transformation:
#' * lambda ≠ 0: `(x + delta)^lambda − 1) / lambda`
#' * lambda = 0: `log(x + delta)`
#'
#' @param x Numeric vector (length ≥ 2), no missing values.
#' @param lambda `NULL` to estimate, or a finite numeric scalar.
#' @param delta `NULL` to compute as `|min(x)| + 1e-6`, or a non-negative
#'   scalar.  Must ensure `x + delta > 0`.
#' @return Named list: `transformed` (numeric vector), `lambda` (scalar),
#'   `delta` (scalar).
#' @noRd
sirBoxCox <- function(x, lambda = NULL, delta = NULL) {
  checkmate::assertNumeric(x, min.len = 2L, any.missing = FALSE, finite = TRUE)

  if (is.null(delta)) {
    delta <- abs(min(x)) + 1e-6
  }
  checkmate::assertNumber(delta, lower = 0, finite = TRUE)

  x_shifted <- x + delta
  if (any(x_shifted <= 0)) {
    cli::cli_abort(
      "{.arg delta} is too small: {.code x + delta} must be strictly positive."
    )
  }

  .bc <- function(xs, lam) {
    if (abs(lam) < 1e-10) log(xs) else (xs^lam - 1) / lam
  }

  if (is.null(lambda)) {
    if (stats::var(x) < .Machine$double.eps) {
      lambda <- 1 # degenerate: all identical, any lambda works
    } else {
      nscores <- qnorm(stats::ppoints(length(x)))
      obj <- function(lam) {
        xt <- sort(.bc(x_shifted, lam))
        v <- suppressWarnings(cor(xt, nscores))
        if (is.na(v)) 0 else -v
      }
      lambda <- stats::optimize(obj, interval = c(-3, 3))$minimum
    }
  }
  checkmate::assertNumber(lambda, finite = TRUE)

  list(transformed = .bc(x_shifted, lambda), lambda = lambda, delta = delta)
}

#' Invert a Box-Cox transformation
#'
#' Recovers the original scale given the transformed values, lambda, and delta
#' returned by `sirBoxCox()`.
#'
#' @param x_transformed Numeric vector of transformed values.
#' @param lambda Scalar lambda used in the forward transform.
#' @param delta Scalar shift used in the forward transform.
#' @return Numeric vector on the original scale.
#' @noRd
sirBoxCoxInverse <- function(x_transformed, lambda, delta) {
  checkmate::assertNumeric(x_transformed, any.missing = FALSE, finite = TRUE)
  checkmate::assertNumber(lambda, finite = TRUE)
  checkmate::assertNumber(delta, finite = TRUE)

  if (abs(lambda) < 1e-10) {
    exp(x_transformed) - delta
  } else {
    base <- lambda * x_transformed + 1
    if (any(base <= 0)) {
      cli::cli_abort(
        "Inverse Box-Cox: {.code lambda * x_transformed + 1} must be positive."
      )
    }
    base^(1 / lambda) - delta
  }
}

# Step 8 -----------------------------------------------------------------------

#' Update the proposal covariance from resampled parameter vectors
#'
#' Computes the empirical covariance of a resampled parameter matrix.  When
#' `boxcox = TRUE`, each column is first Box-Cox transformed via `sirBoxCox()`
#' to reduce skewness before the covariance is calculated; the resulting
#' covariance is in the transformed space and the per-column lambda/delta values
#' are returned so callers can back-transform sampled points if needed.
#'
#' @param resampledMat Numeric matrix with at least 2 rows (resampled parameter
#'   vectors, one per row).  Column names are preserved in the output.
#' @param boxcox Logical.  If `TRUE` (default), apply column-wise Box-Cox
#'   before computing covariance.
#' @param capCorrelation Numeric in \[0, 1\]. Caps absolute correlations in
#'   the updated proposal. Default `0.8`.
#' @return Named list:
#'   \describe{
#'     \item{`covMat`}{Symmetric positive-(semi)definite covariance matrix.}
#'     \item{`boxcoxParams`}{Data frame with columns `param`, `lambda`, `delta`
#'       (one row per column of `resampledMat`), or `NULL` when
#'       `boxcox = FALSE`.}
#'   }
#' @noRd
sirUpdateProposal <- function(
  resampledMat,
  boxcox = TRUE,
  capCorrelation = 0.8
) {
  checkmate::assertMatrix(
    resampledMat,
    mode = "numeric",
    min.rows = 2L,
    min.cols = 1L
  )
  checkmate::assertFlag(boxcox)
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)

  param_names <- colnames(resampledMat)
  n_col <- ncol(resampledMat)

  if (boxcox) {
    bc_list <- lapply(seq_len(n_col), function(j) sirBoxCox(resampledMat[, j]))

    trans_mat <- matrix(
      unlist(lapply(bc_list, `[[`, "transformed")),
      nrow = nrow(resampledMat),
      ncol = n_col,
      dimnames = list(NULL, param_names)
    )

    bc_params <- data.frame(
      param = if (is.null(param_names)) seq_len(n_col) else param_names,
      lambda = vapply(bc_list, `[[`, numeric(1L), "lambda"),
      delta = vapply(bc_list, `[[`, numeric(1L), "delta"),
      stringsAsFactors = FALSE
    )
  } else {
    trans_mat <- resampledMat
    bc_params <- NULL
  }

  cov_mat <- cov(trans_mat)
  dimnames(cov_mat) <- list(param_names, param_names)
  cov_mat <- .sirCapCovCorrelation(cov_mat, capCorrelation = capCorrelation)
  cov_mat <- .sirEnsurePosDef(cov_mat)

  list(covMat = cov_mat, boxcoxParams = bc_params)
}

# Step 9 -----------------------------------------------------------------------

# Transform a named mu vector into the sampling space using stored Box-Cox params.
.sirBcTransformMu <- function(mu, bc_state) {
  if (is.null(bc_state)) {
    return(mu)
  }
  bc_mu <- vapply(
    seq_along(mu),
    function(j) {
      nm <- names(mu)[j]
      row <- bc_state[bc_state$param == nm, ]
      if (nrow(row) != 1L) {
        cli::cli_abort(
          "Missing Box-Cox parameters for SIR parameter {.val {nm}}."
        )
      }
      xs <- mu[j] + row$delta
      if (!is.finite(xs) || xs <= 0) {
        cli::cli_abort(
          "Box-Cox shift for SIR parameter {.val {nm}} is not positive."
        )
      }
      lam <- row$lambda
      if (abs(lam) < 1e-10) log(xs) else (xs^lam - 1) / lam
    },
    numeric(1L)
  )
  setNames(bc_mu, names(mu))
}

.sirBcInverseMatrix <- function(mat, bc_state) {
  if (is.null(bc_state)) {
    return(mat)
  }
  out <- mat
  for (j in seq_len(ncol(mat))) {
    nm <- colnames(mat)[j]
    row <- bc_state[bc_state$param == nm, ]
    if (nrow(row) != 1L) {
      cli::cli_abort(
        "Missing Box-Cox parameters for SIR parameter {.val {nm}}."
      )
    }
    out[, j] <- vapply(
      mat[, j],
      function(v) {
        tryCatch(
          sirBoxCoxInverse(v, row$lambda, row$delta),
          error = function(e) NA_real_
        )
      },
      numeric(1L)
    )
  }
  out
}

# Returns data frame for all lower-triangle omega elements, ordered column-major
# (matching sirSampleOmegaSigma's internal lower.tri ordering).
# Columns: colName, neta1, neta2, est, se, isDiag
.sirOmegaProposalInfo <- function(
  fit,
  n_sub = NULL,
  omegaFallback = c("wishart"),
  omegaDf = NULL
) {
  omegaFallback <- match.arg(omegaFallback)
  if (is.null(n_sub)) {
    n_sub <- .sirNSubjects(fit)
  }
  df <- if (is.null(omegaDf)) max(n_sub - 1L, 1L) else omegaDf
  checkmate::assertNumber(df, lower = 1, finite = TRUE)
  ini_df <- fit$iniDf
  omega_est <- fit$omega

  omega_rows <- ini_df[!is.na(ini_df$neta1), , drop = FALSE]
  omega_fixed <- if ("fix" %in% names(omega_rows)) {
    !is.na(omega_rows$fix) & omega_rows$fix
  } else {
    rep(FALSE, nrow(omega_rows))
  }
  lower_rows <- omega_rows[
    !omega_fixed & omega_rows$neta1 >= omega_rows$neta2,
    ,
    drop = FALSE
  ]
  lower_rows <- lower_rows[
    order(lower_rows$neta2, lower_rows$neta1),
    ,
    drop = FALSE
  ]

  diag_rows <- omega_rows[omega_rows$neta1 == omega_rows$neta2, , drop = FALSE]
  idx2name <- setNames(diag_rows$name, as.character(diag_rows$neta1))

  if (nrow(lower_rows) == 0L) {
    return(data.frame(
      colName = character(0),
      neta1 = integer(0),
      neta2 = integer(0),
      est = numeric(0),
      se = numeric(0),
      isDiag = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(
    rbind,
    lapply(seq_len(nrow(lower_rows)), function(i) {
      r <- lower_rows$neta1[i]
      c_idx <- lower_rows$neta2[i]
      is_diag <- r == c_idx
      col_nm <- if (is_diag) {
        idx2name[[as.character(r)]]
      } else {
        paste(
          idx2name[[as.character(r)]],
          idx2name[[as.character(c_idx)]],
          sep = ":"
        )
      }
      omega_rc <- unname(omega_est[r, c_idx])
      se <- if (is_diag) {
        sqrt(2 * omega_rc^2 / df)
      } else {
        sqrt((omega_est[r, r] * omega_est[c_idx, c_idx] + omega_rc^2) / df)
      }
      data.frame(
        colName = col_nm,
        neta1 = r,
        neta2 = c_idx,
        est = omega_rc,
        se = se,
        isDiag = is_diag,
        stringsAsFactors = FALSE
      )
    })
  )
}

# Returns sigma-type theta info: thetas in fit$theta but absent from fit$cov.
# NULL when every theta has a covariance estimate.
.sirSigmaInfo <- function(fit, sigmaFallbackRse = 30) {
  checkmate::assertNumber(sigmaFallbackRse, lower = 0, finite = TRUE)
  ini_df <- fit$iniDf
  all_theta <- ini_df$name[!is.na(ini_df$ntheta) & !ini_df$fix]
  sigma_names <- setdiff(all_theta, rownames(fit$cov))
  if (length(sigma_names) == 0L) {
    return(NULL)
  }

  pf <- tryCatch(fit$parFixedDf, error = function(e) NULL)
  do.call(
    rbind,
    lapply(sigma_names, function(nm) {
      est <- fit$theta[[nm]]
      se <- if (!is.null(pf) && nm %in% rownames(pf) && !is.na(pf[nm, "SE"])) {
        pf[nm, "SE"]
      } else {
        abs(est) * sigmaFallbackRse / 100
      }
      data.frame(colName = nm, est = est, se = se, stringsAsFactors = FALSE)
    })
  )
}

# Reconstruct the full symmetric omega matrix from a named numeric vector.
# omega_info: from .sirOmegaProposalInfo; base_omega: fit$omega for defaults.
.sirReconstructOmega <- function(omega_info, vals, base_omega) {
  mat <- base_omega
  for (i in seq_len(nrow(omega_info))) {
    cn <- omega_info$colName[i]
    if (cn %in% names(vals)) {
      r <- omega_info$neta1[i]
      c_idx <- omega_info$neta2[i]
      mat[r, c_idx] <- vals[[cn]]
      mat[c_idx, r] <- vals[[cn]]
    }
  }
  mat
}

# Derive per-subject count for chi-squared omega SE approximation.
.sirNSubjects <- function(fit) {
  n <- tryCatch(as.integer(fit$nsub), error = function(e) NULL)
  if (!is.null(n) && length(n) == 1L && !is.na(n)) {
    return(n)
  }
  dat <- tryCatch(fit$origData, error = function(e) NULL)
  if (!is.null(dat) && "ID" %in% names(dat)) {
    return(length(unique(dat$ID)))
  }
  cli::cli_abort("Cannot determine subject count from {.arg fit}.")
}

.sirInitialProposal <- function(
  fit,
  mu,
  proposalCov,
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8,
  omegaFallback = c("wishart"),
  sigmaFallbackRse = 30,
  omegaDf = NULL
) {
  omegaFallback <- match.arg(omegaFallback)
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, min.len = 1L)
  checkmate::assertMatrix(proposalCov, mode = "numeric")
  checkmate::assertNumber(thetaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(omegaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(sigmaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)

  theta_names <- rownames(fit$cov)
  n_sub <- .sirNSubjects(fit)
  omega_info <- .sirOmegaProposalInfo(
    fit,
    n_sub = n_sub,
    omegaFallback = omegaFallback,
    omegaDf = omegaDf
  )
  sigma_info <- .sirSigmaInfo(fit, sigmaFallbackRse = sigmaFallbackRse)
  sigma_names <- if (is.null(sigma_info)) character(0) else sigma_info$colName
  full_names <- c(theta_names, omega_info$colName, sigma_names)

  mu_full <- setNames(numeric(length(full_names)), full_names)
  for (nm in full_names) {
    if (nm %in% names(mu)) {
      mu_full[[nm]] <- mu[[nm]]
    } else if (nm %in% names(fit$theta)) {
      mu_full[[nm]] <- fit$theta[[nm]]
    } else if (nm %in% omega_info$colName) {
      mu_full[[nm]] <- omega_info$est[match(nm, omega_info$colName)]
    } else {
      cli::cli_abort("Cannot initialize SIR parameter {.val {nm}}.")
    }
  }

  n_theta <- length(theta_names)
  n_full <- length(full_names)
  if (nrow(proposalCov) == n_full && ncol(proposalCov) == n_full) {
    cov_full <- proposalCov[full_names, full_names, drop = FALSE]
  } else if (nrow(proposalCov) == n_theta && ncol(proposalCov) == n_theta) {
    theta_cov <- proposalCov[theta_names, theta_names, drop = FALSE]
    if (thetaInflation != 1) {
      theta_sd <- sqrt(diag(theta_cov)) * sqrt(thetaInflation)
      theta_cor <- cov2cor(theta_cov)
      theta_cov <- outer(theta_sd, theta_sd) * theta_cor
      dimnames(theta_cov) <- list(theta_names, theta_names)
    }
    omega_cov <- diag(omega_info$se^2 * omegaInflation, nrow = nrow(omega_info))
    sigma_cov <- if (length(sigma_names) > 0L) {
      diag(sigma_info$se^2 * sigmaInflation, nrow = length(sigma_names))
    } else {
      matrix(numeric(0), nrow = 0L, ncol = 0L)
    }

    cov_full <- matrix(
      0,
      nrow = n_full,
      ncol = n_full,
      dimnames = list(full_names, full_names)
    )
    cov_full[theta_names, theta_names] <- theta_cov
    cov_full[omega_info$colName, omega_info$colName] <- omega_cov
    if (length(sigma_names) > 0L) {
      cov_full[sigma_names, sigma_names] <- sigma_cov
    }
  } else {
    cli::cli_abort(
      "{.arg proposalCov} must match either theta parameters or the full SIR parameter vector."
    )
  }

  cov_full <- .sirCapCovCorrelation(cov_full, capCorrelation = capCorrelation)
  cov_full <- .sirEnsurePosDef(cov_full)

  list(
    mu = mu_full,
    covMat = cov_full,
    thetaNames = theta_names,
    omegaInfo = omega_info,
    sigmaInfo = sigma_info,
    paramNames = full_names
  )
}

.sirParamBounds <- function(fit, paramNames, thetaNames, sigmaInfo) {
  ini_df <- fit$iniDf
  lower <- setNames(rep(-Inf, length(paramNames)), paramNames)
  upper <- setNames(rep(Inf, length(paramNames)), paramNames)
  bounded <- c(
    thetaNames,
    if (!is.null(sigmaInfo)) sigmaInfo$colName else character(0)
  )
  bounded <- intersect(bounded, ini_df$name)
  if (length(bounded) > 0L) {
    idx <- match(bounded, ini_df$name)
    lo <- ini_df$lower[idx]
    hi <- ini_df$upper[idx]
    lower[bounded] <- ifelse(is.na(lo), -Inf, lo)
    upper[bounded] <- ifelse(is.na(hi), Inf, hi)
  }
  list(lower = lower, upper = upper)
}

.sirOmegaPd <- function(omega_info, vals, base_omega) {
  omega_mat <- .sirReconstructOmega(omega_info, vals, base_omega)
  tryCatch(
    {
      chol(omega_mat)
      TRUE
    },
    error = function(e) FALSE
  )
}

.sirSampleFullProposal <- function(
  mu,
  covMat,
  n,
  lower,
  upper,
  omegaInfo,
  baseOmega,
  thetaNames,
  sigmaNames = character(0),
  boxcoxState = NULL,
  maxAttemptFactor = 10L
) {
  p <- length(mu)
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, min.len = 1L)
  checkmate::assertMatrix(covMat, mode = "numeric", nrows = p, ncols = p)
  checkmate::assertCount(n, positive = TRUE)

  n <- as.integer(n)
  max_attempts <- maxAttemptFactor * n
  param_names <- names(mu)
  sample_mat <- matrix(
    NA_real_,
    nrow = n,
    ncol = p,
    dimnames = list(NULL, param_names)
  )
  original_mat <- sample_mat

  n_filled <- 0L
  n_attempted <- 0L
  inverse_rejected <- 0L
  theta_rejected <- 0L
  omega_rejected <- 0L
  sigma_rejected <- 0L

  while (n_filled < n && n_attempted < max_attempts) {
    batch_n <- min(n - n_filled, max_attempts - n_attempted)
    draws <- mvtnorm::rmvnorm(batch_n, mean = mu, sigma = covMat)
    colnames(draws) <- param_names
    n_attempted <- n_attempted + batch_n

    original <- .sirBcInverseMatrix(draws, boxcoxState)
    has_na <- apply(is.na(original), 1L, any)
    theta_out <- rep(FALSE, nrow(original))
    if (length(thetaNames) > 0L) {
      theta_mat <- original[, thetaNames, drop = FALSE]
      theta_out <- apply(
        sweep(theta_mat, 2L, lower[thetaNames], "<") |
          sweep(theta_mat, 2L, upper[thetaNames], ">"),
        1L,
        any
      )
    }
    sigma_out <- rep(FALSE, nrow(original))
    if (length(sigmaNames) > 0L) {
      sigma_mat <- original[, sigmaNames, drop = FALSE]
      sigma_out <- apply(
        sweep(sigma_mat, 2L, lower[sigmaNames], "<") |
          sweep(sigma_mat, 2L, upper[sigmaNames], ">"),
        1L,
        any
      )
    }
    in_bounds <- !has_na & !theta_out & !sigma_out
    omega_pd <- vapply(
      seq_len(nrow(original)),
      function(i) {
        if (has_na[i] || !in_bounds[i]) {
          return(FALSE)
        }
        .sirOmegaPd(omegaInfo, original[i, ], baseOmega)
      },
      logical(1L)
    )

    ok <- !has_na & in_bounds & omega_pd
    inverse_rejected <- inverse_rejected + sum(has_na)
    theta_rejected <- theta_rejected + sum(!has_na & theta_out)
    omega_rejected <- omega_rejected + sum(!has_na & in_bounds & !omega_pd)
    sigma_rejected <- sigma_rejected + sum(!has_na & sigma_out)

    n_take <- min(sum(ok), n - n_filled)
    if (n_take > 0L) {
      idx <- seq(n_filled + 1L, n_filled + n_take)
      keep <- which(ok)[seq_len(n_take)]
      sample_mat[idx, ] <- draws[keep, , drop = FALSE]
      original_mat[idx, ] <- original[keep, , drop = FALSE]
      n_filled <- n_filled + n_take
    }
  }

  if (n_filled < n) {
    cli::cli_warn(c(
      "Only {n_filled} of {n} requested full SIR samples were valid after {max_attempts} draw attempts.",
      "i" = "Consider widening the proposal or parameter bounds."
    ))
    sample_mat <- sample_mat[seq_len(n_filled), , drop = FALSE]
    original_mat <- original_mat[seq_len(n_filled), , drop = FALSE]
  }

  list(
    samples = original_mat,
    samplesForPdf = sample_mat,
    nAttempted = n_attempted,
    inverseRejected = inverse_rejected,
    thetaRejected = theta_rejected,
    omegaRejected = omega_rejected,
    sigmaRejected = sigma_rejected
  )
}

#' Run one SIR iteration
#'
#' Orchestrates sampling, OFV evaluation, weight computation, resampling, and
#' proposal update for a single iteration of the SIR algorithm.  In the first
#' iteration `boxcoxState` should be `NULL`; thereafter pass the `boxcoxState`
#' element returned by the previous call.
#'
#' @param fit nlmixr2 fit object.
#' @param mu Named numeric vector of proposal means in the **original**
#'   parameter scale.  Iteration 1 may pass theta means only; missing
#'   omega/sigma entries are initialized from the fit.
#' @param proposalCov Covariance matrix in the **sampling** scale.  Iteration 1
#'   may pass theta covariance only; later iterations pass the full empirical
#'   proposal covariance returned by the previous call.
#' @param nSamples Positive integer. Number of parameter vectors to draw.
#' @param nResample Positive integer. Target number of vectors to resample.
#' @param iterNum Positive integer. Iteration index (used for file naming).
#' @param requestedSamples Positive integer. User-requested sample count before
#'   PsN-style attempted-sample compensation. Defaults to `nSamples`.
#' @param capResampling Numeric >= 1. Max appearances per sample. Default `1`
#'   samples without replacement; values greater than `1` allow limited
#'   replacement up to the cap.
#' @param recenter Logical. If `TRUE` and any sampled dOFV < 0, shift `mu` to
#'   the best vector. Default `TRUE`.
#' @param boxcox Logical. Apply per-column Box-Cox before updating the
#'   proposal. Default `TRUE`.
#' @param directory Optional path. If supplied, writes
#'   `raw_results_sir_iteration{N}.csv` there.
#' @param workers Passed to `.withWorkerPlan()` for parallel OFV evaluation.
#' @param boxcoxState `NULL` or a data frame (`param`, `lambda`, `delta`) from
#'   a previous call.  When non-`NULL` the proposal is assumed to be in
#'   Box-Cox space and samples are back-transformed before OFV evaluation.
#' @param thetaInflation,omegaInflation,sigmaInflation Non-negative variance
#'   multipliers for theta, fallback omega, and fallback sigma proposal blocks.
#' @param capCorrelation Maximum absolute proposal correlation after covariance
#'   updates.
#' @param omegaFallback Method for omega fallback uncertainty. Currently
#'   `"wishart"`.
#' @param sigmaFallbackRse Percent relative standard error for sigma-like
#'   fallback uncertainty when no SE is available.
#' @param omegaDf Optional degrees of freedom for omega Wishart-style fallback;
#'   defaults to `nsub - 1`.
#' @param isLastIteration Logical. If `TRUE`, update the proposal covariance on
#'   the original scale even when `boxcox = TRUE`.
#' @return Named list with elements `resampledMat`, `newMu`, `newCov`,
#'   `iterSummary`, `boxcoxState`, `rawResults`.
#' @noRd
sirRunIteration <- function(
  fit,
  mu,
  proposalCov,
  nSamples,
  nResample,
  iterNum,
  capResampling = 1,
  recenter = TRUE,
  boxcox = TRUE,
  directory = NULL,
  workers = NULL,
  boxcoxState = NULL,
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8,
  omegaFallback = c("wishart"),
  sigmaFallbackRse = 30,
  omegaDf = NULL,
  requestedSamples = nSamples,
  isLastIteration = FALSE
) {
  omegaFallback <- match.arg(omegaFallback)

  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertNumeric(mu, finite = TRUE, any.missing = FALSE, min.len = 1L)
  checkmate::assertMatrix(proposalCov, mode = "numeric")
  checkmate::assertCount(nSamples, positive = TRUE)
  checkmate::assertCount(nResample, positive = TRUE)
  checkmate::assertCount(iterNum, positive = TRUE)
  checkmate::assertCount(requestedSamples, positive = TRUE)
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)
  checkmate::assertFlag(recenter)
  checkmate::assertFlag(boxcox)
  checkmate::assertFlag(isLastIteration)
  checkmate::assertNumber(thetaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(omegaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(sigmaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)
  checkmate::assertNumber(sigmaFallbackRse, lower = 0, finite = TRUE)
  if (!is.null(directory)) {
    checkmate::assertString(directory)
  }

  proposal <- .sirInitialProposal(
    fit = fit,
    mu = mu,
    proposalCov = proposalCov,
    thetaInflation = thetaInflation,
    omegaInflation = omegaInflation,
    sigmaInflation = sigmaInflation,
    capCorrelation = capCorrelation,
    omegaFallback = omegaFallback,
    sigmaFallbackRse = sigmaFallbackRse,
    omegaDf = omegaDf
  )

  theta_names <- proposal$thetaNames
  omega_info <- proposal$omegaInfo
  sigma_info <- proposal$sigmaInfo
  sigma_names <- if (is.null(sigma_info)) character(0) else sigma_info$colName
  param_names <- proposal$paramNames
  bounds <- .sirParamBounds(fit, param_names, theta_names, sigma_info)

  # ---- 1-3. Sample full proposal vectors and reject invalid draws ----
  bc_mu <- .sirBcTransformMu(proposal$mu, boxcoxState)
  sampled <- .sirSampleFullProposal(
    mu = bc_mu,
    covMat = proposal$covMat,
    n = nSamples,
    lower = bounds$lower,
    upper = bounds$upper,
    omegaInfo = omega_info,
    baseOmega = fit$omega,
    thetaNames = theta_names,
    sigmaNames = sigma_names,
    boxcoxState = boxcoxState
  )

  param_mat <- sampled$samples
  samples_for_weights <- sampled$samplesForPdf
  n_collected <- nrow(param_mat)
  if (n_collected == 0L) {
    cli::cli_abort(c(
      "No valid SIR samples were collected.",
      "i" = "Check proposal covariance, bounds, and omega positive-definiteness diagnostics."
    ))
  }

  # ---- 4. Evaluate OFV ----
  ofv_vals <- sirEvalOFV(fit, param_mat, workers = workers)
  dofv <- ofv_vals - fit$objf

  # ---- 5. Handle failures ----
  n_failed <- sum(is.na(dofv))
  n_success <- n_collected - n_failed
  n_resample_adj <- nResample
  if (n_success == 0L) {
    cli::cli_abort(c(
      "All SIR OFV evaluations failed.",
      "i" = "No valid importance weights can be computed."
    ))
  }
  if (abs(n_success - requestedSamples) / requestedSamples > 0.05) {
    n_resample_adj <- max(
      1L,
      as.integer(floor(nResample * n_success / requestedSamples))
    )
    cli::cli_warn(c(
      "{n_success}/{requestedSamples} requested SIR samples had usable OFV evaluations.",
      "i" = "Adjusting nResample to {n_resample_adj}."
    ))
  }

  # ---- 6. Compute weights in the same full proposal scale used for sampling ----
  weights <- sirCalcWeights(
    samples_for_weights,
    bc_mu,
    proposal$covMat,
    dOFV = dofv
  )

  # ---- 7. Resample ----
  resampled <- sirResample(
    param_mat,
    weights,
    m = n_resample_adj,
    capResampling = capResampling
  )

  # ---- 8. Recenter ----
  new_mu <- proposal$mu
  if (recenter) {
    valid_dofv <- ifelse(is.na(dofv), Inf, dofv)
    best_idx <- which.min(valid_dofv)
    if (isTRUE(valid_dofv[best_idx] < 0)) {
      new_mu <- param_mat[best_idx, ]
      cli::cli_inform(
        "  Iter {iterNum}: recentered mu (dOFV = {round(dofv[best_idx], 4)})."
      )
    }
  } else if (any(dofv < 0, na.rm = TRUE)) {
    cli::cli_warn(
      "At least one SIR sample had negative dOFV, but {.arg recenter} is FALSE."
    )
  }

  # ---- 9. Update full proposal in original scale ----
  updated <- sirUpdateProposal(
    resampled$samples[, param_names, drop = FALSE],
    boxcox = boxcox && !isLastIteration,
    capCorrelation = capCorrelation
  )
  new_cov <- updated$covMat
  new_bc_state <- updated$boxcoxParams # NULL when boxcox = FALSE

  # ---- 10. Raw results data frame ----
  raw_df <- .sirBuildRawResults(
    paramMat = param_mat,
    weights = weights,
    dOFV = dofv,
    resampled = resampled,
    mu = proposal$mu,
    capResampling = capResampling
  )

  if (!is.null(directory)) {
    csv_path <- file.path(
      directory,
      sprintf("raw_results_sir_iteration%d.csv", iterNum)
    )
    utils::write.csv(raw_df, csv_path, row.names = FALSE)
  }

  # ---- Summary ----
  iter_summary <- data.frame(
    iter = iterNum,
    nSamples = requestedSamples,
    nAttempted = nSamples,
    nDrawAttempts = sampled$nAttempted,
    nCollected = n_collected,
    nSuccessful = n_success,
    nFailed = n_failed,
    nResample = nResample,
    nResampled = nrow(resampled$samples),
    minDOFV = if (all(is.na(dofv))) NA_real_ else min(dofv, na.rm = TRUE),
    meanDOFV = if (all(is.na(dofv))) NA_real_ else mean(dofv, na.rm = TRUE),
    nNegativeDOFV = sum(dofv < 0, na.rm = TRUE),
    thetaRejected = sampled$thetaRejected,
    omegaRejected = sampled$omegaRejected,
    sigmaRejected = sampled$sigmaRejected,
    inverseRejected = sampled$inverseRejected
  )

  list(
    resampledMat = resampled$samples,
    newMu = new_mu,
    newCov = new_cov,
    iterSummary = iter_summary,
    boxcoxState = new_bc_state,
    rawResults = raw_df
  )
}

.sirAdjustedAttemptedSamples <- function(
  requestedSamples,
  previousAttempted = NULL,
  previousSuccessful = NULL
) {
  checkmate::assertCount(requestedSamples, positive = TRUE)
  if (is.null(previousAttempted) || is.null(previousSuccessful)) {
    return(as.integer(requestedSamples))
  }
  checkmate::assertCount(previousAttempted, positive = TRUE)
  checkmate::assertCount(previousSuccessful, positive = TRUE)

  if (previousSuccessful < 0.95 * previousAttempted) {
    return(max(
      as.integer(requestedSamples),
      as.integer(ceiling(
        requestedSamples * previousAttempted / previousSuccessful
      ))
    ))
  }
  as.integer(requestedSamples)
}

.sirResolveOutputDir <- function(directory, fitName) {
  if (is.null(directory)) {
    fit_name <- gsub("[^A-Za-z0-9_.]", "_", fitName)
    pattern <- paste0("^", fit_name, "_sir_[0-9]+$")
    existing <- sort(grep(
      pattern,
      list.dirs(".", full.names = FALSE, recursive = FALSE),
      value = TRUE
    ))
    directory <- paste0(fit_name, "_sir_", length(existing) + 1L)
  }

  directory <- normalizePath(directory, mustWork = FALSE)
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  directory
}

.sirOriginalEstimates <- function(fit, paramNames) {
  vals <- setNames(rep(NA_real_, length(paramNames)), paramNames)
  theta_names <- intersect(paramNames, names(fit$theta))
  vals[theta_names] <- fit$theta[theta_names]

  omega_info <- .sirOmegaProposalInfo(fit)
  if (nrow(omega_info) > 0L) {
    omega_names <- intersect(paramNames, omega_info$colName)
    vals[omega_names] <- omega_info$est[match(omega_names, omega_info$colName)]
  }
  vals
}

#' Summarize final SIR resampled parameter vectors
#'
#' `sirSummary()` computes empirical summary statistics from final SIR
#' resampled parameter vectors. It reports the reference estimate from the
#' input fit, empirical standard deviation and RSE, and percentile intervals
#' for each sampled parameter. The empirical covariance and correlation
#' matrices are attached as attributes.
#'
#' @param resampledMat Numeric matrix of final resampled parameter vectors, one
#'   row per vector and one column per parameter.
#' @param fit An nlmixr2 fit object used to obtain reference parameter
#'   estimates.
#' @return A data frame with columns `param`, `estimate`, `sd`, `rse`, `p2.5`,
#'   `p5`, `p25`, `p50`, `p75`, `p95`, and `p97.5`. Attributes `covMatrix`
#'   and `corMatrix` contain empirical covariance and correlation matrices.
#' @export
sirSummary <- function(resampledMat, fit) {
  if (is.data.frame(resampledMat)) {
    resampledMat <- as.matrix(resampledMat)
  }
  checkmate::assertMatrix(resampledMat, mode = "numeric", min.rows = 2L)
  checkmate::assertClass(fit, "nlmixr2FitCore")
  if (is.null(colnames(resampledMat))) {
    cli::cli_abort("{.arg resampledMat} must have parameter column names.")
  }

  param_names <- colnames(resampledMat)
  estimate <- .sirOriginalEstimates(fit, param_names)
  probs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  qmat <- t(apply(
    resampledMat,
    2L,
    stats::quantile,
    probs = probs,
    na.rm = TRUE,
    names = FALSE
  ))
  colnames(qmat) <- c("p2.5", "p5", "p25", "p50", "p75", "p95", "p97.5")
  sd_vals <- apply(resampledMat, 2L, stats::sd, na.rm = TRUE)
  rse <- ifelse(
    is.finite(estimate) & estimate != 0,
    sd_vals / abs(estimate) * 100,
    NA_real_
  )
  data.frame(
    param = param_names,
    estimate = unname(estimate),
    sd = unname(sd_vals),
    rse = unname(rse),
    qmat,
    row.names = NULL,
    check.names = FALSE
  ) |>
    structure(
      covMatrix = stats::cov(resampledMat),
      corMatrix = stats::cor(resampledMat)
    )
}

.sirSummarizeResamples <- function(resampledMat, fit) {
  sirSummary(resampledMat, fit)
}

.sirWriteIterationSummary <- function(iterSummary, directory) {
  out <- iterSummary
  out$requested_sample_resample_ratio <- out$nSamples / out$nResample
  out$actual_sample_resample_ratio <- out$nSuccessful / out$nResampled
  utils::write.csv(
    out,
    file.path(directory, "summary_iterations.csv"),
    row.names = FALSE
  )
  invisible(out)
}

.sirWriteRejectionSummary <- function(iterSummary, directory) {
  s <- iterSummary[nrow(iterSummary), , drop = FALSE]
  lines <- c(
    "SIR sample rejection summary",
    sprintf("Iteration: %s", s$iter),
    sprintf("Requested samples: %s", s$nSamples),
    sprintf("Attempted samples: %s", s$nAttempted),
    sprintf("Raw draw attempts: %s", s$nDrawAttempts),
    sprintf("Collected samples: %s", s$nCollected),
    sprintf("Successful OFV evaluations: %s", s$nSuccessful),
    sprintf("Failed OFV evaluations: %s", s$nFailed),
    sprintf("Inverse Box-Cox rejections: %s", s$inverseRejected),
    sprintf("Theta bound rejections: %s", s$thetaRejected),
    sprintf("Omega positive-definiteness rejections: %s", s$omegaRejected),
    sprintf("Sigma bound rejections: %s", s$sigmaRejected)
  )
  writeLines(lines, file.path(directory, "sample_rejection_summary.txt"))
  invisible(lines)
}

.sirSaveState <- function(directory, state) {
  saveRDS(state, file.path(directory, "sir_state.rds"))
  invisible(state)
}

#' Run sampling importance resampling for an nlmixr2 fit
#'
#' `runSIR()` runs one or more sampling importance resampling iterations using
#' an nlmixr2 fit as the reference model. It follows the PsN SIR workflow where
#' practical for nlmixr2: sample parameter vectors from a proposal covariance,
#' evaluate them with population parameters fixed, compute importance weights,
#' resample, and use the empirical resampled covariance as the next proposal.
#'
#' @param fit An nlmixr2 fit object with a covariance matrix.
#' @param nSamples Integer vector. Requested number of samples per iteration.
#' @param nResample Integer vector. Requested number of resamples per
#'   iteration. Must have the same length as `nSamples`.
#' @param thetaInflation,omegaInflation,sigmaInflation Non-negative variance
#'   multipliers used for the first iteration proposal.
#' @param capCorrelation Numeric in `[0, 1]`. Maximum absolute proposal
#'   correlation after covariance construction and updates.
#' @param capResampling Numeric greater than or equal to one. `1` samples
#'   without replacement; larger values allow limited replacement.
#' @param recenter Logical. If `TRUE`, recenter the next proposal on the best
#'   sampled vector when any sample has negative dOFV.
#' @param boxcox Logical. If `TRUE`, use Box-Cox transformed empirical
#'   covariances for non-final iterations.
#' @param workers `NULL`, `"auto"`, `1`, or a positive integer. Controls
#'   parallel OFV evaluation through `future`. `NULL` leaves the current
#'   `future::plan()` unchanged, `1` forces sequential execution, a positive
#'   integer temporarily uses a multisession plan, and `"auto"` uses
#'   `future::availableCores(omit = 1L)`. SIR iterations remain sequential
#'   because each iteration updates the proposal for the next iteration.
#' @param directory Output directory. If `NULL`, a numbered
#'   `<fitName>_sir_<N>` directory is created.
#' @param recover Logical. If `TRUE` and `directory` contains `sir_state.rds`,
#'   resume from the last completed iteration when possible.
#' @param addIterations Logical. If `TRUE`, append the supplied `nSamples` and
#'   `nResample` schedule after an existing completed state in `directory`.
#' @param omegaFallback Method for omega fallback uncertainty. Currently only
#'   `"wishart"`.
#' @param sigmaFallbackRse Percent relative standard error for sigma-like
#'   fallback uncertainty when no SE is available.
#' @param omegaDf Optional degrees of freedom for the omega Wishart-style
#'   fallback; defaults to `nsub - 1`.
#' @param fitName Name used when creating an automatic output directory.
#' @param ... Reserved for future PsN-compatible inputs.
#' @return A data frame of final SIR summary statistics with class
#'   `c("nlmixr2SIR", "data.frame")`. Attributes include
#'   `iterationSummary`, `iterations`, `resampledMat`, `covMatrix`,
#'   `corMatrix`, and `outputDir`.
#' @examples
#' \dontrun{
#' sir <- runSIR(
#'   fit,
#'   nSamples = c(1000, 1000, 1000),
#'   nResample = c(200, 400, 500),
#'   workers = 4
#' )
#' }
#' @export
runSIR <- function(
  fit,
  nSamples = c(1000, 1000, 1000, 2000, 2000),
  nResample = c(200, 400, 500, 1000, 1000),
  thetaInflation = 1,
  omegaInflation = 1,
  sigmaInflation = 1,
  capCorrelation = 0.8,
  capResampling = 1,
  recenter = TRUE,
  boxcox = TRUE,
  workers = NULL,
  directory = NULL,
  recover = TRUE,
  addIterations = FALSE,
  omegaFallback = c("wishart"),
  sigmaFallbackRse = 30,
  omegaDf = NULL,
  fitName = as.character(substitute(fit)),
  ...
) {
  omegaFallback <- match.arg(omegaFallback)
  dots <- list(...)
  if (length(dots) > 0L) {
    cli::cli_abort("Unsupported SIR argument(s): {names(dots)}.")
  }

  checkmate::assertClass(fit, "nlmixr2FitCore")
  checkmate::assertIntegerish(
    nSamples,
    lower = 1,
    any.missing = FALSE,
    min.len = 1L
  )
  checkmate::assertIntegerish(
    nResample,
    lower = 1,
    any.missing = FALSE,
    len = length(nSamples)
  )
  checkmate::assertNumber(thetaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(omegaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(sigmaInflation, lower = 0, finite = TRUE)
  checkmate::assertNumber(capCorrelation, lower = 0, upper = 1, finite = TRUE)
  checkmate::assertNumber(capResampling, lower = 1, finite = TRUE)
  checkmate::assertFlag(recenter)
  checkmate::assertFlag(boxcox)
  checkmate::assertFlag(recover)
  checkmate::assertFlag(addIterations)
  checkmate::assertNumber(sigmaFallbackRse, lower = 0, finite = TRUE)
  .validateWorkers(workers)

  nSamples <- as.integer(nSamples)
  nResample <- as.integer(nResample)
  output_dir <- .sirResolveOutputDir(directory, fitName)
  state_file <- file.path(output_dir, "sir_state.rds")
  saved_state <- if (file.exists(state_file) && (recover || addIterations)) {
    readRDS(state_file)
  } else {
    NULL
  }

  if (!is.null(saved_state) && isTRUE(addIterations)) {
    iter_offset <- saved_state$completedIterations
    iter_index <- seq_along(nSamples)
    iter_numbers <- iter_offset + iter_index
    mu <- saved_state$nextMu
    proposal_cov <- saved_state$nextCov
    boxcox_state <- saved_state$nextBoxcoxState
    iter_results <- saved_state$iterations
    iter_summary <- saved_state$iterationSummary
    prev_attempted <- tail(iter_summary$nAttempted, 1L)
    prev_successful <- tail(iter_summary$nSuccessful, 1L)
  } else if (!is.null(saved_state) && isTRUE(recover)) {
    completed <- saved_state$completedIterations
    if (completed >= length(nSamples) && !is.null(saved_state$result)) {
      cli::cli_inform("Recovered completed SIR run from {.path {output_dir}}.")
      return(saved_state$result)
    }
    iter_numbers <- seq.int(completed + 1L, length(nSamples))
    iter_index <- iter_numbers
    mu <- saved_state$nextMu
    proposal_cov <- saved_state$nextCov
    boxcox_state <- saved_state$nextBoxcoxState
    iter_results <- saved_state$iterations
    iter_summary <- saved_state$iterationSummary
    prev_attempted <- tail(iter_summary$nAttempted, 1L)
    prev_successful <- tail(iter_summary$nSuccessful, 1L)
  } else {
    iter_numbers <- seq_along(nSamples)
    iter_index <- iter_numbers
    mu <- fit$theta[rownames(fit$cov)]
    proposal_cov <- sirGetProposalCov(
      fit,
      thetaInflation = 1,
      omegaInflation = 1,
      sigmaInflation = 1,
      capCorrelation = capCorrelation
    )
    boxcox_state <- NULL
    iter_results <- list()
    iter_summary <- data.frame()
    prev_attempted <- NULL
    prev_successful <- NULL
  }

  for (j in seq_along(iter_numbers)) {
    iter_num <- iter_numbers[[j]]
    schedule_idx <- iter_index[[j]]
    requested_samples <- nSamples[[schedule_idx]]
    attempted_samples <- .sirAdjustedAttemptedSamples(
      requested_samples,
      previousAttempted = prev_attempted,
      previousSuccessful = prev_successful
    )
    is_last <- j == length(iter_numbers)

    cli::cli_inform(
      "Running SIR iteration {iter_num}: {attempted_samples} attempted samples, {nResample[[schedule_idx]]} requested resamples."
    )
    iter_res <- sirRunIteration(
      fit = fit,
      mu = mu,
      proposalCov = proposal_cov,
      nSamples = attempted_samples,
      requestedSamples = requested_samples,
      nResample = nResample[[schedule_idx]],
      iterNum = iter_num,
      capResampling = capResampling,
      recenter = recenter,
      boxcox = boxcox,
      directory = output_dir,
      workers = workers,
      boxcoxState = boxcox_state,
      thetaInflation = thetaInflation,
      omegaInflation = omegaInflation,
      sigmaInflation = sigmaInflation,
      capCorrelation = capCorrelation,
      omegaFallback = omegaFallback,
      sigmaFallbackRse = sigmaFallbackRse,
      omegaDf = omegaDf,
      isLastIteration = is_last
    )

    iter_results[[as.character(iter_num)]] <- iter_res
    iter_summary <- rbind(iter_summary, iter_res$iterSummary)
    .sirWriteIterationSummary(iter_summary, output_dir)
    .sirWriteRejectionSummary(iter_summary, output_dir)

    mu <- iter_res$newMu
    proposal_cov <- iter_res$newCov
    boxcox_state <- iter_res$boxcoxState
    prev_attempted <- iter_res$iterSummary$nAttempted
    prev_successful <- iter_res$iterSummary$nSuccessful

    .sirSaveState(
      output_dir,
      list(
        completedIterations = iter_num,
        nextMu = mu,
        nextCov = proposal_cov,
        nextBoxcoxState = boxcox_state,
        iterations = iter_results,
        iterationSummary = iter_summary,
        result = NULL
      )
    )
  }

  final_iter <- iter_results[[length(iter_results)]]
  summary_df <- sirSummary(final_iter$resampledMat, fit)
  cov_mat <- attr(summary_df, "covMatrix")
  cor_mat <- attr(summary_df, "corMatrix")
  utils::write.csv(
    summary_df,
    file.path(output_dir, "sir_results.csv"),
    row.names = FALSE
  )

  class(summary_df) <- c("nlmixr2SIR", "data.frame")
  attr(summary_df, "iterationSummary") <- iter_summary
  attr(summary_df, "iterations") <- iter_results
  attr(summary_df, "resampledMat") <- final_iter$resampledMat
  attr(summary_df, "covMatrix") <- cov_mat
  attr(summary_df, "corMatrix") <- cor_mat
  attr(summary_df, "outputDir") <- output_dir
  attr(summary_df, "fitName") <- fitName
  attr(summary_df, "call") <- match.call()

  .sirSaveState(
    output_dir,
    list(
      completedIterations = tail(iter_summary$iter, 1L),
      nextMu = mu,
      nextCov = proposal_cov,
      nextBoxcoxState = boxcox_state,
      iterations = iter_results,
      iterationSummary = iter_summary,
      result = summary_df
    )
  )

  summary_df
}

.sirParameterPlotData <- function(x) {
  resampled <- attr(x, "resampledMat")
  if (is.null(resampled)) {
    cli::cli_abort("SIR object is missing {.field resampledMat}.")
  }
  data.frame(
    param = rep(colnames(resampled), each = nrow(resampled)),
    value = as.vector(resampled),
    stringsAsFactors = FALSE
  )
}

.sirRawResults <- function(x) {
  iterations <- attr(x, "iterations")
  if (is.null(iterations) || length(iterations) == 0L) {
    cli::cli_abort("SIR object is missing iteration raw results.")
  }
  out <- lapply(seq_along(iterations), function(i) {
    raw <- iterations[[i]]$rawResults
    raw$iter <- iterations[[i]]$iterSummary$iter
    raw
  })
  do.call(rbind, out)
}

.sirSignifDataFrame <- function(x, digits) {
  out <- as.data.frame(x)
  num_cols <- vapply(out, is.numeric, logical(1L))
  out[num_cols] <- lapply(out[num_cols], signif, digits = digits)
  out
}

#' Print an SIR result
#'
#' Prints the final parameter uncertainty summary and, when available, the
#' per-iteration SIR diagnostics stored on the `nlmixr2SIR` object returned by
#' [runSIR()].
#'
#' @param x An object returned by [runSIR()].
#' @param ... Unused.
#' @param digits Number of significant digits to print.
#' @return Invisibly returns `x`.
#' @export
print.nlmixr2SIR <- function(x, ..., digits = 3) {
  checkmate::assertDataFrame(x)
  checkmate::assertCount(digits, positive = TRUE)

  output_dir <- attr(x, "outputDir", exact = TRUE)
  cli::cli_h1("SIR Summary")
  if (!is.null(output_dir)) {
    cli::cli_alert_info("Output directory: {.path {output_dir}}")
  }

  summary_cols <- intersect(
    c("param", "estimate", "sd", "rse", "p2.5", "p50", "p97.5"),
    names(x)
  )
  cli::cli_h2("Final parameter uncertainty")
  print(
    .sirSignifDataFrame(x[, summary_cols, drop = FALSE], digits = digits),
    row.names = FALSE
  )

  iter_summary <- attr(x, "iterationSummary", exact = TRUE)
  if (!is.null(iter_summary) && nrow(iter_summary) > 0L) {
    iter_cols <- intersect(
      c(
        "iter",
        "nSamples",
        "nAttempted",
        "nSuccessful",
        "nResample",
        "nResampled",
        "nNegativeDOFV",
        "minDOFV"
      ),
      names(iter_summary)
    )
    cli::cli_h2("Iteration diagnostics")
    print(
      .sirSignifDataFrame(
        iter_summary[, iter_cols, drop = FALSE],
        digits = digits
      ),
      row.names = FALSE
    )
  }

  invisible(x)
}

#' Plot an SIR result
#'
#' Plots diagnostics for an `nlmixr2SIR` object. The default plot shows the
#' final resampled parameter distributions with reference estimates marked.
#' Additional diagnostic views show dOFV distributions or resampling
#' probabilities by iteration.
#'
#' @param x An object returned by [runSIR()].
#' @param y Unused; included for S3 compatibility.
#' @param type Plot type: `"parameters"`, `"dofv"`, or `"resampling"`.
#' @param bins Number of histogram bins for parameter and dOFV plots.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @export
#' @importFrom ggplot2 .data
plot.nlmixr2SIR <- function(
  x,
  y,
  ...,
  type = c("parameters", "dofv", "resampling"),
  bins = 30
) {
  type <- match.arg(type)
  checkmate::assertCount(bins, positive = TRUE)

  if (type == "parameters") {
    plot_df <- .sirParameterPlotData(x)
    estimate_df <- data.frame(
      param = x$param,
      estimate = x$estimate,
      stringsAsFactors = FALSE
    )
    return(
      ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$value)) +
        ggplot2::geom_histogram(
          bins = bins,
          fill = "#6BAED6",
          color = "white"
        ) +
        ggplot2::geom_vline(
          data = estimate_df,
          ggplot2::aes(xintercept = .data$estimate),
          color = "#D94801",
          linewidth = 0.6
        ) +
        ggplot2::facet_wrap(stats::as.formula("~ param"), scales = "free") +
        ggplot2::labs(x = "Parameter value", y = "Resampled vectors") +
        ggplot2::theme_bw()
    )
  }

  raw_df <- .sirRawResults(x)
  raw_df <- raw_df[raw_df$sample_id > 0, , drop = FALSE]

  if (type == "dofv") {
    return(
      ggplot2::ggplot(raw_df, ggplot2::aes(x = .data$dOFV)) +
        ggplot2::geom_histogram(
          bins = bins,
          fill = "#74C476",
          color = "white"
        ) +
        ggplot2::facet_wrap(stats::as.formula("~ iter"), scales = "free_y") +
        ggplot2::labs(x = "dOFV", y = "Sampled vectors") +
        ggplot2::theme_bw()
    )
  }

  ggplot2::ggplot(
    raw_df,
    ggplot2::aes(x = .data$sample_id, y = .data$probability_resample)
  ) +
    ggplot2::geom_col(fill = "#9E9AC8") +
    ggplot2::facet_wrap(stats::as.formula("~ iter"), scales = "free_x") +
    ggplot2::labs(x = "Sample id", y = "Probability resample") +
    ggplot2::theme_bw()
}
