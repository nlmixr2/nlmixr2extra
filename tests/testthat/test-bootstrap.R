skip_on_cran()

# =============================================================================
# Unit tests for sampling() — no model fitting required
# =============================================================================

test_that("sampling: returns a data frame with the same columns", {
  d   <- data.frame(ID = rep(1:6, each = 3), y = rnorm(18), TIME = rep(0:2, 6))
  out <- nlmixr2extra:::sampling(d, nsamp = 6L)
  expect_s3_class(out, "data.frame")
  expect_true(all(names(d) %in% names(out)))
})

test_that("sampling: nsamp controls the number of unique IDs returned", {
  d   <- data.frame(ID = rep(1:10, each = 4), y = rnorm(40))
  out <- nlmixr2extra:::sampling(d, nsamp = 4L)
  expect_equal(length(unique(out$ID)), 4L)
})

test_that("sampling: default nsamp equals number of unique subjects in data", {
  d   <- data.frame(ID = rep(1:5, each = 3), y = rnorm(15))
  out <- nlmixr2extra:::sampling(d, uid_colname = "ID")
  expect_equal(length(unique(out$ID)), 5L)
})

test_that("sampling: IDs are renumbered 1..nsamp", {
  d   <- data.frame(ID = rep(1:5, each = 3), y = rnorm(15))
  out <- nlmixr2extra:::sampling(d, nsamp = 5L)
  expect_setequal(unique(out$ID), 1:5)
})

test_that("sampling: total rows equal nsamp times rows-per-subject", {
  d   <- data.frame(ID = rep(1:8, each = 4), y = rnorm(32))
  out <- nlmixr2extra:::sampling(d, nsamp = 8L)
  expect_equal(nrow(out), 8L * 4L)
})

test_that("sampling: nsamp=2 is the minimum accepted value (lower=2 constraint)", {
  d <- data.frame(ID = rep(1:5, each = 2), y = rnorm(10))
  expect_no_error(nlmixr2extra:::sampling(d, nsamp = 2L))
})

test_that("sampling: nsamp=1 throws an error", {
  d <- data.frame(ID = rep(1:5, each = 2), y = rnorm(10))
  expect_error(nlmixr2extra:::sampling(d, nsamp = 1L))
})

test_that("sampling: uid_colname parameter is respected", {
  d   <- data.frame(SUBJ = rep(11:14, each = 3), y = rnorm(12))
  out <- nlmixr2extra:::sampling(d, nsamp = 4L, uid_colname = "SUBJ")
  expect_equal(length(unique(out$SUBJ)), 4L)
})

test_that("sampling: pvalues concentrating weight on one subject runs without error", {
  set.seed(123)
  d     <- data.frame(ID = rep(1:5, each = 3), y = rnorm(15))
  pvals <- c(1, 0, 0, 0, 0)
  expect_no_error(nlmixr2extra:::sampling(d, nsamp = 5L, pvalues = pvals))
})

test_that("sampling: result only contains rows from sampled subjects", {
  d   <- data.frame(ID = rep(1:5, each = 3), y = seq_len(15))
  out <- nlmixr2extra:::sampling(d, nsamp = 5L)
  # y values must come from valid rows of d
  expect_true(all(out$y %in% d$y))
})

test_that("sampling: stratified run completes without error", {
  set.seed(42)
  d <- data.frame(
    ID    = 1:20,
    group = c(rep("A", 12), rep("B", 8)),
    y     = rnorm(20)
  )
  expect_no_error(
    nlmixr2extra:::sampling(d, nsamp = 10L, performStrat = TRUE, stratVar = "group")
  )
})

test_that("sampling: performStrat=TRUE with missing stratVar throws an error", {
  d <- data.frame(ID = 1:10, y = rnorm(10))
  expect_error(nlmixr2extra:::sampling(d, nsamp = 10L, performStrat = TRUE))
})

# Subject-level (not row-level) sampling with replacement,
# variable rows per subject
# ---------------------------------------------------------------------
# Build a dataset where subjects have *different* row counts so any
# row-level sampling would give a different row total than subject-level.
#   Subject 1: 2 rows  (y = 1, 2)
#   Subject 2: 5 rows  (y = 3, 4, 5, 6, 7)
#   Subject 3: 3 rows  (y = 8, 9, 10)

test_that("sampling: samples by subject, not by row (unit is subject)", {
  d <- data.frame(
    ID = c(rep(1L, 2L), rep(2L, 5L), rep(3L, 3L)),
    y  = 1:10
  )
  # Pin subject 2 (5 rows) with full weight so we know exactly what is drawn.
  out <- nlmixr2extra:::sampling(d, nsamp = 3L, pvalues = c(0, 1, 0))
  # All three draws are subject 2 → 3 × 5 = 15 rows
  expect_equal(nrow(out), 15L)
  # Every y value must belong to subject 2 (original y 3..7)
  expect_true(all(out$y %in% 3:7))
})

test_that("sampling: subject sampled multiple times yields multiple copies of all its rows", {
  d <- data.frame(
    ID   = c(rep(1L, 3L), rep(2L, 3L)),
    TIME = rep(c(0, 1, 2), 2L),
    y    = c(10, 11, 12, 20, 21, 22)
  )
  # Force subject 1 (y = 10, 11, 12) to always be drawn.
  out <- nlmixr2extra:::sampling(d, nsamp = 2L, pvalues = c(1, 0))
  # 2 draws of subject 1 × 3 rows = 6 rows
  expect_equal(nrow(out), 6L)
  expect_true(all(out$y %in% c(10, 11, 12)))
  # Two distinct new IDs are assigned
  expect_equal(length(unique(out$ID)), 2L)
})

test_that("sampling: with-replacement allows the same subject to appear more than once", {
  d <- data.frame(
    ID = c(rep(1L, 4L), rep(2L, 4L), rep(3L, 4L)),
    y  = 1:12
  )
  # Give all weight to subject 1; nsamp=3 means three copies of subject 1.
  out <- nlmixr2extra:::sampling(d, nsamp = 3L, pvalues = c(1, 0, 0))
  # 3 copies × 4 rows = 12 rows
  expect_equal(nrow(out), 12L)
  # All y values come from subject 1 (original y 1..4)
  expect_true(all(out$y %in% 1:4))
})

test_that("sampling: variable row counts — subject with 1 row sampled 3× gives 3 rows", {
  # Subjects 1, 2, 3 have 1, 3, 5 rows respectively.
  d <- data.frame(
    ID = c(1L, rep(2L, 3L), rep(3L, 5L)),
    y  = 1:9
  )
  # Pin all weight on subject 1 (1 row): 3 draws × 1 row = 3 rows total.
  out <- nlmixr2extra:::sampling(d, nsamp = 3L, pvalues = c(1, 0, 0))
  expect_equal(nrow(out), 3L)
  expect_true(all(out$y %in% 1L))
})

test_that("sampling: variable row counts — subject with 5 rows sampled 3× gives 15 rows", {
  # Subjects 1, 2, 3 have 1, 3, 5 rows respectively.
  d <- data.frame(
    ID = c(1L, rep(2L, 3L), rep(3L, 5L)),
    y  = 1:9
  )
  # Pin all weight on subject 3 (5 rows): 3 draws × 5 rows = 15 rows total.
  out <- nlmixr2extra:::sampling(d, nsamp = 3L, pvalues = c(0, 0, 1))
  expect_equal(nrow(out), 15L)
  expect_true(all(out$y %in% 5:9))
})

test_that("sampling: each new ID in output maps to exactly one subject's full row block", {
  d <- data.frame(
    ID   = c(rep(1L, 2L), rep(2L, 4L)),
    TIME = c(1, 2, 1, 2, 3, 4),
    y    = c(10, 11, 20, 21, 22, 23)
  )
  # Draw subject 2 twice, subject 1 once (nsamp=3)
  out <- nlmixr2extra:::sampling(d, nsamp = 3L, pvalues = c(0, 1))
  # 3 draws of subject 2 × 4 rows = 12 rows; 3 unique new IDs
  expect_equal(nrow(out), 12L)
  expect_equal(length(unique(out$ID)), 3L)
  # Every new-ID group has exactly 4 rows (subject 2's row count)
  counts <- table(out$ID)
  expect_true(all(counts == 4L))
})

# =============================================================================
# Integration tests (require model fitting)
# =============================================================================

withr::with_tempdir({

  test_that("sampling should return different datasets at each call", {
    a <- digest::digest(nlmixr2extra:::sampling(nlmixr2data::theo_sd))
    b <- digest::digest(nlmixr2extra:::sampling(nlmixr2data::theo_sd))
    expect_false(isTRUE(all.equal(a, b)))
  })

  test_that("resuming the fit should not return the same datasets as before", {
    skip_on_cran()
    one.cmt <- function() {
      ini({
        tka <- 0.45 ; label("Log Ka")
        tcl <- 1 ; label("Log Cl")
        tv <- 3.45 ; label("log V")
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
    skip_if_not_installed("rxode2")
    suppressMessages(suppressWarnings(
      fit <-
        nlmixr(
          one.cmt,
          nlmixr2data::theo_sd,
          est = "focei",
          control = list(print = 0, eval.max = 10),
          table = list(npde = TRUE, cwres = TRUE)
        )
    ))

    fit1 <- suppressMessages(nlmixr2extra:::runBootstrap(fit, nboot = 2, restart = TRUE))
    fit2 <- suppressMessages(nlmixr2extra:::runBootstrap(fit, nboot = 4, restart = FALSE))

    output_dir <- fit1$outputDir   # absolute path set by runBootstrap

    fnameBootDataPattern <- paste0("boot_data", "_", "[0-9]+", ".rds")

    files <- list.files(output_dir, pattern = fnameBootDataPattern, full.names = TRUE)

    fitdata <- lapply(files, readRDS)

    a <- digest::digest(fitdata[[1]])
    b <- digest::digest(fitdata[[3]])
    expect_false(isTRUE(all.equal(a, b)))

    a <- digest::digest(fitdata[[2]])
    b <- digest::digest(fitdata[[4]])
    expect_false(isTRUE(all.equal(a, b)))

    unlink(output_dir, recursive = TRUE, force = TRUE)
  })

  test_that("different confidence levels should result in different bands", {

    one.cmt <- function() {
      ini({
        tka <- 0.45 ; label("Log Ka")
        tcl <- 1 ; label("Log Cl")
        tv <- 3.45 ; label("log V")
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
    skip_if_not_installed("rxode2")
    suppressMessages(suppressWarnings(
      fit <-
        nlmixr2(
          one.cmt,
          nlmixr2data::theo_sd,
          est = "focei",
          control = list(print = 0, eval.max = 10),
          table = list(npde = TRUE, cwres = TRUE)
        )
    ))

    suppressMessages(
      fitlist <- nlmixr2extra:::modelBootstrap(fit, nboot = 4, restart = TRUE)[[1]]
    )
    bootSummary1 <- nlmixr2extra:::getBootstrapSummary(fitlist, ci = 0.95)
    bootSummary2 <- nlmixr2extra:::getBootstrapSummary(fitlist, ci = 0.75)

    expect_true(all(
      bootSummary1$Estimate < bootSummary2$Estimate
    ))

    expect_true(all(
      bootSummary1$parFixedDf$confUpper > bootSummary2$parFixedDf$confUpper
    ))
  })

  test_that("expected columns in fit$parFixedDf object should match", {

    one.cmt <- function() {
      ini({
        tka <- 0.45 ; label("Log Ka")
        tcl <- 1 ; label("Log Cl")
        tv <- 3.45 ; label("log V")
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
    skip_if_not_installed("rxode2")
    suppressMessages(suppressWarnings(
      fit <-
        nlmixr2(
          one.cmt,
          nlmixr2data::theo_sd,
          est = "focei",
          control = list(print = 0, eval.max = 10),
          table = list(npde = TRUE, cwres = TRUE)
        )
    ))

    colsBefore <- colnames(fit$parFixedDf)
    suppressMessages(
      fitlist <- nlmixr2extra:::modelBootstrap(fit, nboot = 4, restart = TRUE)[[1]]
    )

    bootSummary <- nlmixr2extra:::getBootstrapSummary(fitlist, ci = 0.95)

    colsAfter <- colnames(fit$parFixedDf)

    expect_equal(colsAfter, colsBefore)

    lapply(
      list.files("./", pattern = "nlmixr2BootstrapCache_.*"),
      function(x) {
        unlink(x, recursive = TRUE, force = TRUE)
      })
  })

  test_that("saem bootstrap", {

    one.cmt <- function() {
      ini({
        tka <- 0.45 ; label("Log Ka")
        tcl <- 1 ; label("Log Cl")
        tv <- 3.45 ; label("log V")
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

    suppressMessages(suppressWarnings(
      fit <-
        nlmixr(
          one.cmt,
          nlmixr2data::theo_sd,
          est = "saem",
          control = list(print = 0, nBurn = 10, nEm = 20),
          table = list(npde = TRUE, cwres = TRUE)
        )
    ))
    skip_if_not_installed("rxode2")
    suppressMessages(
      expect_error(fit1 <- nlmixr2extra:::runBootstrap(fit, nboot = 2, restart = TRUE), NA)
    )

    output_dir <-
      paste0("nlmixr2BootstrapCache_", "fit", "_", fit$bootstrapMd5)
    unlink(output_dir, recursive = TRUE, force = TRUE)
  })

  # ---------------------------------------------------------------------------
  # workers parameter
  # ---------------------------------------------------------------------------

  test_that("runBootstrap: workers=1 completes without error", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
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
    skip_if_not_installed("rxode2")
    suppressMessages(suppressWarnings(
      fit <- nlmixr2(
        one.cmt,
        nlmixr2data::theo_sd,
        est = "focei",
        control = list(print = 0, eval.max = 10),
        table = list(npde = TRUE, cwres = TRUE)
      )
    ))
    suppressMessages(
      expect_no_error(
        nlmixr2extra:::runBootstrap(
          fit, nboot = 2, restart = TRUE, workers = 1L
        )
      )
    )
    lapply(
      list.files("./", pattern = "^fit_boot_[0-9]+$", full.names = TRUE),
      function(x) unlink(x, recursive = TRUE, force = TRUE)
    )
  })

  test_that("runBootstrap: workers=NULL and workers=1 produce same structure", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
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
    skip_if_not_installed("rxode2")
    suppressMessages(suppressWarnings(
      fit <- nlmixr2(
        one.cmt,
        nlmixr2data::theo_sd,
        est = "focei",
        control = list(print = 0, eval.max = 10)
      )
    ))

    suppressMessages({
      res_seq <- nlmixr2extra:::runBootstrap(
        fit, nboot = 2, restart = TRUE, workers = NULL
      )
      res_w1 <- nlmixr2extra:::runBootstrap(
        fit, nboot = 2, restart = FALSE, workers = 1L
      )
    })

    # Both return a named list with the same top-level keys
    expect_named(res_seq, names(res_w1))

    lapply(
      list.files("./", pattern = "^fit_boot_[0-9]+$", full.names = TRUE),
      function(x) unlink(x, recursive = TRUE, force = TRUE)
    )
  })
})
