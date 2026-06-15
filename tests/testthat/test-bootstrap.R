skip_on_cran()

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

    fit1 <- suppressMessages(nlmixr2extra:::bootstrapFit(fit, nboot = 2, restart = TRUE))
    fit2 <- suppressMessages(nlmixr2extra:::bootstrapFit(fit, nboot = 4, restart = FALSE))

    output_dir <-
      paste0("nlmixr2BootstrapCache_", "fit", "_", fit$bootstrapMd5)

    fnameBootDataPattern <- paste0("boot_data",
                                   "_", "[0-9]+", ".rds",
                                   sep = "")

    files <- list.files(paste0("./", output_dir), pattern = fnameBootDataPattern, full.names=TRUE)

    fitdata <- lapply(files, function(x) {
      readRDS(x)
    })

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
      fitlist <- nlmixr2extra:::modelBootstrap(fit, nboot = 2, restart = TRUE)[[1]]
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
    suppressMessages(
      expect_error(fit1 <- nlmixr2extra:::bootstrapFit(fit, nboot = 2, restart = TRUE), NA)
    )

    output_dir <-
      paste0("nlmixr2BootstrapCache_", "fit", "_", fit$bootstrapMd5)
    unlink(output_dir, recursive = TRUE, force = TRUE)
  })

  test_that("bootstrapFit handles a single estimated population parameter", {
    skip_on_cran()
    # A model with one estimated population parameter (a single theta row in
    # parFixedDf) and a single random effect (a 1x1 omega).  Both cases made
    # the 3-D quantile array in getBootstrapSummary() collapse to a vector when
    # sliced (quants[k, , ]), so bootstrapFit() crashed with
    # "dim(X) must have a positive length" (omega branch) /
    # "incorrect number of dimensions" (parFixedDf branch).
    one.par <- function() {
      ini({
        tcl <- 1 ; label("Log Cl")
        eta.cl ~ 0.3
      })
      model({
        ka <- exp(0.45)
        cl <- exp(tcl + eta.cl)
        v <- exp(3.45)
        addSd <- 0.7
        linCmt() ~ add(addSd)
      })
    }
    suppressMessages(suppressWarnings(
      fit <-
        nlmixr(
          one.par,
          nlmixr2data::theo_sd,
          est = "focei",
          control = list(print = 0, eval.max = 10)
        )
    ))
    # confirm the reprex really is the single-parameter / single-eta case
    expect_equal(nrow(fit$parFixedDf), 1L)
    expect_equal(dim(fit$omega), c(1L, 1L))

    # bootstrapFit() used to error here; it must now run to completion
    suppressMessages(suppressWarnings(
      expect_error(
        fitB <- nlmixr2extra:::bootstrapFit(fit, nboot = 5, restart = TRUE),
        NA
      )
    ))

    # the quantile slices must keep their matrix shape, not collapse to vectors
    bootSummary <- fitB$bootSummary
    expect_equal(dim(bootSummary$parFixedDf$median), c(1L, 2L))
    expect_equal(dim(bootSummary$parFixedDf$confLower), c(1L, 2L))
    expect_equal(dim(bootSummary$parFixedDf$confUpper), c(1L, 2L))
    expect_equal(dim(bootSummary$omega$median), c(1L, 1L))
    expect_equal(dim(bootSummary$omega$confLower), c(1L, 1L))
    expect_equal(dim(bootSummary$omega$confUpper), c(1L, 1L))

    # the bootstrap covariance should have been stored on the fit
    expect_equal(fitB$env$covMethod, "boot5")

    lapply(
      list.files("./", pattern = "nlmixr2BootstrapCache_.*"),
      function(x) unlink(x, recursive = TRUE, force = TRUE)
    )
  })

  test_that("bootstrapFit handles a model with no random effects", {
    skip_on_cran()
    # A fixed-effects-only model: fit$omega is NULL.  The bootstrap omega
    # summary has nothing to aggregate and the combined covariance must be built
    # from the fixed effects alone; previously bootstrapFit() errored
    # ("dim(X) must have a positive length" / "'data' must be of a vector type,
    # was 'NULL'").  Printing the summary (with its empty omega entries) also
    # used to error in print.nlmixr2BoostrapSummary().
    no.eta <- function() {
      ini({
        tcl <- 1 ; label("Log Cl")
        addSd <- 0.7
      })
      model({
        ka <- exp(0.45)
        cl <- exp(tcl)
        v <- exp(3.45)
        linCmt() ~ add(addSd)
      })
    }
    suppressMessages(suppressWarnings(
      fit <-
        nlmixr(
          no.eta,
          nlmixr2data::theo_sd,
          est = "focei",
          control = list(print = 0, eval.max = 10)
        )
    ))
    # confirm the reprex really has no random effects
    expect_true(is.null(fit$omega) || length(fit$omega) == 0L)

    # bootstrapFit() used to error here; it must now run to completion
    suppressMessages(suppressWarnings(
      expect_error(
        fitB <- nlmixr2extra:::bootstrapFit(fit, nboot = 5, restart = TRUE),
        NA
      )
    ))

    # there is no omega to summarize, but the combined covariance of the fixed
    # effects must still be a named matrix so setCov() can use it
    bootSummary <- fitB$bootSummary
    expect_null(bootSummary$omega$mean)
    expect_true(is.matrix(bootSummary$omega$covMatrixCombined))
    expect_equal(dim(bootSummary$omega$covMatrixCombined), c(2L, 2L))
    expect_equal(rownames(bootSummary$omega$covMatrixCombined), c("tcl", "addSd"))
    expect_equal(fitB$env$covMethod, "boot5")

    # the bootstrap summary (with empty omega entries) must print without error
    expect_error(
      suppressMessages(utils::capture.output(print(bootSummary))),
      NA
    )

    lapply(
      list.files("./", pattern = "nlmixr2BootstrapCache_.*"),
      function(x) unlink(x, recursive = TRUE, force = TRUE)
    )
  })
})
