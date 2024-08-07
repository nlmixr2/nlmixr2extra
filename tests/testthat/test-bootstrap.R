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
    skip_if_not(rxode2::.linCmtSensB())
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
    skip_if_not(rxode2::.linCmtSensB())
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
    skip_if_not(rxode2::.linCmtSensB())
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
    skip_if_not(rxode2::.linCmtSensB())
    suppressMessages(
      expect_error(fit1 <- nlmixr2extra:::bootstrapFit(fit, nboot = 2, restart = TRUE), NA)
    )

    output_dir <-
      paste0("nlmixr2BootstrapCache_", "fit", "_", fit$bootstrapMd5)
    unlink(output_dir, recursive = TRUE, force = TRUE)
  })
})
