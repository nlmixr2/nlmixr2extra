

test_that("linearize error models", {

  pk.turnover.emax3 <- function() {
    ini({
      tktr <- log(1)
      tka <- log(1)
      tcl <- log(0.1)
      tv <- log(10)
      ##
      eta.ktr ~ 1
      eta.ka ~ 1
      eta.cl ~ 2
      eta.v ~ 1
      prop.err <- 0.1
      pkadd.err <- 0.1
      ##
      temax <- logit(0.8)
      tec50 <- log(0.5)
      tkout <- log(0.05)
      te0 <- log(100)
      ##
      eta.emax ~ .5
      eta.ec50  ~ .5
      eta.kout ~ .5
      eta.e0 ~ .5
      ##
      pdadd.err <- 10
    })
    model({
      ktr <- exp(tktr + eta.ktr)
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      emax = expit(temax+eta.emax)
      ec50 =  exp(tec50 + eta.ec50)
      kout = exp(tkout + eta.kout)
      e0 = exp(te0 + eta.e0)
      ##
      DCP = center/v
      PD=1-emax*DCP/(ec50+DCP)
      ##
      effect(0) = e0
      kin = e0*kout
      ##
      d/dt(depot) = -ktr * depot
      d/dt(gut) =  ktr * depot -ka * gut
      d/dt(center) =  ka * gut - cl / v * center
      d/dt(effect) = kin*PD -kout*effect
      ##
      cp = center / v
      cp ~ prop(prop.err) + add(pkadd.err)
      effect ~ add(pdadd.err) | pca
    })
  }

  f <- rxode2::rxode2(pk.turnover.emax3)

  expect_equal(f$linearizeError,
               list(rxR2 = c("if (OCMT == 5) {",
                             "    rxR2 <- (pkadd.err)^2 + (OPRED)^2 * (prop.err)^2",
                             "    fct <- prop.err^2/prop.err.l^2",
                             "}",
                             "if (OCMT == 6) {",
                             "    rxR2 <- (pdadd.err)^2",
                             "    fct <- 0",
                             "}"),
                    tipred = "TIPRED <- y",
                    err = c("y1 <- y",
                            "y2 <- y",
                            "y1 ~ add(rxR)",
                            "y2 ~ add(rxR)")))

  f1 <- f %>% model(effect ~ lnorm(pdadd.err) + prop(pdprop.err))

  expect_equal(f1$linearizeError,
               list(rxR2 = c("if (OCMT == 5) {",
                             "    rxR2 <- (pkadd.err)^2 + (OPRED)^2 * (prop.err)^2",
                             "    fct <- prop.err^2/prop.err.l^2",
                             "}",
                             "if (OCMT == 4) {",
                             "    rxR2 <- (pdadd.err)^2 + (exp(OPRED))^2 * (pdprop.err)^2",
                             "    fct <- pdprop.err^2/pdprop.err.l^2",
                             "}"),
                    tipred = c("if (OCMT == 5) {",
                               "    TIPRED <- y",
                               "}", "if (OCMT == 4) {",
                               "    TIPRED <- exp(y)",
                               "}"),
                    err = c("y1 <- y",
                            "y2 <- y",
                            "y1 ~ add(rxR)",
                            "y2 ~ lnorm(rxR) + dv()")))

  f1 <- f %>% model(effect ~ logitNorm(pdadd.err, 10, 20) + prop(pdprop.err) + yeoJohnson(lambda))

  expect_equal(f1$linearizeError,
               list(rxR2 = c("if (OCMT == 5) {",
                             "    rxR2 <- (pkadd.err)^2 + (OPRED)^2 * (prop.err)^2",
                             "    fct <- prop.err^2/prop.err.l^2",
                             "}",
                             "if (OCMT == 4) {",
                             "    rxR2 <- (pdadd.err)^2 + (rxTBSi(OPRED, lambda, 5, 10, 20))^2 * ",
                             "        (pdprop.err)^2",
                             "    fct <- pdprop.err^2/pdprop.err.l^2",
                             "}"),
                    tipred = c("if (OCMT == 5) {",
                               "    TIPRED <- y",
                               "}",
                               "if (OCMT == 4) {",
                               "    TIPRED <- rxTBSi(y, lambda, 5, 10, 20)",
                               "}"),
                    err = c("y1 <- y",
                            "y2 <- y",
                            "y1 ~ add(rxR)",
                            "y2 ~ logitNorm(rxR, 10, 20) + yeoJohnson(lambda) + dv()"))
               )

  f1 <- f %>% model(effect ~ logitNorm(pdadd.err, 10, 20) + prop(pdprop.err) + yeoJohnson(lambda) + comb1())

  expect_equal(f1$linearizeError,
               list(rxR2 = c("if (OCMT == 5) {",
                             "    rxR2 <- (pkadd.err)^2 + (OPRED)^2 * (prop.err)^2",
                             "    fct <- prop.err^2/prop.err.l^2",
                             "}",
                             "if (OCMT == 4) {",
                             "    rxR2 <- ((pdadd.err) + (rxTBSi(OPRED, lambda, 5, 10, 20)) * ",
                             "        (pdprop.err))^2",
                             "    fct <- (pdprop.err^2 * (rxTBSi(OPRED, lambda, 5, 10, 20)) + ",
                             "        pdprop.err * pdadd.err)/(pdprop.err.l^2 * (rxTBSi(OPRED, ",
                             "        lambda, 5, 10, 20)) + pdprop.err.l * pdadd.err.l)",
                             "}"),
                    tipred = c("if (OCMT == 5) {",
                               "    TIPRED <- y",
                               "}",
                               "if (OCMT == 4) {",
                               "    TIPRED <- rxTBSi(y, lambda, 5, 10, 20)",
                               "}"),
                    err = c("y1 <- y",
                            "y2 <- y",
                            "y1 ~ add(rxR)",
                            "y2 ~ logitNorm(rxR, 10, 20) + yeoJohnson(lambda) + dv()"))
               )



})

test_that("Linearize add err model ", {
    one.cmpt.adderr <- function() {
        ini({
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            tka <- log(1.56) #  Ka
            eta.cl ~ 0.3
            eta.v ~ 0.1
            eta.ka ~ 0.6
            add.sd <- 0.7
        })
        model({
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd)
        })
    }

    # rxode2::rxode(one.cmpt.adderr)$linearizeError
    fit <- nlmixr(one.cmpt.adderr, nlmixr2data::theo_md, est = "focei")
    linModGen(fit)
    derv <- getDeriv(fit)

    sum(grepl("O_ETA\\d+", names(derv))) |> expect_equal(ncol(fit$eta) - 1)
    all(derv$D_ResVar == 1) |> expect_equal(TRUE)
    all(derv$D_VAR_ETA_1_1 == 0) |> expect_equal(TRUE)
    all(derv$D_VAR_ETA_1_2 == 0) |> expect_equal(TRUE)
    all(derv$D_VAR_ETA_1_3 == 0) |> expect_equal(TRUE)

    fitLin <- linearize(fit)

    linearizePlot(fit, fitLin)
    isLinearizeMatch(fit, fitLin)$ofv[[1]] |> expect_true()

    # derv2 <- getDeriv(fitLin)
    # plot(derv2$D_ETA1, derv$D_ETA1)
})


test_that("Linearize prop err model ", {
    one.cmpt.properr <- function() { # non-linear base
        ini({
            tka <- log(1.56) # Ka
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            eta.cl ~ 0.3
            eta.v ~ 0.2
            prop.sd <- 0.10
        })
        model({
            ka <- exp(tka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ prop(prop.sd)
        })
    }
    # set.seed(42)
    # ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
    #     rxode2::et(amt = 300, cmt = "depot") |>
    #     rxode2::et(time = c(0.25, 0.5, 1,2,3,6,8,12,16,24))
    # sim <- rxode2::rxSolve(one.cmpt.properr, ev, nSub = 200, addDosing = TRUE)
    # sim$dv <- sim$sim
    # sim$id <- sim$sim.id
    # sim$sim.id <- NULL
    # sim <- sim[,c("id", "time", "amt", "dv", "evid")]

    # saveRDS(sim, "one.cmpt.properr_data.RDS")
    sim <- readRDS(system.file("one.cmpt.properr_data.RDS", package="nlmixr2extra"))
    fit <- nlmixr(one.cmpt.properr, sim, est = "focei",
            control = nlmixr2est::foceiControl(mceta=10))
    # linMod <- linModGen(fit, FALSE)
    # derv <- getDeriv(fit)

    # fitLin <- nlmixr(linMod, derv, est="focei",
    #         control = nlmixr2est::foceiControl(etaMat = fit, mceta=10,
    #         covMethod = "",
    #         calcTables=FALSE,
    #         maxInnerIterations=100, maxOuterIterations=100))


    # fit$scaleInfo$scaleC
    # fitLin$scaleInfo$scaleC
    # INNER ETA
    # OUTER THETA OMEGA
    #

    fitLin <- linearize(fit, relTol = 0.2, mceta = c(-1, 10, 100))
    isLinearizeMatch(fit, fitLin, tol = 0.15)$ofv[[1]] |> expect_true()

})


test_that("Linearize combined2 model ", {
    one.cmpt.combinederr <- function() {
        ini({
            tka <- log(1.56) # Ka
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            eta.cl ~ 0.3
            eta.v ~ 0.1
            add.sd <- 0.7
            prop.sd <- 0.1
        })
        model({
            ka <- exp(tka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd) + prop(prop.sd) + combined2()
        })
    }

    fit <- nlmixr(one.cmpt.combinederr, nlmixr2data::theo_sd, est = "focei")
    # derv <- getDeriv(fit)

    # sum(grepl("O_ETA\\d+", names(derv))) |> expect_equal(ncol(fit$eta) - 1)
    # all(derv$D_ResVar == 1) |> expect_equal(TRUE)

    fitLin <- linearize(fit)
    isLinearizeMatch(fit, fitLin)$ofv[[1]] |> expect_true()
})


test_that("Linearize combined1 model ", {
    one.cmpt.combinederr <- function() {
        ini({
            tka <- log(1.56) # Ka
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            eta.cl ~ 0.3
            eta.v ~ 0.1
            add.sd <- 0.7
            prop.sd <- 0.1
        })
        model({
            ka <- exp(tka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd) + prop(prop.sd) + combined1()
        })
    }

    fit <- nlmixr(one.cmpt.combinederr, nlmixr2data::theo_sd, est = "focei")
    # derv <- getDeriv(fit)

    # sum(grepl("O_ETA\\d+", names(derv))) |> expect_equal(ncol(fit$eta) - 1)
    # all(derv$D_ResVar == 1) |> expect_equal(TRUE)

    fitLin <- linearize(fit)
    isLinearizeMatch(fit, fitLin)$ofv[[1]] |> expect_true()
})



test_that("Linearize multiple endpoints ", {
    skip_on_cran()

    pk.turnover.emax3 <- function() {
    ini({
        tktr <- 0.326787337229061
        tka <- 0.573847838322594
        tcl <- -2.02615620220267
        tv <- 2.06443263542308
        # prop.err <- c(0, 0.145017473272093)
        pkadd.err <- c(0, 0.198298188759467)
        temax <- 4.75044304930452
        tec50 <- 0.139562783557545
        tkout <- -2.93951993228171
        te0 <- 4.56394773529773
        pdadd.err <- c(0, 3.80147740823949)
        eta.ktr ~ 0.899432614108315
        eta.ka ~ 1.20111175373864
        eta.cl ~ 0.0756455519562788
        eta.v ~ 0.0547982604779906
        eta.emax ~ 0.547762672352332
        eta.ec50 ~ 0.204313071604677
        eta.kout ~ 0.0228239781890516
        eta.e0 ~ 0.0107988390114995
    })
    model({
        ktr <- exp(tktr + eta.ktr)
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        emax = expit(temax + eta.emax, 0, 1)
        ec50 = exp(tec50 + eta.ec50)
        kout = exp(tkout + eta.kout)
        e0 = exp(te0 + eta.e0)
        DCP = center/v
        PD = 1 - emax * DCP/(ec50 + DCP)
        effect(0) = e0
        kin = e0 * kout
        d/dt(depot) = -ktr * depot
        d/dt(gut) = ktr * depot - ka * gut
        d/dt(center) = ka * gut - cl/v * center
        d/dt(effect) = kin * PD - kout * effect
        cp = center/v
        # cp ~ prop(prop.err) + add(pkadd.err)
        cp ~  add(pkadd.err)
        effect ~ add(pdadd.err) | pca
    })
}

    # pk.turnover.emax3 <- pk.turnover.emax3()
    # y <- linModGen(pk.turnover.emax3)
    # set.seed(999)
    # ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
    #     rxode2::et(amt = 100, cmt = "depot") |>
    #     rxode2::et(time =  c(0, 0.5, 1, 1.5, 2, 3, 6, 9, 12, 24, 36, 48, 72, 96, 120) , cmt="cp")|>
    #     rxode2::et(time = c(0, 24, 36, 48, 72, 96, 120, 144), cmt="pca")
    # sim <- rxode2::rxSolve(pk.turnover.emax3, ev, nSub = 50, addDosing = TRUE)
    # plot(sim, "sim")
    # sim$dv <- sim$sim
    # sim$id <- sim$sim.id
    # sim$sim.id <- NULL
    # sim$dvid <- ifelse(sim$CMT == 5, "cp", "pca")

    # sim <- sim[,c("id", "time", "amt", "dv", "dvid", "evid")]
    # # sim <- sim[,c("id", "time", "amt", "dv",  "evid")]

    # fit <- nlmixr(pk.turnover.emax3, sim, est = "focei",
    #         control = nlmixr2est::foceiControl(mceta=10,
    #         calcTables=TRUE))
    # # saveRDS(fit, "warfarin.RDS")
    fit <- readRDS(system.file("warfarin.RDS", package="nlmixr2extra"))

    # derv <- getDeriv(fit)
    # # derv$CMT <- NULL

    # lmod <- linModGen(fit, TRUE)
    # fitLin <- evalLinModel(fit, lmod, derv, 1000)
    # isLinearizeMatch(fit, fitLin)
    # linearizePlot(fit, fitLin)
    # fitLin <- nlmixr(lmod, derv, est = "focei")

    fitLin <- linearize(fit, mceta = 10, focei = FALSE)
    isLinearizeMatch(fit, fitLin)$ofv[[1]] |> expect_true()
    linearizePlot(fit, fitLin)

    # fit$dataMergeInner$nlmixrLlikObs[1:20]
    # fitLin$dataMergeInner$nlmixrLlikObs[1:20]
    # derv2 <- getDeriv(fitLin) # FIXME consider droping names starts with D_ETA ...etc for avoid conflict with new derv
    # plot(derv2$D_ETA1, derv$D_ETA1)
})


test_that("linearization covariance", {
    one.cmpt.adderr <- function() {
        ini({
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            tka <- log(1.56) #  Ka
            eta.cl ~ 0.3
            eta.v ~ 0.1
            eta.ka ~ 0.6
            add.sd <- 0.7
        })
        model({
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd)
        })
    }
    fit <- nlmixr(one.cmpt.adderr, nlmixr2data::theo_md, est = "focei")
    derv <- getDeriv(fit)
    linMod <- linModGen(fit, FALSE)
    fitLin <- evalLinModel(fit, linMod, derv, 0, covMethod = "r,s")

    fitLin

})

test_that("linearize correlated eta ", {
    one.cmpt.adderr <- function() {
        ini({
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            tka <- log(1.56) #  Ka
            eta.cl + eta.v ~ c(0.3,
                                0.1, 0.1)
            eta.ka ~ 0.6
            add.sd <- 0.7
        })
        model({
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd)
        })
    }
    fit <- nlmixr(one.cmpt.adderr, nlmixr2data::theo_md, est = "focei")
    fitLin <- linearize(fit)
    isLinearizeMatch(fit, fitLin)$ofv[[1]] |> expect_true()
})

test_that("Adding covariates to lin models", {


    one.cmpt.adderr <- function() {
        ini({
            tcl <- log(2.7) # Cl
            tv <- 30 # V
            tka <- log(1.56) #  Ka
            eta.cl ~ 0.3
            eta.v ~ 0.1
            eta.ka ~ 0.6
            add.sd <- 0.7
        })
        model({
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            v <- tv*exp(eta.v)
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd)
        })
    }
    one.cmpt.adderr.cov <- function() {
        ini({
            tcl <- log(2.7) # Cl
            tv <- 30 # V
            tka <- log(1.56) #  Ka
            WTVTheta <- 1.5
            eta.cl ~ 0.3
            eta.v ~ 0.1
            eta.ka ~ 0.6
            add.sd <- 0.7
        })
        model({
            WTVCOV = (WT/70)^WTVTheta
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            v <- tv*exp(eta.v)*WTVCOV
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd)
        })
    }

    set.seed(42)
    ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
        rxode2::et(amt = 200, cmt = "depot") |>
        rxode2::et(time = c(0.25, 0.5, 1,2,3,6,8,12,16,24)) |>
        rxode2::et(id = 1:200)
    theo_sd <- rxode2::rxSolve(one.cmpt.adderr.cov, ev, nSub = 200, addDosing = TRUE,
        iCov=data.frame(id=1:200, WT=rnorm(200, 70, 10)))
    theo_sd$dv <- theo_sd$sim
    theo_sd <- theo_sd[,c("id", "time", "amt", "dv", "evid", "WT")]

    nlfitNoCov <- nlmixr(one.cmpt.adderr, theo_sd, est = "focei")
    fitLinNoCov <- linearize(nlfitNoCov)
    fitLinCov <- addCovariate(fitLinNoCov, eta.v~WT/70, effect = "power")
      
    fitLinCov <- nlmixr(fitLinCov, nlme::getData(fitLinNoCov), est = "focei")

    nlfitCov <- nlmixr(one.cmpt.adderr.cov, nlme::getData(nlfitNoCov), est = "focei")

    nlfitCov$parFixed # FIXME why 0.5 and not 1.5?
    fitLinCov$parFixed

    expect_true(fitLinCov$objDf$OBJF < fitLinNoCov$objDf$OBJF)
    expect_true(nlfitCov$objDf$OBJF < nlfitNoCov$objDf$OBJF)
    expect_equal(nlfitCov$objDf$OBJF, fitLinCov$objDf$OBJF, tolerance = 0.1)
})
