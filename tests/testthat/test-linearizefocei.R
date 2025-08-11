

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
    fit <- nlmixr(one.cmpt.adderr, nlmixr2data::theo_sd, est = "focei")
    derv <- getDeriv(fit)

    all(c("O_eta.cl", "O_eta.v", "O_eta.ka") %in% names(derv)) %>% expect_true()
    all(derv$D_ResVar == 1) %>% expect_equal(TRUE)
    all(derv$D_VAR_eta.cl == 0) %>% expect_equal(TRUE)
    all(derv$D_VAR_eta.v == 0) %>% expect_equal(TRUE)
    all(derv$D_VAR_eta.ka == 0) %>% expect_equal(TRUE)

    suppressWarnings(
        fitLin <- linearize(fit)
    )

    linearizePlot(fit, fitLin) %>% expect_no_error()
    expect_true(
        all(sapply(isLinearizeMatch(fit, fitLin), function(x){x[[1]]}))
    )

})



test_that("Linearize prop err model ", {
    skip_on_cran()

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
    set.seed(42)
    ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
        rxode2::et(amt = 300, cmt = "depot") |>
        rxode2::et(time = c(0.125, 0.5, 1,2,3,6,8,12,16,24))
    sim <- rxode2::rxSolve(one.cmpt.properr, ev, nSub = 200, addDosing = TRUE)
    sim$dv <- sim$sim
    sim$id <- sim$sim.id
    sim$sim.id <- NULL
    sim <- sim[,c("id", "time", "amt", "dv", "evid")]

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

    suppressWarnings(
        fitLin <- linearize(fit, relTol = 0.2, mceta = c(-1, 10))
    )

    isLinearizeMatch(fit, fitLin, tol = 0.15)$ofv[[1]] %>% expect_true()

})


test_that("Linerize pheno prop err", {
    one.cmpt.prop.iv <- function() {
        ini({
            tcl <- log(0.01) # Cl
            tv <- log(0.9) # V
            eta.cl ~ 0.1
            eta.v ~ 0.1
            prop.sd <- 0.1
        })
        model({
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d / dt(center) <- - cl / v * center
            cp <- center / v
            cp <- cp
            cp ~ prop(prop.sd) 
        })
    }
    fit <- nlmixr(one.cmpt.prop.iv, nlmixr2data::pheno_sd, est = "focei")

    suppressWarnings(
        fitLin <- linearize(fit)
    )
    
    # increase error to 10% for few outliers
    expect_true( 
        all(sapply(isLinearizeMatch(fit, fitLin, 0.1), function(x){x[[1]]}))
    )

})

test_that("Linearize combined2 model ", {
    one.cmpt.combinederr <- function() {
        ini({
            tka <- log(1.56) # Ka
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            eta.cl ~ 0.3
            eta.v ~ 0.1
            add.sd <- 1
            prop.sd <- 0.3
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

    fit <- nlmixr(one.cmpt.combinederr, nlmixr2data::theo_sd, est = "focei") # 163.56
    suppressWarnings(
        fitLin <- linearize(fit) # saem doesn't work well
    )

    expect_true( 
        all(sapply(isLinearizeMatch(fit, fitLin, 0.1), function(x){x[[1]]}))
    )
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
    suppressWarnings(
        fitLin <- linearize(fit)
    )

    expect_true( 
        all(sapply(isLinearizeMatch(fit, fitLin, 0.1), function(x){x[[1]]}))
    )
})



test_that("Linearize multiple endpoints ", {
    skip_on_cran()

    pk.turnover.emax3 <- function() {
    ini({
        tktr <- 0.326787337229061
        tka <- 0.573847838322594
        tcl <- -2.02615620220267
        tv <- 2.06443263542308
        prop.err <- c(0, 0.145017473272093)
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
        cp ~ prop(prop.err) + add(pkadd.err) | cp
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

    fit <- nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, est = "focei",
            control = nlmixr2est::foceiControl(mceta=10,
            calcTables=TRUE))
    # # saveRDS(fit, "warfarin.RDS")
    # fit <- readRDS(system.file("warfarin.RDS", package="nlmixr2extra"))

    suppressWarnings(
        # this will switch to FOCE even after successful evaluation. 
        # The match was 4% for OFV, but mismatch for omega/eta
        fitLin <- linearize(fit, mceta = 10)
    )

    expect_true(isLinearizeMatch(fit, fitLin)$ofv[[1]])

    linearizePlot(fit, fitLin)
})


test_that("Linearize multiple endpoints SAEM", {
    skip_on_cran()
    
     

    pk.turnover.emax3 <- function() {
    ini({
        tktr <- 0.326787337229061
        tka <- 0.573847838322594
        tcl <- -2.02615620220267
        tv <- 2.06443263542308
        prop.err <- c(0, 0.145017473272093)
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
        # eta.e0 ~ 0.0107988390114995
    })
    model({
        ktr <- exp(tktr + eta.ktr)
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        emax = expit(temax + eta.emax, 0, 1)
        ec50 = exp(tec50 + eta.ec50)
        kout = exp(tkout + eta.kout)
        e0 = exp(te0)
        DCP = center/v
        PD = 1 - emax * DCP/(ec50 + DCP)
        effect(0) = e0
        kin = e0 * kout
        d/dt(depot) = -ktr * depot
        d/dt(gut) = ktr * depot - ka * gut
        d/dt(center) = ka * gut - cl/v * center
        d/dt(effect) = kin * PD - kout * effect
        cp = center/v
        cp ~ prop(prop.err) + add(pkadd.err) | cp
        effect ~ add(pdadd.err) | pca
    })
}

    fit <- nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, est = "saem") 

    # pdf("fitplot.pdf")
    # plot(fit)
    # dev.off()
    
    suppressWarnings(
        fitLin <- linearize(fit, mceta = 10)
    )

    expect_true(isLinearizeMatch(fit, fitLin)$ofv[[1]]) # FIXME this fails despite good fit

    linearizePlot(fitLin)
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
    suppressWarnings(
        fitLin <- linearize(fit)
    )
    isLinearizeMatch(fit, fitLin)$ofv[[1]] %>% expect_true()
    
    expect_true( 
        all(sapply(isLinearizeMatch(fit, fitLin, 0.1), function(x){x[[2]]}))
    )

})

test_that("covariate parse", {
    dat <- nlmixr2data::warfarin
    
   expect_error(parseCovExpr(eta.v~wt/70+BMI, dat, effect = "power"), 
    "Not all covariates are present"
   )

    parseCovExpr(eta.v~wt/70+sex, dat, effect = "power")
    x <- parseCovExpr(eta.v~wt/70+sex, dat, effect = "hockyStick")

    expect_error(parseCovExpr(eta.v~wt/70+sex/70, dat, effect = "power"))
    expect_error(parseCovExpr(eta.v~wt/med+sex, dat, effect = "power"), "Divide only by")

    expect_true(parseCovExpr(eta.v~WT/median, nlmixr2data::theo_sd, effect = "power")$normFactor == 70.5)

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
    one.cmpt.adderr.cov.all <- function() {
        ini({
            tcl <- log(2.7) # Cl
            tv <- 30 # V
            tka <- log(1.56) #  Ka
            WTVTheta <- 1.5
            AGECLTheta <- 2
            eta.cl ~ 0.3
            eta.v ~ 0.1
            eta.ka ~ 0.6
            add.sd <- 0.7
        })
        model({
            WTVCOV = (WT/70)^WTVTheta
            AGECLCOV = (AGE/30)^AGECLTheta
            ka <- exp(tka + eta.ka)
            cl <- tcl*exp(eta.cl)*AGECLCOV
            v <- tv*exp(eta.v)*WTVCOV
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd)
        })
    }
    
    one.cmpt.adderr.cov.wtv <- function() {
        ini({
            tcl <- log(2.7) # Cl
            tv <- 30 # V
            tka <- log(1.56) #  Ka
            WTVTheta <- 1.5
            # AGECLTheta <- 2
            eta.cl ~ 0.3
            eta.v ~ 0.1
            eta.ka ~ 0.6
            add.sd <- 0.7
        })
        model({
            WTVCOV = (WT/70)^WTVTheta
            # AGECLCOV = (AGE/30)^AGECLTheta
            ka <- exp(tka + eta.ka)
            cl <- tcl*exp(eta.cl)
            v <- tv*exp(eta.v)*WTVCOV
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd)
        })
    }
    
    
    one.cmpt.adderr.cov.agecl <- function() {
        ini({
            tcl <- log(2.7) # Cl
            tv <- 30 # V
            tka <- log(1.56) #  Ka
            # WTVTheta <- 1.5
            AGECLTheta <- 2
            eta.cl ~ 0.3
            eta.v ~ 0.1
            eta.ka ~ 0.6
            add.sd <- 0.7
        })
        model({
            # WTVCOV = (WT/70)^WTVTheta
            AGECLCOV = (AGE/30)^AGECLTheta
            ka <- exp(tka + eta.ka)
            cl <- tcl*exp(eta.cl)*AGECLCOV
            v <- tv*exp(eta.v)
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
        iCov=data.frame(id=1:200, WT=rnorm(200, 70, 10), AGE=rnorm(200, 30, 10)))
    theo_sd$dv <- theo_sd$sim
    theo_sd <- theo_sd[,c("id", "time", "amt", "dv", "evid", "WT", "AGE")]

    nlfitNoCov <- nlmixr(one.cmpt.adderr, theo_sd, est = "focei")
    fitLinNoCov <- linearize(nlfitNoCov)
    expect_no_error(
        x <- addCovariate(fitLinNoCov, eta.v~WT/70+AGE/median, effect = "power") %>% 
            addCovariate(eta.cl~WT/80) 
        
    )
    expect_error(
        addCovariate(fitLinNoCov, eta.v~WT/70, effect = "power") %>% 
            addCovariate(eta.v~WT/70), 
            "Duplicated names found"
    )
    expect_error(
        addCovariate(fitLinNoCov, eta.v~AGPR/70, effect = "power") %>% 
            addCovariate(eta.v~WT/70)
    )
    # only Age CL 
    fitLinCov.age <- addCovariate(fitLinNoCov, eta.cl~AGE/30, effect = "power")
    fitLinCov.age <- nlmixr(fitLinCov.age, nlme::getData(fitLinNoCov), est = "focei")
    nlfitCov.age <- nlmixr(one.cmpt.adderr.cov.agecl, nlme::getData(nlfitNoCov), est = "focei")
    
    
    # only Wt Vd 
    
    fitLinCov.wt <- addCovariate(fitLinNoCov, eta.v~WT/70, effect = "power") 

    fitLinCov.wt <- nlmixr(fitLinCov.wt, nlme::getData(fitLinNoCov), est = "focei")

    nlfitCov.wt <- nlmixr(one.cmpt.adderr.cov.wtv, nlme::getData(nlfitNoCov), est = "focei")
    
    # correct 
    fitLinCov <- addCovariate(fitLinNoCov, eta.v~WT/70, effect = "power") |> 
                  addCovariate(eta.cl~AGE/30, effect = "power")

    fitLinCov <- nlmixr(fitLinCov, nlme::getData(fitLinNoCov), est = "focei")

    nlfitCov <- nlmixr(one.cmpt.adderr.cov, nlme::getData(nlfitNoCov), est = "focei")

    nlfitCov$parFixed # why 0.5 and not 1.5? normalization
    fitLinCov$parFixed

    nlfitCov.wt$objDf[,c(1,3)]
    nlfitCov.age$objDf[,c(1,3)]
    nlfitCov$objDf[,c(1,3)]
    
    fitLinCov.wt$objDf[,c(1,3)]
    fitLinCov.age$objDf[,c(1,3)]
    fitLinCov$objDf[,c(1,3)]
    
    expect_true(sum(nlfitCov$time) > sum(fitLinCov$time))
    expect_true(fitLinCov$objDf$OBJF < fitLinNoCov$objDf$OBJF)
    expect_true(nlfitCov$objDf$OBJF < nlfitNoCov$objDf$OBJF)
    expect_equal(nlfitCov$objDf$OBJF, fitLinCov$objDf$OBJF, tolerance = 0.1)
})

test_that("linModGen from any object", {
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

    # call 
    linModGen(one.cmpt.adderr) %>% expect_no_error()
    linModGen(one.cmpt.adderr, derivFct= TRUE) %>% expect_error()

    # rxUi 
    linModGen(one.cmpt.adderr()) %>% expect_no_error()

    # fit => in previous models

})


test_that("linearize mavoglorant lnorm", {
    pbpk <- function(){
      ini({
        ##theta=exp(c(1.1, .3, 2, 7.6, .003, .3))
        lKbBR = 1.1
        lKbMU = 0.3
        lKbAD = 2
        lCLint = 7.6
        lKbBO = 0.03
        lKbRB = 0.3
        eta.LClint ~ 4
        lnormsd <- 10
      })
      model({
        KbBR = exp(lKbBR)
        KbMU = exp(lKbMU)
        KbAD = exp(lKbAD)
        CLint= exp(lCLint + eta.LClint)
        KbBO = exp(lKbBO)
        KbRB = exp(lKbRB)

        ## Regional blood flows
        CO  = (187.00*WT^0.81)*60/1000;         # Cardiac output (L/h) from White et al (1968)
        QHT = 4.0 *CO/100;
        QBR = 12.0*CO/100;
        QMU = 17.0*CO/100;
        QAD = 5.0 *CO/100;
        QSK = 5.0 *CO/100;
        QSP = 3.0 *CO/100;
        QPA = 1.0 *CO/100;
        QLI = 25.5*CO/100;
        QST = 1.0 *CO/100;
        QGU = 14.0*CO/100;
        QHA = QLI - (QSP + QPA + QST + QGU); # Hepatic artery blood flow
        QBO = 5.0 *CO/100;
        QKI = 19.0*CO/100;
        QRB = CO - (QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI);
        QLU = QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI + QRB;

        ## Organs' volumes = organs' weights / organs' density
        VLU = (0.76 *WT/100)/1.051;
        VHT = (0.47 *WT/100)/1.030;
        VBR = (2.00 *WT/100)/1.036;
        VMU = (40.00*WT/100)/1.041;
        VAD = (21.42*WT/100)/0.916;
        VSK = (3.71 *WT/100)/1.116;
        VSP = (0.26 *WT/100)/1.054;
        VPA = (0.14 *WT/100)/1.045;
        VLI = (2.57 *WT/100)/1.040;
        VST = (0.21 *WT/100)/1.050;
        VGU = (1.44 *WT/100)/1.043;
        VBO = (14.29*WT/100)/1.990;
        VKI = (0.44 *WT/100)/1.050;
        VAB = (2.81 *WT/100)/1.040;
        VVB = (5.62 *WT/100)/1.040;
        VRB = (3.86 *WT/100)/1.040;

        ## Fixed parameters
        BP = 0.61;      # Blood:plasma partition coefficient
        fup = 0.028;    # Fraction unbound in plasma
        fub = fup/BP;   # Fraction unbound in blood

        KbLU = exp(0.8334);
        KbHT = exp(1.1205);
        KbSK = exp(-.5238);
        KbSP = exp(0.3224);
        KbPA = exp(0.3224);
        KbLI = exp(1.7604);
        KbST = exp(0.3224);
        KbGU = exp(1.2026);
        KbKI = exp(1.3171);

        ##-----------------------------------------
        S15 = VVB*BP/1000;
        C15 = Venous_Blood/S15

        ##-----------------------------------------
        d/dt(Lungs) = QLU*(Venous_Blood/VVB - Lungs/KbLU/VLU);
        d/dt(Heart) = QHT*(Arterial_Blood/VAB - Heart/KbHT/VHT);
        d/dt(Brain) = QBR*(Arterial_Blood/VAB - Brain/KbBR/VBR);
        d/dt(Muscles) = QMU*(Arterial_Blood/VAB - Muscles/KbMU/VMU);
        d/dt(Adipose) = QAD*(Arterial_Blood/VAB - Adipose/KbAD/VAD);
        d/dt(Skin) = QSK*(Arterial_Blood/VAB - Skin/KbSK/VSK);
        d/dt(Spleen) = QSP*(Arterial_Blood/VAB - Spleen/KbSP/VSP);
        d/dt(Pancreas) = QPA*(Arterial_Blood/VAB - Pancreas/KbPA/VPA);
        d/dt(Liver) = QHA*Arterial_Blood/VAB + QSP*Spleen/KbSP/VSP + QPA*Pancreas/KbPA/VPA + QST*Stomach/KbST/VST + QGU*Gut/KbGU/VGU - CLint*fub*Liver/KbLI/VLI - QLI*Liver/KbLI/VLI;
        d/dt(Stomach) = QST*(Arterial_Blood/VAB - Stomach/KbST/VST);
        d/dt(Gut) = QGU*(Arterial_Blood/VAB - Gut/KbGU/VGU);
        d/dt(Bones) = QBO*(Arterial_Blood/VAB - Bones/KbBO/VBO);
        d/dt(Kidneys) = QKI*(Arterial_Blood/VAB - Kidneys/KbKI/VKI);
        d/dt(Arterial_Blood) = QLU*(Lungs/KbLU/VLU - Arterial_Blood/VAB);
        d/dt(Venous_Blood) = QHT*Heart/KbHT/VHT + QBR*Brain/KbBR/VBR + QMU*Muscles/KbMU/VMU + QAD*Adipose/KbAD/VAD + QSK*Skin/KbSK/VSK + QLI*Liver/KbLI/VLI + QBO*Bones/KbBO/VBO + QKI*Kidneys/KbKI/VKI + QRB*Rest_of_Body/KbRB/VRB - QLU*Venous_Blood/VVB;
        d/dt(Rest_of_Body) = QRB*(Arterial_Blood/VAB - Rest_of_Body/KbRB/VRB);

        C15 ~ lnorm(lnormsd)
      })
    }

    dat <- nlmixr2data::mavoglurant
    dat$occ = unlist(with(dat, tapply(EVID, ID, function(x) cumsum(x>0))))
    dat = subset(dat, occ==1)
    dat = subset(dat, ID<812) ## First 20
    dat = subset(dat, EVID>0 | DV>0)
    dat$CMT[dat$CMT == 0]  <- 1;
    dat$CMT[dat$EVID == 1]  <- "Venous_Blood" ## Compartment dosed to is Venous Blood
    dat$CMT[dat$EVID != 1]  <- "C15" ## Observing C15

    fit <- nlmixr(pbpk, dat, est="focei", control=list(print=20))

    suppressWarnings(
        fitLin <- linearize(fit, mceta = c(-1, 10))
    )


    expect_true( 
        all(sapply(isLinearizeMatch(fit, fitLin, 0.1), function(x){x[[2]]}))
    )

})

