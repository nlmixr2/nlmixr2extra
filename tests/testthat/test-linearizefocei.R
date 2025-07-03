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
               list(rxR2 = c("if (CMT == 5) {",
                             "    rxR2 <- (pkadd.err)^2 + (OPRED)^2 * (prop.err)^2",
                             "}",
                             "if (CMT == 6) {",
                             "    rxR2 <- (pdadd.err)^2",
                             "}"),
                    tipred = "TIPRED <- y",
                    err = c("y1 <- y",
                            "y2 <- y",
                            "y1 ~ add(rxR) + var()",
                            "y2 ~ add(rxR) + var()")))

  f1 <- f %>% model(effect ~ lnorm(pdadd.err) + prop(pdprop.err))

  expect_equal(f1$linearizeError,
               list(rxR2 = c("if (CMT == 5) {",
                             "    rxR2 <- (pkadd.err)^2 + (OPRED)^2 * (prop.err)^2",
                             "}",
                             "if (CMT == 4) {",
                             "    rxR2 <- (pdadd.err)^2 + (exp(OPRED))^2 * (pdprop.err)^2",
                             "}"),
                    tipred = c("if (CMT == 5) {",
                               "    TIPRED <- y",
                               "}", "if (CMT == 4) {",
                               "    TIPRED <- exp(y)",
                               "}"),
                    err = c("y1 <- y",
                            "y2 <- y",
                            "y1 ~ add(rxR) + var()",
                            "y2 ~ lnorm(rxR) + var() + dv()")))

  f1 <- f %>% model(effect ~ logitNorm(pdadd.err, 10, 20) + prop(pdprop.err) + yeoJohnson(lambda))

  expect_equal(f1$linearizeError,
               list(rxR2 = c("if (CMT == 5) {",
                             "    rxR2 <- (pkadd.err)^2 + (OPRED)^2 * (prop.err)^2",
                             "}",
                             "if (CMT == 4) {",
                             "    rxR2 <- (pdadd.err)^2 + (rxTBSi(OPRED, lambda, 5, 10, 20))^2 * ",
                             "        (pdprop.err)^2",
                             "}"),
                    tipred = c("if (CMT == 5) {",
                               "    TIPRED <- y",
                               "}",
                               "if (CMT == 4) {",
                               "    TIPRED <- rxTBSi(y, lambda, 5, 10, 20)",
                               "}"),
                    err = c("y1 <- y",
                            "y2 <- y",
                            "y1 ~ add(rxR) + var()",
                            "y2 ~ logitNorm(rxR, 10, 20) + var() + yeoJohnson(lambda) + dv()"))
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

    fit <- nlmixr(one.cmpt.adderr, nlmixr2data::theo_md, est = "focei")
    derv <- getDerv(fit)

    sum(grepl("O_ETA\\d+", names(derv))) |> expect_equal(ncol(fit$eta) - 1)
    all(derv$D_ResVar == 1) |> expect_equal(TRUE)
    all(derv$D_VAR_ETA_1_1 == 0) |> expect_equal(TRUE)
    all(derv$D_VAR_ETA_1_2 == 0) |> expect_equal(TRUE)
    all(derv$D_VAR_ETA_1_3 == 0) |> expect_equal(TRUE)

    lmod <- linModGen(fit)

    fitLin <- nlmixr(lmod, derv, est = "focei")

    diffObjF <- abs((fit$objDf$OBJF - fitLin$objDf$OBJF)/ fitLin$objDf$OBJF)
    expect_true(diffObjF < 0.1)

    all.equal(fit$omega, fitLin$omega, tolerance = 0.1) |> expect_true()
})


test_that("Linearize prop err model ", {
    one.cmpt.properr <- function() { # non-linear base
        ini({
            tka <- log(1.56) # Ka
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            eta.cl ~ 0.3
            eta.v ~ 0.1
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
    ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
        rxode2::et(amt = 320, cmt = "depot")
    sim <- rxSolve(one.cmpt.properr, ev, nSub = 100, addDosing = TRUE)
    sim$DV <- sim$sim
    sim$ID <- sim$sim.id
    sim$sim.id <- NULL

    fit <- nlmixr(one.cmpt.properr, sim, est = "focei")
    derv <- getDerv(fit)

    sum(grepl("O_ETA\\d+", names(derv))) |> expect_equal(ncol(fit$eta) - 1)
    all(derv$D_ResVar == 1) |> expect_equal(TRUE)

    lmod <- linModGen(fit)

# TODO try remove the var() or
# TODO try squaring the term
lmod <- function() {
ini({
    eta.cl ~ 0.3
    eta.v ~ 0.1
    # eta.ka ~ 0.6
    prop.sd <- 0.1 # actual prop.sd
  })
  model({

    base1 = D_ETA1*(-O_ETA1 + eta.cl)
    base2 = D_ETA2*(- O_ETA2 + eta.v)
    # base3 = D_ETA3*(-O_ETA3 + eta.ka)

    BASE_TERMS = base1 + base2 #+ base3

    IPRED = BASE_TERMS + OPRED

    ERR1 = D_VAR_ETA_1_1*(-O_ETA1 + eta.cl)
    ERR2 = D_VAR_ETA_1_2*(-O_ETA2 + eta.v)
    # ERR3 = D_EPSETA_1_3*(-O_ETA3 + eta.ka)

    BASE_ERROR1 = (ERR1 + ERR2)+(prop.sd*OPRED)**2 # endpoint1
    BASE_ERROR2 = (ERR3 + ERR4)+(prop.sd*OPRED)**2 # endpoint2

    R2 = (BASE_ERROR1+BASE_ERROR2)

    y = IPRED
    y ~ add(R2) + var()
  })
}


    fitLin <- nlmixr(lmod, derv, est = "focei",
        control = nlmixr2::foceiControl(etaMat = as.matrix(fit$eta[-1]), mceta = -1))

    fitLin <- nlmixr(lmod, derv, est = "focei",
        control = nlmixr2::foceiControl(etaMat = as.matrix(fit$eta[-1]), mceta = 100))

    fitLin$parFixedDf
    diffObjF <- abs((fit$objDf$OBJF - fitLin$objDf$OBJF)/ fitLin$objDf$OBJF)
    diffObjF
    expect_true(diffObjF < 0.1)

    all.equal(fit$omega, fitLin$omega, tolerance = 0.1) |> expect_true()
    all.equal(fit$eta, fitLin$eta, tolerance = 0.1) |> expect_true()
})


test_that("Linearize combined model ", {
    one.cmpt.combinederr <- function() {
        ini({
            tka <- log(1.56) # Ka
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            eta.ka ~ 0.6
            eta.cl ~ 0.3
            eta.v ~ 0.1
            add.sd <- 0.7
            prop.sd <- 0.1
        })
        model({
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ add(add.sd) + prop(prop.sd)
        })
    }

    # ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
    #     rxode2::et(amt = 10000, cmt = "depot")
    # sim <- rxSolve(one.cmpt.combinederr, ev, nSub = 100, addDosing = TRUE)
    # sim$DV <- sim$sim
    # sim$ID <- sim$sim.id
    # sim$sim.id <- NULL

    fit <- nlmixr(one.cmpt.combinederr, nlmixr2data::theo_sd, est = "focei")
    derv <- getDerv(fit)

    sum(grepl("O_ETA\\d+", names(derv))) |> expect_equal(ncol(fit$eta) - 1)
    all(derv$D_ResVar == 1) |> expect_equal(TRUE)


})

test_that("linearize wrapper", {
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

    fitLin <- linearize(fit)

    diffObjF <- abs((fit$objDf$OBJF - fitLin$objDf$OBJF)/ fitLin$objDf$OBJF)
    expect_true(diffObjF < 0.1)

    all.equal(fit$omega, fitLin$omega, tolerance = 0.1) |> expect_true()

})


## test_that("Linearize multiple endpoints ", {
## pk.turnover.emax3 <- function() {
##     ini({
##     tktr <- log(1)
##     tka <- log(1)
##     tcl <- log(0.1)
##     tv <- log(10)
##     ##
##     eta.ktr ~ 1
##     eta.ka ~ 1
##     eta.cl ~ 2
##     eta.v ~ 1
##     prop.err <- 0.1
##     pkadd.err <- 0.1
##     ##
##     temax <- logit(0.8)
##     tec50 <- log(0.5)
##     tkout <- log(0.05)
##     te0 <- log(100)
##     ##

##     eta.emax ~ .5
##     eta.ec50  ~ .5
##     eta.kout ~ .5
##     eta.e0 ~ .5
##     ##
##     pdadd.err <- 10
##     })
##     model({
##     ktr <- exp(tktr + eta.ktr)
##     ka <- exp(tka + eta.ka)
##     cl <- exp(tcl + eta.cl)
##     v <- exp(tv + eta.v)
##     emax = expit(temax+eta.emax)
##     ec50 =  exp(tec50 + eta.ec50)
##     kout = exp(tkout + eta.kout)
##     e0 = exp(te0 + eta.e0)
##     ##
##     DCP = center/v
##     PD=1-emax*DCP/(ec50+DCP)
##     ##
##     effect(0) = e0
##     kin = e0*kout
##     ##
##     d/dt(depot) = -ktr * depot
##     d/dt(gut) =  ktr * depot -ka * gut
##     d/dt(center) =  ka * gut - cl / v * center
##     d/dt(effect) = kin*PD -kout*effect
##     ##
##     cp = center / v
##     cp ~ prop(prop.err) + add(pkadd.err)
##     effect ~ add(pdadd.err) | pca
##     })
## }

## fit <- nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, est = "focei")

## derv <- getDerv(fit)

## sum(grepl("O_ETA\\d+", names(derv))) |> expect_equal(ncol(fit$eta) - 1)
## all(derv$D_ResVar == 1) |> expect_equal(TRUE)

## lmod <- linModGen(fit)

## lmod <- function(){
##     ini({
##         eta.ktr ~ 1
##         eta.ka ~ 1
##         eta.cl ~ 2
##         eta.v ~ 1

##         eta.emax ~ .5
##         eta.ec50  ~ .5
##         eta.kout ~ .5
##         eta.e0 ~ .5

##         prop.err <- 0.1
##         pkadd.err <- 0.1

##         pdadd.err <- 10
##     })

##     model({

##     })
## )

## })
