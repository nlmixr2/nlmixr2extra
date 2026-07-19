test_that("add missing etas", {
  
  mod <- function() {
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
          # eta.ktr ~ 0.899432614108315
          # eta.ka ~ 1.20111175373864
          # eta.cl ~ 0.0756455519562788
          # eta.v ~ 0.0547982604779906
          # eta.emax ~ 0.547762672352332
          # eta.ec50 ~ 0.204313071604677
          # eta.kout ~ 0.0228239781890516
          # eta.e0 ~ 0.0107988390114995
      })
      model({
          ktr <- exp(tktr)
          ka <- exp(tka)
          cl <- exp(tcl)
          v <- exp(tv)
          emax = expit(temax, 0, 1)
          ec50 = exp(tec50)
          kout = exp(tkout)
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
  # no etas
  modupdate <- addAllEtas(mod)
  expect_true(length(modupdate$eta) == 8)
  
  # one etas 
  mod <- function() {
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
          # eta.ktr ~ 0.899432614108315
          # eta.ka ~ 1.20111175373864
          # eta.cl ~ 0.0756455519562788
          # eta.v ~ 0.0547982604779906
          eta.emax ~ 0.547762672352332
          eta.ec50 ~ 0.204313071604677
          # eta.kout ~ 0.0228239781890516
          # eta.e0 ~ 0.0107988390114995
      })
      model({
          ktr <- exp(tktr)
          ka <- exp(tka)
          cl <- exp(tcl)
          v <- exp(tv)
          emax = expit(temax + eta.emax, 0, 1)
          ec50 = exp(tec50 + eta.ec50)
          kout = exp(tkout)
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
  mod <- mod()
  modupdate <- addAllEtas(mod)
  expect_true(length(modupdate$eta) == 8)
  
  # all etas
  mod <- function() {
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
  
  modupdate <- addAllEtas(mod)
  expect_true(length(modupdate$eta) == 8)
})

test_that("iiv combinations", {
  x <- iivCombn(c("a", "b", "c"))
  expect_equal(length(x), 17)
})


test_that("linearized eta search", {
  skip_on_cran()
  one.cmpt.adderr <- function() {
    ini({
      tcl <- log(2.7) # Cl
      tv <- log(30) # V
      tka <- log(1.56) #  Ka
      eta.cl + eta.v ~ sd(cor(0.3, 0.99, 0.5))
      eta.ka ~ 0
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
  set.seed(42)
  ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
      rxode2::et(amt = 350, cmt = "depot") |>
      rxode2::et(time = c(0.25, 0.5, 1,2,3,6,8,12,16,24))
  sim <- rxode2::rxSolve(one.cmpt.adderr, ev, nSub = 200, addDosing = TRUE)
  # plot(sim)
  sim$dv <- sim$sim
  sim$id <- sim$sim.id
  sim$sim.id <- NULL
  sim <- sim[,c("id", "time", "amt", "dv", "evid")]
  

  one.cmpt.adderr <- function() {
    ini({
        tcl <- log(2.7) # Cl
        tv <- log(30) # V
        tka <- log(1.56) #  Ka
        # eta.cl ~ 0.5
        # eta.v ~  0.4
        # eta.ka ~ 0.6
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
  # one.cmpt.adderr <- addAllEtas(one.cmpt.adderr)
  fit <- nlmixr(one.cmpt.adderr, sim, est = "focei")
  
  suppressWarnings(
    fitLin <- linearize(fit, addEtas = TRUE, focei = TRUE)
  )
  isLinearizeMatch(fitLin, 0.2)
  linearizePlot(fitLin)
  
  res <- iivSearch(fitLin)
  expect_true(inherits(res, "linIIVSearch"))
  res$summary[order(res$summary$BIC),]
  
  resLast <- rerunTopN(res)
  resLast$summary  |> expect_no_error()

})
