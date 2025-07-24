test_that("linearized eta search", {
    one.cmpt.adderr <- function() {
  ini({
    tcl <- log(2.7) # Cl
    tv <- log(30) # V
    tka <- log(1.56) #  Ka
    eta.cl + eta.v ~ c(0.5, 0.1, 0.3)
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
ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
    rxode2::et(amt = 350, cmt = "depot") |>
    rxode2::et(time = c(0.25, 0.5, 1,2,3,6,8,12,16,24))
sim <- rxode2::rxSolve(one.cmpt.adderr, ev, nSub = 200, addDosing = TRUE)
plot(sim)
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
fit <- addFixedEtas(one.cmpt.adderr)
expect_false(hasUnFixedEta(fit))
fit <- nlmixr(fit, sim, est = "focei")

fitLin <- linearize(fit)
res <- iivSearch(fitLin)

})