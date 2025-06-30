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
    all(derv$D_EPSETA_1_1 == 0) |> expect_equal(TRUE)
    all(derv$D_EPSETA_1_2 == 0) |> expect_equal(TRUE)
    all(derv$D_EPSETA_1_3 == 0) |> expect_equal(TRUE)
})


test_that("Linearize prop err model ", {
    one.cmpt.properr <- function() { # non-linear base
        ini({
            tka <- log(1.56) # Ka
            tcl <- log(2.7) # Cl
            tv <- log(30) # V
            eta.ka ~ 0.6
            eta.cl ~ 0.3
            eta.v ~ 0.1
            prop.sd <- 0.10
        })
        model({
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            v <- exp(tv + eta.v)
            d / dt(depot) <- -ka * depot
            d / dt(center) <- ka * depot - cl / v * center
            cp <- center / v
            cp ~ prop(prop.sd)
        })
    }
    ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
        rxode2::et(amt = 10000, cmt = "depot")
    sim <- rxSolve(one.cmpt.adderr, ev, nSub = 100, addDosing = TRUE)
    sim$DV <- sim$sim
    sim$ID <- sim$sim.id
    sim$sim.id <- NULL

    fit <- nlmixr(one.cmpt.properr, sim, est = "focei")
    derv <- getDerv(fit)
    
    sum(grepl("O_ETA\\d+", names(derv))) |> expect_equal(ncol(fit$eta) - 1)
    all(derv$D_ResVar == 1) |> expect_equal(TRUE)

})