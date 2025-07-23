test_that("linearized eta search", {
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
    


    fitLin <- linearize(fit)
    
    # start: model with no eta ==> stop if has eta
    # add random effects with small var
    fitLin <- linearizeVar(model)
    nlmixr2lib::addEta()

    iivSearch(fitLin)
})