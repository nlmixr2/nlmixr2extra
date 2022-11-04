wd <- getwd()
setwd(devtools::package_file())

one.compartment <- function() {
  ini({
    tka <- 0.45 # Log Ka
    tcl <- 1 # Log Cl
    tv <- 3.45    # Log V
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  # and a model block with the error sppecification and model specification
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(depot) = -ka * depot
    d/dt(center) = ka * depot - cl / v * center
    cp = center / v
    cp ~ add(add.sd)
  })
}

 
theoFitOde <-  nlmixr2est::nlmixr(one.compartment, nlmixr2data::theo_sd, est="focei")

if (file.exists("data/theoFitOde.rda")) unlink("data/theoFitOde.rda")

save(theoFitOde, file="data/theoFitOde.rda", compress="bzip2", version=2, ascii=FALSE)

setwd(wd)
