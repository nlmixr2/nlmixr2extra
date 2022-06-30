


test_that("Testing out a .getThetaNameif a parameter is present", {
  ## basic example code
  ## The basic model consiss of an ini block that has initial estimates
  one.cmt <- function() {
    ini({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- log(c(0, 2.7, 100)) # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45; label("log V")
      ## the label("Label name") works with all models
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
  
  ###
  ###
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
    # and a model block with the error specification and model specification
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - exp(tcl + eta.cl) / v * center
      cp = center / v
      cp ~ add(add.sd)
    })
  }
  
  
  tainted  <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45    # Log V
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    # and a model block with the error specification and model specification
    model({
      ka <- exp(tka + eta.ka)
      v <- exp(tv + eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - exp(tcl + eta.cl) / v * center
      cp = center / v
      cp ~ add(add.sd)
    })
  }
  
  
  
  ui1 <- nlmixr2(one.compartment)
  ui2 <- nlmixr2(one.cmt)
  ui3 <- nlmixr2(tainted)
  
  
  
  expect_equal(.getThetaName(ui1,'ka'),'tka') 
  expect_equal(.getThetaName(ui1,'tka'),'tka') 
  expect_equal(.getThetaName(ui1,'v'),'tv')
  expect_equal(.getThetaName(ui2,'ka'),'tka') 
  expect_equal(.getThetaName(ui2,'tka'),'tka') 
  expect_error(.getThetaName(ui2,'gamma'),'gamma')
  expect_equal(.getThetaName(ui3,'tcl'),'tcl') 
  expect_error(.getThetaName(ui3,'cl'),'cl')
  
  
  ##.addCovariate
  
  expect_equal(as.character(.addCovariate(ui2,'ka','WT',norm=FALSE))[[1]],'ka <- exp(tka + eta.ka + cov_WT_ka * WT)')
  expect_error(as.character(.addCovariate(ui3,'cl','WT',norm=FALSE))[[1]],'ka <- exp(tka + eta.ka + cov_WT_ka * WT)')
  expect_equal(as.character(.addCovariate(ui1,'cl','WT',norm=FALSE))[[2]],'cl <- exp(tcl + eta.cl + cov_WT_cl * WT)')
  expect_equal(as.character(.addCovariate(ui2,'ka','normalized_WT',norm=TRUE,addfactor=TRUE,factor=1))[[1]],'ka <- exp(tka + eta.ka + cov_WT_ka * normalized_WT * 1)')
  
  
})



