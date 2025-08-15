# .libPaths(c("~fidlema3/src/nlmixr2/lib", .libPaths()))

library(clustermq)
library(tidyverse)

devtools::load_all("~/Documents/nlmixr2extra/")
# library(nlmixr2data)
library(nlmixr2)

testfun <- function(x, fix, focei){
  print(x)
  
  .libPaths(c("~/Documents/nlmixr2extra/lib", "~fidlema3/src/nlmixr2/lib", .libPaths()))
  devtools::load_all("~/Documents/nlmixr2extra")
  library(nlmixr2, lib.loc = "~fidlema3/src/nlmixr2/lib")
  
  # library(nlmixr2extra, lib.loc = "~/Documents/nlmixr2extra/lib")
  
  
  one.cmpt.adderr <- function() {
    ini({
      tcl <- log(2.7) # Cl
      tv <- log(30) # V
      tka <- log(1.56) #  Ka
      eta.cl + eta.v ~ sd(cor(0.3, 0.99, 0.5))
      # eta.ka ~ 0.1
      add.sd <- 0.7
      prop.sd <- 0.2
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka*depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd) + combined2()
    })
  }
  
  set.seed(NULL)
  ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
      rxode2::et(amt = 300, cmt = "depot") |>
      # rxode2::et(amt = 200) |>
      rxode2::et(time = c(0.25, 0.5, 1,2,3,6,8,12,16,24))
  sim <- rxode2::rxSolve(one.cmpt.adderr, ev, nSub = 200, addDosing = TRUE)
  plot(sim)
  sim$dv <- sim$sim
  sim$id <- sim$sim.id
  sim$sim.id <- NULL
  sim <- sim[,c("id", "time", "amt", "dv", "evid")]
  
  
  ofit <- nlmixr2est::nlmixr(one.cmpt.adderr,sim, est = "focei")
  
  
  one.cmpt.adderr <- function() {
    ini({
      tcl <- log(2.7) # Cl
      tv <- log(30) # V
      tka <- log(1.56) #  Ka
      # eta.cl ~ 0.5
      # eta.v ~  0.4
      # eta.ka ~ 0.6
      add.sd <- 0.7
      prop.sd <- 0.2
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka*depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) + prop(prop.sd) + combined2()
    })
  }
  
  # fit <- nlmixr2extra::addAllEtas(one.cmpt.adderr, fix = fix)
  fit <- nlmixr2est::nlmixr(one.cmpt.adderr, sim, est = "focei")
  
  fitLin <- linearize(fit, focei=focei, addEtas = TRUE)
  res <- iivSearch(fitLin)
  
  res$summary |> dplyr::arrange(OBJF) 
  res$summary |> dplyr::arrange(BIC) 
  
  passed <- res$summary #[res$summary$covMethod=="r,s",] # FIXME sometimes it will break, sometimes not
  passed |> dplyr::arrange(BIC)
  # expect_true("etaTcl-etaTv-etaTcl,etaTv" %in% passed$search[order(passed$BIC)][1:3])
  
  # top1 <- c(top1,  which("etaTcl-etaTv-etaTcl,etaTv" == passed$search[order(passed$BIC)]))
  
  res$finalFit$parFixed
  
  res$res[[5]]$fit$iniDf
  res$res[[5]]$fit$parFixed
  res$res[[5]]$search
  
  res$finalFit # FIXME
  
  list(rank = which("etaTcl+etaTv+etaTcl~etaTv" == passed$search[order(passed$BIC)]), 
       matched = isLinearizeMatch(fitLin)$ofv[[1]])
  
}
onetst <- testfun(1,TRUE, TRUE)
options(clustermq.ssh.timeout = 30) # in seconds
options(
  clustermq.scheduler = "lsf",
  clustermq.template = "./lsf.tmpl"
)


res_fixed <- Q(testfun,x=1:100,n_jobs=100, 
               const= list(fix = TRUE, focei=TRUE), 
               verbose = TRUE, fail_on_error = TRUE, log_worker = TRUE, timeout = Inf)


res_fixed_foce <- Q(testfun,x=1:100,n_jobs=100, 
               const= list(fix = TRUE, focei=FALSE), 
               verbose = TRUE, fail_on_error = TRUE, log_worker = TRUE, timeout=Inf)
# saveRDS(res_fixed, "sim_res_top_fixed.RDS")


res_unfixed <- Q(testfun,x=1:100,n_jobs=100, const= list(fix = FALSE))
# saveRDS(res_unfixed, "sim_res_top_unfixed.RDS")

# top1 <- readRDS("sim_res_top_fixed.RDS")
# top1 <- readRDS("sim_res_top_unfixed.RDS")

sum(sapply(res_fixed, \(x)x[["matched"]]))
sum(sapply(res_unfixed, \(x)x[["matched"]]))

sum(top1>=3)
ggplot(data.frame(rank = as.factor(unlist(lapply(res_fixed, \(x) x[["rank"]]))))) + 
  geom_bar(aes(x = rank, y =  (..count..)/sum(..count..))) + 
  labs(y = "Frequency", title = "FOCE+I Linearization", x = "True BSV Rank") + theme_classic()


ggplot(data.frame(rank = as.factor(unlist(lapply(res_fixed_foce, \(x) x[["rank"]]))))) + 
  geom_bar(aes(x = rank, y =  (..count..)/sum(..count..))) + 
  labs(y = "Frequency", title = "FOCE Linearization", x = "True BSV Rank") + theme_classic()


data_cumulative <- data.frame(
    iteration = c(1:100, 1:100),
    rank = c(unlist(lapply(res_fixed_foce, \(x) x[["rank"]])), unlist(lapply(res_fixed, \(x) x[["rank"]]))),
    method = sort(rep(c("FOCEI", "FOCE"), times = 100))
  ) |> 
  arrange(method, rank) |>
  group_by(method) |>
  mutate(cumulative_fraction = cumsum(rank) / sum(rank))  

ggplot(data_cumulative, aes(fill = method, x = rank)) + 
  geom_bar(aes(y =  (..count..)/sum(..count..)))  + 
  labs(y = "Frequency", title = "True BSV Rank", x = "True BSV Rank") + 
  scale_fill_manual(values = c("FOCEI" = "#D73027", "FOCE" = "#4575B4")) +
  facet_grid(~method)  + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  scale_x_continuous(breaks = 1:10)

ggplot(data_cumulative, aes(x = rank, y = cumulative_fraction, color = method)) +
  geom_line(size = 1) +
  labs(
    title = "Cumulative Rank",
    x = "Rank",
    y = "Cumulative Fraction"
  ) +
  theme_minimal() +
  theme(legend.position = "none") + 
  scale_x_continuous(breaks = 1:10) +
  geom_hline(yintercept=0.80, linetype = 3) + 
  scale_color_manual(values = c("FOCEI" = "#D73027", "FOCE" = "#4575B4")) +
  theme(legend.title = element_blank())


###

ggplot(data.frame(rank = as.factor(sapply(res_unfixed, \(x) x[["rank"]])))) + 
  geom_bar(aes(x = rank, y =  (..count..)/sum(..count..))) + 
  labs(y = "freq", title = "UnFixed")

# Case linear vs non

one.cmpt.adderr <- function() {
  ini({
    tcl <- log(2.7) # Cl
    tv <- log(30) # V
    tka <- log(1.56) #  Ka
    eta.cl + eta.v ~ sd(cor(0.3, 0.99, 0.5))
    # eta.ka ~ 0.1
    add.sd <- 0.7
    prop.sd <- 0.2
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka*depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd) + prop(prop.sd) + combined2()
  })
}

set.seed(42)
ev <- rxode2::et(amountUnits = "mg", timeUnits = "hours") |>
    rxode2::et(amt = 300, cmt = "depot") |>
    # rxode2::et(amt = 200) |>
    rxode2::et(time = c(0.25, 0.5, 1,2,3,6,8,12,16,24))
sim <- rxode2::rxSolve(one.cmpt.adderr, ev, nSub = 200, addDosing = TRUE)
plot(sim)
sim$dv <- sim$sim
sim$id <- sim$sim.id
sim$sim.id <- NULL
sim <- sim[,c("id", "time", "amt", "dv", "evid")]


ofit <- nlmixr2est::nlmixr(one.cmpt.adderr,sim, est = "focei")


one.cmpt.adderr <- function() {
  ini({
    tcl <- log(2.7) # Cl
    tv <- log(30) # V
    tka <- log(1.56) #  Ka
    # eta.cl ~ 0.5
    # eta.v ~  0.4
    # eta.ka ~ 0.6
    prop.sd <- 0.2
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl)
    v <- exp(tv)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka*depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd) + prop(prop.sd) + combined2()
  })
}

fit <- nlmixr2est::nlmixr(one.cmpt.adderr, sim, est = "focei")

fitLin <- linearize(fit, focei= T, addEtas = TRUE)

fitLin$env$originalUi 
fitLin$env$origData


res <- iivSearch(fitLin)
res$summary |> arrange(BIC)  |> select(4, 1, 3, 5, 6) |> 
    mutate(covMethod = ifelse(covMethod == "r,s", "success", "fail")) |> 
    rename(covariance = covMethod) |> gt::gt()




x <- rerunTopN(res, 7)
x$summary |> arrange(O.BIC) |> select(7, 1,3) |> gt()

# - [x] additive 
# - [x] combined2 



#### Test Residual Model Search 

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

fit.add <- nlmixr(one.cmpt.adderr, nlmixr2data::theo_sd, est = "focei")

fit.prop <-  one.cmpt.adderr |> model(cp ~ prop(prop.sd)) |> 
              nlmixr(nlmixr2data::theo_sd, est = "focei")

fit.combined2 <- one.cmpt.adderr |> model(cp ~ add(add.sd) + prop(prop.sd) + combined2()) |> 
                nlmixr(nlmixr2data::theo_sd, est = "focei")

fit.combined1 <- one.cmpt.adderr |> model(cp ~ add(add.sd) + prop(prop.sd) + combined1()) |> 
                nlmixr(nlmixr2data::theo_sd, est = "focei")
        

# deriv <- getDeriv(fit.add)

linFit <- linearize(fit.add, focei = FALSE)



linFit.prop <- linFit |>  
              model(rxR2 <- (prop.sd^2*OPRED^2)) |> ini(prop.sd <- 0.1) |> 
              nlmixr(getData(linFit), est = "focei")
  
linFit.combined2 <- linFit |> 
                  model(rxR2 <- (prop.sd^2*OPRED^2 + add.sd^2)) |> ini(prop.sd <- 0.1) |>
                  nlmixr(getData(linFit), est = "focei")

linFit.combined1 <- linFit |> 
                    model(rxR2 <- (prop.sd*OPRED + add.sd)^2) |> ini(prop.sd <- 0.1) |>
                    nlmixr(getData(linFit), est = "focei")


c(linFit$objDf$OBJF,  linFit.prop$objDf$OBJF,  linFit.combined1$objDf$OBJF, linFit.combined2$objDf$OBJF)
c(fit.add$objDf$OBJF, fit.prop$objDf$OBJF,  fit.combined1$objDf$OBJF, fit.combined2$objDf$OBJF)
                    