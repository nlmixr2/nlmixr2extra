test_that("precondition tests", {

  one.compartment <- function() {
    ini({
      tka <- 0.45 ; label("Log Ka")
      tcl <- 1 ; label("Log Cl")
      tv <- 3.45 ; label("Log V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
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

  fit2 <-
    suppressMessages(suppressWarnings(
      nlmixr(
        one.compartment, nlmixr2data::theo_sd,
        est = "focei",
        control = list(print = 0, eval.max = 200)
      )
    ))

  df1 <- fit2$parFixedDf
  cov1 <- fit2$cov

  ## Simply re-evaluate with no estimation (including inner estimation)
  suppressWarnings(preconditionFit(fit2, estType = "none"))

  df2 <- fit2$parFixedDf
  cov2 <- fit2$cov

  ## In this case there isn't a theta/omega estimate so these should be the same
  expect_equal(df1$Estimate, df2$Estimate)
  expect_equal(df1$`Back-transformed`, df2$`Back-transformed`)
  expect_equal(df1$`BSV(CV%)`, df2$`BSV(CV%)`)
  expect_equal(df1$`Shrink(SD)%`, df2$`Shrink(SD)%`)

  expect_false(isTRUE(all.equal(df1$SE, df2$SE)))
  expect_false(isTRUE(all.equal(df1$`%RSE`, df2$`%RSE`)))
  expect_false(isTRUE(all.equal(df1$`CI Lower`, df2$`CI Lower`)))
  expect_false(isTRUE(all.equal(df1$`%RSE`, df2$`%RSE`)))
  expect_false(isTRUE(all.equal(cov1, cov2)))

  skip_if_not(any(names(fit2$covList) == "r,s"))

  setCov(fit2, "r,s")

  df3 <- fit2$parFixedDf
  cov3 <- fit2$cov

  expect_equal(df1$Estimate, df3$Estimate)
  expect_equal(df1$`Back-transformed`, df3$`Back-transformed`)
  expect_equal(df1$`BSV(CV%)`, df3$`BSV(CV%)`)
  expect_equal(df1$`Shrink(SD)%`, df3$`Shrink(SD)%`)

  expect_equal(df1$SE, df3$SE)
  expect_equal(df1$`%RSE`, df3$`%RSE`)
  expect_equal(df1$`CI Lower`, df3$`CI Lower`)
  expect_equal(df1$`%RSE`, df3$`%RSE`)
  expect_equal(cov1, cov3)

  setCov(fit2, "precondition")
  df4 <- fit2$parFixedDf
  cov4 <- fit2$cov

  expect_equal(df2$Estimate, df4$Estimate)
  expect_equal(df2$`Back-transformed`, df4$`Back-transformed`)
  expect_equal(df2$`BSV(CV%)`, df4$`BSV(CV%)`)
  expect_equal(df2$`Shrink(SD)%`, df4$`Shrink(SD)%`)

  expect_equal(df2$SE, df4$SE)
  expect_equal(df2$`%RSE`, df4$`%RSE`)
  expect_equal(df2$`CI Lower`, df4$`CI Lower`)
  expect_equal(df2$`%RSE`, df4$`%RSE`)
  expect_equal(cov2, cov4)
})

test_that(".preCondModExtra handles dotted parameter names (#124)", {
  # symengine cannot parse an identifier containing a `.`, so building these
  # lines symbolically made preconditionFit() fail for a conventional residual
  # name like add.sd.  Build them by string assembly instead.
  pre <- matrix(c(1.5, -2.25,
                  0,    3.125), 2, 2, byrow = TRUE)

  expect_equal(
    .preCondModExtra(pre, c("tka", "add.sd")),
    paste0("tka=(1.5)*nlmixr2Pre_tka+(-2.25)*nlmixr2Pre_add.sd\n",
           "add.sd=(3.125)*nlmixr2Pre_add.sd")
  )

  # a negative coefficient must not produce `+-1.5`
  expect_false(grepl("+-", .preCondModExtra(pre, c("a", "b")), fixed = TRUE))

  # the generated lines must parse as R (and so as a model block)
  expect_silent(str2lang(paste0("{", .preCondModExtra(pre, c("tka", "add.sd")), "}")))

  # coefficients round-trip through text at full double precision
  expect_identical(as.numeric(.preCondNum(1/3)), 1/3)

  # an all-zero row still yields a valid line
  expect_equal(.preCondModExtra(matrix(0, 1, 1), "a"), "a=0")
})

test_that(".preCondExpand widens the preconditioner past the theta block (#124)", {
  # fit$R spans only the population parameters while fit$cov also carries the
  # omega elements, so the transform must be the identity off the theta block.
  pre <- matrix(c(2, 1,
                  0, 3), 2, 2, byrow = TRUE)
  covNames <- c("nlmixr2Pre_tka", "nlmixr2Pre_add.sd", "om.eta.ka")

  a <- .preCondExpand(pre, covNames, c("tka", "add.sd"))
  expect_equal(dim(a), c(3L, 3L))
  expect_equal(a[1:2, 1:2], pre)
  expect_equal(a[3, ], c(0, 0, 1))
  expect_equal(a[, 3], c(0, 0, 1))

  # order of the covariance is followed, not assumed
  covNames2 <- c("om.eta.ka", "nlmixr2Pre_add.sd", "nlmixr2Pre_tka")
  a2 <- .preCondExpand(pre, covNames2, c("tka", "add.sd"))
  expect_equal(a2[3, 3], pre[1, 1])
  expect_equal(a2[3, 2], pre[1, 2])
  expect_equal(a2[1, 1], 1)

  expect_error(.preCondExpand(pre, c("om.eta.ka", "nope"), c("tka", "add.sd")),
               "could not find the preconditioned parameters")
})
