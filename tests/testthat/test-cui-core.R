# caches default to the working directory; keep them out of the package
withr::local_dir(withr::local_tempdir())

.cur <- loadNamespace("nlmixr2extra")

# CUI-MET (Zhang et al. 2025) Table 2: modeled marginal probabilities, the
# utility mean (UM) and the two weighted means (UWM) for three examples
.cuimetTable2 <- list(
  ex1 = list(
    tox = c(0.038, 0.058, 0.106, 0.232, 0.533),
    eff = c(0.022, 0.119, 0.343, 0.568, 0.681),
    tol = c(0.147, 0.206, 0.280, 0.368, 0.466),
    um = c(0.377, 0.422, 0.506, 0.568, 0.538),
    uwm1 = c(0.248, 0.310, 0.434, 0.548, 0.574),
    uwm2 = c(0.367, 0.409, 0.477, 0.528, 0.509),
    w2 = c(tox = 0.3, eff = 0.2, tol = 0.5),
    best = c(um = "4", uwm1 = "5", uwm2 = "4")
  ),
  ex2 = list(
    tox = c(0.095, 0.149, 0.224, 0.323, 0.442),
    eff = c(0.414, 0.568, 0.618, 0.642, 0.657),
    tol = c(0.285, 0.653, 0.760, 0.804, 0.828),
    um = c(0.534, 0.691, 0.718, 0.708, 0.681),
    uwm1 = c(0.473, 0.650, 0.692, 0.698, 0.688),
    uwm2 = c(0.584, 0.698, 0.710, 0.688, 0.652),
    w2 = c(tox = 0.4, eff = 0.4, tol = 0.2),
    best = c(um = "3", uwm1 = "4", uwm2 = "3")
  ),
  ex3 = list(
    tox = c(0.105, 0.152, 0.233, 0.374, 0.583),
    eff = c(0.223, 0.387, 0.470, 0.442, 0.312),
    tol = c(0.705, 0.607, 0.500, 0.393, 0.295),
    um = c(0.607, 0.614, 0.579, 0.487, 0.341),
    uwm1 = c(0.502, 0.545, 0.538, 0.464, 0.328),
    uwm2 = c(0.684, 0.660, 0.601, 0.496, 0.347),
    w2 = c(tox = 0.4, eff = 0.2, tol = 0.4),
    best = c(um = "2", uwm1 = "2", uwm2 = "1")
  )
)

.cuimetEndpoints <- list(cuiEndpoint(tox, direction = "lower"),
                         cuiEndpoint(eff), cuiEndpoint(tol))

.cuimetValues <- function(ex) {
  data.frame(regimen = rep(1:5, 3),
             endpoint = rep(c("tox", "eff", "tol"), each = 5),
             value = c(ex$tox, ex$eff, ex$tol))
}

test_that("CUI-MET Table 2 UM and UWM are reproduced", {
  # inputs and outputs are rounded to 3 decimals in the paper
  for (.nm in names(.cuimetTable2)) {
    .ex <- .cuimetTable2[[.nm]]
    .um <- cuiFromValues(.cuimetValues(.ex), endpoints = .cuimetEndpoints)
    expect_equal(.um$summary$est, .ex$um, tolerance = 1.1e-3, ignore_attr = TRUE)
    expect_equal(.um$summary$regimen[.um$summary$best], .ex$best[["um"]])
    .w1 <- update(.um, weights = c(tox = 0.2, eff = 0.5, tol = 0.3))
    expect_equal(.w1$summary$est, .ex$uwm1, tolerance = 1.1e-3, ignore_attr = TRUE)
    expect_equal(.w1$summary$regimen[.w1$summary$best], .ex$best[["uwm1"]])
    .w2 <- update(.um, weights = .ex$w2)
    expect_equal(.w2$summary$est, .ex$uwm2, tolerance = 1.1e-3, ignore_attr = TRUE)
    expect_equal(.w2$summary$regimen[.w2$summary$best], .ex$best[["uwm2"]])
  }
})

test_that("weights are normalized and validated", {
  .v <- .cuimetValues(.cuimetTable2$ex1)
  .a <- cuiFromValues(.v, .cuimetEndpoints, weights = c(tox = 2, eff = 5, tol = 3))
  .b <- cuiFromValues(.v, .cuimetEndpoints, weights = c(tol = 0.3, tox = 0.2, eff = 0.5))
  expect_equal(.a$summary$est, .b$summary$est)
  expect_equal(unname(.a$weightsNorm), c(0.2, 0.5, 0.3))
  expect_error(cuiFromValues(.v, .cuimetEndpoints, weights = c(tox = 1, eff = 1)), "exactly")
  expect_error(cuiFromValues(.v, .cuimetEndpoints, weights = c(1, 1, 1)), "exactly")
  expect_error(cuiFromValues(.v, .cuimetEndpoints, weights = c(tox = 0, eff = 0, tol = 0)),
               "positive total")
  expect_error(cuiFromValues(.v, .cuimetEndpoints, weights = c(tox = -1, eff = 1, tol = 1)))
  expect_error(cuiFromValues(.v, .cuimetEndpoints, weights = c(tox = NA, eff = 1, tol = 1)))
  # a zero weight drops the endpoint
  .z <- cuiFromValues(.v, .cuimetEndpoints, weights = c(tox = 0, eff = 1, tol = 1))
  expect_equal(.z$summary$est, (.cuimetTable2$ex1$eff + .cuimetTable2$ex1$tol) / 2)
})

test_that("utilities and direction", {
  .lin <- cuiLinear(0, 50)
  expect_equal(.cur$.cuiApplyUtility(.lin, "higher", c(-10, 0, 25, 50, 80), "x"),
               c(0, 0, 0.5, 1, 1))
  # direction is applied after clamping
  expect_equal(.cur$.cuiApplyUtility(.lin, "lower", c(-10, 0, 25, 50, 80), "x"),
               c(1, 1, 0.5, 0, 0))
  expect_equal(.cur$.cuiApplyUtility(cuiStep(-1), "higher", c(-2, -1, 0), "x"), c(0, 1, 1))
  expect_equal(.cur$.cuiApplyUtility(cuiStep(-1, "below"), "higher", c(-2, -1, 0), "x"), c(1, 1, 0))
  expect_equal(.cur$.cuiApplyUtility(cuiStep(-1, "below"), "lower", c(-2, -1, 0), "x"), c(0, 0, 1))
  expect_error(.cur$.cuiApplyUtility(cuiIdentity(), "higher", c(0.5, 1.2), "x"), "outside")
  expect_error(cuiLinear(1, 1), "less than")
  expect_error(cuiEndpoint(x, utility = function(x) x), "utility")
  expect_equal(format(cuiStep(2, "below")), "step(<= 2)")
})

test_that("cuiEndpoint takes a symbol, a string or .name", {
  expect_equal(cuiEndpoint(resp)$name, "resp")
  expect_equal(cuiEndpoint("resp")$name, "resp")
  .nm <- "resp"
  expect_equal(cuiEndpoint(.name = .nm)$name, "resp")
  # a bare symbol is the name itself, never a variable's value
  expect_equal(cuiEndpoint(.nm)$name, ".nm")
  .wrap <- function(n) cuiEndpoint(.name = n, direction = "lower")
  expect_equal(.wrap("tox")$name, "tox")
  expect_equal(.wrap("tox")$direction, "lower")
  expect_error(cuiEndpoint(), "required")
  expect_error(cuiEndpoint(resp, .name = "resp"), "only one")
  expect_error(cuiEndpoint(a + b), "name or a single string")
  expect_error(cuiEndpoint(c("a", "b")), "name or a single string")
})

test_that("geometric combine handles zero utilities and zero weights", {
  .v <- data.frame(regimen = rep(c("a", "b", "c"), 2),
                   endpoint = rep(c("e1", "e2"), each = 3),
                   value = c(0.5, 0, 0.8, 0.5, 0.9, 0))
  .g <- cuiFromValues(.v, combine = "geometric")
  expect_equal(.g$summary$est, c(0.5, 0, 0))
  expect_false(anyNA(.g$summary$est))
  # zero-weight endpoint with zero utility does not take part
  .g2 <- update(.g, weights = c(e1 = 1, e2 = 0))
  expect_equal(.g2$summary$est, c(0.5, 0, 0.8))
  .g3 <- update(.g, weights = c(e1 = 1, e2 = 3))
  expect_equal(.g3$summary$est[1], 0.5)
  expect_equal(update(.g, combine = "additive")$summary$est, c(0.5, 0.45, 0.4))
})

# a small dataset with uncertainty draws
.drawValues <- function() {
  .d <- expand.grid(regimen = c("lo", "mid", "hi"), draw = 0:4,
                    endpoint = c("e1", "e2"), stringsAsFactors = FALSE)
  .base <- c(lo = 0.2, mid = 0.5, hi = 0.8)
  .shift <- c(0, 0.3, -0.3, 0.1, -0.1)
  .d$value <- ifelse(.d$endpoint == "e1",
                     .base[.d$regimen] + 0.1 * .shift[.d$draw + 1],
                     1 - .base[.d$regimen])
  .d
}

test_that("point estimate, intervals and selection frequencies", {
  .r <- cuiFromValues(.drawValues(), regimens = cuiRegimens(dose = 1:3, label = c("lo", "mid", "hi")),
                      weights = c(e1 = 2, e2 = 1))
  .s <- .r$summary
  expect_equal(.s$regimen, c("lo", "mid", "hi"))
  expect_equal(.s$est, c(2 * 0.2 + 0.8, 2 * 0.5 + 0.5, 2 * 0.8 + 0.2) / 3)
  expect_equal(.s$regimen[.s$best], "hi")
  expect_equal(.s$selFreq, c(0, 0, 1))
  expect_true(all(.s$lower <= .s$est & .s$est <= .s$upper))
  # the point estimate does not depend on which draws exist
  .only0 <- .drawValues()
  .only0 <- .only0[.only0$draw == 0, ]
  .r0 <- cuiFromValues(.only0, weights = c(e1 = 2, e2 = 1))
  expect_equal(.r0$summary$est[match(.s$regimen, .r0$summary$regimen)], .s$est)
  expect_true(all(is.na(.r0$summary$selFreq)))
  expect_true(all(is.na(.r0$summary$mostSelected)))
  # with equal weights every regimen ties in every draw
  .eq <- update(.r, weights = c(e1 = 1, e2 = 1))
  expect_equal(sum(.eq$summary$selFreq), 1)
})

test_that("ties split a draw equally", {
  .v <- data.frame(regimen = rep(c("a", "b", "c"), 2), draw = rep(c(0L, 1L), each = 3),
                   endpoint = "e", value = c(0.5, 0.5, 0.1, 0.7, 0.7, 0.2))
  .r <- cuiFromValues(.v)
  expect_equal(.r$summary$selFreq, c(0.5, 0.5, 0))
  expect_equal(.r$summary$best, c(TRUE, TRUE, FALSE))
  expect_equal(.r$summary$mostSelected, c(TRUE, TRUE, FALSE))
})

test_that("eligibility restricts the recommendation and recomputes selection", {
  .r <- cuiFromValues(.drawValues(), weights = c(e1 = 2, e2 = 1))
  .e <- update(.r, eligible = c("lo", "mid"))
  expect_equal(.e$summary$regimen[.e$summary$best], "mid")
  # recomputed among eligible regimens, not filtered from the unrestricted ones
  expect_equal(.e$summary$selFreq[.e$summary$regimen == "mid"], 1)
  expect_true(is.na(.e$summary$selFreq[.e$summary$regimen == "hi"]))
  .l <- update(.r, eligible = c(lo = TRUE, mid = TRUE, hi = FALSE))
  expect_equal(.l$summary, .e$summary)
  expect_error(update(.r, eligible = "nope"), "unknown")
  expect_error(update(.r, eligible = c(lo = TRUE)), "every regimen")
  expect_error(update(.r, eligible = character(0)), "no eligible")
})

test_that("assessments with eligibility per assessment", {
  .v <- rbind(cbind(.drawValues(), assessment = "m3"),
              cbind(.drawValues(), assessment = "m6"))
  .r <- cuiFromValues(.v, weights = c(e1 = 2, e2 = 1),
                      eligible = data.frame(regimen = c("lo", "mid", "hi", "lo", "mid"),
                                            assessment = c("m3", "m3", "m3", "m6", "m6")))
  .s <- .r$summary
  expect_equal(.s$regimen[.s$best & .s$assessment == "m3"], "hi")
  expect_equal(.s$regimen[.s$best & .s$assessment == "m6"], "mid")
  expect_equal(.r$assessments, c("m3", "m6"))
})

test_that("incomplete or duplicated grids are errors", {
  .v <- .drawValues()
  expect_error(cuiFromValues(.v[-1, ]), "incomplete")
  expect_error(cuiFromValues(rbind(.v, .v[1, ])), "more than one")
  expect_error(cuiFromValues(.v[.v$draw > 0, ]), "draw 0")
  .na <- .v
  .na$value[3] <- NA
  expect_error(cuiFromValues(.na), "missing")
  expect_error(cuiFromValues(.v, endpoints = list(cuiEndpoint(e1))), "do not match")
  expect_error(cuiFromValues(.v, regimens = cuiRegimens(dose = 1:2, label = c("lo", "mid"))),
               "not in 'regimens'")
  # an assessment missing a draw is incomplete
  .a <- rbind(cbind(.v, assessment = "m3"), cbind(.v[.v$draw != 4, ], assessment = "m6"))
  expect_error(cuiFromValues(.a), "incomplete")
})

test_that("failed draws are excluded from the selection denominator", {
  .v <- .drawValues()
  .v <- .v[.v$draw != 4, ]
  .draws <- data.frame(draw = 0:4, status = c(rep("ok", 4), "failed"),
                       source = "test", message = "", stringsAsFactors = FALSE)
  .r <- .cur$.cuiNew(.v, endpoints = NULL, regimens = NULL, weights = c(e1 = 2, e2 = 1),
                combine = "additive", eligible = NULL, level = 0.95,
                tieTol = sqrt(.Machine$double.eps), draws = .draws, settings = list())
  expect_equal(sum(.r$summary$selFreq), 1)
  expect_output(print(.r), "3 ok, 1 failed")
  .bad <- .draws
  .bad$status[2] <- "failed"
  expect_error(.cur$.cuiNew(.v, NULL, NULL, NULL, "additive", NULL, 0.95, 1e-8, .bad, list()),
               "must not have values")
})

test_that("update() equals building from scratch and keeps unchanged settings", {
  .v <- .drawValues()
  .r <- cuiFromValues(.v, weights = c(e1 = 2, e2 = 1), level = 0.8)
  .u <- update(.r, weights = c(e1 = 1, e2 = 3))
  .d <- cuiFromValues(.v, weights = c(e1 = 1, e2 = 3), level = 0.8)
  expect_equal(.u$summary, .d$summary)
  expect_equal(.u$level, 0.8)
  expect_equal(update(.u, weights = NULL)$summary, cuiFromValues(.v, level = 0.8)$summary)
})

test_that("regimens, printing and data frames", {
  .reg <- cuiRegimens(dose = c(150, 200), ii = c(24, 12), units = "mg")
  expect_equal(.reg$label, c("150 q24", "200 q12"))
  expect_s3_class(.reg, "cuiRegimens")
  expect_error(cuiRegimens(dose = c(1, 1)), "unique")
  .r <- cuiFromValues(.drawValues())
  expect_output(print(.r), "Recommended")
  expect_output(print(.r), "Most frequently selected")
  expect_equal(nrow(as.data.frame(.r)), 3)
  expect_equal(nrow(as.data.frame(.r, type = "cui")), 15)
  expect_equal(nrow(as.data.frame(.r, type = "utility")), 30)
  expect_output(print(cuiEndpoint(tox, direction = "lower")), "lower is better")
  expect_output(print(cuiLinear(0, 1)), "linear")
  # numeric regimens are ordered by dose
  .n <- cuiFromValues(data.frame(regimen = c(10, 2, 5), endpoint = "e", value = c(0.1, 0.2, 0.3)))
  expect_equal(.n$summary$regimen, c("2", "5", "10"))
  # the result carries no environments, so it serializes compactly
  expect_false(any(vapply(unlist(unclass(.r), recursive = TRUE), is.environment, logical(1))))
})

# draws where the winner changes, so misaligned draws or regimens show up
.varyValues <- function() {
  .cui <- rbind(c(0.5, 0.6, 0.4),  # draw 0: mid
                c(0.7, 0.6, 0.4),  # draw 1: lo
                c(0.5, 0.6, 0.4),  # draw 2: mid
                c(0.5, 0.6, 0.9),  # draw 3: hi
                c(0.5, 0.8, 0.4))  # draw 4: mid
  data.frame(regimen = rep(c("lo", "mid", "hi"), 5), draw = rep(0:4, each = 3),
             endpoint = "e", value = as.vector(t(.cui)))
}

test_that("selection follows the winner of each draw", {
  .reg <- cuiRegimens(dose = 1:3, label = c("lo", "mid", "hi"))
  .r <- cuiFromValues(.varyValues(), regimens = .reg)
  expect_equal(.r$summary$selFreq, c(0.25, 0.5, 0.25))
  expect_equal(.r$summary$mean, c(0.55, 0.65, 0.525))
  # regimens given in another order than the values give the same answer
  .rev <- cuiFromValues(.varyValues()[15:1, ],
                        regimens = cuiRegimens(dose = 3:1, label = c("hi", "mid", "lo")))
  .m <- match(.r$summary$regimen, .rev$summary$regimen)
  expect_equal(.rev$summary$selFreq[.m], .r$summary$selFreq)
  expect_equal(.rev$summary$lower[.m], .r$summary$lower)
  # per-assessment selection and intervals
  .two <- rbind(cbind(.varyValues(), assessment = "a1"),
                cbind(.varyValues()[c(3:1, 6:4, 9:7, 12:10, 15:13), ], assessment = "a2"))
  .two$value[.two$assessment == "a2"] <- rev(.two$value[.two$assessment == "a2"])
  .t <- cuiFromValues(.two, regimens = .reg)
  .s1 <- .t$summary[.t$summary$assessment == "a1", ]
  expect_equal(.s1$selFreq, c(0.25, 0.5, 0.25))
  expect_equal(sum(.t$summary$selFreq[.t$summary$assessment == "a2"]), 1)
})

test_that("an eligible data frame honors its eligible column and the object's own format", {
  .r <- cuiFromValues(.drawValues(), weights = c(e1 = 2, e2 = 1))
  .df <- data.frame(regimen = c("lo", "mid", "hi"), assessment = "all",
                    eligible = c(TRUE, TRUE, FALSE))
  .e <- update(.r, eligible = .df)
  expect_equal(.e$eligible$eligible, c(TRUE, TRUE, FALSE))
  expect_equal(.e$summary$regimen[.e$summary$best], "mid")
  # feeding back the result's own eligibility table keeps it
  expect_equal(update(.e, eligible = .e$eligible)$summary, .e$summary)
  expect_error(update(.r, eligible = data.frame(regimen = "lo", assessment = "All")),
               "unknown assessment")
  expect_error(update(.r, eligible = transform(.df, eligible = c(1, 1, 0))), "TRUE/FALSE")
})

test_that("update() equals building from scratch for every setting", {
  .v <- .varyValues()
  .r <- cuiFromValues(.v)
  expect_equal(update(.r, eligible = c("lo", "hi"))$summary,
               cuiFromValues(.v, eligible = c("lo", "hi"))$summary)
  expect_equal(update(.r, combine = "geometric")$summary,
               cuiFromValues(.v, combine = "geometric")$summary)
  expect_equal(update(.r, tieTol = 0.15)$summary, cuiFromValues(.v, tieTol = 0.15)$summary)
  expect_equal(update(.r, level = 0.5)$summary, cuiFromValues(.v, level = 0.5)$summary)
})

test_that("numeric regimens match cuiRegimens() labels and inputs are validated", {
  .v <- data.frame(regimen = rep(c(2.5, 10), 2), draw = rep(0:1, each = 2),
                   endpoint = "e", value = c(0.1, 0.2, 0.3, 0.1))
  .r <- cuiFromValues(.v, regimens = cuiRegimens(dose = c(2.5, 10)))
  expect_equal(.r$summary$regimen, c("2.5", "10"))
  expect_equal(.r$summary$selFreq, c(1, 0))
  .big <- data.frame(regimen = c(1e5, 1.5e5), endpoint = "e", value = c(0.1, 0.2))
  expect_equal(cuiFromValues(.big, regimens = cuiRegimens(dose = c(1e5, 1.5e5)))$summary$est,
               c(0.1, 0.2))
  expect_error(cuiFromValues(data.frame(regimen = 1:2, endpoint = "e",
                                        value = factor(c("90", "10")))), "numeric")
  expect_error(cuiFromValues(transform(.v, draw = draw + 0.5)), "whole numbers")
  # a 'draws' column is not partially matched as 'draw'
  .p <- data.frame(regimen = c(1, 2), endpoint = "e", value = c(0.1, 0.2), draws = c(3, 4))
  expect_equal(unique(cuiFromValues(.p)$cui$draw), 0L)
})
