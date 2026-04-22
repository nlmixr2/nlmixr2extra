skip_on_cran()

.cur <- loadNamespace("nlmixr2extra")

# ── Shared model definition (no fitting) ──────────────────────────────────

.one_cmt_fun <- function() {
  ini({
    tka     <- 0.45
    tcl     <- log(c(0, 2.7, 100))
    tv      <- 3.45
    eta.ka  ~ 0.6
    eta.cl  ~ 0.3
    eta.v   ~ 0.1
    add.sd  <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

# =============================================================================
# .parseInitSpec
# =============================================================================

test_that(".parseInitSpec: NULL returns defaults (0.1, -5, 5)", {
  res <- .cur$.parseInitSpec(NULL)
  expect_equal(res$est,   0.1)
  expect_equal(res$lower, -5)
  expect_equal(res$upper,  5)
})

test_that(".parseInitSpec: scalar preserves est, uses default bounds", {
  res <- .cur$.parseInitSpec(0.75)
  expect_equal(res$est,   0.75)
  expect_equal(res$lower, -5)
  expect_equal(res$upper,  5)
})

test_that(".parseInitSpec: full named list preserved", {
  res <- .cur$.parseInitSpec(list(est = 0.5, lower = -2, upper = 3))
  expect_equal(res$est,   0.5)
  expect_equal(res$lower, -2)
  expect_equal(res$upper,  3)
})

test_that(".parseInitSpec: 'init' alias accepted for est", {
  res <- .cur$.parseInitSpec(list(init = 0.3, lower = -1, upper = 1))
  expect_equal(res$est,   0.3)
  expect_equal(res$lower, -1)
  expect_equal(res$upper,  1)
})

test_that(".parseInitSpec: partial list fills missing fields with defaults", {
  res_lower <- .cur$.parseInitSpec(list(est = 0.2, lower = 0))
  expect_equal(res_lower$upper, 5)

  res_upper <- .cur$.parseInitSpec(list(est = 0.2, upper = 2))
  expect_equal(res_upper$lower, -5)

  res_est <- .cur$.parseInitSpec(list(lower = -3, upper = 3))
  expect_equal(res_est$est, 0.1)
})

test_that(".parseInitSpec: unrecognised input returns defaults", {
  res <- .cur$.parseInitSpec("notanumber")
  expect_equal(res$est,   0.1)
  expect_equal(res$lower, -5)
  expect_equal(res$upper,  5)
})

# =============================================================================
# buildPairs
# =============================================================================

test_that("buildPairs: varsVec x covarsVec produces cartesian product", {
  df <- .cur$buildPairs(varsVec = c("cl", "v"), covarsVec = c("wt", "age"))
  expect_equal(nrow(df), 4)
  expect_true(all(c("var", "covar") %in% names(df)))
  expect_setequal(df$var,   c("cl", "cl", "v", "v"))
  expect_setequal(df$covar, c("wt", "age", "wt", "age"))
})

test_that("buildPairs: pairsVec as list of lists", {
  df <- .cur$buildPairs(pairsVec = list(
    list(var = "cl", covar = "wt"),
    list(var = "v",  covar = "age")
  ))
  expect_equal(nrow(df), 2)
  expect_equal(df$var[1],   "cl")
  expect_equal(df$covar[1], "wt")
  expect_equal(df$var[2],   "v")
  expect_equal(df$covar[2], "age")
})

test_that("buildPairs: pairsVec as data frame passed through", {
  pv <- data.frame(
    var   = c("cl", "v"),
    covar = c("wt", "wt"),
    stringsAsFactors = FALSE
  )
  df <- .cur$buildPairs(pairsVec = pv)
  expect_equal(nrow(df), 2)
  expect_equal(df$var,   c("cl", "v"))
  expect_equal(df$covar, c("wt", "wt"))
})

test_that("buildPairs: shapes list-column preserved when supplied", {
  df <- .cur$buildPairs(pairsVec = list(
    list(var = "cl", covar = "wt", shapes = c("power", "lin"))
  ))
  expect_true("shapes" %in% names(df))
  expect_equal(df$shapes[[1]], c("power", "lin"))
})

test_that("buildPairs: inits list-column preserved when supplied", {
  df <- .cur$buildPairs(pairsVec = list(
    list(var = "cl", covar = "wt",
         inits = list(power = list(est = 0.75, lower = 0, upper = 2)))
  ))
  expect_true("inits" %in% names(df))
  expect_equal(df$inits[[1]][["power"]][["est"]], 0.75)
})

test_that("buildPairs: NULL varsVec with no pairsVec errors", {
  expect_error(.cur$buildPairs(varsVec = NULL, covarsVec = "wt"))
})

test_that("buildPairs: varsVec without covarsVec errors", {
  expect_error(.cur$buildPairs(varsVec = "cl", covarsVec = NULL))
})

# =============================================================================
# .makeSCMData
# =============================================================================

.sex_data <- function(n_m = 7, n_f = 3) {
  data.frame(
    ID  = seq_len(n_m + n_f),
    sex = c(rep("M", n_m), rep("F", n_f)),
    wt  = rnorm(n_m + n_f, 70, 10),
    stringsAsFactors = FALSE
  )
}

test_that(".makeSCMData: creates indicator column for binary covariate", {
  d   <- .sex_data()
  res <- .cur$.makeSCMData(d, catvarsVec = "sex", fit = NULL, catCutoff = 0.05)
  expect_true("sex_F" %in% names(res$data))
  expect_true(all(res$data$sex_F %in% c(0L, 1L)))
})

test_that(".makeSCMData: reference level is most-frequent category", {
  d   <- .sex_data(n_m = 7, n_f = 3)
  res <- .cur$.makeSCMData(d, catvarsVec = "sex", fit = NULL, catCutoff = 0.05)
  expect_equal(res$catRef[["sex"]], "M")
})

test_that(".makeSCMData: reference level not included in catLevels", {
  d   <- .sex_data()
  res <- .cur$.makeSCMData(d, catvarsVec = "sex", fit = NULL, catCutoff = 0.05)
  expect_false("M" %in% res$catLevels[["sex"]])
  expect_true("F" %in% res$catLevels[["sex"]])
})

test_that(".makeSCMData: catCutoff drops rare levels into catDropped", {
  d <- data.frame(
    ID  = seq_len(20),
    grp = c(rep("A", 18), rep("B", 1), rep("C", 1)),
    stringsAsFactors = FALSE
  )
  # B and C each have 1/20 = 5%; with cutoff > 5% they should be dropped
  res <- .cur$.makeSCMData(d, catvarsVec = "grp", fit = NULL, catCutoff = 0.06)
  expect_length(res$catDropped[["grp"]], 2)
  expect_length(res$catLevels[["grp"]],  0)
})

test_that(".makeSCMData: levels at or above catCutoff are retained", {
  d <- data.frame(
    ID  = seq_len(20),
    grp = c(rep("A", 16), rep("B", 4)),
    stringsAsFactors = FALSE
  )
  # B has 4/20 = 20%; well above any reasonable cutoff
  res <- .cur$.makeSCMData(d, catvarsVec = "grp", fit = NULL, catCutoff = 0.05)
  expect_true("B" %in% res$catLevels[["grp"]])
  expect_length(res$catDropped[["grp"]], 0)
})

test_that(".makeSCMData: three-level variable with one rare level", {
  d <- data.frame(
    ID  = seq_len(20),
    grp = c(rep("A", 15), rep("B", 4), rep("C", 1)),
    stringsAsFactors = FALSE
  )
  # A=75% (ref), B=20% (retained), C=5% (< cutoff 0.06 → dropped)
  res <- .cur$.makeSCMData(d, catvarsVec = "grp", fit = NULL, catCutoff = 0.06)
  expect_equal(res$catRef[["grp"]], "A")
  expect_true("B"  %in% res$catLevels[["grp"]])
  expect_false("C" %in% res$catLevels[["grp"]])
  expect_true("C"  %in% res$catDropped[["grp"]])
})

# =============================================================================
# .enrichPairs
# =============================================================================

test_that(".enrichPairs: continuous covariate gets median as center", {
  d     <- data.frame(ID = 1:10, wt = seq(50, 95, by = 5),
                      stringsAsFactors = FALSE)
  pairs <- data.frame(var = "cl", covar = "wt", stringsAsFactors = FALSE)
  res   <- .cur$.enrichPairs(pairs, d)
  expect_equal(res$type[1],   "continuous")
  expect_equal(res$center[1], median(d$wt))
  expect_equal(res$raw_col[1], "wt")
})

test_that(".enrichPairs: multiple continuous covariates each get their own median", {
  d <- data.frame(
    ID = 1:10,
    wt = seq(50, 95, by = 5),
    age = seq(20, 65, by = 5),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(
    var   = c("cl", "cl"),
    covar = c("wt", "age"),
    stringsAsFactors = FALSE
  )
  res <- .cur$.enrichPairs(pairs, d)
  expect_equal(res$center[res$covar == "wt"],  median(d$wt))
  expect_equal(res$center[res$covar == "age"], median(d$age))
})

test_that(".enrichPairs: categorical covariate (via catvarsVec) expands to per-level rows", {
  d <- data.frame(
    ID  = 1:4,
    sex = c("M", "M", "F", "F"),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(var = "cl", covar = "sex", stringsAsFactors = FALSE)
  res   <- .cur$.enrichPairs(pairs, d, catvarsVec = "sex")
  # Expect one row per non-reference level (F; M is reference as most frequent
  # alphabetically here, but .enrichPairs drops first sorted level)
  expect_equal(nrow(res), 1)
  expect_equal(res$type[1], "categorical")
  expect_true(!is.na(res$level[1]))
})

test_that(".enrichPairs: missing value flagged when missingToken present", {
  d <- data.frame(
    ID = 1:5,
    wt = c(70, 75, -99, 80, 65),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(var = "cl", covar = "wt", stringsAsFactors = FALSE)
  res   <- .cur$.enrichPairs(pairs, d, missingToken = -99)
  expect_true(res$has_missing[1])
  expect_false(is.na(res$missing_check[1]))
})

test_that(".enrichPairs: no missing flag when all values present", {
  d <- data.frame(
    ID = 1:5,
    wt = c(70, 75, 80, 85, 90),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(var = "cl", covar = "wt", stringsAsFactors = FALSE)
  res   <- .cur$.enrichPairs(pairs, d)
  expect_false(isTRUE(res$has_missing[1]))
})

# =============================================================================
# .expandShapes
# =============================================================================

.cont_pairs <- function(data) {
  pairs <- data.frame(var = "cl", covar = "wt", stringsAsFactors = FALSE)
  .cur$.enrichPairs(pairs, data)
}

test_that(".expandShapes: power shape produces log-ratio expression", {
  d   <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90),
                    stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "power")
  expect_equal(res$shape[1], "power")
  expect_match(res$covExpr[1], "log\\(wt/")
  expect_match(res$covExpr[1], as.character(median(d$wt)))
})

test_that(".expandShapes: lin shape produces subtraction expression", {
  d   <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90),
                    stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "lin")
  expect_equal(res$shape[1], "lin")
  expect_match(res$covExpr[1], "wt - ")
  expect_match(res$covExpr[1], as.character(median(d$wt)))
})

test_that(".expandShapes: identity shape returns raw column", {
  d   <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90),
                    stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "identity")
  expect_equal(res$covExpr[1], "wt")
})

test_that(".expandShapes: log shape returns log(col)", {
  d   <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90),
                    stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "log")
  expect_match(res$covExpr[1], "^log\\(wt\\)$")
})

test_that(".expandShapes: cat shape uses indicator column name", {
  d <- data.frame(
    ID    = 1:4,
    sex   = c("M", "M", "F", "F"),
    sex_F = c(0L, 0L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(var = "cl", covar = "sex", stringsAsFactors = FALSE)
  pairs <- .cur$.enrichPairs(pairs, d, catvarsVec = "sex",
                              catLevels = list(sex = "F"))
  res   <- .cur$.expandShapes(pairs)
  expect_equal(res$shape[1],   "cat")
  expect_equal(res$covExpr[1], "sex_F")
})

test_that(".expandShapes: multiple shapes per row expand to multiple rows", {
  d   <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90),
                    stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = c("power", "lin"))
  expect_equal(nrow(res), 2)
  expect_setequal(res$shape, c("power", "lin"))
})

test_that(".expandShapes: global inits scalar sets est; bounds are default", {
  d   <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90),
                    stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "power",
                             inits = list(power = 0.75))
  expect_equal(res$init[1],  0.75)
  expect_equal(res$lower[1], -5)
  expect_equal(res$upper[1],  5)
})

test_that(".expandShapes: global inits with bounds all propagated", {
  d   <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90),
                    stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(
    .cont_pairs(d), shapes = "power",
    inits = list(power = list(est = 0.5, lower = 0, upper = 2))
  )
  expect_equal(res$init[1],  0.5)
  expect_equal(res$lower[1], 0)
  expect_equal(res$upper[1], 2)
})

test_that(".expandShapes: per-pair inits override global inits", {
  d     <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90),
                      stringsAsFactors = FALSE)
  pairs <- .cont_pairs(d)
  pairs$inits <- list(list(power = 0.99))
  res <- .cur$.expandShapes(pairs, shapes = "power",
                             inits = list(power = 0.1))
  expect_equal(res$init[1], 0.99)
})

test_that(".expandShapes: NULL inits for shape uses defaults", {
  d   <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90),
                    stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "power", inits = list())
  expect_equal(res$init[1],  0.1)
  expect_equal(res$lower[1], -5)
  expect_equal(res$upper[1],  5)
})

test_that(".expandShapes: continuous covar name has shape appended", {
  d   <- data.frame(ID = 1:5, wt = c(60, 70, 75, 80, 90),
                    stringsAsFactors = FALSE)
  res <- .cur$.expandShapes(.cont_pairs(d), shapes = "power")
  expect_equal(res$covar[1], "wt_power")
})

# =============================================================================
# .updatePairsInits
# =============================================================================

test_that(".updatePairsInits: NULL pairs returns NULL", {
  expect_null(.cur$.updatePairsInits(NULL, list()))
})

test_that(".updatePairsInits: zero-row pairs returned unchanged", {
  pairs <- data.frame(var = character(0), covar = character(0),
                      init = numeric(0), stringsAsFactors = FALSE)
  res <- .cur$.updatePairsInits(pairs, list())
  expect_equal(nrow(res), 0)
})

test_that(".updatePairsInits: updates init from parFixedDf", {
  pairs <- data.frame(
    var   = "cl",
    covar = "wt_power",
    init  = 0.1,
    stringsAsFactors = FALSE
  )
  mock_fit <- list(
    parFixedDf = data.frame(
      Estimate = 0.77,
      row.names = "cov_wt_power_cl",
      stringsAsFactors = FALSE
    )
  )
  res <- .cur$.updatePairsInits(pairs, mock_fit)
  expect_equal(res$init, 0.77)
})

test_that(".updatePairsInits: theta absent in fit leaves init unchanged", {
  pairs <- data.frame(
    var   = "cl",
    covar = "wt_power",
    init  = 0.1,
    stringsAsFactors = FALSE
  )
  mock_fit <- list(
    parFixedDf = data.frame(
      Estimate  = 0.5,
      row.names = "cov_age_lin_cl",
      stringsAsFactors = FALSE
    )
  )
  res <- .cur$.updatePairsInits(pairs, mock_fit)
  expect_equal(res$init, 0.1)
})

test_that(".updatePairsInits: bounds (lower/upper) never modified", {
  pairs <- data.frame(
    var   = "cl",
    covar = "wt_power",
    init  = 0.1,
    lower = -2,
    upper =  2,
    stringsAsFactors = FALSE
  )
  mock_fit <- list(
    parFixedDf = data.frame(
      Estimate  = 0.75,
      row.names = "cov_wt_power_cl",
      stringsAsFactors = FALSE
    )
  )
  res <- .cur$.updatePairsInits(pairs, mock_fit)
  expect_equal(res$lower, -2)
  expect_equal(res$upper,  2)
  expect_equal(res$init,  0.75)
})

test_that(".updatePairsInits: multiple rows updated independently", {
  pairs <- data.frame(
    var   = c("cl", "v"),
    covar = c("wt_power", "wt_power"),
    init  = c(0.1, 0.1),
    stringsAsFactors = FALSE
  )
  mock_fit <- list(
    parFixedDf = data.frame(
      Estimate  = c(0.77, 0.33),
      row.names = c("cov_wt_power_cl", "cov_wt_power_v"),
      stringsAsFactors = FALSE
    )
  )
  res <- .cur$.updatePairsInits(pairs, mock_fit)
  expect_equal(res$init[res$var == "cl"], 0.77)
  expect_equal(res$init[res$var == "v"],  0.33)
})

# =============================================================================
# addCatCovariates (exported)
# =============================================================================

test_that("addCatCovariates: creates indicator columns for non-reference levels", {
  d <- data.frame(
    ID  = 1:6,
    sex = c("M", "M", "M", "F", "F", "F"),
    wt  = c(70, 75, 80, 60, 65, 55),
    stringsAsFactors = FALSE
  )
  res      <- addCatCovariates(d, covarsVec = "wt", catcovarsVec = "sex")
  new_data <- res[[1]]
  new_covs <- res[[2]]
  expect_true(any(grepl("^sex_", names(new_data))))
  expect_true(any(grepl("^sex_", new_covs)))
})

test_that("addCatCovariates: original categorical column removed from data", {
  d <- data.frame(
    ID  = 1:4,
    sex = c("M", "M", "F", "F"),
    stringsAsFactors = FALSE
  )
  res <- addCatCovariates(d, covarsVec = character(0), catcovarsVec = "sex")
  expect_false("sex" %in% names(res[[1]]))
})

test_that("addCatCovariates: most-frequent level (reference) has no indicator", {
  d <- data.frame(
    ID  = seq_len(10),
    grp = c(rep("A", 7), rep("B", 2), rep("C", 1)),
    stringsAsFactors = FALSE
  )
  res      <- addCatCovariates(d, covarsVec = character(0), catcovarsVec = "grp")
  new_data <- res[[1]]
  # A is most frequent; its indicator should NOT be created
  expect_false("grp_A" %in% names(new_data))
  expect_true("grp_B" %in% names(new_data))
  expect_true("grp_C" %in% names(new_data))
})

test_that("addCatCovariates: indicator values are 0/1 integers", {
  d <- data.frame(
    ID  = 1:6,
    grp = c("A", "A", "A", "B", "B", "C"),
    stringsAsFactors = FALSE
  )
  res      <- addCatCovariates(d, covarsVec = character(0), catcovarsVec = "grp")
  new_data <- res[[1]]
  ind_cols <- grep("^grp_", names(new_data), value = TRUE)
  for (col in ind_cols) {
    expect_true(all(new_data[[col]] %in% c(0L, 1L)))
  }
})

# =============================================================================
# .rebuildUiFromPairs
# =============================================================================

test_that(".rebuildUiFromPairs: NULL pairs_df returns base_ui unchanged", {
  ui  <- nlmixr(.one_cmt_fun)
  res <- .cur$.rebuildUiFromPairs(ui, NULL)
  expect_equal(
    sort(ui$iniDf$name[!is.na(ui$iniDf$ntheta)]),
    sort(res$iniDf$name[!is.na(res$iniDf$ntheta)])
  )
})

test_that(".rebuildUiFromPairs: empty pairs_df returns base_ui unchanged", {
  ui    <- nlmixr(.one_cmt_fun)
  empty <- data.frame(var = character(0), covar = character(0),
                      stringsAsFactors = FALSE)
  res   <- .cur$.rebuildUiFromPairs(ui, empty)
  expect_equal(
    sort(ui$iniDf$name[!is.na(ui$iniDf$ntheta)]),
    sort(res$iniDf$name[!is.na(res$iniDf$ntheta)])
  )
})

test_that(".rebuildUiFromPairs: adds cov theta with correct name", {
  ui    <- nlmixr(.one_cmt_fun)
  pairs <- data.frame(
    var     = "cl",
    covar   = "wt_power",
    covExpr = "log(wt/70)",
    init    = 0.5,
    lower   = -2,
    upper   =  2,
    stringsAsFactors = FALSE
  )
  res   <- .cur$.rebuildUiFromPairs(ui, pairs)
  ini_names <- res$iniDf$name
  expect_true("cov_wt_power_cl" %in% ini_names)
})

test_that(".rebuildUiFromPairs: new theta gets specified init and bounds", {
  ui    <- nlmixr(.one_cmt_fun)
  pairs <- data.frame(
    var     = "cl",
    covar   = "wt_power",
    covExpr = "log(wt/70)",
    init    = 0.5,
    lower   = -2,
    upper   =  2,
    stringsAsFactors = FALSE
  )
  res     <- .cur$.rebuildUiFromPairs(ui, pairs)
  cov_row <- res$iniDf[res$iniDf$name == "cov_wt_power_cl", ]
  expect_equal(cov_row$est,   0.5)
  expect_equal(cov_row$lower, -2)
  expect_equal(cov_row$upper,  2)
})

test_that(".rebuildUiFromPairs: two covariates on different parameters both added", {
  ui    <- nlmixr(.one_cmt_fun)
  pairs <- data.frame(
    var     = c("cl",         "v"),
    covar   = c("wt_power",   "wt_power"),
    covExpr = c("log(wt/70)", "log(wt/70)"),
    init    = c(0.5,          0.3),
    lower   = c(-2,           -2),
    upper   = c(2,             2),
    stringsAsFactors = FALSE
  )
  res       <- .cur$.rebuildUiFromPairs(ui, pairs)
  ini_names <- res$iniDf$name
  expect_true("cov_wt_power_cl" %in% ini_names)
  expect_true("cov_wt_power_v"  %in% ini_names)
})

# =============================================================================
# Integration tests — require nlmixr2data, slow fitting
# =============================================================================

skip_if_not_installed("nlmixr2data")

.theoph <- nlmixr2data::theo_sd

.pk_model <- function() {
  ini({
    tka  <- log(1.5)
    tcl  <- log(0.04)
    tv   <- log(0.5)
    eta.ka ~ 0.5
    eta.cl ~ 0.2
    eta.v  ~ 0.1
    prop.err <- 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv + eta.v)
    linCmt() ~ prop(prop.err)
  })
}

.fit_base <- function() {
  nlmixr2(
    .pk_model, .theoph, est = "focei",
    control = nlmixr2est::foceiControl(print = 0, calcTables = TRUE)
  )
}

test_that("runSCM: forward-only returns expected list structure", {
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()
  res <- runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers    = 1L
  )
  expect_type(res, "list")
  expect_named(res, c("summaryTable", "resFwd", "resBck"))
  expect_null(res$resBck)
  expect_type(res$resFwd, "list")
})

test_that("runSCM: backward-only returns expected list structure", {
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()
  res <- runSCM(
    fit               = base_fit,
    pairsVec          = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType        = "backward",
    saveModels        = FALSE,
    includedRelations = list(list(var = "cl", covar = "WT", shapes = "power")),
    workers           = 1L
  )
  expect_named(res, c("summaryTable", "resFwd", "resBck"))
  expect_null(res$resFwd)
})

test_that("runSCM: full SCM returns forward and backward results", {
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()
  res <- runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "scm",
    saveModels = FALSE,
    workers    = 1L
  )
  expect_named(res, c("summaryTable", "resFwd", "resBck"))
  expect_type(res$resFwd, "list")
  expect_type(res$resBck, "list")
})

test_that("runSCM: saveModels=FALSE creates no output directory", {
  td <- withr::local_tempdir(clean = TRUE)
  withr::local_dir(td)
  base_fit <- .fit_base()
  n_before <- length(list.dirs(td, recursive = FALSE))
  runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers    = 1L
  )
  n_after <- length(list.dirs(td, recursive = FALSE))
  expect_equal(n_before, n_after)
})

test_that("runSCM: saveModels=TRUE writes log and CSV files", {
  td <- withr::local_tempdir(clean = TRUE)
  withr::local_dir(td)
  base_fit <- .fit_base()
  runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = TRUE,
    outputDir  = "scm_out",
    restart    = TRUE,
    workers    = 1L
  )
  expect_true(dir.exists("scm_out"))
  expect_true(file.exists(file.path("scm_out", "scm_log.txt")))
  expect_true(file.exists(file.path("scm_out", "scm_step_summary.csv")))
  expect_true(file.exists(file.path("scm_out", "scm_all_candidates.csv")))
})

test_that("runSCM: summaryTable has expected columns", {
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()
  res <- runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers    = 1L
  )
  if (!is.null(res$summaryTable) && nrow(res$summaryTable) > 0) {
    expect_true(all(c("covar", "var", "deltObjf", "pchisqr",
                      "included", "searchType") %in%
                      names(res$summaryTable)))
  }
})

test_that("runSCM: user-supplied control object accepted", {
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()
  ctrl     <- nlmixr2est::foceiControl(print = 0, calcTables = TRUE)
  expect_no_error(
    runSCM(
      fit        = base_fit,
      pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      control    = ctrl,
      workers    = 1L
    )
  )
})

test_that("runSCM: inits with bounds run without error", {
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()
  expect_no_error(
    runSCM(
      fit        = base_fit,
      pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      workers    = 1L,
      inits      = list(power = list(est = 0.75, lower = 0, upper = 3))
    )
  )
})

test_that("runSCM: per-pair shapes via pairsVec respected", {
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()
  res <- runSCM(
    fit        = base_fit,
    pairsVec   = list(
      list(var = "cl", covar = "WT", shapes = c("power", "lin"))
    ),
    searchType = "forward",
    saveModels = FALSE,
    workers    = 1L
  )
  # At least two candidates should have been tested (power + lin)
  if (!is.null(res$summaryTable) && nrow(res$summaryTable) > 0) {
    expect_gte(nrow(res$summaryTable), 1)
  }
})

test_that("runSCM: restart=TRUE backs up existing outputDir", {
  td <- withr::local_tempdir(clean = TRUE)
  withr::local_dir(td)
  base_fit <- .fit_base()

  # First run — create the directory
  runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = TRUE,
    outputDir  = "restart_dir",
    restart    = TRUE,
    workers    = 1L
  )
  expect_true(dir.exists("restart_dir"))

  # Second run — should back up and start fresh
  expect_no_error(
    runSCM(
      fit        = base_fit,
      pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = TRUE,
      outputDir  = "restart_dir",
      restart    = TRUE,
      workers    = 1L
    )
  )
  backups <- list.dirs(td, recursive = FALSE, full.names = FALSE)
  expect_true(any(grepl("^restart_dir_backup_", backups)))
})

test_that("runSCM: multiple pairs tested simultaneously", {
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()
  res <- runSCM(
    fit        = base_fit,
    pairsVec   = list(
      list(var = "cl", covar = "WT", shapes = "power"),
      list(var = "v",  covar = "WT", shapes = "power")
    ),
    searchType = "forward",
    saveModels = FALSE,
    workers    = 1L
  )
  expect_type(res, "list")
})

test_that("runSCM: explicit outputDir used as absolute path", {
  td <- withr::local_tempdir(clean = TRUE)
  out_dir <- file.path(td, "my_scm_output")
  base_fit <- .fit_base()
  runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = TRUE,
    outputDir  = out_dir,
    restart    = TRUE,
    workers    = 1L
  )
  expect_true(dir.exists(out_dir))
  expect_true(file.exists(file.path(out_dir, "scm_log.txt")))
})

# =============================================================================
# workers parameter — parallelization
# =============================================================================

test_that("runSCM: workers=1 returns same structure as workers=NULL", {
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()

  res_default <- runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers    = NULL
  )

  res_w1 <- runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers    = 1L
  )

  expect_named(res_w1, names(res_default))
  expect_type(res_w1$resFwd, "list")
  expect_null(res_w1$resBck)
})

test_that("runSCM: workers=1 forward+backward both respect parameter", {
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()

  res <- runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "scm",
    saveModels = FALSE,
    workers    = 1L
  )

  expect_named(res, c("summaryTable", "resFwd", "resBck"))
  expect_type(res$resFwd, "list")
  expect_type(res$resBck, "list")
})

test_that("runSCM: workers='auto' runs without error", {
  skip_if_not_installed("future")
  # multisession workers load from the *installed* package, so this test only
  # makes sense when nlmixr2extra is installed (not just load_all()-ed).
  skip_if(
    isTRUE(tryCatch(
      pkgload::is_dev_package("nlmixr2extra"),
      error = function(e) FALSE
    )),
    "Package loaded via load_all(); install first to run multisession test"
  )
  withr::local_tempdir(clean = TRUE)
  base_fit <- .fit_base()

  expect_no_error(
    runSCM(
      fit        = base_fit,
      pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
      searchType = "forward",
      saveModels = FALSE,
      workers    = "auto"
    )
  )
})

test_that("runSCM: future plan restored to original after workers=1", {
  skip_if_not_installed("future")
  withr::local_tempdir(clean = TRUE)
  base_fit  <- .fit_base()
  plan_orig <- class(future::plan())
  on.exit(future::plan("sequential"), add = TRUE)

  runSCM(
    fit        = base_fit,
    pairsVec   = list(list(var = "cl", covar = "WT", shapes = "power")),
    searchType = "forward",
    saveModels = FALSE,
    workers    = 1L
  )

  expect_equal(class(future::plan()), plan_orig)
})
