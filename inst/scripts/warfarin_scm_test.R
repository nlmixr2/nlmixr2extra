# =============================================================================
# Warfarin PK covariate search (SCM) — integration test
#
# Data   : nlmixr2data::warfarin, filtered to dvid == "cp" (PK only)
# Model  : 1-CMT oral, proportional residual error
# Covars : wt (continuous), sex (categorical)
# Method : full SCM (forward p<0.05 then backward p<0.01)
#
# Design note — no data pre-transformation:
#   Covariate centering and all other transformations are generated inside the
#   model body by covarSearchAuto().  The data passed to the SCM must be the
#   original, unmodified dataset (except for the dvid filter below).
#
#   For continuous covariates the generated expression is, e.g.:
#     cov_wt_power_v * log(wt / 70.5)     # power shape, center = median(wt)
#     cov_wt_lin_v   * (wt - 70.5)        # lin   shape
#
#   For categorical covariates the generated expression is, e.g.:
#     cov_sex_1_v * ifelse(sex == 1, 1, 0)   # sex coded 0/1 (numeric)
#     cov_sex_male_v * ifelse(sex == "male", 1, 0)   # sex coded as character
#
#   The data argument must be passed explicitly to covarSearchAuto() so that
#   covariate columns absent from the base model (e.g. sex) are available when
#   candidate models are fitted.
# =============================================================================

library(nlmixr2)
library(nlmixr2data)

# Load the development version of nlmixr2extra.
# Run this script with the working directory set to anywhere inside the
# nlmixr2extra package tree (e.g. the package root itself).  devtools will
# walk up from the cwd until it finds DESCRIPTION.
devtools::load_all(quiet = FALSE)

# Safety check: confirm the in-development version was loaded (not an old
# installed copy that lacks the 'data' / 'saveModels' parameters).
stopifnot(
  "'data' parameter missing — load_all may have failed or loaded wrong package" =
    "data"       %in% names(formals(covarSearchAuto)),
  "'saveModels' parameter missing — load_all may have failed or loaded wrong package" =
    "saveModels" %in% names(formals(covarSearchAuto))
)

# ── 1. Prepare data ───────────────────────────────────────────────────────────

# Keep PK observations and their corresponding dose records (evid 0/1, dvid cp)
# Drops the anticoagulation PD endpoint (dvid == "pca").
# No covariate transformations here — the SCM generates them in the model body.
pkdata <- warfarin[warfarin$dvid == "cp", ]

d_subj <- pkdata[!duplicated(pkdata$id), ]
cat(sprintf(
  "Data ready: %d subjects, %d rows\n",
  length(unique(pkdata$id)), nrow(pkdata)
))
cat(sprintf(
  "  wt  : median %.1f  range %.1f – %.1f  (raw, centering done in model)\n",
  median(d_subj$wt), min(d_subj$wt), max(d_subj$wt)
))
cat(sprintf("  sex : %s\n\n",
  paste(names(table(d_subj$sex)), table(d_subj$sex), sep = "=", collapse = ", ")
))

# ── 2. Base warfarin 1-CMT model ─────────────────────────────────────────────
# The base model contains no covariate terms; covarSearchAuto adds them.

warf_pk <- function() {
  ini({
    tka  <- log(1.15)   # log absorption rate constant (h^-1)
    tcl  <- log(0.135)  # log clearance (L/h)
    tv   <- log(7.0)    # log volume of distribution (L)
    eta.ka ~ 0.40
    eta.cl ~ 0.25
    eta.v  ~ 0.10
    prop.err <- 0.10
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    linCmt() ~ prop(prop.err)
  })
}

# ── 3. Fit base model ─────────────────────────────────────────────────────────

cat("Fitting base model (SAEM)...\n")
fit_base <- nlmixr2(
  warf_pk, pkdata,
  est     = "saem", table = tableControl(cwres = T)
)

cat("\n--- Base model summary ---\n")
print(fit_base)

cat("\nOmega (BSV variances):\n")
print(round(diag(fit_base$omega), 4))

cat(sprintf(
  "\nBase model OFV: %.3f   AIC: %.3f   BIC: %.3f\n\n",
  fit_base$objf, fit_base$AIC, fit_base$BIC
))

# ── 4. Run full SCM ───────────────────────────────────────────────────────────
# covarSearchAuto writes a cache directory; run inside a temporary directory.
#
# pairsVec specifies pairs and shapes explicitly:
#   wt~v : power shape -> log(wt / median(wt))  in model body
#          lin   shape -> (wt - median(wt))      in model body
#   sex~v: categorical, auto-detected because "sex" is in catvarsVec;
#          generates ifelse(sex == <non-ref level>, 1, 0) in model body
#
# data = pkdata must be passed explicitly: the base model does not reference
# wt or sex, so nlme::getData(fit_base) may drop those columns.

cat("Starting SCM covariate search (forward p<0.05, backward p<0.01)...\n")
cat("Pairs: wt~v (power + lin shapes), sex~v (categorical)\n\n")

rxode2::.rxWithWd(tempdir(), {
  scm <- covarSearchAuto(
    fit        = fit_base,
    data       = pkdata,          # full data including covariate columns
    pairsVec   = list(
      list(var = "v", covar = "wt",  shapes = c("power", "lin")),
      list(var = "v", covar = "sex") # catvarsVec triggers categorical expansion
    ),
    catvarsVec = c("sex"),
    pVal       = list(fwd = 0.05, bck = 0.01),
    searchType = "scm",
    saveModels = TRUE,
    verbose    = TRUE,
    restart    = TRUE
  )
})

# ── 5. Results ────────────────────────────────────────────────────────────────

cat("\n\n=======================================================\n")
cat("  SCM RESULTS\n")
cat("=======================================================\n\n")

cat("--- Full summary table ---\n")
print(scm$summaryTable)

# Forward search detail
cat("\n--- Forward search steps ---\n")
fwd_tbl <- scm$resFwd[[2]]
if (!is.null(fwd_tbl) && nrow(fwd_tbl) > 0) {
  cols <- intersect(
    c("step", "covar", "var", "shape", "deltObjf", "pchisqr",
      "bsvReduction", "covarEffect", "included"),
    colnames(fwd_tbl)
  )
  print(fwd_tbl[, cols])
} else {
  cat("(no forward steps recorded)\n")
}

# Backward search detail
cat("\n--- Backward search steps ---\n")
bck_tbl <- scm$resBck[[2]]
if (!is.null(bck_tbl) && nrow(bck_tbl) > 0) {
  cols <- intersect(
    c("step", "covar", "var", "shape", "deltObjf", "pchisqr",
      "bsvReduction", "covarEffect", "included"),
    colnames(bck_tbl)
  )
  print(bck_tbl[, cols])
} else {
  cat("(no backward steps recorded)\n")
}

# Final model
fit_final <- scm$resFwd[[1]]
cat("\n--- Final model parameters ---\n")
print(fit_final$parFixedDf)

cat("\n--- Final model omega (BSV variances) ---\n")
print(round(diag(fit_final$omega), 4))

cat(sprintf(
  "\nFinal model OFV: %.3f   AIC: %.3f   BIC: %.3f\n",
  fit_final$objf, fit_final$AIC, fit_final$BIC
))
cat(sprintf(
  "ΔOFV vs base:    %.3f\n",
  fit_base$objf - fit_final$objf
))

# Covariates accepted (included == "yes" in forward step)
if (!is.null(fwd_tbl) && "included" %in% colnames(fwd_tbl)) {
  accepted <- fwd_tbl[fwd_tbl$included == "yes",
                      c("covar", "var", "shape", "deltObjf", "pchisqr", "bsvReduction")]
  cat("\n--- Accepted covariates (forward, p < 0.05) ---\n")
  if (nrow(accepted) > 0) print(accepted) else cat("None\n")
}
