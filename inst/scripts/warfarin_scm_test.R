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
#   model body by runSCM().  The data passed to the SCM must be the
#   original, unmodified dataset (except for the dvid filter below).
#
#   For continuous covariates the generated expression multiplies a theta
#   parameter by a shape transform, e.g. log(wt / 70.5) for power or
#   (wt - 70.5) for linear, where 70.5 is the observed median of wt.
#
#   For categorical covariates runSCM() automatically creates 0/1
#   indicator columns (one per non-reference level) in a temporary dataset,
#   then generates expressions using those columns directly.  For example,
#   a binary sex column produces sex_1 or sex_male multiplied by its theta.
#   No ifelse() is emitted; the indicator column must exist in the data.
#
#   The data argument must be passed explicitly to runSCM() so that
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
  "load_all failed: 'data' parameter missing from runSCM" =
    "data"       %in% names(formals(runSCM)),
  "load_all failed: 'saveModels' parameter missing from runSCM" =
    "saveModels" %in% names(formals(runSCM)),
  "load_all failed: 'catCutoff' parameter missing from runSCM" =
    "catCutoff"  %in% names(formals(runSCM))
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
# The base model contains no covariate terms; runSCM adds them.

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
  est     = "saem",
  control = saemControl(print = 100, nBurn = 500, nEm = 500),
  table   = tableControl(cwres = TRUE)
)

cat("\n--- Base model summary ---\n")
print(fit_base)

cat("\nOmega (BSV variances):\n")
print(round(diag(fit_base$omega), 4))

cat(sprintf(
  "\nBase model OFV: %.3f   AIC: %.3f   BIC: %.3f\n\n",
  fit_base$objf, fit_base$AIC, fit_base$BIC
))

# ── Parallel capability check (used throughout sections 4 and 7) ─────────────
# Multisession workers load from the *installed* package.  If this script is
# run via devtools::load_all() the installed copy may be outdated, causing
# internal functions added in the current dev session to be missing in workers.
# Parallel is therefore only enabled when the package is installed and current.
.is_dev <- isTRUE(tryCatch(
  pkgload::is_dev_package("nlmixr2extra"),
  error = function(e) FALSE
))
has_parallel <- !.is_dev &&
  requireNamespace("future",       quietly = TRUE) &&
  requireNamespace("future.apply", quietly = TRUE) &&
  requireNamespace("progressr",    quietly = TRUE)

if (.is_dev && (requireNamespace("future.apply", quietly = TRUE))) {
  message(
    "NOTE: parallel SCM/bootstrap disabled — package loaded via load_all().\n",
    "      Run devtools::install() then re-source this script to enable parallel."
  )
}

# ── 4. Run full SCM ───────────────────────────────────────────────────────────
# runSCM writes a cache directory; run inside a temporary directory.
#
# pairsVec specifies pairs and shapes explicitly:
#   wt~v : power shape -> log(wt / median(wt))  in model body
#          lin   shape -> (wt - median(wt))      in model body
#   sex~v: categorical; runSCM() creates sex_<level> indicator
#          columns automatically and uses them directly (no ifelse)
#
# data = pkdata must be passed explicitly: the base model does not reference
# wt or sex, so nlme::getData(fit_base) may drop those columns.

cat("Starting SCM covariate search (forward p<0.05, backward p<0.01)...\n")
cat("Pairs: wt~v (power + lin shapes), sex~v (categorical)\n\n")

.scm_pairs <- list(
  list(var = "v",  covar = "wt",  shapes = c("power")),
  list(var = "v",  covar = "sex"),
  list(var = "cl", covar = "wt",  shapes = c("power")),
  list(var = "cl", covar = "sex")
)

# ── Sequential SCM (workers = 1) ─────────────────────────────────────────────
cat("Running SCM sequential (workers = 1)...\n")
t_scm_seq <- system.time(
  scm_seq <- runSCM(
    fit        = fit_base,
    data       = pkdata,
    pairsVec   = .scm_pairs,
    catvarsVec = c("sex"),
    pVal       = list(fwd = 0.05, bck = 0.01),
    searchType = "scm",
    saveModels = TRUE,
    verbose    = TRUE,
    restart    = TRUE,
    print      = 100,
    workers    = 1L,
    confirm    = FALSE
  )
)
cat(sprintf("  Sequential SCM elapsed: %.1f s\n\n", t_scm_seq["elapsed"]))

# ── Parallel SCM ─────────────────────────────────────────────────────────────
if (has_parallel) {
  n_workers_scm <- min(4L, future::availableCores(omit = 1L))
  cat(sprintf("Running SCM parallel (%d workers)...\n", n_workers_scm))
  t_scm_par <- system.time(
    scm_par <- runSCM(
      fit        = fit_base,
      data       = pkdata,
      pairsVec   = .scm_pairs,
      catvarsVec = c("sex"),
      pVal       = list(fwd = 0.05, bck = 0.01),
      searchType = "scm",
      saveModels = TRUE,
      verbose    = TRUE,
      restart    = TRUE,
      print      = 100,
      workers    = n_workers_scm,
      confirm    = FALSE
    )
  )
  cat(sprintf("  Parallel SCM elapsed:    %.1f s\n\n", t_scm_par["elapsed"]))
} else {
  t_scm_par   <- NULL
  n_workers_scm <- NA_integer_
  scm_par     <- NULL
  cat("Skipping parallel SCM: future/future.apply/progressr not installed.\n\n")
}

# Use the sequential result for downstream sections.
scm <- scm_seq

# ── 5. Verify categorical covariate format ────────────────────────────────────
# Confirm that the SCM used indicator columns (sex_<level>) rather than
# ifelse() expressions in the generated model body.

cat("\n--- Categorical covariate format check ---\n")
final_ui <- scm$resFwd[[1]]$ui
model_lines <- as.character(nlmixr2est::.saemDropMuRefFromModel(final_ui))

# Check: no ifelse() in any model line involving sex
ifelse_in_model <- any(
  grepl("ifelse", model_lines, fixed = TRUE) &
    grepl("sex", model_lines, fixed = TRUE)
)
if (ifelse_in_model) {
  warning("UNEXPECTED: ifelse() found in sex covariate expression — ",
          "indicator column approach may not have been applied.")
} else {
  cat("OK: no ifelse() in sex covariate lines\n")
}

# Check: indicator-column expression present (sex_<something> * cov_sex_...)
indicator_present <- any(
  grepl("sex_[^=]+ \\* cov_sex_", model_lines, perl = TRUE) |
    grepl("sex_[^=]+=", model_lines, perl = TRUE)
)
if (indicator_present) {
  cat("OK: indicator column expression found in model body\n")
} else {
  cat("NOTE: sex covariate may not have been accepted into the final model ",
      "(check forward search table)\n")
}

# ── 6. Results ────────────────────────────────────────────────────────────────

cat("\n\n=======================================================\n")
cat("  SCM RESULTS (section 6)\n")
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
  keep_cols <- c("covar", "var", "shape", "deltObjf", "pchisqr", "bsvReduction")
  accepted  <- fwd_tbl[fwd_tbl$included == "yes", keep_cols]
  cat("\n--- Accepted covariates (forward, p < 0.05) ---\n")
  if (nrow(accepted) > 0) print(accepted) else cat("None\n")
}

# =============================================================================
# Section 7 — Bootstrap parameter uncertainty on the final SCM model
#
# The final model from the backward search is the natural target for bootstrap.
# Here we use fit_final (forward-only result) so the script works end-to-end
# even when the backward search drops all covariates.
#
# workers = NULL  -> respects the current future::plan() (default: sequential)
# workers = 4     -> spawns 4 multisession workers for the duration of the run
# workers = "auto"-> uses future::availableCores(omit = 1) workers
#
# Output directories follow the SCM naming convention: <fitName>_boot_<N>.
# Each directory contains:
#   boot_data_*.rds          per-replicate resampled datasets
#   modelsEnsemble_*.rds     per-replicate model parameter lists
#   fitEnsemble_*.rds        per-replicate full fit objects
#   bootstrap_seed.rds       full RNG state before any resampling
#   bootstrap_results.csv   per-replicate parameter estimates (one row each)
#   bootstrap_summary.txt   text summary: seed, original and bootstrap CIs
#
# restart = TRUE  creates a new numbered directory (e.g. warf_final_seq_boot_2)
# restart = FALSE (default) resumes from the most recently numbered directory
# =============================================================================

cat("\n\n=======================================================\n")
cat("  BOOTSTRAP (section 7)\n")
cat("=======================================================\n\n")

# ── 7a. Sequential bootstrap on the final (post-SCM) model ────────────────────
# runBootstrap() returns a named list AND (by default, updateFit = TRUE)
# attaches results to the fit object as a single new $bootstrap slot without
# modifying any pre-existing slot:
#
#   Return value / fit$bootstrap:
#     $seed        — full RNG state before resampling (integer vector)
#     $results     — per-replicate parameter estimates (data.frame)
#     $summary     — original parFixedDf augmented with Bootstrap CI columns
#     $model       — fitName string
#     $outputDir   — absolute path to the output directory
#     $timestamp   — POSIXct when the function returned
#
#   fit$bootstrap additionally contains:
#     $omegaSummary, $covMatrix, $corMatrix, $bias, $nboot, $ci
#
# nboot = 50 is enough to verify the workflow; use >= 200 for real uncertainty.

cat("Running sequential bootstrap (workers = 1)...\n")
t_boot_seq <- system.time(
  boot_final_seq <- runBootstrap(
    fit      = fit_final,
    nboot    = 50,
    restart  = TRUE,
    fitName  = "warf_final_seq",
    workers  = 1L
  )
)
cat(sprintf("  Sequential bootstrap elapsed: %.1f s\n\n", t_boot_seq["elapsed"]))

cat("\n--- Bootstrap CI columns (from return value) ---\n")
boot_cols <- grep("^Bootstrap", colnames(boot_final_seq$summary), value = TRUE)
print(boot_final_seq$summary[, boot_cols])

cat("\n--- Same results accessible via fit_final$bootstrap$summary ---\n")
print(fit_final$bootstrap$summary[, boot_cols])

cat("\n--- Per-replicate results (first 6 rows) ---\n")
print(head(boot_final_seq$results))

cat(sprintf("\nOutput folder: %s\n", boot_final_seq$outputDir))

# ── 7b. Parallel bootstrap ────────────────────────────────────────────────────
if (has_parallel) {
  n_workers_boot <- min(4L, future::availableCores(omit = 1L))
  cat(sprintf(
    "\nRunning parallel bootstrap (%d workers)...\n",
    n_workers_boot
  ))

  t_boot_par <- system.time(
    boot_final_par <- runBootstrap(
      fit      = fit_final,
      nboot    = 50,
      restart  = TRUE,
      fitName  = "warf_final_par",
      workers  = n_workers_boot
    )
  )
  cat(sprintf("  Parallel bootstrap elapsed: %.1f s\n\n", t_boot_par["elapsed"]))

  boot_cols_par <- grep(
    "^Bootstrap", colnames(boot_final_par$summary), value = TRUE
  )
  cat("\n--- Parallel bootstrap: Bootstrap CI columns ---\n")
  print(boot_final_par$summary[, boot_cols_par])

} else {
  t_boot_par    <- NULL
  n_workers_boot <- NA_integer_
  cat("\nSkipping parallel bootstrap: future/future.apply/progressr missing.\n")
  cat("Install: install.packages(c('future', 'future.apply', 'progressr'))\n")
}

# ── 7c. Bootstrap on the base model (no covariates) ───────────────────────────
cat("\n--- Base model original parameter estimates ---\n")
orig_cols <- intersect(
  c("Estimate", "SE", "%RSE", "Back-transformed"),
  colnames(fit_base$parFixedDf)
)
print(fit_base$parFixedDf[, orig_cols])

cat("\nRunning 20-replicate bootstrap on the base model (no covariates)...\n")
boot_base <- runBootstrap(
  fit     = fit_base,
  nboot   = 20,
  restart = TRUE,
  fitName = "warf_base_boot"
)

cat("\n--- Bootstrap summary (fit_base unchanged; results in boot_base) ---\n")
boot_cols_base <- grep(
  "^Bootstrap", colnames(boot_base$summary), value = TRUE
)
print(boot_base$summary[, boot_cols_base])

# =============================================================================
# Section 8 — Timing summary
# =============================================================================

cat("\n\n=======================================================\n")
cat("  TIMING SUMMARY (section 8)\n")
cat("=======================================================\n\n")

.fmt_time <- function(t) {
  if (is.null(t)) return("skipped")
  sprintf("%.1f s", t["elapsed"])
}

.speedup <- function(t_seq, t_par) {
  if (is.null(t_par)) return("—")
  sprintf("%.2fx", t_seq["elapsed"] / t_par["elapsed"])
}

timing_df <- data.frame(
  Task       = c(
    sprintf("SCM (1 worker)"),
    sprintf("SCM (%d workers)", n_workers_scm),
    sprintf("Bootstrap 50 (1 worker)"),
    sprintf("Bootstrap 50 (%d workers)", if (has_parallel) n_workers_boot else NA_integer_)
  ),
  Elapsed    = c(
    .fmt_time(t_scm_seq),
    .fmt_time(t_scm_par),
    .fmt_time(t_boot_seq),
    .fmt_time(t_boot_par)
  ),
  Speedup    = c(
    "1.00x",
    .speedup(t_scm_seq, t_scm_par),
    "1.00x",
    .speedup(t_boot_seq, t_boot_par)
  ),
  stringsAsFactors = FALSE
)

print(timing_df, row.names = FALSE)
cat("\n(Speedup = sequential elapsed / parallel elapsed; >1 means parallel was faster)\n")
