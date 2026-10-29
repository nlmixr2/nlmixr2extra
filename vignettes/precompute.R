# Populate (or refresh) vignettes/cache/ for the articles that cache their
# fits and results with nlmixr2save's `:=`.
#
#   cd vignettes && Rscript precompute.R            # compute only what is missing
#   cd vignettes && Rscript precompute.R --clean    # clear cache/ and recompute all
#   cd vignettes && Rscript precompute.R articles/cui-met.Rmd   # selected articles
#
# Each article renders in its own fresh R process: building many rxode2
# models in one long-lived session eventually fails.  The package is loaded
# from the source tree (../) when pkgload is installed, so the articles see
# the development version; otherwise the installed nlmixr2extra is used.
# Commit vignettes/cache/<prefix>* alongside the .Rmd.

vignettes <- c(
  "reporting-helpers.Rmd",
  "articles/likelihood-profiling.Rmd",
  "articles/model-linearization.Rmd",
  "articles/multistart.Rmd",
  "articles/nlmixrFormula.Rmd",
  "articles/precondition.Rmd",
  "articles/clinical-utility-index.Rmd",
  "articles/cui-met.Rmd"
)

args <- commandArgs(trailingOnly = TRUE)
if (!file.exists("precompute.R")) stop("run from the vignettes/ directory", call. = FALSE)
if ("--clean" %in% args) {
  unlink(list.files("cache", full.names = TRUE))
  message("cleared cache/")
}
selected <- setdiff(args, "--clean")
if (length(selected) > 0) {
  unknown <- setdiff(selected, vignettes)
  if (length(unknown) > 0) {
    stop("not a registered article: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  vignettes <- selected
}
dir.create("cache", showWarnings = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")
for (v in vignettes) {
  message("rendering ", v)
  code <- sprintf(paste0(
    "if (requireNamespace('pkgload', quietly = TRUE)) pkgload::load_all('..', quiet = TRUE); ",
    "rmarkdown::render('%s', output_dir = tempdir(), quiet = TRUE)"
  ), v)
  status <- system2(rscript, c("-e", shQuote(code)))
  if (status != 0) stop("rendering ", v, " failed", call. = FALSE)
}
message("done; commit cache/ with the articles")
