.nlmixr2extra <- loadNamespace("nlmixr2extra")

## Newer 'lotri' adds a `prior` column to the ini data frame for prior
## distributions, so `$iniDf` gains a column.  Rows that are built by
## hand here have to work with an iniDf that has the column and with one
## that does not, since which one you get depends on the installed
## version of 'lotri'.  These build both shapes by hand so the test does
## not depend on that version.

.iniDfNoPrior <- function() {
  data.frame(ntheta=c(1L, NA_integer_),
             neta1=c(NA_real_, 1),
             neta2=c(NA_real_, 1),
             name=c("tka", "eta.ka"),
             lower=c(-Inf, -Inf),
             est=c(0.45, 0.6),
             upper=c(Inf, Inf),
             fix=c(FALSE, FALSE),
             label=NA_character_,
             backTransform=NA_character_,
             condition=c(NA_character_, "id"),
             err=NA_character_,
             stringsAsFactors=FALSE)
}

.iniDfWithPrior <- function() {
  .df <- .iniDfNoPrior()
  .df$prior <- c("dnorm(0, 10)", NA_character_)
  .df[, c("ntheta", "neta1", "neta2", "name", "lower", "est", "upper",
          "fix", "label", "backTransform", "condition", "prior", "err")]
}

test_that(".iniDfMatchColumns() lets hand built rows rbind with either shape", {

  .row <- data.frame(ntheta=2L, neta1=NA_real_, neta2=NA_real_,
                     name="cov_wt_cl", lower=-Inf, est=0, upper=Inf,
                     fix=FALSE, label=NA_character_,
                     backTransform=NA_character_, condition=NA_character_,
                     err=NA_character_, stringsAsFactors=FALSE)

  ## a target without `prior` is left alone
  .noPrior <- .nlmixr2extra$.iniDfMatchColumns(.row, .iniDfNoPrior())
  expect_equal(names(.noPrior), names(.iniDfNoPrior()))
  expect_error(rbind(.iniDfNoPrior(), .noPrior), NA)

  ## a target with `prior` gains the column, filled with NA
  .withPrior <- .nlmixr2extra$.iniDfMatchColumns(.row, .iniDfWithPrior())
  expect_equal(names(.withPrior), names(.iniDfWithPrior()))
  expect_true(is.na(.withPrior$prior))
  expect_error(rbind(.iniDfWithPrior(), .withPrior), NA)
  expect_equal(nrow(rbind(.iniDfWithPrior(), .withPrior)), 3L)

  ## an extra column on the row is dropped rather than breaking rbind
  .extra <- .row
  .extra$prior <- "dnorm(0, 1)"
  expect_equal(names(.nlmixr2extra$.iniDfMatchColumns(.extra, .iniDfNoPrior())),
               names(.iniDfNoPrior()))
  expect_error(rbind(.iniDfNoPrior(),
                     .nlmixr2extra$.iniDfMatchColumns(.extra, .iniDfNoPrior())), NA)
})
