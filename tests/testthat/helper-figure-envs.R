# Helpers for checking what a saved figure carries with it.  serialize()
# writes every environment a ggplot references in full, except namespaces,
# package environments and the global, base and empty environments (written
# as references), so those are the environments that matter here.

# Value bound to `nm` in `env`, or NULL when it is a missing argument (the
# empty symbol, which cannot be assigned to a variable and then used).  This
# forces a promise bound to `nm`.
.figureEnvValue <- function(env, nm) {
  .l <- tryCatch(mget(nm, envir = env), error = function(e) list(NULL))
  if (is.symbol(.l[[1]]) && !nzchar(as.character(.l[[1]]))) {
    return(NULL)
  }
  .l[[1]]
}

# Add to `acc$envs` (and `acc$todo`) the environments that serialize() meets
# while writing `x`, apart from `x` itself.  serialize() asks `refhook` about
# each environment it would write in full (and each external pointer and weak
# reference, which it writes as usual when the hook returns NULL); answering
# with a string writes only that string, so each environment is visited on its
# own and nothing is evaluated.  This is how the environments of unevaluated
# promises (which serialize() writes with their expressions) are found
# without forcing them.
.figureFindEnvs <- function(x, acc) {
  serialize(x, NULL, refhook = function(e) {
    if (!is.environment(e) || identical(e, x)) {
      return(NULL)
    }
    for (.e in acc$envs) {
      if (identical(.e, e)) {
        return("seen")
      }
    }
    acc$envs <- c(acc$envs, list(e))
    acc$todo <- c(acc$todo, list(e))
    "new"
  })
  invisible()
}

# Every environment that serialize() writes out in full for `x`
.figureEnvs <- function(x) {
  .acc <- new.env(parent = emptyenv())
  .acc$envs <- list()
  .acc$todo <- list()
  .figureFindEnvs(x, .acc)
  while (length(.acc$todo) > 0L) {
    .e <- .acc$todo[[1L]]
    .acc$todo <- .acc$todo[-1L]
    .figureFindEnvs(.e, .acc)
  }
  .acc$envs
}

# Data frames (more than one row), fits and plots bound in the environments a
# figure references, as "name <class>" strings.  A figure should hold its
# data only in `$data`; one-row data frames are ggplot2's own (like the
# intercept and slope of geom_abline()).  Reading a binding forces a promise
# bound there, so the bindings are read from a copy of the figure after all
# of its environments have been found.
.figureHeldData <- function(fig) {
  fig <- unserialize(serialize(fig, NULL))
  .ret <- character(0)
  for (.e in .figureEnvs(fig)) {
    for (.nm in setdiff(ls(.e, all.names = TRUE), "...")) {
      .v <- .figureEnvValue(.e, .nm)
      if ((is.data.frame(.v) && nrow(.v) > 1L) || inherits(.v, c("ggplot", "gglist"))) {
        .ret <- c(.ret, paste0(.nm, " <", class(.v)[1], ">"))
      }
    }
  }
  .ret
}
