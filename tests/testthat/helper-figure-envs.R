# Helpers for checking what a saved figure carries with it.  serialize()
# writes every environment a ggplot references in full, except namespaces,
# package environments and the global, base and empty environments (written
# as references), so those are the environments that matter here.

# Value bound to `nm` in `env`, or NULL when it is a missing argument (the
# empty symbol, which cannot be assigned to a variable and then used)
.figureEnvValue <- function(env, nm) {
  .l <- tryCatch(mget(nm, envir = env), error = function(e) list(NULL))
  if (is.symbol(.l[[1]]) && !nzchar(as.character(.l[[1]]))) {
    return(NULL)
  }
  .l[[1]]
}

# Add to `acc$envs` every environment that serialize() writes out in full and
# that is reachable from `x` through environment bindings and enclosures,
# closure environments, list elements and attributes (which is how quosures,
# formulas, ggproto objects and S7 properties hold them)
.figureCollectEnvs <- function(x, acc) {
  if (is.environment(x)) {
    if (
      identical(x, globalenv()) ||
        identical(x, baseenv()) ||
        identical(x, emptyenv()) ||
        isNamespace(x) ||
        nzchar(environmentName(x))
    ) {
      return(invisible())
    }
    for (.e in acc$envs) {
      if (identical(.e, x)) {
        return(invisible())
      }
    }
    acc$envs <- c(acc$envs, list(x))
    for (.nm in setdiff(ls(x, all.names = TRUE), "...")) {
      .figureCollectEnvs(.figureEnvValue(x, .nm), acc)
    }
    .figureCollectEnvs(parent.env(x), acc)
  } else if (is.function(x)) {
    .figureCollectEnvs(environment(x), acc)
  } else if (is.list(x) || is.pairlist(x)) {
    for (.el in as.list(x)) {
      .figureCollectEnvs(.el, acc)
    }
  }
  for (.a in attributes(x)) {
    .figureCollectEnvs(.a, acc)
  }
  invisible()
}

# Data frames (more than one row), fits and plots bound in the environments a
# figure references, as "name <class>" strings.  A figure should hold its
# data only in `$data`; one-row data frames are ggplot2's own (like the
# intercept and slope of geom_abline()).
.figureHeldData <- function(fig) {
  .acc <- new.env(parent = emptyenv())
  .acc$envs <- list()
  .figureCollectEnvs(fig, .acc)
  .ret <- character(0)
  for (.e in .acc$envs) {
    for (.nm in setdiff(ls(.e, all.names = TRUE), "...")) {
      .v <- .figureEnvValue(.e, .nm)
      if ((is.data.frame(.v) && nrow(.v) > 1L) || inherits(.v, c("ggplot", "gglist"))) {
        .ret <- c(.ret, paste0(.nm, " <", class(.v)[1], ">"))
      }
    }
  }
  .ret
}
