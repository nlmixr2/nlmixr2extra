#' @importFrom knitr knit_print
#' @export
knitr::knit_print

#' Extract the equations from an nlmixr2/rxode2 model to produce a 'LaTeX'
#' equation.
#'
#' @param x The model to extract equations from
#' @param ... Ignored
#' @param output The type of output to request (currently, just "equations")
#' @export
knit_print.nlmixr2FitCore <- function(x,  ..., output = "equations") {
  output <- match.arg(output)
  if ("equations" %in% output) {
    ret <-
      paste0(
        "\\begin{align*}\n",
        paste(
          extractEqHelper(x = x, ...),
          collapse = " \\\\\n"
        ),
        "\n\\end{align*}\n"
      )
  }
  knitr::asis_output(ret)
}

#' @rdname knit_print.nlmixr2FitCore
#' @export
knit_print.rxUi <- function(x, ...) {
  knit_print.nlmixr2FitCore(x, ...)
}

extractEqHelper <- function(x, ..., inModel) {
  UseMethod("extractEqHelper")
}

# Ignore one (or more) items in a sequence.  Typically used for {} and ()
extractEqHelperSeqDrop <- function(x, ..., inModel, dropIdx = 1) {
  ret <- character()
  retNames <- names(x)
  for (idx in setdiff(seq_along(x), dropIdx)) {
    ret <- c(ret, extractEqHelper(x[[idx]], name=retNames[[idx]], ..., inModel = inModel))
  }
  ret
}

# Generate something with a left and right hand side
extractEqHelperLhsRhs <- function(x, ..., inModel, prefix = "", middle, suffix = "", lhsForce = NULL, rhsForce = NULL) {
  stopifnot(length(prefix) == 1)
  stopifnot(length(middle) == 1)
  stopifnot(length(suffix) == 1)
  if (length(x) != 3) {
    stop("extractEqHelperLhsRhs requires length of 3, please report a bug") # nocov
  }
  stopifnot(length(x) == 3)
  if (!is.null(lhsForce)) {
    lhs <- lhsForce
  } else {
    lhs <- extractEqHelper(x[[2]], ..., inModel = inModel)
  }
  if (!is.null(rhsForce)) {
    rhs <- rhsForce
  } else {
    rhs <- extractEqHelper(x[[3]], ..., inModel = inModel)
  }
  stopifnot(length(lhs) == 1)
  stopifnot(length(rhs) == 1)
  paste0(prefix, lhs, middle, rhs, suffix)
}

#' @export
extractEqHelper.nlmixr2FitCore <- function(x, ..., inModel) {
  extractEqHelper(as.function(x), ..., inModel=FALSE)
}

#' @export
extractEqHelper.rxUi <- function(x, ..., inModel, name) {
  extractEqHelper(as.function(x), ..., inModel=FALSE)
}

#' @export
extractEqHelper.function <- function(x, ..., inModel, name) {
  extractEqHelper(methods::functionBody(x), ..., inModel=FALSE)
}

#' @export
"extractEqHelper.{" <- function(x, ..., inModel, name) {
  extractEqHelperSeqDrop(x, ..., inModel = inModel)
}

#' @export
"extractEqHelper.(" <- function(x, ..., inModel, dropParen = FALSE, name) {
  stopifnot(length(x) == 2)
  if (dropParen) {
    ret <- extractEqHelper(x[[2]], ..., inModel = inModel)
  } else {
    # This inherently respects inModel = FALSE because sprintf() with
    # character(0) input returns character(0)
    ret <-
      sprintf(
        "\\left(%s\\right)",
        extractEqHelper(x[[2]], ..., inModel = inModel)
      )
  }
  ret
}

matchesTemplate <- function(x, template) {
  ret <- class(x) == class(template)
  if (ret) {
    if (length(x) == length(template)) {
      if (length(x) > 1) {
        for (idx in seq_along(x)) {
          ret <- ret && matchesTemplate(x[[idx]], template[[idx]])
        }
      } else if (is.name(x)) {
        if (as.character(template) != "name") {
          # Check for a value if the name is not "name"
          ret <- x == template
        }
      } else {
        # Do nothing for numeric, character, etc.
      }
    } else {
      ret <- FALSE
    }
  }
  ret
}

isOdeAssign <- function(x) {
  matchesTemplate(x[[2]], str2lang("d/dt(name)"))
}

# @param alignment gives the text to provide before the equal sign to align on
#   the sign
# This function is named `extractEqHelperAssign` instead of `extractEqHelper.<-`
# to fix an R CMD check issue.  It is used via the `extractEqHelper.default`
# function.
extractEqHelperAssign <- function(x, ..., inModel, alignment = "&", name) {
  if (inModel) {
    lhsForce <- NULL
    rhsForce <- NULL
    if (isOdeAssign(x)) {
      cmtName <- as.character(x[[2]][[3]][[2]])
      if (nchar(cmtName) > 1) {
        # \: forces a medium space in LaTeX equations, do it for state names
        # with more than one character.
        spacesep <- " \\: "
      } else {
        spacesep <- ""
      }
      lhsForce <- sprintf("\\frac{d%s%s}{dt}", spacesep, cmtName)
      rhsForce <- extractEqHelper(x[[3]], ..., inModel = inModel)
    }
    middleForce <- sprintf(" %s = ", alignment)
    ret <-
      extractEqHelperLhsRhs(
        x, ..., inModel = inModel,
        middle = middleForce,
        lhsForce = lhsForce,
        rhsForce = rhsForce
      )
  } else {
    ret <- character()
  }
  ret
}

latexOpMap <-
  list(
    "<"="<",
    "<="="{\\leq}",
    "=="="{\\equiv}",
    ">="="{\\geq}",
    ">"=">",
    "&"="{\\land}",
    "&&"="{\\land}",
    "|"="{\\lor}",
    "||"="{\\lor}",
    "!="="{\\ne}",
    "!"="{\\lnot}"
  )

#' @export
extractEqHelper.call <- function(x, ..., inModel, name) {
  if (inModel) {
    if (is.name(x[[1]])) {
      x1c <- as.character(x[[1]])
      if (x1c %in% c("exp", "log")) {
        stopifnot(length(x) == 2)
        ret <-
          sprintf(
            "\\%s\\left(%s\\right)",
            x1c,
            extractEqHelper(x[[2]], ..., inModel = inModel)
          )
        stopifnot(length(ret) == 1)
      } else if (x1c %in% c("<", "<=",  "==", ">=", ">", "&", "&&", "|", "||", "!=")) {
        # binary logical operators
        ret <- extractEqHelperLhsRhs(x, ..., inModel = inModel, middle = latexOpMap[[x1c]])
      } else if (x1c %in% "!") {
        # unary logical operator
        stopifnot(length(x) == 2)
        ret <-
          sprintf(
            "%s %s",
            latexOpMap[[x1c]],
            extractEqHelper(x[[2]], ..., inModel = inModel)
          )
      } else if (x1c %in% c("+", "-")) {
        if (length(x) == 2) {
          # Handle unary plus and minus
          ret <- sprintf("%s%s", x1c, extractEqHelper(x[[2]], ..., inModel = inModel))
        } else if (length(x) == 3) {
          # Handle binary plus and minus
          ret <- extractEqHelperLhsRhs(x, ..., inModel = inModel, middle = x1c)
        } else {
          cli::cli_abort("Cannot handle addition or subtraction with length that is not 2 or 3") # nocov
        }
      } else if (x1c == "*") {
        ret <- extractEqHelperLhsRhs(x, ..., inModel = inModel, middle = " {\\times} ")
      } else if (x1c == "/") {
        ret <-
          extractEqHelperLhsRhs(
            x, ..., inModel = inModel,
            prefix = "\\frac{",
            middle = "}{",
            suffix = "}"
          )
      } else if (x1c %in% c("**", "^")) {
        ret <-
          extractEqHelperLhsRhs(
            x, ..., inModel = inModel,
            prefix = "{",
            middle = "}^{",
            suffix = "}"
          )
      } else if (x1c == "~") {
        ret <-
          extractEqHelperLhsRhs(
            x, ...,
            inModel = inModel,
            middle = " & \\sim "
          )
      } else {
        # A generic call that is printed using parentheses
        ret <-
          sprintf(
            "%s(%s)",
            x1c,
            paste(
              extractEqHelperSeqDrop(x, ..., inModel = inModel),
              collapse = ", "
            )
          )
      }
    } else {
      stop("a call that does not start with a name is in the model, please report a bug") # nocov
    }
  } else if (length(x) == 1) {
    stop("a one-length call is not inModel, please report a bug") # nocov
  } else {
    if (is.name(x[[1]])) {
      if (x[[1]] == as.name("model")) {
        ret <- extractEqHelperSeqDrop(x, ..., inModel = TRUE, dropBrace = TRUE)
      } else {
        # When not in the model, calls don't matter unless they are `model()`
        ret <- character()
      }
    } else {
      stop("a call that does not start with a name does not have inModel set, please report a bug") # nocov
    }
  }
  ret
}

#' @export
extractEqHelper.name <- function(x, ..., inModel, underscoreToSubscript = FALSE, name) {
  if (inModel) {
    ret <- as.character(x)
    if (underscoreToSubscript && grepl(x = ret, pattern = "_", fixed = TRUE)) {
      firstSecondPattern <- "^(.*?)_(.*)$"
      # The first underscore becomes a subscript
      firstPart <- gsub(x = ret, pattern = firstSecondPattern, replacement = "\\1")
      # Remaining underscores become commas
      secondPartPrep <- gsub(x = ret, pattern = firstSecondPattern, replacement = "\\2")
      secondPart <- gsub(x = secondPartPrep, pattern = "_", replacement = ", ", fixed = TRUE)
      ret <- sprintf("%s_{%s}", firstPart, secondPart)
    } else {
      # Otherwise underscores are protected from LaTeX interpretation
      ret <- gsub(x = ret, pattern = "_", replacement = "\\_", fixed = TRUE)
    }
    ret <- paste0("{", ret, "}")
  } else {
    ret <- character()
  }
  ret
}

#' @export
extractEqHelper.numeric <- function(x, ..., inModel, name = NULL) {
  if (inModel) {
    ret <- format(x)
    # Convert SI to exponents
    patternSi <- "^([0-9]+(?:\\.[0-9]+)?)e([+-][0-9]+)$"
    if (grepl(x = ret, pattern = patternSi)) {
      base <- gsub(x = ret, pattern = patternSi, replacement = "\\1")
      # as.numeric will remove the leading zeros
      exponent <- as.numeric(gsub(x = ret, pattern = patternSi, replacement = "\\2"))
      ret <- sprintf("%s \\times 10^{%d}", base, exponent)
    }
    if (!is.null(name) && !(name == "")) {
      # Hopefully this is right; it occurs in named arguments to calls
      ret <- sprintf("%s=%s", name, ret)
    }
    # Add braces to protect from combination with the prior command
    ret <- paste0("{", ret, "}")
  } else {
    ret <- character()
  }
  ret
}

# copied from knitr:::escape_latex
escapeLatex <- function (x, newlines = FALSE, spaces = FALSE) {
  x = gsub("\\\\", "\\\\textbackslash", x)
  x = gsub("([#$%&_{}])", "\\\\\\1", x)
  x = gsub("\\\\textbackslash", "\\\\textbackslash{}", x)
  x = gsub("~", "\\\\textasciitilde{}", x)
  x = gsub("\\^", "\\\\textasciicircum{}", x)
  if (newlines)
    x = gsub("(?<!\n)\n(?!\n)", "\\\\\\\\", x, perl = TRUE)
  if (spaces)
    x = gsub("(?<= ) ", "\\\\ ", x, perl = TRUE)
  x
}

#' @export
extractEqHelper.character <- function(x, ..., inModel, name = NULL) {
  if (inModel) {
    ret <- sprintf('\\text{"%s"}', escapeLatex(x))
    if (!is.null(name) && !(name == "")) {
      # Hopefully this is right; it occurs in named arguments to calls
      ret <- sprintf("%s=%s", name, ret)
    }
  } else {
    ret <- character()
  }
  ret
}

braceWrap <- function(x, ..., inModel, indent) {
  if (inherits(x, "{")) {
    innerPart <- extractEqHelperSeqDrop(x, ..., inModel = inModel, indent = indent + 1L)
    ret <-
      c(
        "\\{",
        paste(" &", innerPart),
        "\\}"
      )
  } else {
    ret <- extractEqHelper(x, ..., inModel = inModel)
  }
  ret
}

#' Generate LaTeX for if blocks
#'
#' @inheritParams extractEqHelper
#' @param firstIf tracks if this is the first if at the current indent level
#' @param indent tracks the current indention level
#' @noRd
extractEqHelper.if <- function(x, ..., inModel, alignment, indent = 0L, firstIf = 0L, name) {
  # This is used at the same level as assignment
  stopifnot(length(x) %in% c(3, 4))
  # Generate the text for what is in the braces "{}" for this part of the if
  # block
  conditionalPart <- extractEqHelper(x[[2]], ..., inModel = inModel, alignment = "", indent = indent + 1L, firstIf = firstIf + 1L)
  bracedPart <- braceWrap(x[[3]], ..., inModel = inModel, alignment = "", indent = indent + 1L, firstIf = firstIf + 1L)

  if (length(x) == 4) {
    appendPart <- extractEqHelper(x[[4]], ..., inModel = inModel, alignment = "", indent = indent, firstIf = firstIf + 1L)
    appendPart[1] <-
      sprintf(
        " \\quad & \\mathrm{else} \\: %s",
        gsub(x = appendPart[1], pattern = "&", replacement = "", fixed = TRUE)
      )
  } else {
    appendPart <- NULL
  }

  firstPart <- sprintf("\\mathrm{if} & \\left(%s\\right) %s", conditionalPart, bracedPart[1])
  lastPart <-
    if (length(bracedPart) > 1) {
      bracedPart[length(bracedPart)] <- paste(bracedPart[length(bracedPart)], appendPart[1])
      c(bracedPart[-1], appendPart[-1])
    } else {
      appendPart
    }
  ret <- c(firstPart, lastPart)
  ret
}

#' @export
extractEqHelper.default <- function(x, ..., inModel) {
  if (inherits(x, "<-") | inherits(x, "=")) {
    # The assignment classes go via extractEqHelper.default to fix an R CMD
    # check issue:
    #   Warning:     'extractEqHelper.<-'
    #     The argument of a replacement function which corresponds to the right
    #     hand side must be named 'value'.
    ret <- extractEqHelperAssign(x, ..., inModel = inModel)
  } else {
    stop("cannot handle class, please report a bug: ", class(x)[1]) # nocov
  }
  ret
}
