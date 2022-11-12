#' @importFrom equatiomatic extract_eq
#' @export
equatiomatic::extract_eq

#' Extract the equations from an nlmixr2/rxode2 model to produce a 'LaTeX'
#' equation.
#'
#' @param model The model to extract equations from
#' @param intercept,greek,greek_colors,subscript_colors,var_colors,var_subscript_colors,raw_tex,swap_var_names,swap_subscript_names,ital_vars,label,index_factors,show_distribution,wrap,terms_per_line,operator_location,align_env,use_coefs,coef_digits,fix_signs,font_size,mean_separate,return_variances,se_subscripts Ignored
#' @export
extract_eq.nlmixr2FitCore <- function(model, intercept = "alpha", greek = "beta",
                                      greek_colors = NULL, subscript_colors = NULL,
                                      var_colors = NULL, var_subscript_colors = NULL,
                                      raw_tex = FALSE,
                                      swap_var_names = NULL, swap_subscript_names = NULL,
                                      ital_vars = FALSE, label = NULL,
                                      index_factors = FALSE, show_distribution = FALSE,
                                      wrap = FALSE, terms_per_line = 4,
                                      operator_location = "end", align_env = "aligned",
                                      use_coefs = FALSE, coef_digits = 2,
                                      fix_signs = TRUE, font_size = NULL,
                                      mean_separate, return_variances = FALSE,
                                      se_subscripts = FALSE, ...) {
  stopifnot(identical(intercept, "alpha"))
  stopifnot(identical(greek, "beta"))
  stopifnot(is.null(greek_colors))
  stopifnot(is.null(subscript_colors))
  stopifnot(is.null(var_colors))
  stopifnot(is.null(var_subscript_colors))
  stopifnot(identical(raw_tex, FALSE))
  stopifnot(is.null(swap_var_names))
  stopifnot(is.null(swap_subscript_names))
  stopifnot(identical(ital_vars, FALSE))
  stopifnot(is.null(label))
  stopifnot(identical(index_factors, FALSE))
  stopifnot(identical(show_distribution, FALSE))
  stopifnot(identical(wrap, FALSE))
  stopifnot(identical(terms_per_line, 4))
  stopifnot(identical(operator_location, "end"))
  stopifnot(identical(align_env, "aligned"))
  stopifnot(identical(use_coefs, FALSE))
  stopifnot(identical(coef_digits, 2))
  stopifnot(identical(fix_signs, TRUE))
  stopifnot(is.null(font_size))
  stopifnot(missing(mean_separate))
  stopifnot(identical(return_variances, FALSE))
  stopifnot(identical(se_subscripts, FALSE))
  ret <-
    paste0(
      "\\begin{align*}\n",
      paste(
        extractEqHelper(x = model, ...),
        collapse = " \\\\\n"
      ),
      "\n\\end{align*}\n"
    )
  knitr::asis_output(ret)
}

#' @rdname extract_eq.nlmixr2FitCore
#' @export
extract_eq.rxUi <- function(model, intercept = "alpha", greek = "beta",
                            greek_colors = NULL, subscript_colors = NULL,
                            var_colors = NULL, var_subscript_colors = NULL,
                            raw_tex = FALSE,
                            swap_var_names = NULL, swap_subscript_names = NULL,
                            ital_vars = FALSE, label = NULL,
                            index_factors = FALSE, show_distribution = FALSE,
                            wrap = FALSE, terms_per_line = 4,
                            operator_location = "end", align_env = "aligned",
                            use_coefs = FALSE, coef_digits = 2,
                            fix_signs = TRUE, font_size = NULL,
                            mean_separate, return_variances = FALSE,
                            se_subscripts = FALSE, ...) {
  extract_eq.nlmixr2FitCore(model, intercept = intercept, greek = greek,
                            greek_colors = greek_colors, subscript_colors = subscript_colors,
                            var_colors = var_colors, var_subscript_colors = var_colors,
                            raw_tex = raw_tex,
                            swap_var_names = swap_var_names, swap_subscript_names = swap_subscript_names,
                            ital_vars = ital_vars, label = label,
                            index_factors = index_factors, show_distribution = show_distribution,
                            wrap = wrap, terms_per_line = terms_per_line,
                            operator_location = operator_location, align_env = align_env,
                            use_coefs = use_coefs, coef_digits = coef_digits,
                            fix_signs = fix_signs, font_size = font_size,
                            mean_separate = mean_separate, return_variances = return_variances,
                            se_subscripts = se_subscripts, ...)
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

extractEqHelper.nlmixr2FitCore <- function(x, ..., inModel) {
  extractEqHelper(as.function(x), ..., inModel=FALSE)
}

extractEqHelper.rxUi <- function(x, ..., inModel, name) {
  extractEqHelper(as.function(x), ..., inModel=FALSE)
}

extractEqHelper.function <- function(x, ..., inModel, name) {
  extractEqHelper(methods::functionBody(x), ..., inModel=FALSE)
}

"extractEqHelper.{" <- function(x, ..., inModel, name) {
  extractEqHelperSeqDrop(x, ..., inModel = inModel)
}

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
    "<="="\\leq",
    "=="="\\equiv",
    ">="="\\geq",
    ">"=">",
    "&"="\\land",
    "&&"="\\land",
    "|"="\\lor",
    "||"="\\lor",
    "!="="\\ne",
    "!"="\\lnot"
  )

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
  } else {
    ret <- character()
  }
  ret
}

extractEqHelper.numeric <- function(x, ..., inModel, name = NULL) {
  if (inModel) {
    ret <- format(x)
    # Convert SI to exponents
    patternSi <- "^([0-9]+(?:\\.[0-9]+)?)e([+-][0-9]+)$"
    if (grepl(x = ret, pattern = patternSi)) {
      base <- gsub(x = ret, pattern = patternSi, replacement = "\\1")
      # as.numeric will remove the leading zeros
      exponent <- as.numeric(gsub(x = ret, pattern = patternSi, replacement = "\\2"))
      ret <- sprintf("{%s \\times 10^{%d}}", base, exponent)
    }
    if (!is.null(name) && !(name == "")) {
      # Hopefully this is right; it occurs in named arguments to calls
      ret <- sprintf("%s=%s", name, ret)
    }
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
