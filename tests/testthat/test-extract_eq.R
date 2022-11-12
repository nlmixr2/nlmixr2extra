test_that("extract_eq, simple version", {
  mod <- function() {
    ini({
      lka <- 0.45
      lcl <- 1
      lvc  <- 3.45
      propSd <- 0.5
    })
    model({
      ka <- exp(lka)
      cl <- exp(lcl)
      vc  <- exp(lvc)

      cp <- linCmt()
      cp ~ prop(propSd)
    })
  }
  ui <- rxode2::rxode(mod)
  expect_equal(
    extract_eq(ui),
    knitr::asis_output("\\begin{align*}\nka & = \\exp\\left(lka\\right) \\\\\ncl & = \\exp\\left(lcl\\right) \\\\\nvc & = \\exp\\left(lvc\\right) \\\\\ncp & = linCmt() \\\\\ncp & \\sim prop(propSd)\n\\end{align*}\n")
  )

  suppressMessages(
    fit <- nlmixr2est::nlmixr(mod, data = nlmixr2data::theo_sd, est = "focei", control = nlmixr2::foceiControl(eval.max = 1, print = 0))
  )
  expect_equal(
    extract_eq(fit),
    knitr::asis_output("\\begin{align*}\nka & = \\exp\\left(lka\\right) \\\\\ncl & = \\exp\\left(lcl\\right) \\\\\nvc & = \\exp\\left(lvc\\right) \\\\\ncp & = linCmt() \\\\\ncp & \\sim prop(propSd)\n\\end{align*}\n")
  )
})

test_that("extract_eq, less common models", {
  mod <- function() {
    ini({
      lka <- 0.45
      lcl <- 1
      lvc  <- 3.45
    })
    model({
      ka <- exp(lka) < 1
      cl <- exp(lcl) <= 2
      vc  <- exp(lvc) == 3
      vc4 <- vc >= 4
      vc5 <- vc > 5
      vc6 <- vc & 6
      vc7 <- vc && 7
      vc8 <- vc | 8
      vc9 <- vc || 9
      vc10 <- vc != 10
      vc11 <- !vc
      if (vc > 11) {
        cl <- 12
      } else if (vc > 13) {
        cl <- 14
        cl <- 15
      } else {
        cl <- 16
      }

      cp <- linCmt()
      # ordinal model
      cp ~ c(p0=0, p1=1, p2=2, 3)
    })
  }
  ui <- rxode2::rxode(mod)
  expect_equal(
    extract_eq(ui),
    knitr::asis_output("\\begin{align*}\nka & = \\exp\\left(lka\\right)<1 \\\\\ncl & = \\exp\\left(lcl\\right)\\leq2 \\\\\nvc & = \\exp\\left(lvc\\right)\\equiv3 \\\\\nvc4 & = vc\\geq4 \\\\\nvc5 & = vc>5 \\\\\nvc6 & = vc\\land6 \\\\\nvc7 & = vc\\land7 \\\\\nvc8 & = vc\\lor8 \\\\\nvc9 & = vc\\lor9 \\\\\nvc10 & = vc\\ne10 \\\\\nvc11 & = \\lnot vc \\\\\n\\mathrm{if} & \\left(vc>11\\right) \\{ \\\\\n & cl  = 12 \\\\\n\\}  \\quad & \\mathrm{else} \\: \\mathrm{if}  \\left(vc>13\\right) \\{ \\\\\n & cl  = 14 \\\\\n & cl  = 15 \\\\\n\\}  \\quad & \\mathrm{else} \\: cl  = 16 \\\\\ncp & = linCmt() \\\\\ncp & \\sim c(p0=0, p1=1, p2=2, 3)\n\\end{align*}\n")
  )
})

expect_equal(
  extractEqHelper.function(function() {ini({a <- 1});model({log(foo)+exp(bar)})}, inModel = FALSE),
  "\\log\\left(foo\\right)+\\exp\\left(bar\\right)"
)

test_that("extractEqHelper.(", {
  expect_equal(
    "extractEqHelper.("(str2lang("(foo)"), inModel = TRUE),
    "\\left(foo\\right)"
  )
  expect_equal(
    "extractEqHelper.("(str2lang("(foo)"), inModel = FALSE),
    character()
  )
})

test_that("extractEqHelper.= and extractEqHelper.<-", {
  # ODE, space with multiple-character state ("abc")
  expect_equal(
    `extractEqHelper.=`(str2lang("d/dt(abc) = d"), inModel = TRUE),
    "\\frac{d \\: abc}{dt} & = d"
  )
  # ODE, no space with single-character state ("a")
  expect_equal(
    `extractEqHelper.=`(str2lang("d/dt(a) = b"), inModel = TRUE),
    "\\frac{da}{dt} & = b"
  )
  expect_equal(
    `extractEqHelper.=`(str2lang("a = b"), inModel = TRUE),
    "a & = b"
  )
  expect_equal(
    `extractEqHelper.<-`(str2lang("a <- b"), inModel = TRUE),
    "a & = b"
  )
  expect_equal(
    `extractEqHelper.<-`(str2lang("a <- b"), inModel = FALSE),
    character()
  )
})

test_that("extractEqHelper.call", {
  # named functions are escaped
  expect_equal(
    extractEqHelper.call(str2lang("model({log(foo)+exp(bar)})"), inModel = FALSE),
    "\\log\\left(foo\\right)+\\exp\\left(bar\\right)"
  )
  # binary plus and unary minus
  expect_equal(
    extractEqHelper.call(str2lang("model({-foo+add(bar)})"), inModel = FALSE),
    "-foo+add(bar)"
  )
  # multiplication
  expect_equal(
    extractEqHelper.call(str2lang("model({foo*add(bar)})"), inModel = FALSE),
    "foo {\\times} add(bar)"
  )
  # division
  expect_equal(
    extractEqHelper.call(str2lang("model({foo/add(bar)})"), inModel = FALSE),
    "\\frac{foo}{add(bar)}"
  )
  # exponent with both operators
  expect_equal(
    extractEqHelper.call(str2lang("model({foo**add(bar)})"), inModel = FALSE),
    "{foo}^{add(bar)}"
  )
  expect_equal(
    extractEqHelper.call(str2lang("model({foo^add(bar)})"), inModel = FALSE),
    "{foo}^{add(bar)}"
  )
  expect_equal(
    extractEqHelper.call(str2lang("model({foo~add(bar)})"), inModel = FALSE),
    "foo & \\sim add(bar)"
  )
  # transition to inModel works
  expect_equal(
    extractEqHelper.call(str2lang("model({foo(bar)})"), inModel = FALSE),
    "foo(bar)"
  )
  expect_equal(
    extractEqHelper.call(str2lang("foo(bar)"), inModel = FALSE),
    character()
  )
})

test_that("extractEqHelper.name", {
  expect_equal(
    extractEqHelper.name(as.name("fo_o_bar"), inModel = TRUE, underscoreToSubscript = FALSE),
    "fo\\_o\\_bar"
  )
  expect_equal(
    extractEqHelper.name(as.name("fo_o_bar"), inModel = TRUE, underscoreToSubscript = TRUE),
    "fo_{o, bar}"
  )
  expect_equal(
    extractEqHelper.name(as.name("fo_o"), inModel = TRUE, underscoreToSubscript = TRUE),
    "fo_{o}"
  )
  expect_equal(
    extractEqHelper.name(as.name("foo"), inModel = TRUE, underscoreToSubscript = TRUE),
    "foo"
  )
  expect_equal(
    extractEqHelper.name(as.name("foo"), inModel = FALSE),
    character()
  )
})

test_that("extractEqHelper.numeric", {
  expect_equal(
    extractEqHelper.numeric(5, inModel = TRUE),
    "5"
  )
  expect_equal(
    extractEqHelper.numeric(pi, inModel = TRUE),
    "3.141593"
  )
  # test SI
  expect_equal(
    extractEqHelper.numeric(1.234567890123456789e15, inModel = TRUE),
    "{1.234568 \\times 10^{15}}"
  )
  expect_equal(
    extractEqHelper.numeric(5, inModel = FALSE),
    character()
  )
})

test_that("extractEqHelper.character", {
  # no escape required
  expect_equal(
    extractEqHelper.character("foo", inModel = TRUE),
    '\\text{"foo"}'
  )
  # escape LaTeX
  expect_equal(
    extractEqHelper.character("fo%o", inModel = TRUE),
    '\\text{"fo\\%o"}'
  )
  # not in model outputs empty
  expect_equal(
    extractEqHelper.character("foo", inModel = FALSE),
    character()
  )
})

test_that("extractEqHelper.if", {
  ifOnly <- str2lang("
  if (vc > 11) {
    cl <- 12
  }")
  ifOnlyNoBrace <- str2lang("if (vc > 11) cl <- 12")
  ifOnlyMultiline <- str2lang("
  if (vc > 11) {
    cl <- 12
    cl2 <- 13
  }")
  ifElse <- str2lang("if (vc > 11) {
    cl <- 12
  } else {
    cl <- 16
  }")
  ifElseIf <- str2lang("if (vc > 11) {
    cl <- 12
  } else if (vc > 13) {
    cl <- 14
    cl <- 15
  } else {
    cl <- 16
  }")
  expect_equal(
    extractEqHelper.if(ifOnly, inModel = TRUE),
    c("\\mathrm{if} & \\left(vc>11\\right) \\{", " & cl  = 12", "\\} ")
  )
  expect_equal(
    extractEqHelper.if(ifOnlyNoBrace, inModel = TRUE),
    "\\mathrm{if} & \\left(vc>11\\right) cl  = 12"
  )
  expect_equal(
    extractEqHelper.if(ifOnlyMultiline, inModel = TRUE),
    c("\\mathrm{if} & \\left(vc>11\\right) \\{", " & cl  = 12", " & cl2  = 13", "\\} ")
  )
  expect_equal(
    extractEqHelper.if(ifElse, inModel = TRUE),
    c("\\mathrm{if} & \\left(vc>11\\right) \\{", " & cl  = 12", "\\}  \\quad & \\mathrm{else} \\: cl  = 16")
  )
  expect_equal(
    extractEqHelper.if(ifElseIf, inModel = TRUE),
    c("\\mathrm{if} & \\left(vc>11\\right) \\{", " & cl  = 12", "\\}  \\quad & \\mathrm{else} \\: \\mathrm{if}  \\left(vc>13\\right) \\{", " & cl  = 14", " & cl  = 15", "\\}  \\quad & \\mathrm{else} \\: cl  = 16")
  )
})
