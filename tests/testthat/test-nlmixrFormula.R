test_that(".nlmixrFormulaParser breaks the formula up into the correct bits", {
  # without random effects
  expect_equal(
    .nlmixrFormulaParser(object = y ~ m*x + b),
    list(
      DV=as.name("y"),
      predictor=list(str2lang("m*x + b")),
      ranef=NULL
    )
  )
  # with random effects
  expect_equal(
    .nlmixrFormulaParser(object = y ~ m*x + b ~ (m|c)),
    list(
      DV=as.name("y"),
      predictor=list(str2lang("m*x + b")),
      ranef=
        list(
          list(
            ranefVar=as.name("m"),
            ranefGroup=as.name("c"),
            start=1
          )
        )
    )
  )
})

test_that(".nlmixrFormulaParser gives expected errors for invalid formula", {
  expect_error(
    .nlmixrFormulaParser(object = ~ m*x + b),
    regexp = "formula must be two-sided"
  )

  # a weird amalgam of one-sided and two-sided formula that looks like a
  # two-sided formula with simple parsing
  expect_error(
    .nlmixrFormulaParser(object = ~ m*x + b ~ c),
    regexp = "formula left-hand-side must be a single variable"
  )

  expect_error(
    .nlmixrFormulaParser(object = m*x + b ~ c),
    regexp = "formula left-hand-side must be a single variable, not m * x + b",
    fixed = TRUE
  )
})

test_that(".nlmixrFormulaParserRanef correctly parses random effects", {
  # single random effect
  expect_equal(
    .nlmixrFormulaParserRanef(str2lang("c|id")),
    list(
      list(
        ranefVar=as.name("c"),
        ranefGroup=as.name("id"),
        start=1
      )
    )
  )
  # single random effect (grouped with parentheses)
  expect_equal(
    .nlmixrFormulaParserRanef(str2lang("(c|id)")),
    list(
      list(
        ranefVar=as.name("c"),
        ranefGroup=as.name("id"),
        start=1
      )
    )
  )
  expect_equal(
    .nlmixrFormulaParserRanef(str2lang("(c|id)+(d|id2)")),
    list(
      list(
        ranefVar=as.name("c"),
        ranefGroup=as.name("id"),
        start=1
      ),
      list(
        ranefVar=as.name("d"),
        ranefGroup=as.name("id2"),
        start=1
      )
    )
  )
})

test_that(".nlmixrFormulaParserRanef expected errors", {
  expect_error(.nlmixrFormulaParserRanef("A"))
  expect_error(
    .nlmixrFormulaParserRanef(str2lang("a*b")),
    regexp = "Invalid random effect: a * b",
    fixed = TRUE
  )
})

test_that("nlmixrFormula creates factor parameters correctly", {
  # Base case with no reordering
  d_factor <- data.frame(A=factor(c("A", "B")))
  expect_equal(
    .nlmixrFormulaExpandStartParamFactor(startName="myest", startValue=1, param="A", data=d_factor),
    list(
      ini=
        list(
          str2lang('myest.A.A <- 1'),
          str2lang('myest.A.B <- 0')
        ),
      model=list(str2lang('myest <- myest.A.A + myest.A.B * (A == "B")'))
    )
  )

  # reorder factors based on prevalence
  d_factor <- data.frame(A=factor(c("A", "B", "B")))
  expect_message(
    v1 <- .nlmixrFormulaExpandStartParamFactor(startName="myest", startValue=1, param="A", data=d_factor),
    regexp = "ordering the parameters by factor frequency: myest with parameter A"
  )
  expect_equal(
    v1,
    list(
      ini=
        list(
          str2lang('myest.A.B <- 1'),
          str2lang('myest.A.A <- 0')
        ),
      model=list(str2lang('myest <- myest.A.B + myest.A.A * (A == "A")'))
    )
  )

  # do not reorder when start is the same length as the number of factors
  d_factor <- data.frame(A=factor(c("A", "B", "B")))
  expect_message(
    v1 <- .nlmixrFormulaExpandStartParamFactor(startName="myest", startValue=c(1, 2), param="A", data=d_factor),
    NA
  )
  expect_equal(
    v1,
    list(
      ini=
        list(
          str2lang('myest.A.A <- 1'),
          str2lang('myest.A.B <- 2')
        ),
      model=list(str2lang('myest <- myest.A.A + myest.A.B * (A == "B")'))
    )
  )

  # do not reorder when the factors are ordered
  d_factor <- data.frame(A=ordered(c("A", "B", "B")))
  expect_message(
    v1 <- .nlmixrFormulaExpandStartParamFactor(startName="myest", startValue=1, param="A", data=d_factor),
    NA
  )
  expect_equal(
    v1,
    list(
      ini=
        list(
          str2lang('myest.A.A <- 1'),
          str2lang('myest.A.B <- 0')
        ),
      model=list(str2lang('myest <- myest.A.A + myest.A.B * (A == "B")'))
    )
  )
})

test_that(".nlmixrFormulaSetupIniRandom", {
  expect_equal(
    .nlmixrFormulaSetupIniRandom(NULL),
    str2lang("{}")
  )
  expect_equal(
    .nlmixrFormulaSetupIniRandom(
      list(
        list(ranefVar = "foo", start = 1.1),
        list(ranefVar = "bar", start = 2.2)
      )
    ),
    str2lang("{
             foo ~ 1.1
             bar ~ 2.2
             }"
    )
  )
})

test_that(".nlmixrFormulaExpandStartParamSingle", {
  expect_equal(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, param = NA, data = data.frame(z = factor(c("a", "b")))),
    list(
      ini = list(str2lang("b <- 5")),
      model = list(NULL)
    )
  )
  # `data` has no effect if is.na(param)
  expect_equal(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, param = NA, data = data.frame(z = factor(c("a", "b")))),
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, param = NA)
  )

  expect_error(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, param = "z", data = NULL),
    regexp = "data must be given when parameters are not single fixed effects"
  )
  expect_error(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, param = "z", data = data.frame(z = factor(c("a", "b", NA)))),
    regexp = "NA found in data column: z"
  )
  expect_error(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, param = "z", data = data.frame(z = c("a", "b"))),
    regexp = "Can only handle factors for fixed effect grouping levels"
  )

  expect_equal(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, param = "z", data = data.frame(z = factor(c("a", "b")))),
    .nlmixrFormulaExpandStartParamFactor(startName = "b", startValue = 5, param = "z", data = data.frame(z = factor(c("a", "b"))))
  )
})

test_that(".nlmixrFormulaExpandStartParamFactor", {
  expect_equal(
    .nlmixrFormulaExpandStartParamFactor(startName = "b", startValue = 5, param = "z", data = data.frame(z = factor(c("a", "b")))),
    list(
      ini = list(str2lang("b.z.b <- 5"), str2lang("b.z.a <- 0")),
      model = list(str2lang('b <- b.z.b + b.z.a * (z == "a")'))
    )
  )
})

test_that(".nlmixrFormulaSetupModel", {
  expect_equal(
    .nlmixrFormulaSetupModel(
      start =
        list(
          list(
            ini = list(str2lang("a <- 1")),
            model = list(NULL)
          ),
          list(
            ini = list(str2lang("b <- 2")),
            model = list(NULL)
          )
        ),
      predictor = list(str2lang("a*x + b*y + z")),
      residualModel = ~add(addSd)
    ),
    str2lang("
    {
      value <- a * x + b * y + z
      value ~ add(addSd)
    }
    ")
  )
})
