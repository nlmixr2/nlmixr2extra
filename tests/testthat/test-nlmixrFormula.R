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
  # Base case: no reordering needed because frequencies are tied and the
  # stable sort preserves original level order
  d_factor <- data.frame(A=factor(c("A", "B")))
  expect_equal(
    .nlmixrFormulaExpandStartParamFactor(startName="myest", startValue=1, param="A", data=d_factor),
    list(
      ini=
        list(
          str2lang('myest.A.A <- 1'),
          str2lang('myest.A.B <- 0')
        ),
      rhs='myest.A.A + myest.A.B * (A == "B")'
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
      rhs='myest.A.B + myest.A.A * (A == "A")'
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
      rhs='myest.A.A + myest.A.B * (A == "B")'
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
      rhs='myest.A.A + myest.A.B * (A == "B")'
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

test_that(".nlmixrFormulaExpandStartParamSingle scalar fixed effect", {
  expect_equal(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5),
    list(
      ini = list(str2lang("b <- 5")),
      model = list(NULL)
    )
  )
  # `data` has no effect when there is no covariate
  expect_equal(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, data = data.frame(z = factor(c("a", "b")))),
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5)
  )
  # paramLink other than identity is meaningless without a covariate model
  expect_error(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, link = "log"),
    regexp = "paramLink for 'b' is 'log'",
    fixed = TRUE
  )
})

test_that(".nlmixrFormulaExpandStartParamSingle dispatches on column type", {
  # Missing data raises a clear error
  expect_error(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, param = "z", data = NULL),
    regexp = "data must be given when parameters are not single fixed effects"
  )
  # NA in the covariate column
  expect_error(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, param = "z", data = data.frame(z = factor(c("a", "b", NA)))),
    regexp = "NA found in data column: z"
  )
  # Character columns must be converted to factor first
  expect_error(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = 5, param = "z", data = data.frame(z = c("a", "b"))),
    regexp = "Column 'z' in `data` is character; convert it to a factor",
    fixed = TRUE
  )
  # Unsupported column types are rejected
  expect_error(
    .nlmixrFormulaExpandStartParamSingle(
      startName = "b", startValue = 5, param = "z",
      data = data.frame(z = as.Date("2024-01-01") + 0:1)
    ),
    regexp = "Unsupported column type"
  )

  # Factor dispatch produces a single model line of the assembled rhs
  factorOut <- .nlmixrFormulaExpandStartParamSingle(
    startName = "b", startValue = 5, param = "z",
    data = data.frame(z = factor(c("a", "b")))
  )
  expect_equal(
    factorOut,
    list(
      ini = list(str2lang("b.z.a <- 5"), str2lang("b.z.b <- 0")),
      model = list(str2lang('b <- b.z.a + b.z.b * (z == "b")'))
    )
  )
})

test_that(".nlmixrFormulaExpandStartParamFactor", {
  expect_equal(
    .nlmixrFormulaExpandStartParamFactor(startName = "b", startValue = 5, param = "z", data = data.frame(z = factor(c("a", "b")))),
    list(
      ini = list(str2lang("b.z.a <- 5"), str2lang("b.z.b <- 0")),
      rhs = 'b.z.a + b.z.b * (z == "b")'
    )
  )
})

test_that(".nlmixrFormulaExpandStartParamContinuous", {
  # Length-1 startValue: intercept only, slope defaults to 0
  expect_equal(
    .nlmixrFormulaExpandStartParamContinuous(startName = "b", startValue = 5, param = "w", includeIntercept = TRUE),
    list(
      ini = list(str2lang("pop.b <- 5"), str2lang("cov_w_b <- 0")),
      rhs = "pop.b + cov_w_b * w"
    )
  )
  # Length-2 startValue: c(intercept, slope)
  expect_equal(
    .nlmixrFormulaExpandStartParamContinuous(startName = "b", startValue = c(5, 0.3), param = "w", includeIntercept = TRUE),
    list(
      ini = list(str2lang("pop.b <- 5"), str2lang("cov_w_b <- 0.3")),
      rhs = "pop.b + cov_w_b * w"
    )
  )
  # includeIntercept = FALSE: a sister factor supplies the intercept, so the
  # helper only receives the slope value (single number)
  expect_equal(
    .nlmixrFormulaExpandStartParamContinuous(startName = "b", startValue = 0.3, param = "w", includeIntercept = FALSE),
    list(
      ini = list(str2lang("cov_w_b <- 0.3")),
      rhs = "cov_w_b * w"
    )
  )
  # Length > 2 (with intercept) is a hard error
  expect_error(
    .nlmixrFormulaExpandStartParamContinuous(startName = "b", startValue = c(1, 2, 3), param = "w", includeIntercept = TRUE),
    regexp = "must be 1 (intercept only) or 2 (intercept, slope)",
    fixed = TRUE
  )
  # includeIntercept = FALSE requires exactly one slope value
  expect_error(
    .nlmixrFormulaExpandStartParamContinuous(startName = "b", startValue = c(5, 0.3), param = "w", includeIntercept = FALSE),
    regexp = "must contribute a single",
    fixed = TRUE
  )
})

test_that(".nlmixrFormulaExpandStartParamSingle continuous covariate dispatch", {
  d <- data.frame(w = c(1.0, 2.0, 3.0))
  expect_equal(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = c(5, 0.3), param = "w", data = d),
    list(
      ini = list(str2lang("pop.b <- 5"), str2lang("cov_w_b <- 0.3")),
      model = list(str2lang("b <- pop.b + cov_w_b * w"))
    )
  )
  # Log link wraps the rhs in exp()
  expect_equal(
    .nlmixrFormulaExpandStartParamSingle(startName = "b", startValue = c(5, 0.3), param = "w", link = "log", data = d),
    list(
      ini = list(str2lang("pop.b <- 5"), str2lang("cov_w_b <- 0.3")),
      model = list(str2lang("b <- exp(pop.b + cov_w_b * w)"))
    )
  )
})

test_that(".nlmixrFormulaExpandStartParamSingle mixed factor + continuous", {
  d <- data.frame(
    z = factor(c("a", "b", "a", "b")),
    w = c(1.0, 2.0, 3.0, 4.0)
  )
  # start contract for mixed: factor first (one per level), then continuous
  # slopes. Here z has 2 levels and w is one continuous covariate, so length 3.
  out <- .nlmixrFormulaExpandStartParamSingle(
    startName = "b", startValue = c(5, 1, 0.3), param = c("z", "w"), data = d
  )
  expect_equal(
    out$ini,
    list(
      str2lang("b.z.a <- 5"),
      str2lang("b.z.b <- 1"),
      str2lang("cov_w_b <- 0.3")
    )
  )
  expect_equal(
    out$model,
    list(str2lang('b <- b.z.a + b.z.b * (z == "b") + cov_w_b * w'))
  )
  # Length mismatch produces a clear error
  expect_error(
    .nlmixrFormulaExpandStartParamSingle(
      startName = "b", startValue = c(5, 0.3), param = c("z", "w"), data = d
    ),
    regexp = "start must have length 3"
  )
})

test_that(".nlmixrFormulaExpandStartParamSingle log link wraps mixed rhs", {
  d <- data.frame(
    z = factor(c("a", "b")),
    w = c(1.0, 2.0)
  )
  out <- .nlmixrFormulaExpandStartParamSingle(
    startName = "b", startValue = c(log(5), log(2), 0.3),
    param = c("z", "w"), link = "log", data = d
  )
  expect_equal(
    out$model,
    list(str2lang('b <- exp(b.z.a + b.z.b * (z == "b") + cov_w_b * w)'))
  )
})

test_that(".paramExpand", {
  # NULL → empty named list
  expect_equal(.paramExpand(NULL), stats::setNames(list(), character()))
  # Single covariate
  expect_equal(.paramExpand(b ~ z), list(b = "z"))
  # Multiple covariates on one parameter via `+`
  expect_equal(.paramExpand(b ~ z + w), list(b = c("z", "w")))
  # List of formulas combines by name
  expect_equal(
    .paramExpand(list(b ~ z, b ~ w)),
    list(b = c("z", "w"))
  )
  # Different parameters stay separate
  expect_equal(
    .paramExpand(list(b ~ z, m ~ w)),
    list(b = "z", m = "w")
  )
  # Invalid RHS (interaction) is rejected
  expect_error(
    .paramExpand(b ~ z * w),
    regexp = "Invalid right-hand side in `param`"
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

test_that(".renameOrOverwrite errors when the destination column already exists", {
  d <- data.frame(A=1:3, B=4:6)
  expect_error(
    .renameOrOverwrite(d, newName="B", oldName="A"),
    regexp = "Cannot rename column 'A' to 'B'",
    fixed = TRUE
  )
  # No-op rename (newName == oldName) is allowed even when the column exists
  expect_equal(
    .renameOrOverwrite(d, newName="A", oldName="A"),
    d
  )
  # Renaming when only the source exists succeeds and creates the new column
  expect_equal(
    .renameOrOverwrite(data.frame(A=1:3), newName="B", oldName="A"),
    data.frame(A=1:3, B=1:3)
  )
})

# ----------------------------------------------------------------------------
# Integration tests: exercise the full .nlmixrFormulaBuild() pipeline so we
# verify that the parser, parameter expansion, ini/model assembly, and S3
# entrypoint all wire together. These tests do not fit or simulate; they only
# inspect the assembled ini/model bodies.
# ----------------------------------------------------------------------------

test_that("nlmixrFormula integrates: continuous covariate only", {
  set.seed(1)
  d <- data.frame(x = 1:6, y = NA_real_, w = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1))
  built <- .nlmixrFormulaBuild(
    y ~ m*x + b,
    data = d,
    start = list(m = 3, b = c(5, 0.2), addSd = 1),
    param = list(b ~ w)
  )
  expect_equal(
    built$ini,
    str2lang("{
      m <- 3
      pop.b <- 5
      cov_w_b <- 0.2
      addSd <- 1
    }")
  )
  expect_equal(
    built$model,
    str2lang("{
      b <- pop.b + cov_w_b * w
      value <- m * x + b
      value ~ add(addSd)
    }")
  )
})

test_that("nlmixrFormula integrates: factor + continuous covariate together", {
  d <- data.frame(
    x = 1:6, y = NA_real_,
    w = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1),
    z = factor(c("a", "b", "a", "b", "a", "b"))
  )
  built <- .nlmixrFormulaBuild(
    y ~ m*x + b,
    data = d,
    start = list(m = 3, b = c(5, 1, 0.2), addSd = 1),
    param = list(b ~ z + w)
  )
  expect_equal(
    built$ini,
    str2lang("{
      m <- 3
      b.z.a <- 5
      b.z.b <- 1
      cov_w_b <- 0.2
      addSd <- 1
    }")
  )
  expect_equal(
    built$model,
    str2lang('{
      b <- b.z.a + b.z.b * (z == "b") + cov_w_b * w
      value <- m * x + b
      value ~ add(addSd)
    }')
  )
})

test_that("nlmixrFormula integrates: log link wraps the assembled rhs", {
  d <- data.frame(x = 1:6, y = NA_real_, w = c(0.1, 0.3, 0.5, 0.7, 0.9, 1.1))
  built <- .nlmixrFormulaBuild(
    y ~ m*x + b,
    data = d,
    start = list(m = 3, b = c(log(5), 0.2), addSd = 1),
    param = list(b ~ w),
    paramLink = c(b = "log")
  )
  expect_equal(
    built$model,
    str2lang("{
      b <- exp(pop.b + cov_w_b * w)
      value <- m * x + b
      value ~ add(addSd)
    }")
  )
})

test_that("nlmixrFormula integrates: rich residual model passes through", {
  d <- data.frame(x = 1:6, y = NA_real_)
  built <- .nlmixrFormulaBuild(
    y ~ m*x + b,
    data = d,
    start = list(m = 3, b = 5, addSd = 1, propSd = 0.1),
    residualModel = ~ add(addSd) + prop(propSd)
  )
  # Both sigma names appear in ini and the residual line preserves the user's
  # residual model expression verbatim.
  expect_equal(
    built$ini,
    str2lang("{
      m <- 3
      b <- 5
      addSd <- 1
      propSd <- 0.1
    }")
  )
  # The residual line is the last statement in the model body.
  expect_equal(
    built$model[[length(built$model)]],
    str2lang("value ~ add(addSd) + prop(propSd)")
  )
})

test_that("nlmixrFormula errors when ID column already exists with different values", {
  d <- data.frame(
    id = rep(c("X", "Y"), each = 3),
    ID = rep(c("Z", "W"), each = 3),
    x  = 1:6,
    y  = NA_real_
  )
  expect_error(
    .nlmixrFormulaBuild(
      y ~ m*x + b ~ (bRe|id),
      data = d,
      start = list(m = 3, b = 5, addSd = 1)
    ),
    regexp = "Cannot rename column 'id' to 'ID'",
    fixed = TRUE
  )
})

test_that("nlmixrFormula rejects multiple grouping variables", {
  d <- data.frame(id = 1:6, occ = 1:6, x = 1:6, y = NA_real_)
  expect_error(
    .nlmixrFormulaBuild(
      y ~ m*x + b ~ (mRe|id) + (bRe|occ),
      data = d,
      start = list(m = 3, b = 5, addSd = 1)
    ),
    regexp = "Only one random-effect grouping variable is supported"
  )
})

test_that("nlmixr2.formula S3 dispatch produces an equivalent build", {
  skip_on_cran()
  d <- data.frame(x = 1:6, y = c(1, 2, 3, 4, 5, 6))
  fit_s3 <- nlmixr2(
    y ~ m*x + b,
    data = d,
    start = list(m = 1, b = 1, addSd = 1),
    est = "rxSolve"
  )
  fit_direct <- nlmixrFormula(
    y ~ m*x + b,
    data = d,
    start = list(m = 1, b = 1, addSd = 1),
    est = "rxSolve"
  )
  # class() carries rxode2-attached attributes; strip them before comparing.
  expect_equal(as.character(class(fit_s3)), as.character(class(fit_direct)))
})
