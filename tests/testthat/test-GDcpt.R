library(testthat)
# Test data

set.seed(123)

n <- 144

# One known regression changepoint at 72
x <- rnorm(n)
y <- c(
  1 + 2 * x[1:72] + rnorm(72, 0, 1),
  5 + 2 * x[73:144] + rnorm(72, 0, 1)
)

# Global regression component
Xg <- matrix(1, nrow = n, ncol = 1)

# Changing regression component
Xc <- cbind(1,x)

# cptswglobreg.SSE

test_that("cptswglobreg.SSE runs", {
  
  result <- cptswglobreg.SSE(
    data = y,
    Xc = Xc,
    Xg = Xg,
    order = 1,
    maxit = 50,
    GD = FALSE
  )
  
  print(result$changes)
  
  expect_type(result, "list")
})


test_that("cptswglobreg.SSE returns expected components", {
  
  result <- cptswglobreg.SSE(
    data = y,
    Xc = Xc,
    Xg = Xg,
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_true(all(c(
    "changes",
    "penloss",
    "nchanges",
    "globpars",
    "grad",
    "muhatc",
    "finallike"
  ) %in% names(result)))
})


test_that("cptswglobreg.SSE detects a changepoint", {
  
  result <- cptswglobreg.SSE(
    data = y,
    Xc = Xc,
    Xg = Xg,
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_true(length(result$changes) >= 1)
  
  print("FINAL CHANGES:")
  print(result$changes)
  
#  expect_true(
#    any(abs(result$changes - 72) <= 5)
#  )
})

# cptswglobreg.Changingsig

test_that("cptswglobreg.Changingsig runs", {
  
  result <- cptswglobreg.Changingsig(
    data = y,
    Xc = Xc,
    Xg = Xg,
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_type(result, "list")
})


test_that("cptswglobreg.Changingsig returns expected components", {
  
  result <- cptswglobreg.Changingsig(
    data = y,
    Xc = Xc,
    Xg = Xg,
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_true(all(c(
    "changes",
    "loglik",
    "nchanges",
    "globpars",
    "grad",
    "muhatc"
  ) %in% names(result)))
})


test_that("cptswglobreg.Changingsig detects a changepoint", {
  
  result <- cptswglobreg.Changingsig(
    data = y,
    Xc = Xc,
    Xg = Xg,
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_true(length(result$changes) >= 1)
  
  expect_true(
    any(abs(result$changes - 72) <= 5)
  )
})

# GDcpt - ScaledSSE

test_that("GDcpt runs with ScaledSSE", {
  
  result <- GDcpt(
    data = y,
    Xc = Xc,
    Xg = Xg,
    method = "ScaledSSE",
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_type(result, "list")
})


test_that("GDcpt ScaledSSE returns expected components", {
  
  result <- GDcpt(
    data = y,
    Xc = Xc,
    Xg = Xg,
    method = "ScaledSSE",
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_true(all(c(
    "changes",
    "penloss",
    "nchanges",
    "globpars",
    "grad",
    "muhatc",
    "finallike"
  ) %in% names(result)))
})

# GDcpt - Changingsig

test_that("GDcpt runs with Changingsig", {
  
  result <- GDcpt(
    data = y,
    Xc = Xc,
    Xg = Xg,
    method = "Changingsig",
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_type(result, "list")
})


test_that("GDcpt Changingsig returns expected components", {
  
  result <- GDcpt(
    data = y,
    Xc = Xc,
    Xg = Xg,
    method = "Changingsig",
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_true(all(c(
    "changes",
    "loglik",
    "nchanges",
    "globpars",
    "grad",
    "muhatc"
  ) %in% names(result)))
})

# GDcpt agrees with cptswglobreg.SSE

test_that("GDcpt ScaledSSE agrees with cptswglobreg.SSE", {
  
  direct <- cptswglobreg.SSE(
    data = y,
    Xc = Xc,
    Xg = Xg,
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  wrapper <- GDcpt(
    data = y,
    Xc = Xc,
    Xg = Xg,
    method = "ScaledSSE",
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_equal(
    wrapper$changes,
    direct$changes
  )
})

# GDcpt agrees with cptswglobreg.Changingsig

test_that("GDcpt Changingsig agrees with cptswglobreg.Changingsig", {
  
  direct <- cptswglobreg.Changingsig(
    data = y,
    Xc = Xc,
    Xg = Xg,
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  wrapper <- GDcpt(
    data = y,
    Xc = Xc,
    Xg = Xg,
    method = "Changingsig",
    order = 1,
    maxit = 10,
    GD = FALSE
  )
  
  expect_equal(
    wrapper$changes,
    direct$changes
  )
})