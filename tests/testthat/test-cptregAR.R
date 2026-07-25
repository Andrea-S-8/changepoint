library(testthat)
library(changepoint)

context("cpt.regAR tests")

#############################
# Array input regression AR #
#############################

test_that("cpt.regAR detects a regression changepoint", {
  
  set.seed(123)
  
  x <- c(
    rnorm(50,0,1),
    rnorm(50,5,1)
  )
  
  data_matrix <- cbind(
    x,
    rep(1,length(x))
  )
  
  result <- cpt.regAR(
    data=data_matrix,
    penalty="MBIC",
    method="AMOC",
    dist="Normal",
    class=TRUE
  )
  
  expect_s4_class(result,"cpt.reg")
  
})


#############################
# Create regression AR data #
#############################

set.seed(123)

n <- 100

# Mean shift at 50
x <- c(
  rnorm(50, 0, 1),
  rnorm(50, 5, 1)
)

# Regression format:
# response | regressors | intercept
trend <- (1:n)/n

data <- cbind(
  x,
  trend,
  rep(1,n)
)


########################
# AMOC test
########################

test_that("cpt.regAR AMOC detects changepoint", {
  
  result <- cpt.regAR(
    data = data,
    penalty = "MBIC",
    method = "AMOC",
    dist = "Normal",
    class = TRUE
  )
  
  expect_s4_class(result, "cpt.reg")
  
  # First changepoint should be at 50
  expect_true(
    result@cpts[1] %in% 49:51
  )
  
})


########################
# PELT test
########################

test_that("cpt.regAR PELT runs successfully", {
  
  result <- cpt.regAR(
    data = data,
    penalty = "MBIC",
    method = "PELT",
    dist = "Normal",
    class = TRUE
  )
  
  expect_s4_class(result, "cpt.reg")
  
  expect_true(
    length(result@cpts) >= 1
  )
  
})


########################
# Output contains cpts
########################

test_that("cpt.regAR returns changepoints", {
  
  result <- cpt.regAR(
    data = data,
    penalty = "MBIC",
    method = "AMOC",
    dist = "Normal",
    class = TRUE
  )
  
  expect_s4_class(result, "cpt.reg")
  
  expect_true(
    length(result@cpts) >= 1
  )
  
})


########################
# No changepoint regression
########################

test_that("cpt.regAR runs on no changepoint regression", {
  
  nochangedata <- cbind(
    0.2*(1:n)/n + rnorm(n),
    (1:n)/n,
    rep(1,n)
  )
  
  
  result <- cpt.regAR(
    data = nochangedata,
    penalty = "MBIC",
    method = "AMOC",
    dist = "Normal",
    class = TRUE
  )
  
  
  expect_s4_class(result,"cpt.reg")
  
})


########################
# Missing data test
########################

test_that("cpt.regAR rejects NA values", {
  
  baddata <- data
  
  baddata[10,1] <- NA
  
  
  expect_error(
    cpt.regAR(
      data = baddata,
      penalty = "MBIC",
      method = "AMOC",
      dist = "Normal"
    )
  )
  
})


########################
# Non numeric test
########################

test_that("cpt.regAR rejects non numeric data", {
  
  baddata <- data
  
  baddata[,1] <- as.character(baddata[,1])
  
  
  expect_error(
    cpt.regAR(
      data = baddata,
      penalty = "MBIC",
      method = "AMOC",
      dist = "Normal"
    )
  )
  
})


########################
# Short data test
########################

test_that("cpt.regAR rejects very short data", {
  
  shortdata <- matrix(
    c(1,1),
    ncol = 2
  )
  
  
  expect_error(
    cpt.regAR(
      data = shortdata,
      penalty = "MBIC",
      method = "AMOC",
      dist = "Normal"
    )
  )
  
})


########################
# class = FALSE
########################

test_that("cpt.regAR returns changepoints when class=FALSE", {
  
  result <- cpt.regAR(
    data = data,
    penalty = "MBIC",
    method = "AMOC",
    dist = "Normal",
    class = FALSE
  )
  
  expect_type(result, "list")
  
  expect_true("cpts" %in% names(result))
  expect_true("pen.value" %in% names(result))
  
  expect_true(is.numeric(result$cpts))
  expect_true(is.numeric(result$pen.value))  
})



########################
# Invalid method
########################

test_that("cpt.regAR rejects invalid method", {
  
  expect_error(
    
    cpt.regAR(
      data = data,
      penalty = "MBIC",
      method = "TEST",
      dist = "Normal"
    )
    
  )
  
})


########################
# Invalid penalty
########################

test_that("cpt.regAR rejects invalid penalty", {
  
  expect_error(
    
    cpt.regAR(
      data = data,
      penalty = "TEST",
      method = "AMOC",
      dist = "Normal"
    )
    
  )
  
})


########################
# Invalid distribution
########################

test_that("cpt.regAR converts invalid distribution to Normal", {
  
  result <- cpt.regAR(
    data = data,
    penalty = "MBIC",
    method = "AMOC",
    dist = "Gamma",
    class = TRUE
  )
  
  
  expect_s4_class(result,"cpt.reg")
  
})