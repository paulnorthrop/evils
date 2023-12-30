# Check that log1pdx2() gives the correct values

# Values of x to check:
xvals <- c(-1.1, -1, "minL1" = -0.79149064, "-eps" = 0.01,
           "zero" = 0, "eps" = 0.01, 1, 1.1)
# Check all values in one call
res1 <- suppressWarnings(log1pdx3(xvals))
res2 <- suppressWarnings(
  2 * log1p(xvals) / xvals ^ 3 - 2 / (xvals ^ 2 * (1 + xvals)) -
    1 / (xvals * (1 + xvals) ^ 2)
  )
# x = 0. In the limit as x becomes 0, log3pdx() = 2/3
res2[names(xvals) == "zero"] <- 2 / 3
test_that("log1pdx3() and direct calculation agree for various x", {
  testthat::expect_equal(res1, res2)
})

# Check x = -0.01 and x = 0.01 individually

# x = -(the default value of) eps
# res1 uses the explicit 5-term approximation and res2 uses DPQ::logcf()
# res3 uses direct calculation
# The returned values should be very close, but not identical
x <- -0.01
res1 <- log1pdx3(x)
res2 <- log1pdx3(x - .Machine$double.neg.eps)
res3 <- 2 * log1p(x) / x ^ 3 - 2 / (x ^ 2 * (1 + x)) - 1 / (x * (1 + x) ^ 2)
test_that(paste("log1pdx3() and log1pdx3(x) using logcf() agree for x = ", x), {
  testthat::expect_equal(res1, res2)
})
test_that(paste("log1pdx3() and direct calculation agree for x = ", x), {
  testthat::expect_equal(res1, res3)
})

# x = (the default value of) eps
# res1 uses the explicit 5-term approximation and res2 uses DPQ::logcf()
# res3 uses direct calculation
# The returned values should be very close, but not identical
x <- 0.01
res1 <- log1pdx3(x)
res2 <- log1pdx3(x + .Machine$double.neg.eps)
res3 <- 2 * log1p(x) / x ^ 3 - 2 / (x ^ 2 * (1 + x)) - 1 / (x * (1 + x) ^ 2)
test_that(paste("log1pdx3() and log1pdx3(x) using logcf() agree for x = ", x), {
  testthat::expect_equal(res1, res2)
})
test_that(paste("log1pdx3() and direct calculation agree for x = ", x), {
  testthat::expect_equal(res1, res3)
})

