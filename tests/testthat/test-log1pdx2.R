# Check that log1pdx2() gives the correct values

# Values of x to check:
xvals <- c(-1.1, -1, "minL1" = -0.79149064, "-eps" = 0.01,
           "zero" = 0, "eps" = 0.01, 1, 1.1)
# Check all values in one call
res1 <- suppressWarnings(log1pdx2(xvals))
res2 <- suppressWarnings(
  1 / (xvals * (1 + xvals)) - log1p(xvals) / xvals ^ 2
  )
# x = 0. In the limit as x becomes 0, log2pdx() = -1/2
res2[names(xvals) == "zero"] <- -1 / 2
test_that("log1pdx2() and 1/(x*(1 + x)) - log1p(x)/x^2 agree for various x", {
  testthat::expect_equal(res1, res2)
})

# Check x = -0.01 and x = 0.01 individually

# x = -(the default value of) eps
# res1 uses the explicit 5-term approximation and res2 uses DPQ::logcf()
# res3 uses 1 / (x * (1 + x)) - log1px(x) / x ^ 2
# The returned values should be very close, but not identical
x <- -0.01
res1 <- log1pdx2(x)
res2 <- log1pdx2(x - .Machine$double.neg.eps)
res3 <- 1 / (x * (1 + x)) - log1p(x) / x ^ 2
test_that(paste("log1pdx2() and log1pdx2(x) using logcf() agree for x = ", x), {
  testthat::expect_equal(res1, res2)
})
test_that(paste("log1pdx2() and 1/(x*(1 + x))-log1p(x)/x^2 agree for x = ", x), {
  testthat::expect_equal(res1, res3)
})

# x = (the default value of) eps
# res1 uses the explicit 5-term approximation and res2 uses DPQ::logcf()
# res3 uses 1 / (x * (1 + x)) - log1p(x) / x ^ 2
# The returned values should be very close, but not identical
x <- 0.01
res1 <- log1pdx2(x)
res2 <- log1pdx2(x + .Machine$double.neg.eps)
res3 <- 1 / (x * (1 + x)) - log1p(x) / x ^ 2
test_that(paste("log1pdx2() and log1pdx2(x) using logcf() agree for x = ", x), {
  testthat::expect_equal(res1, res2)
})
test_that(paste("log1pdx2() and 1/(x*(1 + x))-log1p(x)/x^2 agree for x = ", x), {
  testthat::expect_equal(res1, res3)
})
