# Check that log1pdx() gives the correct values

# Values of x to check:
xvals <- c(-1.1, -1, "minL1" = -0.79149064, "-eps" = 0.01,
           "zero" = 0, "eps" = 0.01, 1, 1.1, "NA" = NA, "Inf" = Inf)
# Check all values in one call
res1 <- suppressWarnings(log1pdx(xvals))
res2 <- suppressWarnings(log1p(xvals) / xvals)
# x = 0. In the limit as x becomes 0, log1pdx(x) = 1
res2[names(xvals) == "zero"] <- 1
# x = Inf. In the limit as x becomes Inf, log1pdx(x) = 0
res2[names(xvals) == "Inf"] <- 0
test_that("log1pdx() and log1p(x) / x agree for various x", {
  testthat::expect_equal(res1, res2)
})

# Check x = -0.01 and x = 0.01 individually

# x = -(the default value of) eps
# res1 uses the explicit 5-term approximation and res2 uses DPQ::logcf()
# res3 uses log1px(x) / x
# The returned values should be very close, but not identical
x <- -0.01
res1 <- log1pdx(x)
res2 <- log1pdx(x - .Machine$double.neg.eps)
res3 <- log1p(x) / x
test_that(paste("log1pdx() and log1pdx(x) using logcf() agree for x = ", x), {
  testthat::expect_equal(res1, res2)
})
test_that(paste("log1pdx() and log1p(x) / x agree for x = ", x), {
  testthat::expect_equal(res1, res3)
})

# x = (the default value of) eps
# res1 uses the explicit 5-term approximation and res2 uses DPQ::logcf()
# res3 uses log1px(x) / x
# The returned values should be very close, but not identical
x <- 0.01
res1 <- log1pdx(x)
res2 <- log1pdx(x + .Machine$double.neg.eps)
res3 <- log1p(x) / x
test_that(paste("log1pdx() and log1pdx(x) using logcf() agree for x = ", x), {
  testthat::expect_equal(res1, res2)
})
test_that(paste("log1pdx() and log1p(x) / x agree for x = ", x), {
  testthat::expect_equal(res1, res3)
})
