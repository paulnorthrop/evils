# Check that log1pdx() gives the correct values

# Values of x to check:
# -1.1, 1, minL1 = -0.79149064, -eps = 0.01, 0, eps = 0.01, 1, 1.1

# x < -1: expect warning and NAN result
x <- -1.1
test_that(paste("log1pdx() warning for x = ", x), {
  testthat::expect_warning(log1pdx(x))
})
res <- suppressWarnings(log1pdx(x))
test_that(paste("log1pdx() returns NaN for x = ", x), {
  testthat::expect_identical(res, NaN)
})

# x = -1: expect Inf
x <- -1
res <- log1pdx(x)
test_that(paste("log1pdx() returns NaN for x = ", x), {
  testthat::expect_identical(res, Inf)
})

# x = (the default value of) minL1
# res1 uses DPQ::logcf() and res2 uses log1px(x) / x
# The returned values should be very close, but not identical
x <- -0.79149064
res1 <- log1pdx(x)
res2 <- log1p(x) / x
test_that(paste("log1pdx() and log1p(x) / x agree for x = ", x), {
  testthat::expect_equal(res1, res2)
})

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

# x = 0. In the limit as x becomes 0, log1pdx() = 1
x <- 0
res <- log1pdx(x)
test_that(paste("log1pdx() and log1pdx(x) using logcf() agree for x = ", x), {
  testthat::expect_identical(res, 1)
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

# x = 1
# res1 uses DPQ::logcf() and res2 uses log1px(x) / x
# The returned values should be very close, but not identical
x <- 1
res1 <- log1pdx(x)
res2 <- log1p(x) / x
test_that(paste("log1pdx() and log1p(x) / x agree for x = ", x), {
  testthat::expect_equal(res1, res2)
})

# x > 1
# res1 and res2 use log1px(x) / x
# The returned values should be identical
x <- 1.1
res1 <- log1pdx(x)
res2 <- log1p(x) / x
test_that(paste("log1pdx() and log1p(x) / x agree for x = ", x), {
  testthat::expect_identical(res1, res2)
})
