# Check that d2log1pdx() gives the correct values

# Values of x to check:
xvals <- c(-1.1, -1, "minL1" = -0.79149064, "-eps" = 0.01,
           "zero" = 0, "eps" = 0.01, 1, 1.1, "NA" = NA, "Inf" = Inf)
# Check all values in one call
res1 <- suppressWarnings(d2log1pdx(xvals))
res2 <- suppressWarnings(
  2 * log1p(xvals) / xvals ^ 3 - 2 / (xvals ^ 2 * (1 + xvals)) -
    1 / (xvals * (1 + xvals) ^ 2)
  )
# x = 0. In the limit as x becomes 0, log3pdx() = 2/3
res2[names(xvals) == "zero"] <- 2 / 3
# x = Inf. In the limit as x becomes Inf, dlog1pdx(x) = 0
res2[names(xvals) == "Inf"] <- 0
test_that("d2log1pdx() and direct calculation agree for various x", {
  testthat::expect_equal(res1, res2)
})

# Check x = -0.01 and x = 0.01 individually

# x = -(the default value of) eps
# res1 uses the explicit 5-term approximation and res2 uses DPQ::logcf()
# res3 uses direct calculation
# The returned values should be very close, but not identical
x <- -0.01
res1 <- d2log1pdx(x)
res2 <- d2log1pdx(x - .Machine$double.neg.eps)
res3 <- 2 * log1p(x) / x ^ 3 - 2 / (x ^ 2 * (1 + x)) - 1 / (x * (1 + x) ^ 2)
test_that(paste("d2log1pdx() and d2log1pdx(x) using logcf() agree for x = ", x), {
  testthat::expect_equal(res1, res2)
})
test_that(paste("d2log1pdx() and direct calculation agree for x = ", x), {
  testthat::expect_equal(res1, res3)
})

# x = (the default value of) eps
# res1 uses the explicit 5-term approximation and res2 uses DPQ::logcf()
# res3 uses direct calculation
# The returned values should be very close, but not identical
x <- 0.01
res1 <- d2log1pdx(x)
res2 <- d2log1pdx(x + .Machine$double.neg.eps)
res3 <- 2 * log1p(x) / x ^ 3 - 2 / (x ^ 2 * (1 + x)) - 1 / (x * (1 + x) ^ 2)
test_that(paste("d2log1pdx() and d2log1pdx(x) using logcf() agree for x = ", x), {
  testthat::expect_equal(res1, res2)
})
test_that(paste("d2log1pdx() and direct calculation agree for x = ", x), {
  testthat::expect_equal(res1, res3)
})

