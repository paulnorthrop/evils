# Check that log1pxOverxDirect(x) and log1pxOverxApprox(x) agree.
# For log1pxOverxApprox(x) the series expansion used by sumR::infiniteSum()
# converges for x in (-1, 1).
xvec <- c(-0.99, seq(from = -0.9, to = 0.9, by = 0.1), 0.99)
for (x in xvec) {
  res1 <- log1pxOverxDirect(x)
  res2 <- log1pxOverxApprox(x)
  test_that(paste("log1pxOverDirect() = log1pxApprox(), x = ", x), {
    testthat::expect_equal(res1, res2)
  })
}

# Check that these functions return errors when they should
xvec <- c(-1.1, -1, 1, 1.1)
for (x in xvec) {
  test_that(paste("log1pxOverApprox() errors, x = ", x), {
    testthat::expect_error(log1pxOverxApprox(x))
  })
}
test_that(paste("log1pxOverDirect() errors, x = ", x), {
  testthat::expect_error(log1pxOverxDirect(-1.1))
})

# Check log1pxOver(-1) is Inf
test_that("log1pxOverDirect(-1) is Inf", {
  testthat::expect_equal(log1pxOverxDirect(-1), Inf)
})

# Does log1pxOverx(), which calls log1pxOverDirect() or log1pxOverApprox()
# depending on whether abs(x) < tol (Approx) or not (Direct) give the correct
# result?

res1 <- log1pxOverx(1e-2, tol = 1e-3)
res2 <- log1pxOverx(1e-2, tol = 1e-1)
test_that("log1pxOverx correct, x = 1e-2", {
  testthat::expect_equal(res1, res2)
})
