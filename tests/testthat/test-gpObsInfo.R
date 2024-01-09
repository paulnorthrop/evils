# Check that the calculation of the GEV observed information agrees with
# numerical derivatives

set.seed(28122023)
x <- rGP(2)

xi <- c(0, 0.01, -0.01, 0.5, -0.5)
sigma <- c(0.5, 1, 1.5, 2, 2.5)

myTol <- 1e-5 # > testthat_tolerance()

testFunction <- function(i, x) {
  shape <- xi[i]
  scale <- sigma[i]

  res1 <- gpObsInfo(x, scale = scale, shape = shape, sum = TRUE)
  fn <- function(par, data, sum) {
    scale <- par[1]
    shape <- par[2]
    val <- gpLoglik(x = data, scale = scale, shape = shape, sum = sum)
    return(val)
  }
  res2 <- numDeriv::hessian(func = fn, x = c(scale, shape), data = x,
                            sum = TRUE)

  # Note: we need to negate res2 to obtain the observed information
  test_that(paste0("gpObsInfo() vs stats::numHess(), shape = ", shape), {
    testthat::expect_equal(res1, -res2, ignore_attr = "dimnames",
                           tolerance = myTol)
  })
  return(invisible())
}

lapply(1:length(xi), testFunction, x = x)

# Check gevScore() returns NaN when scale <= 0 and 0 when x is out of bounds

x <- 0:4
sigma <- -1:3
xi <- c(0.2, 0, -1e-6, -2/3, -1)
res1 <- gpObsInfo(x, scale = sigma, shape = xi)
res2 <- array(NaN, dim = c(length(x), 2, 2))
res2[3, , ] <- gpObsInfo(x[3], scale = sigma[3], shape = xi[3])
# 1 + xi * x / sigma = 0 in this case
res2[4, , ] <- 0
# 1 + xi * x / sigma < 0 in this case
res2[5, , ] <- 0

test_that("gpObsInfo(): NaN when scale <= 0, 0 when x out of bounds", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
