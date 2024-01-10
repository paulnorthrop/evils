# Check that the calculation of the GEV observed information agrees with
# numerical derivatives

# Check that package numDeriv is available
if (!requireNamespace("tools", quietly = TRUE)) {
  numDerivAvailable <- FALSE
} else {
  numDerivAvailable <- TRUE
}

set.seed(28122023)
x <- rGEV(2)

xi <- c(0, 0.01, -0.01, 0.5, -0.5)
sigma <- c(0.5, 1, 1.5, 2, 2.5)
mu <- -2:2

myTol <- 1e-5 # > testthat_tolerance()

testFunction <- function(i, x) {
  shape <- xi[i]
  scale <- sigma[i]
  loc <- mu[i]

  res1 <- gevObsInfo(x, loc = loc, scale = scale, shape = shape, sum = TRUE)
  fn <- function(par, data, sum) {
    loc <- par[1]
    scale <- par[2]
    shape <- par[3]
    val <- gevLoglik(x = data, loc = loc, scale = scale, shape = shape,
                     sum = sum)
    return(val)
  }
  if (numDerivAvailable) {
    res2 <- numDeriv::hessian(func = fn, x = c(loc, scale, shape), data = x,
                              sum = TRUE)
  } else {
    res2 <- -res1
  }

  # Note: we need to negate res2 to obtain the observed information
  test_that(paste0("gevObsInfo() vs stats::numHess(), shape = ", shape), {
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
res1 <- gevObsInfo(x, scale = sigma, shape = xi)
res2 <- array(NaN, dim = c(length(x), 3, 3))
res2[3, , ] <- gevObsInfo(x[3], scale = sigma[3], shape = xi[3])
# 1 + xi (x - mu) / sigma = 0 in this case
res2[4, , ] <- 0
# 1 + xi (x - mu) / sigma < 0 in this case
res2[5, , ] <- 0

test_that("gevObsInfo(): NaN when scale <= 0, 0 when x out of bounds", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
