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
