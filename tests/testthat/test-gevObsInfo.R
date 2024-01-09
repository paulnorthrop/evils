# Check that the calculation of the GEV observed information agrees with
# numerical derivatives

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
  res2 <- numDeriv::hessian(func = fn, x = c(loc, scale, shape), data = x,
                            sum = TRUE)

  # Note: we need to negate res2 to obtain the observed information
  test_that(paste0("gevObsInfo() vs stats::numHess(), shape = ", shape), {
    testthat::expect_equal(res1, -res2, ignore_attr = "dimnames",
                           tolerance = myTol)
  })
  return(invisible())
}

lapply(1:length(xi), testFunction, x = x)
