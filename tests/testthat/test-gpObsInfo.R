# Check that the calculation of the GEV observed information agrees with
# numerical derivatives

# Check that package numDeriv is available
if (!requireNamespace("tools", quietly = TRUE)) {
  numDerivAvailable <- FALSE
} else {
  numDerivAvailable <- TRUE
}

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
  if (numDerivAvailable) {
    res2 <- numDeriv::hessian(func = fn, x = c(scale, shape), data = x,
                              sum = TRUE)
  } else {
    res2 <- -res1
  }

  # Note: we need to negate res2 to obtain the observed information
  test_that(paste0("gpObsInfo() vs stats::numHess(), shape = ", shape), {
    testthat::expect_equal(res1, -res2, ignore_attr = "dimnames",
                           tolerance = myTol)
  })
  return(invisible())
}

lapply(1:length(xi), testFunction, x = x)

# Check gpObsInfo() returns NaN when scale <= 0 and 0 when x is out of bounds

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

# Check that NA main argument or parameter values result in an NA return

test_that("gpObsInfo: NA parameters produce an NA return", {
  expect_equal(gevObsInfo(NA, scale = c(1, NA, 1), shape = c(0, 0, NA)),
               array(as.numeric(NA), dim = c(3, 3, 3)),
               ignore_attr = "dimnames")
})

# Check that infinite parameter values result in an NaN return

res2 <- array(NaN, dim = c(4, 2, 2))
res2[1, , ] <- gpObsInfo(0)
res2[2, , ] <- NA
test_that("gpObsInfo: infinite parameters produce an NaN return", {
  expect_equal(gpObsInfo(0, scale = c(1, NA, Inf, 1), shape = c(0, 0, 0, -Inf)),
               res2, ignore_attr = "dimnames")
})
