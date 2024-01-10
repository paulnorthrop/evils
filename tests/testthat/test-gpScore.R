# Check that the calculation of the GP score at shape = 0 agrees with
# the theoretical expression

set.seed(28122023)
x <- rGP(10)
sigma <- c(1, 2.5)
xi <- 0
res1 <- gpScore(x, scale = sigma, shape = xi)

scoreScale <- (x / sigma - 1) / sigma
scoreShape <- x * (0.5 * x / sigma - 1) / sigma
res2 <- cbind("scale" = scoreScale, "shape" = scoreShape)

test_that("gpScore(shape = 0) agrees with the theoretical expression", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# Likewise, for sum = TRUE

res1 <- gpScore(x, scale = sigma, shape = xi, sum = TRUE)
res2 <- colSums(res2)
test_that("gpScore(shape = 0, sum = TRUE) = theoretical expression", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# Check that gpScore() returns 0 when scale <= 0

x <- 0:4
sigma <- -1:3
xi <- c(0.2, 0, -1e-6, -2/3, -1)
res1 <- gpScore(x, scale = sigma, shape = xi)
res2 <- matrix(NaN, 5, 2)
res2[3, ] <- gpScore(x[3], scale = sigma[3], shape = xi[3])
# 1 + xi (x - mu) / sigma = 0 in this case
res2[4, ] <- 0
# 1 + xi (x - mu) / sigma < 0 in this case
res2[5, ] <- 0

test_that("gpScore(): NaN when scale <= 0, 0 when x out of bounds", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# Check that the calculation of the GP score agrees with numerical derivatives

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

testFunction <- function(i, x) {
  shape <- xi[i]
  scale <- sigma[i]

  res1 <- gpScore(x, scale = scale, shape = shape, sum = TRUE)
  fn <- function(par, data, sum) {
    scale <- par[1]
    shape <- par[2]
    val <- gpLoglik(x = data, scale = scale, shape = shape, sum = sum)
    return(val)
  }
  if (numDerivAvailable) {
    res2 <- numDeriv::grad(func = fn, x = c(scale, shape), data = x,
                           sum = TRUE)
  } else {
    res2 <- res1
  }

  # Note: we need to negate res2 to obtain the observed information
  test_that(paste0("gpScore() vs stats::numHess(), shape = ", shape), {
    testthat::expect_equal(res1, res2, ignore_attr = "names")
  })
  return(invisible())
}

lapply(1:length(xi), testFunction, x = x)

# Check that NA parameter values result in an NA return

test_that("gpScore: NA parameters produce an NA return", {
  expect_equal(gpScore(0, scale = c(1, NA, 1), shape = c(0, 0, NA)),
               rbind(gpScore(0), matrix(NA, 2, 2)), ignore_attr = "dimnames")
})

# Check that infinite parameter values result in an NaN return

test_that("gpScore: infinite parameters produce an NaN return", {
  expect_equal(gpScore(0, scale = c(1, NA, Inf, 1), shape = c(0, 0, 0, -Inf)),
               rbind(gpScore(0), rep(NA, 2), matrix(NaN, 2, 2)))
})
