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

test_that("gpScore() returns 0 when scale <= 0", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
