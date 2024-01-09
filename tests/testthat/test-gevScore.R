# Check that the calculation of the GEV score at shape = 0 agrees with
# the theoretical expression

set.seed(28122023)
x <- rGEV(10)
mu <- -1:8
sigma <- c(1, 2.5)
xi <- 0

res1 <- gevScore(x, loc = mu, scale = sigma, shape = xi)

w <- x - mu
ws <- w / sigma
scoreLoc <- (1 - exp(-ws)) / sigma
scoreScale <- (w * scoreLoc - 1) / sigma
scoreShape <- ws ^ 2 * (1 - exp(-ws)) / 2 - ws
res2 <- cbind("loc" = scoreLoc, "scale" = scoreScale, "shape" = scoreShape)

test_that("gevScore(shape = 0) agrees with the theoretical expression", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# Likewise, for sum = TRUE

res1 <- gevScore(x, loc = mu, scale = sigma, shape = xi, sum = TRUE)
res2 <- colSums(res2)
test_that("gevScore(shape = 0, sum = TRUE) = theoretical expression", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# Check gevScore() returns NaN when scale <= 0 and 0 when x is out of bounds

x <- 0:4
sigma <- -1:3
xi <- c(0.2, 0, -1e-6, -2/3, -1)
res1 <- gevScore(x, scale = sigma, shape = xi)
res2 <- matrix(NaN, 5, 3)
res2[3, ] <- gevScore(x[3], scale = sigma[3], shape = xi[3])
# 1 + xi (x - mu) / sigma = 0 in this case
res2[4, ] <- 0
# 1 + xi (x - mu) / sigma < 0 in this case
res2[5, ] <- 0

test_that("gevScore(): NaN when scale <= 0, 0 when x out of bounds", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
