# Check that the calculation of the GEV score at shape = 0 agrees with
# the theoretical expression

set.seed(28122023)
x <- rGEV(10)
mu <- 0
sigma <- 1
xi <- 0
w <- x - mu
ws <- w / sigma
res1 <- gevScore(x, shape = xi)

scoreLoc <- (1 - exp(-ws)) / sigma
scoreScale <- (w * scoreLoc - 1) / sigma
scoreShape <- ws ^ 2 * (1 - exp(-ws)) / 2 - ws
res2 <- cbind("loc" = scoreLoc, "scale" = scoreScale, "shape" = scoreShape)

test_that("gevScore() agrees with the theoretical expression", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# Check that gevScore() returns 0 when scale <= 0
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

test_that("gevScore() returns 0 when scale <= 0", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
