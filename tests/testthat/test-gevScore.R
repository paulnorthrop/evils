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

test_that("gevScore agrees with the theoretical expression", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
