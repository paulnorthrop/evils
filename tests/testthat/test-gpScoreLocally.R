# Check that the calculation of the GEV score agrees with nieve::dGEV()

# We ignore attributes in the dGP() cases below because nieve::dGP() returns
# the score vector as an attribute called "gradient"

# In nieve::dGP() a Taylor series approximation is used if |shape| < 2e-4.
# (log-)density: order 2 expansion
# score: order 1 expansion

# For the comparisons of the score, we increase the tolerance a little below
# from the testthat_tolerance() of 1.490116e-08 because for shape = 1e-6 and
# shape = -1e-6 evils and nieve differ by a little more than this, owing to the
# differences in the approximations for |shape| close to 0.

theTolerance <- 1e-7

set.seed(28122023)
x <- rGP(10)

xi <- c(0, 2e-4, -2e-4, 0.5, -0.5, 1e-6, -1e-6)
sigma <- 1:7

testFunction <- function(i, x) {
  shape <- xi[i]
  scale <- sigma[i]

  res1score <- gpScore(x, scale = scale, shape = shape)
  res1dGP <- dGP(x, scale = scale, shape = shape, log = TRUE)
  res2dGP <- nieve::dGPD2(x, scale = scale, shape = shape, log = TRUE,
                          deriv = TRUE)
  res2score <- attr(res2dGP, "gradient")
  desc1 <- paste0("GP density vs nieve::dGP(), shape = ", shape)
  test_that(desc1, {
    testthat::expect_equal(res1dGP, res2dGP, ignore_attr = TRUE)
  })
  desc2 <- paste0("GP score vs nieve::dGP(), shape = ", shape)
  test_that(desc2, {
    testthat::expect_equal(res1score, res2score, ignore_attr = TRUE,
                           tolerance = theTolerance)
  })
  return(invisible())
}

lapply(1:length(xi), testFunction, x = x)
