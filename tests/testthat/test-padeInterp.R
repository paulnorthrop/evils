# Check that padeInterp() does interpolate

## f(x) = log(1 + x) / x

# Orders of the Padé numerator and denominator.
L <- M <- 3
# Taylor series coefficients of expansion about x = 0, where f(0) = 1
j <- seq_len(L + M + 1) - 1
A <- (-1) ^ j / (j + 1)

# Pade approximation
p1 <- padeInterp(L, M, A)

# Pade approximation adjusted to interpolate over [-0.1, 0.1]
f <- function(x) log1p(x) / x
xint <- c(-0.1, 0.1)
p2 <- padeInterp(L, M, A, f, xint)

# Evaluate the approximations at xint and 0
xvals <- c(xint[1], 0, xint[2])
res1 <- predict(p2, xvals)

test_that("padeInterp() works for log(1+x)/x", {
  testthat::expect_equal(res1, c(f(xint[1]), 1, f(xint[2])))
})

# Check that padeInterp() gives the same main coefficients as Pade::Pade()

p2short <- lapply(p2, function(x) x[1:(L+1)])

test_that("padeInterp() agrees with Pade::Pade()", {
  testthat::expect_equal(p2short, p1, ignore_attr = "class")
})

# Check that padeInterp() errors when it should

test_that("padeInterp() errors when L is not an integer", {
  expect_error(padeInterp(L + 0.1, M, A))
})
test_that("padeInterp() errors when A is not longer enough", {
  expect_error(padeInterp(L, M, A[-1]))
})
