# Some basic checks of the Box-Cox transformation function BC()

test_that("BC(x), x < 0 errors", {
  expect_error(BC(-1, lambda = 0))
})

test_that("BC(x = NA) returns NA", {
  expect_equal(BC(c(1, NA), lambda = 0), c(log(1), NA))
})

test_that("BC(lambda = NA) returns NA", {
  expect_equal(BC(x = 2, lambda = c(0, NA)), c(log(2), NA))
})

# Tests for abs(lambda) < 1e-6, the default value of eps

lambda <- -1e-7

test_that("BC(x = Inf, lambda < 0) returns -1/lambda", {
  expect_equal(BC(x = Inf, lambda = -1e-7), -1 / lambda)
})

test_that("BC(x = 0, lambda < 0) returns -Inf", {
  expect_equal(BC(x = 0, lambda = -1e-7), -Inf)
})

test_that("BC(x = (Inf, 0), lambda < 0) returns (1/2, -Inf)", {
  expect_equal(BC(x = c(Inf, 0), lambda = -1e-7), c(-1 / lambda, -Inf))
})
