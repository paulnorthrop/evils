# Basic checks of gpExpInfo in simple cases

test_that("gpExpInfo() is correct for scale = 1, shape = 0", {
  testthat::expect_equal(gpExpInfo(scale = 1, shape = 0),
                         matrix(c(1, 1, 1, 2), 2, 2), ignore_attr = TRUE)
})

test_that("gpExpInfo() is correct for scale = 1/4, shape = 0", {
  testthat::expect_equal(gpExpInfo(scale = 1 / 4, shape = 0),
                         matrix(c(16, 4, 4, 2), 2, 2), ignore_attr = TRUE)
})

# Check that gevExpInfo() throws an error when shape <= -0.5

test_that("gpExpInfo() errors for shape = -1/2", {
  testthat::expect_error(gpExpInfo(shape = -1/2))
})
test_that("gExpInfo() errors for shape = -1", {
  testthat::expect_error(gpExpInfo(shape = -1))
})
