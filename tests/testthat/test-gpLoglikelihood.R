# Check that gpLoglik() agrees with gpLoglikDirect()

set.seed(28122023)
x <- rGP(10)

# shape = 0

pars <- c(2, 0)
res1 <- gpLoglik(x, scale = pars[1], shape = pars[2])
res2 <- gpLoglikDirect(pars, excesses = x, individual = TRUE)
test_that("gpLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
res1 <- gpLoglik(x, scale = pars[1], shape = pars[2],
                  sum = TRUE)
res2 <- gpLoglikDirect(pars, excesses = x, individual = FALSE)
test_that("gpLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# shape = 0.5

pars <- c(1, 0.5)
res1 <- gpLoglik(x, scale = pars[1], shape = pars[2])
res2 <- gpLoglikDirect(pars, excesses = x, individual = TRUE)
test_that("gpLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
res1 <- gpLoglik(x, scale = pars[1], shape = pars[2],
                  sum = TRUE)
res2 <- gpLoglikDirect(pars, excesses = x, individual = FALSE)
test_that("gpLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# shape = 1e-6

pars <- c(3, 1e-6)
res1 <- gpLoglik(x, scale = pars[1], shape = pars[2])
res2 <- gpLoglikDirect(pars, excesses = x, individual = TRUE)
test_that("gpLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
res1 <- gpLoglik(x, scale = pars[1], shape = pars[2],
                  sum = TRUE)
res2 <- gpLoglikDirect(pars, excesses = x, individual = FALSE)
test_that("gpLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
