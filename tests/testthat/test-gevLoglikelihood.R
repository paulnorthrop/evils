# Check that gevLoglik() agrees with gevLoglikDirect()

set.seed(28122023)
x <- rGEV(10)

# shape = 0

pars <- c(0.5, 2, 0)
res1 <- gevLoglik(x, loc = pars[1], scale = pars[2], shape = pars[3])
res2 <- gevLoglikDirect(pars, maxima = x, individual = TRUE)
test_that("gevLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
res1 <- gevLoglik(x, loc = pars[1], scale = pars[2], shape = pars[3],
                  sum = TRUE)
res2 <- gevLoglikDirect(pars, maxima = x, individual = FALSE)
test_that("gevLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# shape = 0.5

pars <- c(-0.5, 1, 0.5)
res1 <- gevLoglik(x, loc = pars[1], scale = pars[2], shape = pars[3])
res2 <- gevLoglikDirect(pars, maxima = x, individual = TRUE)
test_that("gevLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
res1 <- gevLoglik(x, loc = pars[1], scale = pars[2], shape = pars[3],
                  sum = TRUE)
res2 <- gevLoglikDirect(pars, maxima = x, individual = FALSE)
test_that("gevLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})

# shape = 1e-6

pars <- c(1, 3, 1e-6)
res1 <- gevLoglik(x, loc = pars[1], scale = pars[2], shape = pars[3])
res2 <- gevLoglikDirect(pars, maxima = x, individual = TRUE)
test_that("gevLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
res1 <- gevLoglik(x, loc = pars[1], scale = pars[2], shape = pars[3],
                  sum = TRUE)
res2 <- gevLoglikDirect(pars, maxima = x, individual = FALSE)
test_that("gevLoglik vs direct, shape = 0 , sum = FALSE", {
  testthat::expect_equal(res1, res2, ignore_attr = TRUE)
})
