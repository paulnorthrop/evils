# We check the functions pGP, qGP and rGP.

# Set a tolerance for the comparison of the simulated values
my_tol <- 1e-5 # testthat_tolerance()

# 1. Check that calling qGP with probabilities p and then calling pGP with
#    the results gets us back to the initial probabilities.

pqgp_test_fn <- function(x, p) {
  scale <- x[1]
  shape <- x[2]
  qs <- qGP(p = p, scale = scale, shape = shape)
  ps <- pGP(qs, scale = scale, shape = shape)
  return(list(p = p, ps = ps))
}

test_function <- function(x, test_string) {
  testthat::test_that(test_string, {
    testthat::expect_equal(x$p, x$ps, tolerance = my_tol)
  })
}

ep <- 1e-10
scale_check <- 2
shape_check <- c(-1, -0.5, -0.1, -ep, 0, ep, 0.1, 0.5, 1)
par_vals <- cbind(scale_check, shape_check)
p_vals <- c(0.01, 0.1, 0.5, 0.9, 0.99)
for (i in 1:nrow(par_vals)) {
  test_string <- paste("GP shape = ", par_vals[i, 2])
  x <- pqgp_test_fn(x = par_vals[i, ], p = p_vals)
  test_function(x, test_string)
}

# 2. Check that calling rGP and then pGP with the results gets us back to
#    the random U(0,1) variates simulated by stats::runif.

seed <- 28082017

rqgp_test_fn <- function(x) {
  scale <- x[, 1]
  shape <- x[, 2]
  n <- length(scale)
  set.seed(seed)
  qs <- rGP(n = n, scale = scale, shape = shape)
  set.seed(seed)
  us <- stats::runif(length(scale))
  ps <- pGP(qs, scale = scale, shape = shape)
  return(list(us = us, ps = ps))
}

test_function <- function(x, test_string) {
  testthat::test_that(test_string, {
    testthat::expect_equal(x$us, x$ps, tolerance = my_tol)
  })
}

ep <- 1e-10
scale_check <- 2
shape_check <- c(-1, -0.5, -0.1, -ep, 0, ep, 0.1, 0.5, 1)
par_vals <- cbind(scale_check, shape_check)
test_string <- "rGP and pGP"
x <- rqgp_test_fn(x = par_vals)
test_function(x, test_string)
