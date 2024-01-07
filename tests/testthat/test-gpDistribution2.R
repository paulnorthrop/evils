## Example shape parameters

# Positive
xi1 <- 0.1

# Zero
xi2 <- 0

# Negative
xi3 <- -1e-7
up <- -1 / xi3

## Example input vectors
# For testing pdf, log_pdf and cdf
xvec <- c(0, Inf, NA)
x1 <- xvec
x2 <- xvec
x3 <- c(up, xvec)
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("rGP, length of output", {
  expect_length(rGP(n = 1, shape = xi1), 1)
  expect_length(rGP(n = 100, shape = xi1), 100)
  expect_length(rGP(n = 0, shape = xi1), 0)
  expect_error(rGP(n = -2, shape = xi1))

  expect_length(rGP(n = 1, shape = xi2), 1)
  expect_length(rGP(n = 100, shape = xi2), 100)
  expect_length(rGP(n = 0, shape = xi2), 0)
  expect_error(rGP(n = -2, shape = xi2))

  expect_length(rGP(n = 1, shape = xi3), 1)
  expect_length(rGP(n = 100, shape = xi3), 100)
  expect_length(rGP(n = 0, shape = xi3), 0)
  expect_error(rGP(n = -2, shape = xi3))
})

test_that("dGP, values and length of output", {
  p <- pvec[2:4]
  qq <- qGP(p, shape = xi1)
  expect_equal(dGP(x1, shape = xi1), c(1, 0, NA))
  expect_equal(dGP(qq, shape = xi1), (1 - p) ^ (1 + xi1))
  expect_length(dGP(seq_len(0), shape = xi1), 0)
  expect_length(dGP(seq_len(1), shape = xi1), 1)
  expect_length(dGP(seq_len(10), shape = xi1), 10)

  qq <- qGP(p, shape = xi2)
  expect_equal(dGP(x2, shape = xi2), c(1, 0, NA))
  expect_equal(dGP(qq, shape = xi2), (1 - p) ^ (1 + xi2))
  expect_length(dGP(seq_len(0), shape = xi2), 0)
  expect_length(dGP(seq_len(1), shape = xi2), 1)
  expect_length(dGP(seq_len(10), shape = xi2), 10)

  qq <- qGP(p, shape = xi3)
  expect_equal(dGP(x3, shape = xi3), c(0, 1, 0, NA))
  expect_equal(dGP(qq, shape = xi3), (1 - p) ^ (1 + xi3))
  expect_length(dGP(seq_len(0), shape = xi3), 0)
  expect_length(dGP(seq_len(1), shape = xi3), 1)
  expect_length(dGP(seq_len(10), shape = xi3), 10)
})

test_that("pGP, values and length of output", {
  expect_equal(pGP(x1, shape = xi1), c(0, 1, NA))
  expect_length(pGP(seq_len(0), shape = xi1), 0)
  expect_length(pGP(seq_len(1), shape = xi1), 1)
  expect_length(pGP(seq_len(10), shape = xi1), 10)

  expect_equal(pGP(x2, shape = xi2), c(0, 1, NA))
  expect_length(pGP(seq_len(0), shape = xi2), 0)
  expect_length(pGP(seq_len(1), shape = xi2), 1)
  expect_length(pGP(seq_len(10), shape = xi2), 10)

  expect_equal(pGP(x3, shape = xi3), c(1, 0, 1, NA))
  expect_length(pGP(seq_len(0), shape = xi3), 0)
  expect_length(pGP(seq_len(1), shape = xi3), 1)
  expect_length(pGP(seq_len(10), shape = xi3), 10)
})

test_that("qGP, values and length of output", {
  q1 <- ((1 - pvec[2:4]) ^ (-xi1) - 1) / xi1
  expect_equal(qGP(pvec, shape = xi1), c(0, q1, Inf, NA))
  expect_length(qGP(seq_len(0), shape = xi1), 0)
  expect_length(qGP(c(0, 1), shape = xi1), 2)
  expect_length(qGP(seq_len(10) / 10, shape = xi1), 10)

  q2 <- -log(1 - pvec[2:4])
  expect_equal(qGP(pvec, shape = xi2), c(0, q2, Inf, NA))
  expect_length(qGP(seq_len(0), shape = xi2), 0)
  expect_length(qGP(c(0, 1), shape = xi2), 2)
  expect_length(qGP(seq_len(10) / 10, shape = xi2), 10)

  q3 <- ((1 - pvec[2:4]) ^ (-xi3) - 1) / xi3
  expect_equal(qGP(pvec, shape = xi3), c(0, q3, up, NA))
  expect_length(qGP(seq_len(0), shape = xi3), 0)
  expect_length(qGP(c(0, 1), shape = xi3), 2)
  expect_length(qGP(seq_len(10) / 10, shape = xi3), 10)
})

test_that("pGP and qGP are consistent", {
  expect_equal(pGP(qGP(pvec, shape = xi1), shape = xi1), pvec)
  expect_equal(pGP(qGP(pvec, shape = xi2), shape = xi2), pvec)
  expect_equal(pGP(qGP(pvec, shape = xi3), shape = xi3), pvec)
})

# Add the test that caused a problem in distributions3

## Example distributions

# Positive shape, finite lower end point
xi1 <- 0.1
g1 <- list(scale = 1, xi = xi1)
class(g1) <- c("GP", "distribution")

# Zero shape
xi2 <- 0
g2 <- list(scale = 1, xi = xi2)
class(g2) <- c("GP", "distribution")

# Negative shape, finite upper end point
xi3 <- -1e-7
g3 <- list(scale = 1, xi = xi3)
class(g3) <- c("GP", "distribution")
up <- -1 / xi3

log_GP_pdf <- function(d, x) {
  dGP(x = x, scale = d$scale, shape = d$xi, log = TRUE)
}

## Example input vectors

# For testing pdf, log_pdf and cdf
xvec <- c(0, Inf, NA)
x1 <- xvec
x2 <- xvec
x3 <- c(up, xvec)

test_that("log_pdf.GP works correctly", {
  expect_equal(log_GP_pdf(g1, x1), c(0, -Inf, NA))
  expect_length(log_GP_pdf(g1, seq_len(0)), 0)
  expect_length(log_GP_pdf(g1, seq_len(1)), 1)
  expect_length(log_GP_pdf(g1, seq_len(10)), 10)

  expect_equal(log_GP_pdf(g2, x2), c(0, -Inf, NA))
  expect_length(log_GP_pdf(g2, seq_len(0)), 0)
  expect_length(log_GP_pdf(g2, seq_len(1)), 1)
  expect_length(log_GP_pdf(g2, seq_len(10)), 10)

  expect_equal(log_GP_pdf(g3, x3), c(-Inf, 0, -Inf, NA))
  expect_length(log_GP_pdf(g3, seq_len(0)), 0)
  expect_length(log_GP_pdf(g3, seq_len(1)), 1)
  expect_length(log_GP_pdf(g3, seq_len(10)), 10)
})

