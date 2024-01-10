## Example shape parameters

# Positive
xi1 <- 0.1
low <- -1 / xi1

# Zero
xi2 <- 0

# Negative
xi3 <- -1e-7
up <- -1 / xi3

## Example input vectors
# For testing pdf, log_pdf and cdf
xvec <- c(-Inf, 0, Inf, NA)
x1 <- c(low, xvec)
x2 <- xvec
x3 <- c(up, xvec)
# For testing quantile
pvec <- c(0, 0.25, 0.5, 0.75, 1, NA)

test_that("rGEV, length of output", {
  expect_length(rGEV(n = 1, shape = xi1), 1)
  expect_length(rGEV(n = 100, shape = xi1), 100)
  expect_length(rGEV(n = 0, shape = xi1), 0)
  expect_error(rGEV(n = -2, shape = xi1))

  expect_length(rGEV(n = 1, shape = xi2), 1)
  expect_length(rGEV(n = 100, shape = xi2), 100)
  expect_length(rGEV(n = 0, shape = xi2), 0)
  expect_error(rGEV(n = -2, shape = xi2))

  expect_length(rGEV(n = 1, shape = xi3), 1)
  expect_length(rGEV(n = 100, shape = xi3), 100)
  expect_length(rGEV(n = 0, shape = xi3), 0)
  expect_error(rGEV(n = -2, shape = xi3))
})

test_that("dGEV, values and length of output", {
  p <- pvec[2:4]
  qq <- qGEV(p, shape = xi1)
  expect_equal(dGEV(x1, shape = xi1), c(0, 0, exp(-1), 0, NA))
  expect_equal(dGEV(qq, shape = xi1), (-log(p)) ^ (1 + xi1) * p)
  expect_length(dGEV(seq_len(0), shape = xi1), 0)
  expect_length(dGEV(seq_len(1), shape = xi1), 1)
  expect_length(dGEV(seq_len(10), shape = xi1), 10)

  qq <- qGEV(p, shape = xi2)
  expect_equal(dGEV(x1, shape = xi2), c(0, 0, exp(-1), 0, NA))
  expect_equal(dGEV(qq, shape = xi2), (-log(p)) ^ (1 + xi2) * p)
  expect_length(dGEV(seq_len(0), shape = xi2), 0)
  expect_length(dGEV(seq_len(1), shape = xi2), 1)
  expect_length(dGEV(seq_len(10), shape = xi2), 10)

  qq <- qGEV(p, shape = xi3)
  expect_equal(dGEV(x1, shape = xi3), c(0, 0, exp(-1), 0, NA))
  expect_equal(dGEV(qq, shape = xi3), (-log(p)) ^ (1 + xi3) * p)
  expect_length(dGEV(seq_len(0), shape = xi3), 0)
  expect_length(dGEV(seq_len(1), shape = xi3), 1)
  expect_length(dGEV(seq_len(10), shape = xi3), 10)
})

test_that("pGEV, values and length of output", {
  expect_equal(pGEV(x1, shape = xi1), c(0, 0, exp(-1), 1, NA))
  expect_length(pGEV(seq_len(0), shape = xi1), 0)
  expect_length(pGEV(seq_len(1), shape = xi1), 1)
  expect_length(pGEV(seq_len(10), shape = xi1), 10)

  expect_equal(pGEV(x2, shape = xi2), c(0, exp(-1), 1, NA))
  expect_length(pGEV(seq_len(0), shape = xi2), 0)
  expect_length(pGEV(seq_len(1), shape = xi2), 1)
  expect_length(pGEV(seq_len(10), shape = xi2), 10)

  expect_equal(pGEV(x3, shape = xi3), c(1, 0, exp(-1), 1, NA))
  expect_length(pGEV(seq_len(0), shape = xi3), 0)
  expect_length(pGEV(seq_len(1), shape = xi3), 1)
  expect_length(pGEV(seq_len(10), shape = xi3), 10)
})

test_that("qGEV, values and length of output", {
  q1 <- ((-log(pvec[2:4])) ^ (-xi1) - 1) / xi1
  expect_equal(qGEV(pvec, shape = xi1), c(low, q1, Inf, NA))
  expect_length(qGEV(seq_len(0), shape = xi1), 0)
  expect_length(qGEV(c(0, 1), shape = xi1), 2)
  expect_length(qGEV(seq_len(10) / 10, shape = xi1), 10)

  q2 <- -log(-log(pvec[2:4]))
  expect_equal(qGEV(pvec, shape = xi2), c(-Inf, q2, Inf, NA))
  expect_length(qGEV(seq_len(0), shape = xi2), 0)
  expect_length(qGEV(c(0, 1), shape = xi2), 2)
  expect_length(qGEV(seq_len(10) / 10, shape = xi2), 10)

  q3 <- ((-log(pvec[2:4])) ^ (-xi3) - 1) / xi3
  expect_equal(qGEV(pvec, shape = xi3), c(-Inf, q3, up, NA))
  expect_length(qGEV(seq_len(0), shape = xi3), 0)
  expect_length(qGEV(c(0, 1), shape = xi3), 2)
  expect_length(qGEV(seq_len(10) / 10, shape = xi3), 10)
})

test_that("pGEV and qGEV are consistent", {
  expect_equal(pGEV(qGEV(pvec, shape = xi1), shape = xi1), pvec)
  expect_equal(pGEV(qGEV(pvec, shape = xi2), shape = xi2), pvec)
  expect_equal(pGEV(qGEV(pvec, shape = xi3), shape = xi3), pvec)
})

# Repeat the last tests for lower.tail = FALSE and/or log.p = TRUE

test_that("pGEV, qGEV consistent, lower.tail = FALSE and/or log.p = TRUE", {
  expect_equal(pGEV(
    qGEV(pvec, shape = xi1, lower.tail = FALSE),
    shape = xi1, lower.tail = FALSE),
    pvec)
  expect_equal(pGEV(
    qGEV(log(pvec), shape = xi2, log.p = TRUE),
    shape = xi2, log.p = TRUE),
    log(pvec))
  expect_equal(pGEV(
    qGEV(log(pvec), shape = xi3, lower.tail = FALSE, log.p = TRUE),
    shape = xi3, lower.tail = FALSE, log.p = TRUE),
    log(pvec))
})

# Check that NA parameter values result in an NA return

test_that("dGEV: NA parameters produce an NA return", {
  expect_equal(dGEV(0, loc = c(NA, 0, 0), scale = c(1, NA, 1),
                  shape = c(0, 0, NA)), as.numeric(c(NA, NA, NA)))
})

test_that("pGEV: NA parameters produce an NA return", {
  expect_equal(pGEV(0, loc = c(NA, 0, 0), scale = c(1, NA, 1),
                    shape = c(0, 0, NA)), as.numeric(c(NA, NA, NA)))
})

test_that("qGEV: NA parameters produce an NA return", {
  expect_equal(qGEV(0, loc = c(NA, 0, 0), scale = c(1, NA, 1),
                    shape = c(0, 0, NA)), as.numeric(c(NA, NA, NA)))
})

# Check that infinite parameter values result in an NaN return

test_that("dGEV: infinite parameters produce an NaN return", {
  expect_equal(dGEV(0, loc = c(-2, Inf, 0, 0), scale = c(1, 1, Inf, 1),
                    shape = c(0, 0, 0, -Inf)),
               as.numeric(c(dGEV(0, loc = -2), NaN, NaN, NaN)))
})

test_that("pGEV: infinite parameters produce an NaN return", {
  expect_equal(pGEV(0, loc = c(0, -Inf, 0, 0), scale = c(1, 1, -Inf, 1),
                    shape = c(0.1, 0, 0, Inf)),
               as.numeric(c(pGEV(0, shape = 0.1), NaN, NaN, NaN)))
})

test_that("qGEV: infinite parameters, or scale <= 0, produces an NaN return", {
  expect_equal(qGEV(0.6, loc = c(0, 1, 0, 0), scale = c(1, 0, Inf, 1),
                    shape = c(-1e-6, 0, 0, Inf)),
               as.numeric(c(qGEV(0.6, shape = -1e-6), NaN, NaN, NaN)))
})

# Check that qGEV throws an error if p < 0 or p > 1

test_that("qGEV, error if p < 0 or p > 1", {
  expect_error(qGEV(p = -0.1))
  expect_error(qGEV(p = 1.1))
})

#
