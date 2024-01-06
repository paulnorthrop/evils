#' Accurate `log(1+x)/x` computation (and similar)
#'
#' Compute the function \eqn{h(x) = \log(1+x)/x} accurately also for small \eqn{x}, that is,
#' \eqn{|x| \ll 1}{|x| << 1}.  Similarly, compute the first and second
#' derivatives \eqn{h'(x) = x^{-1}(1+x)^{-1} - \log(1+x)/x^2} and
#' \eqn{h''(x) = 2 \log(1+x)/x^3 - 2 x^{-2}(1+x)^{-1} - x^{-1}(1+x)^{-2}}.
#'
#' @param x A numeric vector with values \eqn{x > -1}.
#' @param minL1 A negative cutoff \eqn{m_l}. For \eqn{x \in (m_l, 1)} the
#'   computations are based on the series expansion 4.1.29 on page 68 of
#'   Abramowitz and Stegun (1972).
#' @param eps2 A non-negative cutoff \eqn{\epsilon_2}. For
#'   \eqn{|x| > \epsilon_2}, [`logcf`][`DPQ::logcf`] is used in the computation
#'   of \eqn{\log(1+x)/x} using the series expansion. Otherwise, that is,
#'   \eqn{|x| \leq \epsilon_2}, only a few terms of the expansion are used.
#' @param tol_logcf A non-negative number indicating the tolerance (maximal
#'   relative error) for the [`logcf`][`DPQ::logcf`] function.
#' @param trace.lcf A logical scalar. Used in
#'   [`logcf`][`DPQ::logcf`]`(.., trace = trace.lcf)`.
#' @param logCF The function to be used as [`logcf`][`DPQ::logcf`]. The default
#'   chooses the pure \R `logcfR()` when `x` is not numeric, and chooses the
#'   C-based `logcf()` when `is.numeric(x)` is `TRUE`.
#' @details For \eqn{x} in \eqn{(m_l, 1)} the computations are based on the
#'   series expansion
#' \deqn{\log(1+x) = 2 \left( t+\frac{t^3}{3}+\frac{t^5}{5}+\cdots \right),}
#'   where \eqn{t = x/(2+x)} which is formula 4.1.29 on page 68 of
#'   Abramowitz and Stegun (1972).
#'
#' Different computations are used in 3 different ranges of \eqn{x}.
#' Let \eqn{y = t^2}.
#'
#' For `log1pdx`:
#'
#'  * \eqn{x < m_l} or \eqn{x > 1}:
#'    `log1pdx(x):=` [`log1p`] `/x`;
#'  * \eqn{|x| \leq \epsilon_2}:
#'    \eqn{\frac{2}{2+x}\left((((\frac19y+\frac17)y+\frac15)y+\frac13)y+1\right)};
#'  * \eqn{x \in (m_l, 1)} and \eqn{|x| > \epsilon_2}:
#'    \eqn{\frac{2}{2+x}} [`logcf`]\eqn{(y, 1, 2)}.
#'
#' For `dlog1pdx`:
#'
#'  * \eqn{x < m_l} or \eqn{x > 1}:
#'    `dlog1pdx(x):= 1/(x*(1 + x)) - ` [`log1p`] `/x^2`;
#'  * \eqn{|x| \leq \epsilon_2}:
#'    \eqn{-\frac{1}{(2+x)^2}} \eqn{\left(\frac{2+x}{1+x} + 2 t
#'    ((((\frac{1}{11}y+\frac19)y+\frac17)y+\frac15)y+\frac13\right)};
#'  * \eqn{x \in (m_l, 1)} and \eqn{|x| > \epsilon_2}:
#'    \eqn{-\frac{1}{(2+x)^2}} \eqn{\Big(\frac{2+x}{1+x} + 2 t}
#'    [`logcf`]\eqn{(y, 3, 2)\Big)}.
#'
#' For `d2log1pdx`:
#'
#'  * \eqn{x < m_l} or \eqn{x > 1}:
#'    `dlog1pdx(x):= 2` [`log1p`] `/x^3 - 2/(x^2*(1+x)) - 1/(x*(1 + x)^2)`;
#'  * \eqn{|x| \leq \epsilon_2}:
#'    \eqn{\frac{4}{(2+x)^3}} \eqn{\left(\frac{(2+x)^2}{4(1+x)^2} +
#'    ((((\frac{1}{11}y+\frac19)y+\frac17)y+\frac15)y+\frac13\right)};
#'  * \eqn{x \in (m_l, 1)} and \eqn{|x| > \epsilon_2}:
#'    \eqn{\frac{4}{(2+x)^3}} \eqn{\Big(\frac{(2+x)^2}{4(1+x)^2} + }
#'    [`logcf`]\eqn{(y, 3, 2)\Big)}.
#'
#' @return A numeric vector (with the same attributes as `x`).
#' @author Paul Northrop created these functions by modifying the code in
#'   [`log1pmx`][`DPQ::log1pmx`].
#' @references Abramowitz, M. and Stegun, I. A. (1972) Handbook of Mathematical
#' Functions. New York: Dover.
#' [https://en.wikipedia.org/wiki/Abramowitz_and_Stegun](https://en.wikipedia.org/wiki/Abramowitz_and_Stegun)
#' provides links to the full text, which is in public domain.
#' @examples
#' # In the limit as x tends to 0, log(1+x)/x = 1
#' log1pdx(0)
#'
#' # log1p(x) / x is fine unless x is 0 or extremely close to 0, e.g. 1e-324,
#' # when it returns NaN
#' x <- seq(from = -1, to = 2, by = 0.01)
#' y1 <- log1pdx(x)
#' y2 <- log1p(x) / x
#' y <- cbind(y1, y2)
#' matplot(x, y, type = "l")
#'
#' # In the limit as x tends to 0, 1/(x*(1+x)) - log(1+x)/x^2 tends to -1/2
#' dlog1pdx(0)
#'
#' # 1 / (x * (1 + x)) - log1p(x) / x ^ 2 is fine unless x is 0 (NaN returned)
#' # or extremely close to 0, e.g. 1e-15 returns 0.625, 1e-16 returns 0
#' # when it returns NaN
#' x <- seq(from = -1, to = 2, by = 0.01)
#' y1 <- dlog1pdx(x)
#' y2 <- 1 / (x * (1 + x)) - log1p(x) / x ^ 2
#' y <- cbind(y1, y2)
#' matplot(x, y, type = "l")
#'
#' # In the limit as x tends to 0, 2log(1+x)/x^3-2/(x^2*(1+x))-1/(x(1+x)^2)
#' # tends to 2/3
#' d2log1pdx(0)
#'
#' # 2log(1+x)/x^3-2/(x^2*(1+x))1/(x(1+x)^2) becomes unstable for
#' # abs(x) < 1e-6
#' x <- seq(from = -1, to = 2, by = 0.01)
#' y1 <- d2log1pdx(x)
#' y2 <- 2 * log1p(x) / x ^ 3 - 2 / (x ^ 2 * (1 + x)) - 1 / (x * (1 + x) ^ 2)
#' y <- cbind(y1, y2)
#' matplot(x, y, type = "l")
#' @name log1pdx
NULL
## NULL

#' @rdname log1pdx
#' @export
log1pdx <- function(x, minL1 = -0.79149064, eps2 = 0.01, tol_logcf = 1e-14,
                    trace.lcf = FALSE,
                    logCF = if (is.numeric(x)) DPQ::logcf else DPQ::logcfR.) {

  # Check input values of eps2 and minL1
  stopifnot(is.numeric(eps2), eps2 >= 0, is.numeric(minL1),
            -1 <= minL1, minL1 < 0)

  # Create return vector of the correct length
  r <- x
  # Deal with cases NA and Inf
  if (any(xIsNA <- is.na(x))) {
    r[xIsNA] <- NA
  }
  if (any(xIsInf <- is.infinite(x) & x > 0)) {
    r[xIsInf] <- 0
  }
  # Otherwise, there are 3 cases:
  #
  # 1. x > 1 or x < minL1
  # 2. -eps2 <= x <= eps2
  # 3. x in [minL1, -eps2) or x in (eps2, 1]
  #
  # 1.
  if (any(c1 <- (x > 1 | x < minL1) & !xIsNA & !xIsInf)) {
    r[c1] <- log1p(x[c1]) / x[c1]
  }
  # Not 1.
  if (any(c2 <- !c1 & !xIsNA & !xIsInf)) {
    x <- x[c2]
    term <- x / (2 + x)
    term2 <- 2 / (2 + x)
    y <- term * term
    r[c2] <- term2 * {
      A <- x
      # 2.
      if (any(isSml <- abs(x) <= eps2)) {
        i <- which(isSml)
        y. <- y[i]
        # Becomes y-like (e.g. mpfr)
        z <- if (is.numeric(x)) 1 else 1 + 0 * y.
        A[i] <- ((((z / 9 * y. + z / 7) * y. + z / 5) * y. + z / 3) * y. + 1)
      }
      # 3.
      if (length(iLrg <- which(!isSml))) {
        y. <- y[iLrg]
        A[iLrg] <- logCF(y., 1, 2, tol_logcf, trace = trace.lcf)
      }
      A
    }
  }
  return(r)
}

#' @rdname log1pdx
#' @export
dlog1pdx <- function(x, minL1 = -0.79149064, eps2 = 0.01, tol_logcf = 1e-14,
                     trace.lcf = FALSE,
                     logCF = if (is.numeric(x)) DPQ::logcf else DPQ::logcfR.) {

  # Check input values of eps2 and minL1
  stopifnot(is.numeric(eps2), eps2 >= 0, is.numeric(minL1),
            -1 <= minL1, minL1 < 0)

  # Create return vector of the correct length
  r <- x
  # Deal with cases NA and Inf
  if (any(xIsNA <- is.na(x))) {
    r[xIsNA] <- NA
  }
  if (any(xIsInf <- is.infinite(x) & x > 0)) {
    r[xIsInf] <- 0
  }
  # Otherwise, there are 3 cases:
  #
  # 1. x > 1 or x < minL1
  # 2. -eps2 <= x <= eps2
  # 3. x in [minL1, -eps2) or x in (eps2, 1]
  #
  # 1.
  if (any(c1 <- (x > 1 | x < minL1) & !xIsNA & !xIsInf)) {
    r[c1] <- log1p(x[c1]) / x[c1] ^ 2 - 1 / (x[c1] * (1 + x[c1]))
  }
  # Not 1.
  if (any(c2 <- !c1 & !xIsNA & !xIsInf)) {
    x <- x[c2]
    term <- x / (2 + x)
    term2 <- 1 / (2 + x) ^ 2
    y <- term * term
    r[c2] <- term2 * {
      A <- x
      # 2.
      if (any(isSml <- abs(x) <= eps2)) {
        i <- which(isSml)
        y. <- y[i]
        # Becomes y-like (e.g. mpfr)
        z <- if (is.numeric(x)) 1 else 1 + 0 * y.
        A[i] <- ((((z / 11 * y. + z / 9) * y. + z / 7) * y. + z / 5) * y. +
                   1 / 3)
      }
      # 3.
      if (length(iLrg <- which(!isSml))) {
        y. <- y[iLrg]
        A[iLrg] <- logCF(y., 3, 2, tol_logcf, trace = trace.lcf)
      }
      2 * term * A + (2 + x) / (1 + x)
    }
  }
  return(-r)
}

#' @rdname log1pdx
#' @export
d2log1pdx <- function(x, minL1 = -0.79149064, eps2 = 0.01, tol_logcf = 1e-14,
                     trace.lcf = FALSE,
                     logCF = if (is.numeric(x)) DPQ::logcf else DPQ::logcfR.) {

  # Check input values of eps2 and minL1
  stopifnot(is.numeric(eps2), eps2 >= 0, is.numeric(minL1),
            -1 <= minL1, minL1 < 0)

  # Create return vector of the correct length
  r <- x
  # Deal with cases NA and Inf
  if (any(xIsNA <- is.na(x))) {
    r[xIsNA] <- NA
  }
  if (any(xIsInf <- is.infinite(x) & x > 0)) {
    r[xIsInf] <- 0
  }
  # Otherwise, there are 3 cases:
  #
  # 1. x > 1 or x < minL1
  # 2. -eps2 <= x <= eps2
  # 3. x in [minL1, -eps2) or x in (eps2, 1]
  #
  # 1.
  if (any(c1 <- (x > 1 | x < minL1) & !xIsNA & !xIsInf)) {
    r[c1] <- -2 * log1p(x[c1]) / x[c1] ^ 3  + 2 / (x[c1] ^ 2 * (1 + x[c1])) +
      1 / (x[c1] * (1 + x[c1]) ^ 2)
  }
  # Not 1.
  if (any(c2 <- !c1 & !xIsNA & !xIsInf)) {
    x <- x[c2]
    term <- x / (2 + x)
    term2 <- -4 / (2 + x) ^ 3
    y <- term * term
    r[c2] <- term2 * {
      A <- x
      # 2.
      if (any(isSml <- abs(x) <= eps2)) {
        i <- which(isSml)
        y. <- y[i]
        # Becomes y-like (e.g. mpfr)
        z <- if (is.numeric(x)) 1 else 1 + 0 * y.
        A[i] <- ((((z / 11 * y. + z / 9) * y. + z / 7) * y. + z / 5) * y. +
                   1 / 3)
      }
      # 3.
      if (length(iLrg <- which(!isSml))) {
        y. <- y[iLrg]
        A[iLrg] <- logCF(y., 3, 2, tol_logcf, trace = trace.lcf)
      }
      A + (2 + x) ^ 2 / (4 * (1 + x) ^ 2)
    }
  }
  return(-r)
}
