#' The Generalised Pareto distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the 2-parameter Generalised Pareto (GP) distribution. Care is taken to avoid
#' the numerical problems that can occur when the shape parameter is very close
#' to zero.
#'
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param n A numeric scalar. Number of observations. If
#'   `length(n) > 1`, the length is taken to be the number required.
#' @param scale A numeric vector of **positive** values. Values of the
#'   GP scale parameter \eqn{\sigma}.
#' @param shape A numeric vector. Values of the GP shape parameter \eqn{\xi}.
#' @param log,log.p A logical scalar. If `TRUE` then probabilities, or
#'   probability density values, are given as \eqn{\log(p)}.
#' @param lower.tail A logical scalar. If `TRUE`, probabilities are
#'   \eqn{P(X \leq x)}, otherwise, \eqn{P(X > x)}.
#' @param ... Further arguments to be passed to [`log1pdx`], which evaluates
#'   terms of the form \eqn{\log(1+x)/x} in the GP log-density and
#'   log-distribution functions.
#' @param eps A numeric scalar. For values of \eqn{\xi} in `shape` that lie in
#'   `(-eps, eps)` an approximation to the GP quantile function is used
#'   instead of a direct calculation. See **Details**.
#' @details The distribution function of a GP distribution with parameters
#'  `scale` = \eqn{\sigma (> 0)} and `shape` = \eqn{\xi} is
#'  \deqn{F(x) = P(X \leq x) = 1 - \left[ 1+\xi x / \sigma \right]_+^{-1/\xi},}
#'  where \eqn{x_+ = \max(x, 0)}. If \eqn{\xi = 0} the distribution function is
#'  defined as the limit as \eqn{\xi} tends to zero.
#'  The support of the distribution depends on \eqn{\xi}: it is
#'  \eqn{x \geq 0} for \eqn{\xi \geq 0}; and \eqn{0 \leq x \leq -\sigma / \xi}
#'  for \eqn{\xi < 0}.  Note that if \eqn{\xi < -1} the GP density function
#'  becomes infinite as \eqn{x} approaches \eqn{-\sigma/\xi}.
#'
#'  If `lower.tail = TRUE` then if `p = 0` (`p = 1`) then
#'  the lower (upper) limit of the distribution is returned.
#'  The upper limit is `Inf` if `shape` is non-negative.
#'  Similarly, but reversed, if `lower.tail = FALSE`.
#'
#'  In `dGP` and `pGP`, calculations are performed on a log-scale, which
#'  involves the evaluation of \eqn{\log(1+z)/z}, where
#'  \eqn{z = \xi x / \sigma}. Direct naive calculation using
#'  `log(1+z)/z` is unstable for `z` close to `0`. Use of [`log1p(z)`]`/z`
#'  is much better, but cannot handle the cases where `x` is equal to `0` or is
#'  extremely close to `0`. `dGP` and `pGP` avoid these issues using a series
#'  approximation, implemented by [`log1pdx`].
#'
#'  The GP quantile function is \eqn{-\sigma BC(-\log p, -\xi)}, where
#'  \eqn{BC(x, \lambda) = (x^\lambda - 1)/\lambda} and \eqn{p} is the
#'  probability of the required quantile. If \eqn{\lambda \in} `(-eps,eps)`
#'  then \eqn{BC(x, \lambda)} is approximated by
#'  \eqn{\log x (1+\lambda \log x / 2 + (\lambda \log x)^2 / 6)}.
#' @return `dGP` gives the density function, `pGP` gives the
#'   distribution function, `qGP` gives the quantile function,
#'   and `rGP` generates random deviates.
#'
#'   The length of the result is determined by `n` for `rGP`,
#'   and is the maximum of the lengths of the numerical arguments for the
#'   other functions.
#'
#'   The numerical arguments other than `n` are recycled to the length
#'   of the result.
#'
#'   `NaN` is returned for any component of `scale` that is non-positive, with
#'   no warning given.
#' @examples
#' # example code
#'
#' @name gpDistribution
NULL
## NULL

#' @rdname gpDistribution
#' @export
dGP <- function(x, scale = 1, shape = 0, log = FALSE, ...) {
  if (length(x) == 0) {
    return(numeric(0))
  }
  # Recycle the vector input x, scale and shape, if necessary
  maxLen <- max(length(x), length(scale), length(shape))
  x <- rep_len(x, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  # Return NA if x, scale or shape is NA
  if (any(xIsNA <- !complete.cases(x, scale, shape))) {
    x[xIsNA] <- NA
  }
  # The density is undefined if scale <= 0
  if (any(invalidScale <- scale <= 0 & !xIsNA)) {
    x[invalidScale] <- NaN
  }
  # The density is 0 if 1 + shape * x / scale <= 0
  # The density is 0 if 1 + shape * x / scale <= 0 and/or if x is
  # +/- infinity
  zw <- shape * x / scale
  outOfBounds <- 1 + zw <= 0 | is.infinite(x)
  if (any(zerod <- outOfBounds & !invalidScale & !xIsNA)) {
    x[zerod] <- -Inf
  }
  # Otherwise, the density is positive
  if (any(posd <- !zerod & !invalidScale & !xIsNA)) {
    m1 <- (shape[posd] + 1) * x[posd] / scale[posd]
    logterm <- log1pdx(zw[posd], ...)
    x[posd] <- -log(scale[posd]) - m1 * logterm
  }
  if (!log) {
    x <- exp(x)
  }
  return(x)
}

#' @rdname gpDistribution
#' @export
pGP <- function(q, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE,
                ...) {
  if (length(q) == 0) {
    return(numeric(0))
  }
  # Recycle the vector input q, scale and shape, if necessary
  maxLen <- max(length(q), length(scale), length(shape))
  q <- rep_len(q, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  # Return NA if q, scale or shape is NA
  if (any(qIsNA <- !complete.cases(q, scale, shape))) {
    q[qIsNA] <- NA
  }
  # The cdf is undefined if scale <= 0
  if (any(invalidScale <- scale <= 0 & !qIsNA)) {
    q[invalidScale] <- NaN
  }
  # The cdf is 0 if q <= 0 and is 1 if shape < 0 and 1+shape*q/scale <= 0
  # It is also 1 if q is +Inf
  if (any(cdf0 <- q <= 0 & !invalidScale & !qIsNA)) {
    q[cdf0] <- 0
  }
  zw <- shape * q / scale
  cdf1cond <- (1 + zw <= 0 & shape < 0) | (is.infinite(q) & q > 0)
  if (any(cdf1 <- cdf1cond & !invalidScale & !cdf0)) {
    q[cdf1] <- -Inf
  }
  # Otherwise, the cdf is in (0, 1)
  if (any(cdfp <- !cdf0 & !cdf1 & !invalidScale & !qIsNA)) {
    m2 <- q[cdfp] / scale[cdfp]
    q[cdfp] <- -m2 * log1pdx(zw[cdfp], ...)
  }
  if (lower.tail) {
    if (log.p) {
      p <- DPQ::log1mexp(-q)
    } else {
      p <- 1 - exp(q)
    }
  } else {
    if (!log.p) {
      p <- exp(q)
    }
  }
  return(p)
}

#' @rdname gpDistribution
#' @export
qGP <- function(p, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE,
                eps = 1e-6) {

  if (length(p) == 0) {
    return(numeric(0))
  }
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("invalid p: p must be in [0,1].")
  }
  # Recycle the vector input q, scale and shape, if necessary
  maxLen <- max(length(p), length(scale), length(shape))
  p <- rep_len(p, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  # The quantile is undefined if scale <= 0
  if (any(invalidScale <- scale <= 0 & !is.na(scale))) {
    p[invalidScale] <- NaN
  }
  # Quantiles are scale [(-log(p))^(-shape) - 1] / shape
  #      which is scale BoxCox(-log(p), -shape)
  p[!invalidScale] <- -scale[!invalidScale] *
    BC(x = 1 - p[!invalidScale], lambda = -shape[!invalidScale], eps = eps)
  return(p)
}

#' @rdname gpDistribution
#' @export
rGP <- function (n, scale = 1, shape = 0, eps = 1e-6) {
  maxLen <- ifelse(length(n) > 1, length(n), n)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  return(qGP(stats::runif(n), scale = scale, shape = shape, eps = eps))
}
