#' The Generalised Extreme Value distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the Generalised Extreme Value (GEV) distribution. Care is taken to avoid the
#' numerical problems that can occur when the shape parameter is very close to zero.
#'
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param n A numeric scalar. Number of observations. If
#'   `length(n) > 1`, the length is taken to be the number required.
#' @param loc A numeric vector. Values of the GEV location parameter \eqn{\mu}.
#' @param scale A numeric vector of **positive** values. Values of the
#'   GEV scale parameter \eqn{\sigma}.
#' @param shape A numeric vector. Values of the GEV shape parameter \eqn{\xi}.
#' @param log,log.p A logical scalar. If `TRUE` then probabilities, or
#'   probability density values, are given as \eqn{\log(p)}.
#' @param lower.tail A logical scalar. If `TRUE`, probabilities are
#'   \eqn{P(X \leq x)}, otherwise, \eqn{P(X > x)}.
#' @param ... Further arguments to be passed to [`log1pdx`], which evaluates
#'   terms of the form \eqn{\log(1+x)/x} in the GEV log-density and
#'   log-distribution functions.
#' @param eps A numeric scalar. For values of \eqn{\xi} in `shape` that lie in
#'   `(-eps, eps)` an approximation to the GEV quantile function is used
#'   instead of a direct calculation. See **Details**.
#' @details The distribution function of a GEV distribution with parameters
#'  `loc` = \eqn{\mu}, `scale` = \eqn{\sigma (> 0)} and
#'  `shape` = \eqn{\xi} is
#'   \deqn{F(x) = P(X \leq x) = \exp\left\{ -\left[ 1+\xi\left(\frac{x-\mu}{\sigma}\right)
#'   \right]_+^{-1/\xi} \right\},}
#'  where \eqn{x_+ = \max(x, 0)}. If \eqn{\xi = 0} the distribution function is
#'  defined as the limit as \eqn{\xi} tends to zero. The support of the
#'  distribution depends on \eqn{\xi}: it is
#'  \eqn{x \leq \mu - \sigma / \xi} for \eqn{\xi < 0};
#'  \eqn{x \geq \mu - \sigma / \xi} for \eqn{\xi > 0};
#'  and \eqn{x} is unbounded for \eqn{\xi = 0}.
#'  Note that if \eqn{\xi < -1} the GEV density function becomes infinite
#'  as \eqn{x} approaches \eqn{\mu -\sigma / \xi} from below.
#'
#'  If `lower.tail = TRUE` then if `p = 0` (`p = 1`) then
#'  the lower (upper) limit of the distribution is returned, which is
#'  `-Inf` or `Inf` in some cases.  Similarly, but reversed,
#'  if `lower.tail = FALSE`.
#'
#'  In `dGEV` and `pGEV`, calculations are performed on a log-scale, which
#'  involves the evaluation of \eqn{\log(1+z)/z}, where
#'  \eqn{z = \xi (x - \mu) / \sigma}. Direct naive calculation using
#'  `log(1+z)/z` is unstable for `z` close to `0`. Use of [`log1p`]`(z)/z`
#'  is much better, but cannot handle the cases where `x` is equal to `0` or is
#'  extremely close to `0`. `dGEV` and `pGEV` avoid these issues using a series
#'  approximation, implemented by [`log1pdx`].
#'
#'  The GEV quantile function is \eqn{\mu - \sigma BC(-\log p, -\xi)}, where
#'  \eqn{BC(x, \lambda) = (x^\lambda - 1)/\lambda} and \eqn{p} is the
#'  probability of the required quantile. If \eqn{\lambda \in} `(-eps,eps)`
#'  then \eqn{BC(x, \lambda)} is approximated by
#'  \eqn{\log x (1+\lambda \log x / 2 + (\lambda \log x)^2 / 6)}.
#' @return `dGEV` gives the density function, `pGEV` gives the
#'   distribution function, `qGEV` gives the quantile function,
#'   and `rGEV` generates random deviates.
#'
#'   The length of the result is determined by `n` for `rgev`,
#'   and is the maximum of the lengths of the numerical arguments for the
#'   other functions.
#'
#'   The numerical arguments other than `n` are recycled to the length
#'   of the result.
#'
#'   `NA` is returned for any component of the inputs for which `x`, `q`, `p`,
#'   `loc`, `scale` or `shape` is `NA`.
#'
#'   `NaN` is returned for any component of `scale` that is non-positive and
#'   for any component for which `loc`, `scale` and/or `shape` are either `Inf`
#'   or `-Inf`, with no warning given.
#' @examples
#' # example code
#'
#' @name gevDistribution
NULL
## NULL

#' @rdname gevDistribution
#' @export
dGEV <- function(x, loc = 0, scale = 1, shape = 0, log = FALSE, ...) {
  if (length(x) == 0) {
    return(numeric(0))
  }
  # Recycle the vector input x, loc, scale and shape, if necessary
  maxLen <- max(length(x), length(loc), length(scale), length(shape))
  x <- rep_len(x, maxLen)
  loc <- rep_len(loc, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  # Return NA if x, loc, scale or shape is NA
  if (any(xIsNA <- !stats::complete.cases(x, loc, scale, shape))) {
    x[xIsNA] <- NA
  }
  # Density is undefined if scale <= 0 or any of the parameters are infinite
  isNaN <- scale <= 0 | is.infinite(loc) | is.infinite(scale) |
    is.infinite(shape)
  if (any(isNaN <- isNaN & !xIsNA)) {
    x[isNaN] <- NaN
  }
  # The density is 0 if 1 + shape * (x - loc) / scale <= 0 or if x is
  # +/- infinity
  zw <- shape * (x - loc) / scale
  outOfBounds <- 1 + zw <= 0 | is.infinite(x)
  if (any(zerod <- outOfBounds & !isNaN & !xIsNA)) {
    x[zerod] <- -Inf
  }
  # Otherwise, the density is positive
  if (any(posd <- !zerod & !isNaN & !xIsNA)) {
    m1 <- (shape[posd] + 1) * (x[posd] - loc[posd]) / scale[posd]
    m2 <- (x[posd] - loc[posd]) / scale[posd]
    logterm <- log1pdx(zw[posd], ...)
    x[posd] <- -log(scale[posd]) - m1 * logterm - exp(-m2 * logterm)
  }
  if (!log) {
    x <- exp(x)
  }
  return(x)
}

#' @rdname gevDistribution
#' @export
pGEV <- function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                 log.p = FALSE, ...) {
  if (length(q) == 0) {
    return(numeric(0))
  }
  # Recycle the vector input q, loc, scale and shape, if necessary
  maxLen <- max(length(q), length(loc), length(scale), length(shape))
  q <- rep_len(q, maxLen)
  loc <- rep_len(loc, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  # Return NA if q, loc, scale or shape is NA
  if (any(qIsNA <- !stats::complete.cases(q, loc, scale, shape))) {
    q[qIsNA] <- NA
  }
  # The cdf is undefined if scale <= 0 or any of the parameters are infinite
  isNaN <- scale <= 0 | is.infinite(loc) | is.infinite(scale) |
    is.infinite(shape)
  if (any(isNaN <- isNaN & !qIsNA)) {
    q[isNaN] <- NaN
  }
  # Return 0 if q is -infinity and 1 if q is +infinity
  if (any(qIsInf <- is.infinite(q) & !isNaN & !qIsNA)) {
    q[qIsInf] <- log((1 + sign(q[qIsInf])) / 2)
  }
  # The cdf is 0 (shape > 0) or 1 (shape < 0) if 1+shape*(q-loc)/scale <= 0
  zw <- shape * (q - loc) / scale
  if (any(cdf01 <- 1 + zw <= 0 & !isNaN & !qIsInf & !qIsNA)) {
    q[cdf01] <- log(shape[cdf01] < 0)
  }
  # Otherwise, the cdf is in (0, 1)
  if (any(cdfp <- !cdf01 & !isNaN & !qIsInf & !qIsNA)) {
    m2 <- (q[cdfp] - loc[cdfp]) / scale[cdfp]
    q[cdfp] <- -exp(-m2 * log1pdx(zw[cdfp], ...))
  }
  if (lower.tail) {
    if (!log.p) {
      q <- exp(q)
    }
  } else {
    q <- -expm1(q)
    if (log.p) {
      q <- log(q)
    }
  }
  return(q)
}

#' @rdname gevDistribution
#' @export
qGEV <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                 log.p = FALSE, eps = 1e-6) {
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
  # Recycle the vector input q, loc, scale and shape, if necessary
  maxLen <- max(length(p), length(loc), length(scale), length(shape))
  p <- rep_len(p, maxLen)
  loc <- rep_len(loc, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  # Return NA if p, loc, scale or shape is NA
  if (any(pIsNA <- !stats::complete.cases(p, loc, scale, shape))) {
    p[pIsNA] <- NA
  }
  # Quantile is undefined if scale <= 0 or any of the parameters are infinite
  isNaN <- scale <= 0 | is.infinite(loc) | is.infinite(scale) |
    is.infinite(shape)
  if (any(isNaN <- isNaN & !pIsNA)) {
    p[isNaN] <- NaN
  }
  # The quantile is undefined if scale <= 0
  # Quantiles are loc + scale [(-log(p))^(-shape) - 1] / shape
  #      which is loc - scale BoxCox(-log(p), -shape)
  cond <- !isNaN & !pIsNA
  p[cond] <- loc[cond] - scale[cond] *
    BC(x = -log(p[cond]), lambda = -shape[cond], eps = eps)
  return(p)
}

#' @rdname gevDistribution
#' @export
rGEV <- function (n, loc = 0, scale = 1, shape = 0, eps = 1e-6) {
  maxLen <- ifelse(length(n) > 1, length(n), n)
  loc <- rep_len(loc, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  return(qGEV(stats::runif(n), loc = loc, scale = scale, shape = shape,
              eps = eps))
}
