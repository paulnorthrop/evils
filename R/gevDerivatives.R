#' Generalised Extreme Value log-likelihood derivatives
#'
#' Calculates the first and second derviatives of the log-density of a
#' generalised extreme value (GEV) distribution, including cases where the
#' shape parameter is very close to zero.
#'
#' @param x A numeric vector of observations. Typically, these are
#'   block maxima, that is, the largest observation in a block of contiguous
#'   observations.
#' @param loc A numeric vector. Value(s) of the GEV location parameter \eqn{\mu}.
#' @param scale A numeric vector of **positive** values. Value(s) of the
#'   GEV scale parameter \eqn{\sigma}.
#' @param shape A numeric vector. Value(s) of the GEV shape parameter
#' \eqn{\xi}.
#' @param ... Further arguments to be passed to [`log1pdx`], [`dlog1pdx`] or
#'   [`d2log1pdx`]. See **Details**.
#' @details
#'   The input vectors `x`, `loc`, `scale` and `shape` are not inside these
#'   functions to force them to have the same lengths. Therefore, each of these
#'   input vectors should either have the same length as the longest of these
#'   vectors, or have length 1. These vectors should also satisfy the parameter
#'   constraints \eqn{\sigma > 0} and \eqn{1 + \xi(x - \mu) / \sigma > 0}.
#'   These requirements are explicitly not checked inside these functions,
#'   so the user must take care when calling them.
#'
#'   **First derivatives**. The functions `gev1`, `gev2` and `gev3` are used in
#'   [`gevScore`] to calculate the first derivatives of the GEV log-likelihood
#'   \eqn{\ell} with respect to \eqn{\mu}, \eqn{\sigma} and \eqn{\xi}. `gevi`
#'   provides \eqn{\partial \ell / \partial \theta_i}, where
#'   \eqn{\theta = (\mu, \sigma, \xi)}.
#'
#'   The first derivatives of the log-likelihood with respect to \eqn{\mu} and
#'   \eqn{\sigma} both involve \eqn{\log(1+z)/z}, calculated using [`log1pdx`].
#'   The first derivative of the log-likelihood with respect to \eqn{\xi}
#'   involves \eqn{\log(1+z)/z} and also
#'   \eqn{\log(1+z)/z^2 - z^{-1}(1+z)^{-1}}. The latter is calculated using
#'   [`dlog1pdx`].
#'
#'   **Second derivatives**. The functions `gev11`, `gev12`, `gev13`, `gev22`,
#'   `gev23` and `gev33` are used in [`gevObsInfo`] to calculate the second
#'   derivatives of the GEV log-likelihood \eqn{\ell} with respect to
#'   \eqn{\mu}, \eqn{\sigma} and \eqn{\xi}. `gevij` provides
#'   \eqn{\partial^2 \ell / \partial \theta_i \theta_j}, where
#'   \eqn{\theta = (\mu, \sigma, \xi)}.
#'
#'   All the second derivatives of the log-likelihood involve
#'   \eqn{\log(1+z)/z}, calculated using [`log1pdx`]. Any second derivative
#'   made with respect to \eqn{\xi} also involves
#'   \eqn{\log(1+z)/z^2 - z^{-1}(1+z)^{-1}}, which is calculated using
#'   [`dlog1pdx`]. Differentiating the log-likelihood twice with respect to
#'   \eqn{\xi} produces the function
#'   \eqn{2 \log(1+z)/z^3 - 2 z^{-2}(1+z)^{-1} - z^{-1}(1+z)^{-2}}, which is
#'   calculated using [`d2log1pdx`].
#' @return
#'   A numeric vector of length `length(x)` containing the relevant first or
#'   second derivative of the GEV log-density.
#' @examples
#' ### Simulate some data
#' set.seed(28122023)
#' x <- rGEV(5)
#'
#' # First derivatives
#' gev1(x)
#' gev2(x)
#' gev3(x)
#'
#' # Second derivatives
#' gev11(x)
#' gev12(x)
#' gev13(x)
#' gev22(x)
#' gev23(x)
#' gev33(x)
#' @name gevDerivatives
NULL
## NULL

#' @rdname gevDerivatives
#' @export
gev1 <- function(x, loc = 0, scale = 1, shape = 0, ...) {
  w <- x - loc
  zw <- shape * w / scale
  val <- (shape + 1) / (1 + zw) - exp(-(shape + 1) * w * log1pdx(zw) / scale)
  return(val / scale)
}

#' @rdname gevDerivatives
#' @export
gev2 <- function(x, loc = 0, scale = 1, shape = 0, ...) {
  w <- x - loc
  val <- w * gev1(x, loc, scale, shape) - 1
  return(val / scale)
}

#' @rdname gevDerivatives
#' @export
gev3 <- function(x, loc = 0, scale = 1, shape = 0, ...) {
  w <- x - loc
  zw <- shape * w / scale
  val <- w ^ 2 * dlog1pdx(zw) * (exp(-w * log1pdx(zw) / scale) - 1) /
    scale ^ 2 - w / (scale * (1 + zw))
  return(val)
}

#' @rdname gevDerivatives
#' @export
gev11 <- function(x, loc = 0, scale = 1, shape = 0, ...) {
  w <- x - loc
  zw <- shape * w / scale
  val <- (shape + 1) * (shape / (1 + zw) ^ 2 -
                          exp(-(2 * shape + 1) * w * log1pdx(zw) / scale))
  return(val / scale ^ 2)
}

#' @rdname gevDerivatives
#' @export
gev12 <- function(x, loc = 0, scale = 1, shape = 0, ...) {
  w <- x - loc
  val <- w * gev11(x, loc, scale, shape) - gev1(x, loc, scale, shape)
  return(val / scale)
}

#' @rdname gevDerivatives
#' @export
gev22 <- function(x, loc = 0, scale = 1, shape = 0, ...) {
  w <- x - loc
  val <- (1 - w * gev1(x, loc, scale, shape)) / scale +
    w * gev12(x, loc, scale, shape)
  return(val / scale)
}

#' @rdname gevDerivatives
#' @export
gev13 <- function(x, loc = 0, scale = 1, shape = 0, ...) {
  w <- x - loc
  zw <- shape * w / scale
  hzw <- log1pdx(zw)
  hdashzw <- dlog1pdx(zw)
  zwr <- 1 / (1 + zw)
  term1 <- scale * zwr - w * (1 + shape) * zwr ^ 2
  term2 <- w * hzw + (1 + shape) * w ^ 2 * hdashzw / scale
  term3 <- exp(-w * (1 + shape) * hzw / scale)
  val <- term1 + term2 * term3
  return(val / scale ^ 2)
}

#' @rdname gevDerivatives
#' @export
gev23 <- function(x, loc = 0, scale = 1, shape = 0, ...) {
  w <- x - loc
  val <- w * gev13(x, loc, scale, shape)
  return(val / scale)
}

#' @rdname gevDerivatives
#' @export
gev33 <- function(x, loc = 0, scale = 1, shape = 0, ...) {
  w <- x - loc
  zw <- shape * w / scale
  hzw <- log1pdx(zw)
  hdashzw <- dlog1pdx(zw)
  hdash2zw <- d2log1pdx(zw)
  zwr <- 1 / (1 + zw)
  term1 <- hdash2zw * (exp(-w * hzw / scale) - 1)
  term2 <- w * hdashzw ^ 2 * exp(-w * hzw / scale) / scale
  val <- w ^ 3 * (term1 - term2) / scale + w ^ 2 * zwr ^ 2
  return(val / scale ^ 2)
}
