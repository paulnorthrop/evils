#' Generalised Pareto log-likelihood derivatives
#'
#' Calculates the first and second derviatives of the log-density of a
#' generalised Pareto (GP) distribution, including cases where the shape
#' parameter is very close to zero.
#'
#' @param x A numeric vector of observations. Typically, these are
#'   block maxima, that is, the largest observation in a block of contiguous
#'   observations.
#' @param scale A numeric vector of **positive** values. Value(s) of the
#'   GEV scale parameter \eqn{\sigma}.
#' @param shape A numeric vector. Value(s) of the GEV shape parameter
#'  \eqn{\xi}.
#' @param ... Further arguments to be passed to [`log1pdx`], [`dlog1pdx`] or
#'   [`d2log1pdx`]. See **Details**.
#' @details
#'   The input vectors `x`, `scale` and `shape` are not inside these
#'   functions to force them to have the same lengths. Therefore, each of these
#'   input vectors should either have the same length as the longest of these
#'   vectors, or have length 1. These vectors should also satisfy the parameter
#'   constraints \eqn{\sigma > 0} and \eqn{1 + \xi x / \sigma > 0}.
#'   These requirements are explicitly not checked inside these functions,
#'   so the user must take care when calling them.
#'
#'   **First derivatives**. The functions `gp1` and `gp2` are used in
#'   [`gpScore`] to calculate the first derivatives of the GP log-likelihood
#'   \eqn{\ell} with respect to \eqn{\sigma} and \eqn{\xi}. `gevi`
#'   provides \eqn{\partial \ell / \partial \theta_i}, where
#'   \eqn{\theta = (\sigma, \xi)}.
#'
#'   The first derivative of the log-likelihood with respect to \eqn{\sigma}
#'   involves \eqn{\log(1+z)/z}, calculated using [`log1pdx`].
#'   The first derivative of the log-likelihood with respect to \eqn{\xi}
#'   involves \eqn{\log(1+z)/z} and also
#'   \eqn{\log(1+z)/z^2 - z^{-1}(1+z)^{-1}}. The latter is calculated using
#'   [`dlog1pdx`].
#'
#'   **Second derivatives**. The functions `gp11`, `gp12` and `gp22`are used in
#'   [`gpObsInfo`] to calculate the second derivatives of the GP log-likelihood
#'   \eqn{\ell} with respect to \eqn{\sigma} and \eqn{\xi}. `gevij` provides
#'   \eqn{\partial^2 \ell / \partial \theta_i \theta_j}, where
#'   \eqn{\theta = (\sigma, \xi)}.
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
#' x <- rGP(5)
#'
#' # First derivatives
#' gp1(x)
#' gp2(x)
#'
#' # Second derivatives
#' gp11(x)
#' gp12(x)
#' gp22(x)
#' @name gpDerivatives
NULL
## NULL

#' @rdname gpDerivatives
#' @export
gp1 <- function(x, scale = 1, shape = 0) {
  zw <- shape * x / scale
  val <- x * (shape + 1) / ((1 + zw) * scale) - 1
  return(val / scale)
}

#' @rdname gpDerivatives
#' @export
gp2 <- function(x, scale = 1, shape = 0, ...) {
  zw <- shape * x / scale
  val <- -x ^ 2 * dlog1pdx(zw, ...) / scale - x / (1 + zw)
  return(val / scale)
}

#' @rdname gpDerivatives
#' @export
gp11 <- function(x, scale = 1, shape = 0) {
  zw <- shape * x / scale
  val <- 1 + 2 * x * (shape + 1) / ((1 + zw) * scale) -
    x ^ 2 * shape * (shape  + 1) / ((1 + zw) ^ 2 * scale ^ 2)
  return(val / scale ^ 2)
}

#' @rdname gpDerivatives
#' @export
gp12 <- function(x, scale = 1, shape = 0) {
  zw <- shape * x / scale
  val <- x * (1 / (1 + zw) - x * (shape + 1) / ((1 + zw) ^ 2 * scale))
  return(val / scale ^ 2)
}

#' @rdname gpDerivatives
#' @export
gp22 <- function(x, loc = 0, scale = 1, shape = 0, ...) {
  zw <- shape * x / scale
  val <- x ^ 2 / ((1 + zw) ^ 2) - x ^ 3 * d2log1pdx(zw, ...) / scale
  return(val / scale ^ 2)
}
