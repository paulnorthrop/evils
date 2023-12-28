#' Generalised Extreme Value functions
#'
#' Calculate the log-likelihood, score and observed information for
#' a random sample from a generalised extreme value (GEV) distribution,
#' including cases where the shape parameter is very close to zero.
#'
#' @param x A numeric vector of observations. Typically, these are
#'   block maxima, that is, the largest observation in a block of contiguous
#'   observations.
#' @param loc A numeric vector. Values of the GEV location parameter \eqn{\mu}.
#' @param scale A numeric vector of **positive** values. Values of the
#'   GEV scale parameter \eqn{\sigma}.
#' @param shape A numeric vector. Values of the GEV shape parameter \eqn{\xi}.
#' @param sum A logical scalar. Relevant to `gevLoglik` and
#'   `gevScore`. If `sum = TRUE` then only the sum of
#'   contributions from all observations in `x` is calculated.  If
#'   `sum = FALSE` then individual contributions from each
#'   observation in `x` are calculated.
#' @param ... Further arguments to be passed to [`log1pdx`], which evaluates
#'   terms of the form \eqn{\log(1+x)/x} in the GEV log-likelihood.
#' @details
#'   **Log-likelihood** (`gevLoglik`). The two problematic
#'   terms of the log-likelihood both involve \eqn{\log(1+z)/z},
#'   where \eqn{z = \xi (x - \mu)/ \sigma}. In one part this is exponentiated,
#'   in the other it is not. The function [`log1pdx`] is used to calculate
#'   \eqn{\log(1+z)/z}, with special case taken for cases where `z` is close to
#'   zero.
#'
#'   **Score** (`gevScore`). The first derivatives of the log-likelihood with
#'   respect to \eqn{\mu} and \eqn{\sigma} both involve \eqn{\log(1+z)/z},
#'   calculated using [`log1pdx`]. The first derivative of the log-likelihood
#'   with respect to \eqn{\xi} involves \eqn{\log(1+z)/z} and also
#'   \eqn{\log(1+z)/z^2 - z^{-1}(1+z)^{-1}}. The latter is calculated using
#'   [`log1pdx2`].
#' @return
#'   **Log-likelihood** (`gevLoglik`). If
#'   `sum = TRUE` a scalar, the value of the log-likelihood. If
#'   `sum = FALSE` a vector of length `length{x}`
#'   containing the contributions to the log-likelihood from each of the
#'   observations.
#'
#' **Score** (`gevScore`).  If `sum = TRUE` the value
#'  of the score, a vector of length 3 containing the derivative of the
#'  log-likelihood evaluated at the input parameter values.
#'  If `sum = FALSE` the values of the contributions to the score
#'  from each of the observations, a
#'   `length(x)`\eqn{ \times 3} matrix.
#'   The columns are named `sigma[u]` and `xi`.
#'
#' **Observed information** (`gevInfo`).  The observed information: a
#'   \eqn{2 \times 2} matrix with row and column names
#'   `c(sigma[u], xi)`.
#' @name gevLikelihood
NULL
## NULL

#' Generalised Extreme Value Log-likelihood
#'
#' @examples
#' ### Simulate some data
#'
#' set.seed(17042022)
#' y <- rGenExtremeValue(100, 0, 1, 0)
#'
#' ### Log-likelihood
#'
#' # Calculation using log1pdx()
#' gevLoglik(y, 0, 1, 1e-8, )
#' gevLoglik(y, 0, 1, -1e-8)
#' gevLoglik(y, 0, 1, 0)
#'
#' # Direct calculation, involving (1 / xi) * log1p(xi * (y - mu) / sigma)
#' # Mostly fine, but breaks down eventually
#' gevLoglikDirect(pars = c(0, 1, 1e-309), maxima = y)
#' gevLoglikDirect(pars = c(0, 1, -1e-309), maxima = y)
#'
#' # Score
#' set.seed(28122023)
#' x <- rGEV(10)
#' gevScore(x)
#' gevScore(x, sum = TRUE)
#' nieve::dGEV(x, log = TRUE, deriv = TRUE)
#' @rdname gevLikelihood
#' @export
gevLoglik <- function(x, loc = 0, scale = 1, shape = 0, sum = FALSE, ...) {
  loglik <- dGEV(x, loc, scale, shape, log = TRUE, ...)
  if (sum) {
    loglik <- sum(loglik)
  }
  return(loglik)
}

#' @rdname gevLikelihood
#' @export
gevScore <- function(x, loc = 0, scale = 1, shape = 0, sum = FALSE, ...) {
  # Recycle the vector input x, loc, scale and shape, if necessary
  maxLen <- max(length(x), length(loc), length(scale), length(shape))
  x <- rep_len(x, maxLen)
  loc <- rep_len(loc, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  # Create a 3-column matrix to store the results
  scoreMat <- matrix(NA, maxLen, 3)
  # The density, and, hence the score, is undefined if scale <= 0
  if (any(invalidScale <- scale <= 0)) {
    scoreMat[invalidScale, ] <- NaN
  }
  w <- x - loc
  zw <- shape * w / scale
  # The density is 0 if 1 + shape * (x - loc) / scale <= 0
  # Set the score to 0 for these cases
  zw <- shape * (x - loc) / scale
  if (any(zerod <- 1 + zw <= 0 & !invalidScale)) {
    scoreMat[zerod] <- 0
  }
  # Otherwise, the density is positive
  if (any(posd <- !zerod & !invalidScale)) {
    # Derivative of the log-likelihood with respect to loc = mu
    scoreMat[posd, 1] <- ((shape[posd] + 1) / (1 + zw[posd]) -
                         exp(-(shape[posd] + 1) * w[posd] * log1pdx(zw[posd]) /
                               scale[posd])) / scale[posd]
    # Derivative of the log-likelihood with respect to scale = sigma
    scoreMat[posd, 2] <- (w[posd] * scoreMat[posd, 1] - 1) / scale[posd]
    # Derivative of the log-likelihood with respect to shape = xi
    scoreMat[posd, 3] <- w[posd] ^ 2 * log1pdx2(zw[posd]) *
      (1 - exp(-w[posd] * log1pdx(zw[posd]) / scale[posd])) /
      scale[posd] ^ 2 - w[posd] / (scale[posd] * (1 + zw[posd]))
  }
  if (sum) {
    scoreMat <- colSums(scoreMat)
  }
  colnames(scoreMat) <- c("loc", "scale", "shape")
  return(scoreMat)
}
