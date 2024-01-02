#' Generalised Pareto functions
#'
#' Calculate the log-likelihood, score and observed information for
#' a random sample from a generalised Pareto (GP) distribution,
#' including cases where the shape parameter is very close to zero.
#'
#' @param x A numeric vector of observations. Typically, these are
#'   block maxima, that is, the largest observation in a block of contiguous
#'   observations.
#' @param scale A numeric vector of **positive** values. Values of the
#'   GP scale parameter \eqn{\sigma}.
#' @param shape A numeric vector. Values of the GP shape parameter \eqn{\xi}.
#' @param sum A logical scalar. Relevant to `gpLoglik` and `gpScore`. If
#'   `sum = TRUE` then only the sum of contributions from all combinations of
#'   observations in `x ` and parameter values in `loc`, `scale` and `shape`
#'   is calculated.  If `sum = FALSE` then these individual contributions are
#'   summed to produce the result.
#' @param ... Further arguments to be passed to [`log1pdx`], [`dlog1pdx`] or
#'   [`d2log1pdx`]. See [`gpDerivatives`].
#' @details
#'   **Log-likelihood** (`gpLoglik`). The two problematic
#'   terms of the log-likelihood both involve \eqn{\log(1+z)/z},
#'   where \eqn{z = \xi x / \sigma}. In one part this is exponentiated,
#'   in the other it is not. The function [`log1pdx`] is used to calculate
#'   \eqn{\log(1+z)/z}, with special case taken for cases where `z` is close to
#'   zero.
#'
#'   **Score** (`gpScore`). The first derivative of the log-likelihood
#'   with respect to \eqn{\xi} involves \eqn{\log(1+z)/z}, calculated using
#'   [`log1pdx`].
#'
#'   The functions `gp1` and `gp2` are used in
#'   `gpScore` to calculate the first derivatives of the GP log-likelihood
#'   \eqn{\ell} with respect to \eqn{\sigma} and \eqn{\xi}.
#'   See [`gpDerivatives`].
#'
#'   **Observed information** (`gpObsInfo`). The second derivative of the
#'   log-likelihood respect to \eqn{\xi} involves the function
#'   \eqn{2 \log(1+z)/z^3 - 2 z^{-2}(1+z)^{-1} - z^{-1}(1+z)^{-2}}, which is
#'   calculated using [`d2log1pdx`].
#'
#'   The functions `gp11`, `gp12` and `gp22` are
#'   used in `gpObsInfo` to calculate the second derivatives of the GP
#'   log-likelihood \eqn{\ell} with respect to \eqn{\sigma} and
#'   \eqn{\xi}. `gpij` provides
#'   \eqn{\partial^2 \ell / \partial \theta_i \theta_j}, where
#'   \eqn{\theta = (\sigma, \xi)}.
#'   See [`gpDerivatives`].
#' @return
#'   The number `n`, say, of combinations of the input parameters `x`,
#'   `scale` and `shape` in the returned object is the maximum of the lengths
#'   of these arguments, with these arguments being recycled as necessary.
#'
#'   **Log-likelihood** (`gpLoglik`). If `sum = FALSE` a vector of length `n`
#'   containing individual log-likelihood contributions. If `sum = TRUE` a
#'   scalar, the value of sum of these contributions.
#'
#'   **Score** (`gpScore`). If `sum = FALSE` the values of the contributions
#'   to the score, an `n`\eqn{ \times 2} matrix. The columns are named
#'   `loc`, `scale` and `shape`. If `sum = TRUE`, a vector of length 2
#'   containing the column sums of the matrix returned if `sum = FALSE`.
#'
#'   **Observed information** (`gpObsInfo`).  If `sum = FALSE` the values of
#'   the contributions to the observed information, an
#'   `n` \eqn{ \times 2 \times 2} array. The second and third dimensions are
#'   named columns are named `scale` and `shape`. If `sum = TRUE`, a
#'   \eqn{2 \times 2} matrix giving the observed information matrix, obtained
#'   by applying [`colsums`] to the array returned if `sum = FALSE`.
#' @name gpLikelihood
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
#' gpLoglik(y, 0, 1, 1e-8, )
#' gpLoglik(y, 0, 1, -1e-8)
#' gpLoglik(y, 0, 1, 0)
#'
#' # Direct calculation, involving (1 / xi) * log1p(xi * (y - mu) / sigma)
#' # Mostly fine, but breaks down eventually
#' gpLoglikDirect(pars = c(0, 1, 1e-309), maxima = y)
#' gpLoglikDirect(pars = c(0, 1, -1e-309), maxima = y)
#'
#' # Score
#' set.seed(28122023)
#' x <- rGP(10)
#' gpScore(x)
#' gpScore(x, sum = TRUE)
#' @rdname gpLikelihood
#' @export
gpLoglik <- function(x, scale = 1, shape = 0, sum = FALSE, ...) {
  loglik <- dGP(x, scale, shape, log = TRUE, ...)
  if (sum) {
    loglik <- sum(loglik)
  }
  return(loglik)
}

#' @rdname gpLikelihood
#' @export
gpScore <- function(x, scale = 1, shape = 0, sum = FALSE, ...) {
  # Recycle the vector input x, loc, scale and shape, if necessary
  maxLen <- max(length(x), length(scale), length(shape))
  x <- rep_len(x, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  # Create a 3-column matrix to store the results
  scoreMat <- matrix(NA, maxLen, 2)
  # The density, and, hence the score, is undefined if scale <= 0
  if (any(invalidScale <- scale <= 0)) {
    scoreMat[invalidScale, ] <- NaN
  }
  # The density is 0 if 1 + shape * (x - loc) / scale <= 0
  # Set the score to 0 for these cases
  zx <- shape * x / scale
  if (any(zerod <- 1 + zx <= 0 & !invalidScale)) {
    scoreMat[zerod] <- 0
  }
  # Otherwise, the density is positive
  if (any(pos <- !zerod & !invalidScale)) {
    # Derivative of the log-likelihood with respect to scale = sigma
    scoreMat[pos, 1] <- gp2(x[pos], scale[pos], shape[pos], ...)
    # Derivative of the log-likelihood with respect to shape = xi
    scoreMat[pos, 2] <- gp3(x[pos], scale[pos], shape[pos], ...)
  }
  colnames(scoreMat) <- c("scale", "shape")
  if (sum) {
    scoreMat <- colSums(scoreMat)
  }
  return(scoreMat)
}

#' @rdname gpLikelihood
#' @export
gpObsInfo <- function(x, scale = 1, shape = 0, sum = FALSE, ...) {
  # Recycle the vector input x, loc, scale and shape, if necessary
  maxLen <- max(length(x), length(scale), length(shape))
  x <- rep_len(x, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  # Create an array to store the results
  infoArray <- array(NA, dim = c(length(x), 2, 2))
  # The density, and, hence the information, is undefined if scale <= 0
  if (any(invalidScale <- scale <= 0)) {
    infoArray[invalidScale, ] <- NaN
  }
  # The density is 0 if 1 + shape * (x - loc) / scale <= 0
  # Set the score to 0 for these cases
  zx <- shape * x / scale
  if (any(zerod <- 1 + zx <= 0 & !invalidScale)) {
    infoArray[zerod] <- 0
  }
  # Otherwise, the density is positive
  if (any(pos <- !zerod & !invalidScale)) {
    infoArray[pos, 1, 1] <- gp11(x[pos], scale[pos], shape[pos], ...)
    infoArray[pos, 1, 2] <- infoArray[pos, 2, 1] <-
      gp12(x[pos], scale[pos], shape[pos], ...)
    infoArray[pos, 2, 2] <- gp22(x[pos], scale[pos], shape[pos], ...)
  }
  parNames <- c("scale", "shape")
  dimnames(infoArray) <- list(NULL, parNames, parNames)
  if (sum) {
    infoArray <- colSums(infoArray)
  }
  return(-infoArray)
}
