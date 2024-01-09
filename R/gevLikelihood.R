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
#' @param sum A logical scalar. Relevant to `gevLoglik` and `gevScore`. If
#'   `sum = TRUE` then only the sum of contributions from all combinations of
#'   observations in `x ` and parameter values in `loc`, `scale` and `shape`
#'   is calculated.  If `sum = FALSE` then these individual contributions are
#'   summed to produce the result.
#' @param ... Further arguments to be passed to [`log1pdx`], [`dlog1pdx`] or
#'   [`d2log1pdx`]. See [`gevDerivatives`].
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
#'   [`dlog1pdx`].
#'
#'   The functions `gev1`, `gev2` and `gev3` are used in
#'   `gevScore` to calculate the first derivatives of the GEV log-likelihood
#'   \eqn{\ell} with respect to \eqn{\mu}, \eqn{\sigma} and \eqn{\xi}.
#'   See [`gevDerivatives`].
#'
#'   **Observed information** (`gevObsInfo`). All the second derivatives of the
#'   log-likelihood involve \eqn{\log(1+z)/z}, calculated using [`log1pdx`].
#'   Any second derivative made with respect to \eqn{\xi} also involves
#'   \eqn{\log(1+z)/z^2 - z^{-1}(1+z)^{-1}}, which is calculated using
#'   [`dlog1pdx`]. Differentiating the log-likelihood twice with respect to
#'   \eqn{\xi} produces the function
#'   \eqn{2 \log(1+z)/z^3 - 2 z^{-2}(1+z)^{-1} - z^{-1}(1+z)^{-2}}, which is
#'   calculated using [`d2log1pdx`].
#'
#'   The functions `gev11`, `gev12`, `gev13`, `gev22`, `gev23` and `gev33` are
#'   used in `gevObsInfo` to calculate the second derivatives of the GEV
#'   log-likelihood \eqn{\ell} with respect to \eqn{\mu}, \eqn{\sigma} and
#'   \eqn{\xi}. `gevij` provides
#'   \eqn{\partial^2 \ell / \partial \theta_i \theta_j}, where
#'   \eqn{\theta = (\mu, \sigma, \xi)}.
#'   See [`gevDerivatives`].
#' @return
#'   The number `n`, say, of combinations of the input parameters `x`, `loc`,
#'   `scale` and `shape` in the returned object is the maximum of the lengths
#'   of these arguments, with these arguments being recycled as necessary.
#'
#'   **Log-likelihood** (`gevLoglik`). If `sum = FALSE` a vector of length `n`
#'   containing individual log-likelihood contributions. If `sum = TRUE` a
#'   scalar, the value of sum of these contributions.
#'
#'   **Score** (`gevScore`). If `sum = FALSE` the values of the contributions
#'   to the score, an `n`\eqn{ \times 3} matrix. The columns are named
#'   `loc`, `scale` and `shape`. If `sum = TRUE`, a vector of length 3
#'   containing the column sums of the matrix returned if `sum = FALSE`.
#'
#'   **Observed information** (`gevObsInfo`).  If `sum = FALSE` the values of
#'   the contributions to the observed information, an
#'   `n` \eqn{ \times 3 \times 3} array. The second and third dimensions are
#'   named columns are named `loc`, `scale` and `shape`. If `sum = TRUE`, a
#'   \eqn{3 \times 3} matrix giving the observed information matrix, obtained
#'   by applying [`colSums`] to the array returned if `sum = FALSE`.
#' @examples
#' ### Simulate some data
#' set.seed(28122023)
#' x <- rGEV(10)
#'
#' ### Log-likelihood
#'
#' # Calculation using log1pdx()
#' gevLoglik(x, 0, 1, 1e-8)
#' gevLoglik(x, 0, 1, -1e-8)
#' gevLoglik(x, 0, 1, 0)
#'
#' # Direct calculation, involving (1 / xi) * log1p(xi * (x - mu) / sigma)
#' # Mostly fine, but breaks down eventually
#' gevLoglikDirect(pars = c(0, 1, 1e-309), maxima = x)
#' gevLoglikDirect(pars = c(0, 1, -1e-309), maxima = x)
#'
#' # Score
#'
#' gevScore(x)
#' gevScore(x, sum = TRUE)
#'
#' # Information
#'
#' gpObsInfo(x)
#' gpObsInfo(x, sum = TRUE)
#' @name gevLikelihood
NULL
## NULL

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
  # The density is 0 if 1 + shape * (x - loc) / scale <= 0
  # Set the score to 0 for these cases
  zw <- shape * (x - loc) / scale
  if (any(zerod <- 1 + zw <= 0 & !invalidScale)) {
    scoreMat[zerod] <- 0
  }
  # Otherwise, the density is positive
  if (any(pos <- !zerod & !invalidScale)) {
    # Derivative of the log-likelihood with respect to loc = mu
    scoreMat[pos, 1] <- gev1(x[pos], loc[pos], scale[pos], shape[pos], ...)
    # Derivative of the log-likelihood with respect to scale = sigma
    scoreMat[pos, 2] <- gev2(x[pos], loc[pos], scale[pos], shape[pos], ...)
    # Derivative of the log-likelihood with respect to shape = xi
    scoreMat[pos, 3] <- gev3(x[pos], loc[pos], scale[pos], shape[pos], ...)
  }
  colnames(scoreMat) <- c("loc", "scale", "shape")
  if (sum) {
    scoreMat <- colSums(scoreMat)
  }
  return(scoreMat)
}

#' @rdname gevLikelihood
#' @export
gevObsInfo <- function(x, loc = 0, scale = 1, shape = 0, sum = FALSE, ...) {
  # Recycle the vector input x, loc, scale and shape, if necessary
  maxLen <- max(length(x), length(loc), length(scale), length(shape))
  x <- rep_len(x, maxLen)
  loc <- rep_len(loc, maxLen)
  scale <- rep_len(scale, maxLen)
  shape <- rep_len(shape, maxLen)
  # Create an array to store the results
  infoArray <- array(NA, dim = c(length(x), 3, 3))
  # The density, and, hence the information, is undefined if scale <= 0
  if (any(invalidScale <- scale <= 0)) {
    infoArray[invalidScale, , ] <- NaN
  }
  # The density is 0 if 1 + shape * (x - loc) / scale <= 0
  # Set the score to 0 for these cases
  zw <- shape * (x - loc) / scale
  if (any(zerod <- 1 + zw <= 0 & !invalidScale)) {
    infoArray[zerod, , ] <- 0
  }
  # Otherwise, the density is positive
  if (any(pos <- !zerod & !invalidScale)) {
    infoArray[pos, 1, 1] <-
      gev11(x[pos], loc[pos], scale[pos], shape[pos], ...)
    infoArray[pos, 1, 2] <- infoArray[pos, 2, 1] <-
      gev12(x[pos], loc[pos], scale[pos], shape[pos], ...)
    infoArray[pos, 1, 3] <- infoArray[pos, 3, 1] <-
      gev13(x[pos], loc[pos], scale[pos], shape[pos], ...)
    infoArray[pos, 2, 2] <-
      gev22(x[pos], loc[pos], scale[pos], shape[pos], ...)
    infoArray[pos, 2, 3] <- infoArray[pos, 3, 2] <-
      gev23(x[pos], loc[pos], scale[pos], shape[pos], ...)
    infoArray[pos, 3, 3] <-
      gev33(x[pos], loc[pos], scale[pos], shape[pos], ...)
  }
  parNames <- c("loc", "scale", "shape")
  dimnames(infoArray) <- list(NULL, parNames, parNames)
  if (sum) {
    infoArray <- colSums(infoArray)
  }
  return(-infoArray)
}
