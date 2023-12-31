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
#' @param ... Further arguments to be passed to [`log1pdx`] or [`dlog1pdx`].
#'   See **Details**.
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
#'   \eqn{\ell} with respect to \eqn{\mu}, \eqn{\sigma} and \eqn{\xi}. `gevi` provides
#'   \eqn{\partial \ell / \partial \theta_i}, where
#'   \eqn{\theta = (\mu, \sigma, \xi)}. The input vectors
#'   `x`, `loc`, `scale` and `shape` must have the same lengths and satisfy the
#'   parameter constraints \eqn{\sigma > 0} and
#'   \eqn{1 + \xi(x - \mu) / \sigma > 0}. These conditions are not checked
#'   in `gev1`, `gev2` and `gev3`, so the user must take care when
#'   calling these functions.
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
#'   \eqn{\theta = (\mu, \sigma, \xi)}. The input vectors `x`, `loc`, `scale`
#'   and `shape` must have the same lengths and satisfy the parameter
#'   constraints \eqn{\sigma > 0} and \eqn{1 + \xi(x - \mu) / \sigma > 0}.
#'   These conditions are not checked in `gev11`, `gev12`, `gev13`, `gev22`,
#'   `gev23` and `gev33`, so the user must take care when calling these
#'   functions.
#' @return
#'   For `gevLoglik` and `gevScore` the number `n`, say, of combinations of the
#'   input parameters `x`, `loc`, `scale` and `shape` in the returned object is
#'   the maximum of the lengths of these arguments, with these arguments being
#'   recycled as necessary.
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
#'   **Observed information** (`gevInfo`).  If `sum = FALSE` the values of the
#'   contributions to the observed information, an `n` \eqn{ \times 3 \times 3}
#'   array. The second and third dimensions are named columns are named
#'   `loc`, `scale` and `shape`. If `sum = TRUE`, a \eqn{3 \times 3} matrix
#'   giving the observed information matrix, obtained by applying
#'   [`colsums`] to the array returned if `sum = FALSE`.
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
  if (any(pos <- !zerod & !invalidScale)) {
    # Derivative of the log-likelihood with respect to loc = mu
    scoreMat[pos, 1] <- gev1(x[pos], loc[pos], scale[pos], shape[pos])
    # Derivative of the log-likelihood with respect to scale = sigma
    scoreMat[pos, 2] <- gev2(x[pos], loc[pos], scale[pos], shape[pos])
    # Derivative of the log-likelihood with respect to shape = xi
    scoreMat[pos, 3] <- gev3(x[pos], loc[pos], scale[pos], shape[pos])
  }
  colnames(scoreMat) <- c("loc", "scale", "shape")
  if (sum) {
    scoreMat <- colSums(scoreMat)
  }
  return(scoreMat)
}

#' @rdname gevLikelihood
#' @export
gev1 <- function(x, loc, scale, shape) {
  w <- x - loc
  zw <- shape * w / scale
  val <- (shape + 1) / (1 + zw) - exp(-(shape + 1) * w * log1pdx(zw) / scale)
  return(val / scale)
}

#' @rdname gevLikelihood
#' @export
gev2 <- function(x, loc, scale, shape) {
  w <- x - loc
  val <- w * gev1(x, loc, scale, shape) - 1
  return(val / scale)
}

#' @rdname gevLikelihood
#' @export
gev3 <- function(x, loc, scale, shape) {
  w <- x - loc
  zw <- shape * w / scale
  val <- w ^ 2 * dlog1pdx(zw) * (exp(-w * log1pdx(zw) / scale) - 1) /
    scale ^ 2 - w / (scale * (1 + zw))
  return(val)
}

#' @rdname gevLikelihood
#' @export
gevInfo <- function(x, loc = 0, scale = 1, shape = 0, sum = FALSE, ...) {
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
    infoArray[invalidScale, ] <- NaN
  }
  w <- x - loc
  zw <- shape * w / scale
  # The density is 0 if 1 + shape * (x - loc) / scale <= 0
  # Set the score to 0 for these cases
  zw <- shape * (x - loc) / scale
  if (any(zerod <- 1 + zw <= 0 & !invalidScale)) {
    infoArray[zerod] <- 0
  }
  # Otherwise, the density is positive
  if (any(pos <- !zerod & !invalidScale)) {
    infoArray[pos, 1, 1] <-
      gev11(x[pos], loc[pos], scale[pos], shape[pos])
    infoArray[pos, 1, 2] <- infoArray[pos, 2, 1] <-
      gev12(x[pos], loc[pos], scale[pos], shape[pos])
    infoArray[pos, 1, 3] <- infoArray[pos, 3, 1] <-
      gev13(x[pos], loc[pos], scale[pos], shape[pos])
    infoArray[pos, 2, 2] <-
      gev22(x[pos], loc[pos], scale[pos], shape[pos])
    infoArray[pos, 2, 3] <- infoArray[pos, 3, 2] <-
      gev23(x[pos], loc[pos], scale[pos], shape[pos])
    infoArray[pos, 3, 3] <-
      gev33(x[pos], loc[pos], scale[pos], shape[pos])
  }
  parNames <- c("loc", "scale", "shape")
  dimnames(infoArray) <- list(NULL, parNames, parNames)
  if (sum) {
    infoArray <- colSums(infoArray)
  }
  return(infoArray)
}

#' @rdname gevLikelihood
#' @export
gev11 <- function(x, loc, scale, shape) {
  w <- x - loc
  zw <- shape * w / scale
  val <- (shape + 1) * (shape / (1 + zw) ^ 2 -
                          exp(-(2 * shape + 1) * w * log1pdx(zw) / scale))
  return(val / scale ^ 2)
}

#' @rdname gevLikelihood
#' @export
gev12 <- function(x, loc, scale, shape) {
  w <- x - loc
  val <- w * gev11(x, loc, scale, shape) - gev1(x, loc, scale, shape)
  return(val / scale)
}

#' @rdname gevLikelihood
#' @export
gev22 <- function(x, loc, scale, shape) {
  w <- x - loc
  val <- (1 - w * gev1(x, loc, scale, shape)) / scale +
    w * gev12(x, loc, scale, shape)
  return(val / scale)
}

#' @rdname gevLikelihood
#' @export
gev13 <- function(x, loc, scale, shape) {
  w <- x - loc
  z <- shape / scale
  zw <- shape * w / scale
  hzw <- log1pdx(zw)
  hdashzw <- dlog1pdx(zw)
  zwr <- 1 / (1 + zw)
  term1 <- scale * zwr - w * (1 + scale * z) * zwr ^2
  term2 <- w * hzw + (1 + scale * z) * w ^ 2 * hdashzw / scale
  term3 <- exp(-w * (1 + scale * z) * hzw / scale)
  val <- term1 + term2 * term3
  return(val / scale ^ 2)
}

#' @rdname gevLikelihood
#' @export
gev23 <- function(x, loc, scale, shape) {
  w <- x - loc
  val <- w * gev13(x, loc, scale, shape)
  return(val / scale)
}

#' @rdname gevLikelihood
#' @export
gev33 <- function(x, loc, scale, shape) {
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
