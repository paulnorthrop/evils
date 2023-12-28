#' Generalised Extreme Value functions
#'
#' Calculate the log-likelihood function, score and observed information for
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
#'   where \eqn{z = \xi (y - \mu)/ \sigma} and where \eqn{y} is a
#'   sample maximum. In one part this is exponentiated, in the other it is not.
#'   The function [`log1pdx`] is used to calculate \eqn{\log(1+z)/z}, with
#'   special case taken for cases where `z` is close to zero.
#' @return
#'   **Log-likelihood** (`gevLoglik`). If
#'   `individual = FALSE` the value of the log-likelihood. If
#'   `individual = TRUE` a vector of length `length{x}`
#'   containing the contributions to the log-likelihood from each of the
#'   observations.
#'
#' **Score** (`gevScore`).  If `individual = FALSE` the value
#'  of the score, a vector of length 2 containing the derivative of the
#'  log-likelihood evaluated at the input parameter values.
#'  If `individual = TRUE` the values of the contributions to the score
#'  from each of the observations, a
#'   `length(x)`\eqn{ \times 2}{ x 2} matrix.
#'   The columns are named `sigma[u]` and `xi`.
#'
#' **Observed information** (`gevInfo`).  The observed information: a
#'   \eqn{2 \times 2}{2 x 2} matrix with row and column names
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
#' x <- rGEV(10)
#' gevScore(x)
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
  #
  scoreMat <- matrix(NA, maxLen, 3)
  w <- x - loc
  zw <- shape * w / scale
  # Derivative of the log-likelihood with respect to loc = mu
  scoreMat[, 1] <- ((shape + 1) / (1 + zw) -
                      exp(-(shape + 1) * w * log1pdx(zw) / scale)) / scale
  # Derivative of the log-likelihood with respect to scale = sigma
  scoreMat[, 2] <- (w * scoreMat[, 1] - 1) / scale
  # Derivative of the log-likelihood with respect to shape = xi
  scoreMat[, 3] <- w ^ 2 * log1pdx2(zw) * (1 - exp(-w * log1pdx(zw) / scale)) /
    scale ^ 2 - w / (scale * (1 + zw))
  if (sum) {
    scoreMat <- colSums(scoreMat)
  }
  return(scoreMat)
}
