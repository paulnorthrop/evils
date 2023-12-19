#' Generalised Extreme Value functions
#'
#' Calculate the log-likelihood function, score and observed information for
#' a random sample from a generalised extreme value (GEV) distribution,
#' including cases where the shape parameter is very close to zero.
#'
#' @param pars A numeric parameter vector of length 3 containing the respective
#'   values of the GEV location \eqn{\mu}, scale \eqn{\sigma} and shape
#'   \eqn{\xi} parameters.
#' @param maxima A numeric vector of observations. Typically, these are
#'   block maxima, that is, the largest observation in a block of contiguous
#'   observations.
#' @param individual A logical scalar. Relevant to `gevLoglik` and
#'   `gevScore`. If `individual = FALSE` then only the sum of
#'   contributions from all observations in `maxima` is calculated.  If
#'   `individual = TRUE` then individual contributions from each
#'   observation in `maxima` are calculated.
#' @param tol A positive numeric scalar.  Tolerance used to determine whether
#'   to perform a calculation directly or via a series expansion approximation.
#'   See **Details**.
#' @param epsilon The desired error margin when an approximation is used.
#' @details
#'   **Log-likelihood** (`gevLoglik`). The two problematic
#'   terms of the log-likelihood both involve
#'   \ifelse{html}{log(1+z)/z}{\eqn{\log(1+z)/z}},
#'   where \ifelse{html}{z=\eqn{\xi}\eqn{(y - \mu)} / \eqn{\sigma}}{
#'   \eqn{z} = \eqn{\xi}\eqn{(y - \mu)} / \eqn{\sigma}} and where \eqn{y} is a
#'   sample maximum. In one part this is exponentiated, in the other it is not.
#'   If \eqn{|z| \geq}{|z| >=} `tol` then this is calculated directly,
#'   using `log1p(z)/z`.
#'   If \eqn{|z| <} `tol` then we use [sumR::infiniteSum()]
#'   to approximate the series \ifelse{html}{log(1+z)/z}{\eqn{\log(1+z)/z}}
#'   \ifelse{html}{= 1 - z/2 + z\out{<sup>2</sup>}/3 - z\out{<sup>3</sup>}/4 +
#'   ...}{\eqn{= 1 - z/2 + z^2/3 - z^3/4 + \cdots}}. Before the call to
#'   [sumR::infiniteSum()] the input value of `epsilon`
#'   is adjusted to achieve the desired error margin for the approximation of
#'   the log-likelihood, taking into account the error of approximation from
#'   both terms. If `z = 0` then
#'   \ifelse{html}{log(1+z)/z = 1}{\eqn{\log(1+z)/z} = 1}.
#' @return
#'   **Log-likelihood** (`gevLoglik`). If
#'   `individual = FALSE` the value of the log-likelihood. If
#'   `individual = TRUE` a vector of length \code{length{maxima}}
#'   containing the contributions to the log-likelihood from each of the
#'   observations.
#'
#' **Score** (`gevScore`).  If `individual = FALSE` the value
#'  of the score, a vector of length 2 containing the derivative of the
#'  log-likelihood evaluated at the input parameter values.
#'  If `individual = TRUE` the values of the contributions to the score
#'  from each of the observations, a
#'   `length(maxima)`\eqn{ \times 2}{ x 2} matrix.
#'   The columns are named `sigma[u]` and `xi`.
#'
#' **Observed information** (`gevInfo`).  The observed information: a
#'   \eqn{2 \times 2}{2 x 2} matrix with row and column names
#'   `c(sigma[u], xi)`.
#' @name gev
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
#' # Approximation using sumR::infinitesum()
#' gevLoglik(pars = c(0, 1, 1e-8), maxima = y, individual = TRUE)
#' gevLoglik(pars = c(0, 1, -1e-8), maxima = y, individual = TRUE)
#' gevLoglik(pars = c(0, 1, 0), maxima = y, individual = TRUE)
#'
#' # Direct calculation, involving (1 / xi) * log1p(xi * (y - mu) / sigma)
#' # Mostly fine, but breaks down eventually
#' gevLoglikDirect(pars = c(0, 1, 1e-309), maxima = y)
#' gevLoglikDirect(pars = c(0, 1, -1e-309), maxima = y)
#' @rdname gev
#' @export
gevLoglik <- function(pars, maxima, individual = FALSE, tol = 1e-4,
                      epsilon = 1e-15) {
  if (tol <= 0) {
    stop("'tol' must be positive")
  }
  if (tol > 1) {
    stop("tol must be no larger than 1 and should be close to 0")
  }
  y <- maxima
  # mu
  m <- pars[1]
  # sigma
  s <- pars[2]
  if (s <= 0) {
    stop("The GEV scale parameter must be positive")
  }
  # xi
  x <- pars[3]
  #
  zw <- x * (y - m) / s
  t0 <- 1 + zw
  if (any(t0 <= 0)) {
    stop("The likelihood is 0 for this combination of data and parameters")
  }
  # Adjust epsilon based on the errors of approximation for each of the 2 terms
  # in the log-likelihood.  We work with epsilon / 2, because there are 2
  # terms, and make an adjustment for the exp() in the 2nd term
  mult1 <- (x + 1) * (y - m) / s
  mult2 <- (y - m) / s
  absmult1 <- abs(mult1)
  absmult2 <- abs(mult2)
  delta <- pmin(-log(1 - epsilon * exp(absmult2) / 2) / absmult2,
                log(1 + epsilon * exp(absmult2) / 2) / absmult2,
                epsilon / (2 * absmult1))
  logterm <- mapply(log1pxOverx, x = zw, epsilon = delta,
                    MoreArgs = list(tol = tol))
  gevloglik <- -log(s) - mult1 * logterm - exp(-mult2 * logterm)
  if (!individual) {
    gevloglik <- sum(gevloglik)
  }
  return(gevloglik)
}
