#' Non-homogeneous Poisson Process functions
#'
#' Calculate the log-likelihood function, score and observed information for
#' a random sample from a non-homogeneous Poisson Process (NHPP) model for
#' threshold exceedances, including cases where the shape parameter is very
#' close to zero.
#'
#' @param pars A numeric parameter vector of length 3 containing the respective
#'   values of the GEV location \eqn{\mu}, scale \eqn{\sigma} and shape
#'   \eqn{\xi} parameters of a Generalised Extreme Value (GEV) approximation to
#'   the distribution of \eqn{b}-year block maxima. The value of \eqn{b} is
#'   set using the argument \code{b}.
#' @param data A numeric vector of observations.
#' @param u A numeric scalar. The extreme value threshold applied to the data.
#' @param b A numeric scalar. The value of \eqn{b}.
#' @param nb A numeric scalar. The mean number of observations (including
#'   missing values) per \eqn{b}-year block.
#' @param individual A logical scalar. Relevant to \code{nhppLoglik} and
#'   \code{nhppScore}. If \code{individual = FALSE} then only the sum of
#'   contributions from all observations in \code{maxima} is calculated.  If
#'   \code{individual = TRUE} then individual contributions from each
#'   observation in \code{maxima} are calculated.
#' @param tol A positive numeric scalar.  Tolerance used to determine whether
#'   to perform a calculation directly or via a series expansion approximation.
#'   See \strong{Details}.
#' @param epsilon The desired error margin when an approximation is used.
#' @details
#'   \strong{Log-likelihood} (\code{nhppLoglik}). The two problematic
#'   terms of the log-likelihood both involve the function
#'   \ifelse{html}{log(1+z)/z}{\eqn{\log(1+z)/z}}. In one term
#'   \ifelse{html}{z=\eqn{\xi}\eqn{(u - \mu)} / \eqn{\sigma}}{
#'   \eqn{z} = \eqn{\xi}\eqn{(u - \mu)} / \eqn{\sigma}},
#'   in the other
#'   \ifelse{html}{z=\eqn{\xi}\eqn{(y - \mu)} / \eqn{\sigma}}{
#'   \eqn{z} = \eqn{\xi}\eqn{(y - \mu)} / \eqn{\sigma}},
#'   where \eqn{u} is the threshold and \eqn{y} is a observation that exceeds
#'   the threshold. In the first of these terms this function is exponentiated,
#'   in the other it is not.
#'   If \eqn{|z| \geq}{|z| >=} \code{tol} then this is calculated directly,
#'   using \code{log1p(z)/z}.
#'   If \eqn{|z| <} \code{tol} then we use \code{\link[sumR]{infiniteSum}}
#'   to approximate the series \ifelse{html}{log(1+z)/z}{\eqn{\log(1+z)/z}}
#'   \ifelse{html}{= 1 - z/2 + z\out{<sup>2</sup>}/3 - z\out{<sup>3</sup>}/4 +
#'   ...}{\eqn{= 1 - z/2 + z^2/3 - z^3/4 + \cdots}}. Before the call to
#'   \code{\link[sumR]{infiniteSum}} the input value of \code{epsilon}
#'   is adjusted to achieve the desired error margin for the approximation of
#'   the log-likelihood, taking into account the error of approximation from
#'   both terms. If \code{z = 0} then
#'   \ifelse{html}{log(1+z)/z = 1}{\eqn{\log(1+z)/z} = 1}.
#' @return
#'   \strong{Log-likelihood} (\code{gevLoglik}). If
#'   \code{individual = FALSE} the value of the log-likelihood. If
#'   \code{individual = TRUE} a vector of length \code{length{maxima}}
#'   containing the contributions to the log-likelihood from each of the
#'   observations.
#'
#' \strong{Score} (\code{gevScore}).  If \code{individual = FALSE} the value
#'  of the score, a vector of length 2 containing the derivative of the
#'  log-likelihood evaluated at the input parameter values.
#'  If \code{individual = TRUE} the values of the contributions to the score
#'  from each of the observations, a
#'   \code{length(maxima)}\eqn{ \times 2}{ x 2} matrix.
#'   The columns are named \code{sigma[u]} and \code{xi}.
#'
#' \strong{Observed information} (\code{gevInfo}).  The observed information: a
#'   \eqn{2 \times 2}{2 x 2} matrix with row and column names
#'   \code{c(sigma[u], xi)}.
#' @name nhpp
NULL
## NULL

#' Non-homogeneous Poisson Process Log-likelihood
#'
#' @examples
#' ### Simulate some data
#'
#' set.seed(17042022)
#' y <- rGenExtremeValue(365, 0, 1, 0)
#' u <- quantile(y, probs = 0.9)
#'
#' ### Log-likelihood
#'
#' # Approximation using sumR::infinitesum()
#' nhppLoglik(pars = c(0, 1, 1e-8), data = y, u = u, individual = TRUE)
#' nhppLoglik(pars = c(0, 1, -1e-8), data = y, u = u, individual = TRUE)
#' nhppLoglik(pars = c(0, 1, 0), data = y, u = u, individual = TRUE)
#'
#' # Direct calculation, involving (1 / xi) * log1p(xi * (y - mu) / sigma)
#' # and (1 / xi) * log1p(xi * (u - mu) / sigma)
#' # Mostly fine, but breaks down eventually
#' nhppLoglikDirect(pars = c(0, 1, 1e-323), data = y, u = u)
#' nhppLoglikDirect(pars = c(0, 1, -1e-323), data = y, u = u)
#' @rdname nhpp
#' @export
nhppLoglik <- function(pars, data, u, b = 1, nb = 365, individual = FALSE,
                       tol = 1e-4, epsilon = 1e-15) {
  if (tol <= 0) {
    stop("'tol' must be positive")
  }
  if (tol > 1) {
    stop("tol must be no larger than 1 and should be close to 0")
  }
  y <- data
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
  u <- as.numeric(u)
  zw <- x * (y - m) / s
  t0 <- 1 + zw
  zu <- x * (u - m) / s
  u0 <- 1 + zu
  if (any(t0 <= 0) || any(u0 <= 0)) {
    stop("The likelihood is 0 for this combination of data and parameters")
  }
  # Which observations is y exceed the threshold u?
  isExc <- y > u
  # Adjust epsilon based on the errors of approximation for each of the 2 terms
  # in the log-likelihood.  We work with epsilon / 2, because there are 2
  # terms, and make an adjustment for the exp() in the 2nd term
  mult1 <- (x + 1) * (y - m) / s
  mult2 <- (u - m) / s
  absmult1 <- abs(mult1)
  absmult2 <- abs(mult2)
  delta <- pmin(-log(1 - epsilon * exp(absmult2) / 2) / absmult2,
                log(1 + epsilon * exp(absmult2) / 2) / absmult2,
                epsilon / (2 * absmult1))
  ylogterm <- mapply(log1pxOverx, x = zw, epsilon = delta,
                     MoreArgs = list(tol = tol))
  ulogterm <- mapply(log1pxOverx, x = zu, epsilon = delta,
                     MoreArgs = list(tol = tol))
  nhpploglik <- -isExc * (log(b) + log(s) + mult1 * ylogterm) -
    exp(-mult2 * ulogterm) / nb
  if (!individual) {
    nhpploglik <- sum(nhpploglik)
  }
  return(nhpploglik)
}
