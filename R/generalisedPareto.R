#' Generalised Pareto functions
#'
#' Calculate the log-likelihood function, score and observed information for
#' a random sample from a generalised Pareto (GP) distribution, including cases
#' where the shape parameter is very close to zero.
#'
#' @param pars A numeric parameter vector of length 2 containing the respective
#'   values of the GP scale
#'   \ifelse{html}{\eqn{\sigma}\out{<sub>u</sub>}}{\eqn{\sigma_u}}
#'   and shape \eqn{\xi} parameters.
#' @param excesses A numeric vector containing positive observations.
#'   Typically, these are threshold excesses, that is, amounts by which
#'   exceedances of a threshold exceed that threshold.
#' @param individual A logical scalar. Relevant to `gpLoglik` and
#'   `gpScore`. If `individual = FALSE` then only the sum of
#'   contributions from all observations in `excesses` is calculated.  If
#'   `individual = TRUE` then individual contributions from each
#'   observation in `excesses` are calculated.
#' @param tol A positive numeric scalar.  Tolerance used to determine whether
#'   to perform a calculation directly or via a series expansion approximation.
#'   See **Details**.
#' @param epsilon The desired error margin when an approximation is used.
#' @details
#'   **Log-likelihood** (`gpLoglik`). The problematic part of
#'   the log-likelihood is \ifelse{html}{log(1+z)/z}{\eqn{\log(1+z)/z}},
#'   where \ifelse{html}{z=\eqn{\xi}\eqn{y} / \eqn{\sigma}\out{<sub>u</sub>}}{
#'   \eqn{z} = \eqn{\xi}\eqn{y} / \eqn{\sigma_u}} and where \eqn{y} is a
#'   sample threshold excess.
#'   If \eqn{|z| \geq}{|z| >=} `tol` then this is calculated directly,
#'   using `log1p(z)/z`.
#'   If \eqn{|z| <} `tol` then we use [sumR::infiniteSum()]
#'   to approximate the series \ifelse{html}{log(1+z)/z}{\eqn{\log(1+z)/z}}
#'   \ifelse{html}{= 1 - z/2 + z\out{<sup>2</sup>}/3 - z\out{<sup>3</sup>}/4 +
#'   ...}{\eqn{= 1 - z/2 + z^2/3 - z^3/4 + \cdots}}. Before the call to
#'   [sumR::infiniteSum()] the input value of `epsilon`
#'   is adjusted to achieve the desired error margin for the approximation of
#'   the log-likelihood. If `z = 0` then
#'   \ifelse{html}{log(1+z)/z = 1}{\eqn{\log(1+z)/z} = 1}.
#' @return
#'   **Log-likelihood** (`gpLoglik`). If
#'   `individual = FALSE` the value of the log-likelihood. If
#'   `individual = TRUE` a vector of length \code{length{excesses}}
#'   containing the contributions to the log-likelihood from each of the
#'   observations.
#'
#' **Score** (`gpScore`).  If `individual = FALSE` the value
#'  of the score, a vector of length 2 containing the derivative of the
#'  log-likelihood evaluated at the input parameter values.
#'  If `individual = TRUE` the values of the contributions to the score
#'  from each of the observations, a
#'   `length(excesses)`\eqn{ \times 2}{ x 2} matrix.
#'   The columns are named `sigma[u]` and `xi`.
#'
#' **Observed information** (`gpInfo`).  The observed information: a
#'   \eqn{2 \times 2}{2 x 2} matrix with row and column names
#'   `c(sigma[u], xi)`.
#' @name gp
NULL
## NULL

#' Generalised Pareto Log-likelihood
#'
#' @examples
#' ### Simulate some data
#'
#' set.seed(15042022)
#' y <- rGenPareto(100, 0, 1, 0)
#'
#' ### Log-likelihood
#'
#' # Approximation using sumR::infinitesum()
#' gpLoglik(pars = c(1, 1e-8), excesses = y)
#' gpLoglik(pars = c(1, -1e-8), excesses = y)
#' gpLoglik(pars = c(1, 0), excesses = y)
#'
#' # Direct calculation, involving (1 / xi) * log1p(xi * y / sigmau)
#' # Mostly fine, but breaks down eventually
#' gpLoglikDirect(pars = c(1, 1e-323), excesses = y)
#' gpLoglikDirect(pars = c(1, -1e-323), excesses = y)
#' @rdname gp
#' @export
gpLoglik <- function(pars, excesses, individual = FALSE, tol = 1e-4,
                     epsilon = 1e-15) {
  if (tol <= 0) {
    stop("'tol' must be positive")
  }
  if (tol > 1) {
    stop("tol must be no larger than 1 and should be close to 0")
  }
  y <- excesses
  # sigma_u
  s <- pars[1]
  if (s <= 0) {
    stop("The GP scale parameter must be positive")
  }
  # xi
  x <- pars[2]
  #
  zy <- x * y / s
  t0 <- 1 + zy
  if (any(t0 <= 0)) {
    stop("The likelihood is 0 for this combination of data and parameters")
  }
  # Adjust epsilon based on the multiplier of logterm below
  mult <- (x + 1) * y / s
  logterm <- mapply(log1pxOverx, x = zy, epsilon = epsilon / mult,
                    MoreArgs = list(tol = tol))
  gploglik <- -log(s) - mult * logterm
  if (!individual) {
    gploglik <- sum(gploglik)
  }
  return(gploglik)
}
