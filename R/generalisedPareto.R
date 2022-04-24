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
#' @param individual A logical scalar. Relevant to \code{gpLogLikelihood} and
#'   \code{gpScore}. If \code{individual = FALSE} then only the sum of
#'   contributions from all observations in \code{excesses} is calculated.  If
#'   \code{indidivdual = TRUE} then individual contributions from each
#'   observation in \code{excesses} are calculated.
#' @param tol A positive numeric scalar.  If \eqn{|\xi| <} \code{tol} then some
#'   of the required quantities are approximated using a series expansion for
#'   the quantity.  See \strong{Details}.
#' @param epsilon The desired error margin for an approximation used when
#'    \eqn{|\xi| <} \code{tol}.  This is passed to
#'    \code{\link[sumR]{infiniteSum}} as the argument \code{epsilon}.
#' @details If \eqn{|\xi| \geq}{|\xi| >=} \code{tol} all quantities are
#'   calculated in a direct manner, using R functions that correspond to the
#'   expressions involved.
#'
#'   If \eqn{\xi = 0} then all quantities are calculated directly, using
#'   expressions based on limiting values as \eqn{\xi} tends to zero where
#'   necessary.
#'
#'   If \eqn{|\xi| <} \code{tol}, for those quantities for which this direct
#'   calculation is unreliable when \eqn{|\xi|} is very close to zero,
#'   \code{\link[sumR]{infiniteSum}} is used to approximate their values.  The
#'   theoretical error margin is controlled using the argument \code{epsilon}.
#'   The algorithms are described in the \code{\link[sumR]{infiniteSum}}
#'   documentation.  If \eqn{\xi < 0} then the Sum-to-threshold method is used.
#'   If \eqn{\xi > 0} then the Batches method is used.
#'
#'   The following notes for which quantities we need to take this approach
#'   when \eqn{\xi} is very small.
#'
#'   \strong{Log-likelihood} (\code{gpLogLikelihood}). The term of the
#'   log-likelihood in which the observations appear.
#'
#'   \strong{Score} (\code{gpSCore}).  The second element of the score, that
#'   is, the contributions to the score from the derivatives of the the
#'   log-likelihood with respect to the shape parameter \eqn{\xi}.
#'
#'   \strong{Observed information} (\code{gpInfo}).  The \code{[2, 2]} element
#'   of the matrix, corresponding to the negated second derivative of the
#'   log-likelihood with respect to \eqn{\xi}.
#' @return
#'   \strong{Log-likelihood} (\code{gpLogLikelihood}). If
#'   \code{individual = FALSE} the value of the log-likelihood. If
#'   \code{individual = TRUE} a vector of length \code{length{excesses}}
#'   containing the contributions to the log-likelihood from each of the
#'   observations.
#'
#' \strong{Score} (\code{gpScore}).  If \code{individual = FALSE} the value
#'  of the score, a vector of length 2 containing the derivative of the
#'  log-likelihood evaluated at the input parameter values.
#'  If \code{individual = TRUE} the values of the contributions to the score
#'  from each of the observations, a
#'   \code{length(excesses)}\eqn{ \times 2}{ x 2} matrix.
#'   The columns are named \code{sigma[u]} and \code{xi}.
#'
#' \strong{Observed information} (\code{gpInfo}).  The observed information: a
#'   \eqn{2 \times 2}{2 x 2} matrix with row and column names
#'   \code{c(sigma[u], xi)}.
#'
#' If \code{\link[sumR]{infiniteSum}}, and the the maximum number of iterations
#' \code{maxIter} (= \code{1e+5}) in \code{\link[sumR]{infiniteSum}} has been
#' reached, then the returned object has a attribute that indicates where this
#' sign of a lack of convergence has occurred.
#' @name generalisedPareto
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
#' gpLogLikelihood(pars = c(1, 1e-8), excesses = y)
#' gpLogLikelihood(pars = c(1, -1e-8), excesses = y)
#' gpLogLikelihood(pars = c(1, 0), excesses = y)
#'
#' # Direct calculation, involving (1 + 1 / xi) * log1p(xi * y / sigmau)
#' # Mostly fine, but breaks down eventually
#' gpLogLikDirect(pars = c(1, 1e-309), excesses = y)
#' gpLogLikDirect(pars = c(1, -1e-309), excesses = y)
#' @rdname generalisedPareto
#' @export
gpLogLikelihood <- function(pars, excesses, individual = FALSE, tol = 1e-4,
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
  # Adjust epsilon based on the multiplier of the log term
  mult <- (x + 1) * (y / s)
  logterm <- mapply(log1pxOverx, x = zy, epsilon = epsilon / mult,
                    MoreArgs = list(tol = tol))
  gploglik <- -log(s) - mult * logterm
  if (!individual) {
    gploglik <- sum(gploglik)
  }
  return(gploglik)
}
