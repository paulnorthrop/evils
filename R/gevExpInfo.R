#' GEV Distribution Expected Information
#'
#' Calculates the expected information matrix for the GEV distribution.
#'
#' @param sigma,xi Numeric vectors. Respective values of the GEV parameters
#'   \eqn{\sigma} and \eqn{\xi}. For `gevExpInfo`, `sigma` and `xi` must have
#'   length 1.
#' @param eps A numeric scalar. For values of \eqn{\xi} in `xi` that lie in
#'   `(-eps, eps)` an approximation is used instead of a direct calculation.
#'   See **Details**. If `eps` is a vector then only the first element is used.
#' @details `gevExpInfo` calculates, for single pair of values
#'   \eqn{(\sigma, \xi) = } `(sigma, xi)`, the expected information matrix for a
#'   single observation from a GEV distribution with distribution function
#'   \deqn{F(x) = P(X \leq x) = \exp\left\{ -\left[ 1+\xi\left(\frac{x-\mu}{\sigma}\right)
#'   \right]_+^{-1/\xi} \right\},}
#'   where \eqn{x_+ = \max(x, 0)}.
#'   The GEV expected information is defined only for \eqn{\xi > -0.5} and does
#'   not depend on the value of \eqn{\mu}.
#'
#'   The other functions are vectorised and calculate the individual
#'   contributions to the expected information matrix. For example, `imm`
#'   calculates the expectation \eqn{i_{\mu\mu}} of the negated second
#'   derivative of the GEV log-density with respect to \eqn{\mu}, that is, each
#'   `m` indicates one derivative with respect to \eqn{\mu}. Similarly, `s`
#'   denotes one derivative with respect to \eqn{\sigma} and `x` one derivative
#'   with respect to \eqn{\xi}, so that, for example, `isx` calculates the
#'   expectation \eqn{i_{\sigma\xi}} of the negated GEV log-density after one
#'   taking one derivative with respect to \eqn{\sigma} and one derivative with
#'   respect to \eqn{\xi}. Note that \eqn{i_{\xi\xi}} depends only on
#'   \eqn{\xi}.
#'
#'   The expectation in `imm` can be calculated in a directly way for all
#'   \eqn{\xi > -0.5}. For the other components, direct calculation of the
#'   expectation is unstable when \eqn{\xi} is close to 0. Instead, we use
#'   a quadratic approximation over `(-eps, eps)`, from a Lagrangian
#'   interpolation of the values from the direct calculation for \eqn{\xi = }
#'   `-eps`, 0 and `eps`.
#' @returns `gevExpInfo` returns a 3 by 3 numeric matrix with row and column
#'   names `mu, sigma, xi`. The other functions return a numeric vector of
#'   length equal to the maximum of the lengths of the arguments, excluding
#'   `eps`.
#' @examples
#' # Expected information matrices for ...
#' # ... sigma = 2 and xi = -0.4
#' gevExpInfo(2, -0.4)
#' # ... sigma = 3 and xi = 0.001
#' gevExpInfo(3, 0.001)
#' # ... sigma = 3 and xi = 0
#' gevExpInfo(3, 0)
#' # ... sigma = 1 and xi = 0.1
#' gevExpInfo(1, 0.1)
#'
#' # The individual components of the latter matrix
#' imm(1, 0.1)
#' ims(1, 0.1)
#' imx(1, 0.1)
#' iss(1, 0.1)
#' isx(1, 0.1)
#' ixx(0.1)
#' @name gevExpInfo
NULL
## NULL

#' @rdname gevExpInfo
#' @export
imm <- function(sigma, xi) {
  m <- max(length(sigma), length(xi))
  return(rep_len(pxi(xi), m) / rep_len(sigma ^ 2, m))
}

#' @rdname gevExpInfo
#' @export
iss <- function(sigma, xi, eps = 3e-3) {
  eps <- eps[1]
  val <- gevExpInfoComp(fun = issFn, fun0 = iss0Constant, xi = xi, eps = eps)
  m <- max(length(sigma), length(xi))
  return(rep_len(val, m) / rep_len(sigma ^ 2, m))
}

#' @rdname gevExpInfo
#' @export
ixx <- function(xi, eps = 3e-3) {
  eps <- eps[1]
  val <- gevExpInfoComp(fun = ixxFn, fun0 = ixx0Constant, xi = xi, eps = eps)
  return(val)
}

#' @rdname gevExpInfo
#' @export
ims <- function(sigma, xi, eps = 3e-3) {
  eps <- eps[1]
  val <- gevExpInfoComp(fun = imsFn, fun0 = ims0Constant, xi = xi, eps = eps)
  m <- max(length(sigma), length(xi))
  return(rep_len(val, m) / rep_len(sigma ^ 2, m))
}

#' @rdname gevExpInfo
#' @export
imx <- function(sigma, xi, eps = 3e-3) {
  eps <- eps[1]
  val <- gevExpInfoComp(fun = imxFn, fun0 = imx0Constant, xi = xi, eps = eps)
  m <- max(length(sigma), length(xi))
  return(rep_len(val, m) / rep_len(sigma, m))
}

#' @rdname gevExpInfo
#' @export
isx <- function(sigma, xi, eps = 3e-3) {
  eps <- eps[1]
  val <- gevExpInfoComp(fun = isxFn, fun0 = isx0Constant, xi = xi, eps = eps)
  m <- max(length(sigma), length(xi))
  return(rep_len(val, m) / rep_len(sigma, m))
}

#' @rdname gevExpInfo
#' @export
gevExpInfo <- function(sigma, xi, eps = 3e-3) {
  if (xi <= -0.5) {
    stop("The GEV expected information is undefined for xi <= -0.5")
  }
  # The expected information does not depend on mu
  val <- matrix(NA, 3, 3)
  val[1, 1] <- imm(sigma, xi)
  val[2, 2] <- iss(sigma, xi, eps)
  val[3, 3] <- ixx(xi, eps)
  val[2, 1] <- val[1, 2] <- ims(sigma, xi, eps)
  val[3, 1] <- val[1, 3] <- imx(sigma, xi, eps)
  val[3, 2] <- val[2, 3] <- isx(sigma, xi, eps)
  dimnames(val) <- list(c("mu", "sigma", "xi"), c("mu", "sigma", "xi"))
  return(val)
}
