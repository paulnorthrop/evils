#' Padé approximant coefficients for interpolation
#'
#' Uses power series coefficients \eqn{\{a_n, n = 0, 1, ..., T\}} to calculate
#' Padé \eqn{[L/M]} approximant coefficients, with the option to use one extra
#' pair of coefficients to adjust the Padé approximation to interpolate between
#' two points.
#'
#' @param L An integer. Order of the Padé numerator.
#' @param M An integer. Order of the Padé denominator.
#' @param A A numeric vector of power series coefficients, starting at
#'   \eqn{x^0}. Typically, these are the coefficients of a truncated Taylor
#'   series of a function \eqn{f(x)} about some point.
#' @param f A function that returns the value of \eqn{f(x)}.
#' @param xint A numeric vector of length 2. The values of \eqn{x} between
#'   which interpolation is required.
#' @param x A numeric vector. Values at which to evaluate a Padé approximation
#'   of \eqn{f(x)}.
#' @details If `f` is missing then this function behaves in the same way as
#'   the function `Pade` in the `Pade` package.
#'
#'   If `f` is supplied then the Padé approximant is adjusted by adding a extra
#'   pair of coefficient, one coefficient of \eqn{x^{L+1}} in the numerator and
#'   one coefficient of \eqn{x^{M+1}} in the denominator, so that the
#'   approximation recovers the values of `f(xint[1])` and `f(xint[2])`.
#' @author Paul Northrop created this function by modifying the code in the
#'   `Pade` function in the `Pade` package.
#' @return `padeInterp` returns an object (a list) a class `"pade"` with two
#'   components:
#'
#'   * `Px` Coefficients of the numerator polynomial, starting at \eqn{x^0}.
#'   * `Qx` Coefficients of the denominator polynomial, starting at \eqn{x^0}.
#'
#' `predict.pade` returns a numeric vector of length `length(x)`.
#' @references Adler A (2015). Pade: Padé Approximant Coefficients.
#' \doi{10.5281/zenodo.4270254} R package version 1.0.6,
#' \url{https://CRAN.R-project.org/package=Pade}.
#' @examples
#' ## f(x) = log(1 + x) / x
#'
#' # Orders of the Padé numerator and denominator.
#' L <- M <- 1
#' # Taylor series coefficients of expansion about x = 0, where f(0) = 1
#' j <- seq_len(L + M + 1) - 1
#' A <- (-1) ^ j / (j + 1)
#'
#' # Pade approximation
#' p1 <- padeInterp(L, M, A)
#'
#' # Pade approximation adjusted to interpolate over [-0.1, 0.1]
#' f <- function(x) log1p(x) / x
#' xint <- c(-0.1, 0.1)
#' p2 <- padeInterp(L, M, A, f, xint)
#'
#' # Evaluate the approximations at xint and 0
#' xvals <- c(xint[1], 0, xint[2])
#' predict(p2, xvals)
#' # Check that p2 recovers the values of f at (xint[1], 0, xint[2])
#' f(xint)
#' # The standard approximation does not interpolate
#' predict(p1, xvals)
#' @name padeInterp
NULL
## NULL

#' @rdname padeInterp
#' @export
padeInterp <- function (L, M, A, f, xint) {
  if (floor(L) != L || floor(M) != M) {
    stop("Polynomial orders need to be integers.")
  }
  lPlus1 <- L + 1L
  matSize <- lPlus1 + M
  if (length(A) < matSize) {
    stop("Not enough Taylor series coefficients provided.")
  }
  PQ <- matrix(0, ncol = matSize, nrow = matSize)
  PQ[1:lPlus1, 1:lPlus1] <- -diag(lPlus1)
  for (i in seq_len(M)) {
    PQ[, lPlus1 + i] <- c(rep.int(0, i), head(A, (matSize - i)))
  }
  # If f is not supplied then use Pade::Pade()
  # Otherwise, extend PQ and calculate the extra pair of coefficients
  if (missing(f)) {
    padeCoeff <- solve(PQ, -head(A, matSize))
    numer <- head(padeCoeff, lPlus1)
    denom <- c(1, tail(padeCoeff, M))
  } else {
    mPlus1 <- M + 1L
    x1 <- x[1]
    x2 <- x[2]
    f1 <- f(x1)
    f2 <- f(x2)
    lpow <- seq_len(lPlus1) - 1
    mpow <- seq_len(M)
    coef1 <- c(x1 ^ lpow, -f1 * x1 ^ mpow, x1 ^ lPlus1, -f1 * x1 ^ mPlus1)
    coef2 <- c(x2 ^ lpow, -f2 * x2 ^ mpow, x2 ^ lPlus1, -f2 * x2 ^ mPlus1)
    PQ <- rbind(cbind(PQ, 0, 0), coef1, coef2)
    headA <- c(head(A, matSize), -f1, -f2)
    padeCoeff <- solve(PQ, -headA)
    numer <- c(head(padeCoeff, lPlus1), padeCoeff[lPlus1 + mPlus1])
    Qcoefs <- c((lPlus1 + 1L):(lPlus1 + M), lPlus1 + mPlus1 + 1L)
    denom <- c(1, padeCoeff[Qcoefs])
  }
  val <- list(Px = numer, Qx = denom)
  class(val) <- "pade"
  return(val)
}

#' @rdname padeInterp
#' @export
predict.pade <- function(object, x) {
  nP <- length(object$Px)
  num <- outer(x[1:length(x)], 0:(nP - 1), "^") %*% object$Px
  nQ <- length(object$Qx)
  den <- outer(x[1:length(x)], 0:(nQ - 1), "^") %*% object$Qx
  val <- num / den
  dim(val) <- dim(x)
  return(val)
}
