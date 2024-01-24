#' Accurate `log(1+x)/x` computation (and similar)
#'
#' Compute the function \eqn{h(x) = \log(1+x)/x} accurately also for small \eqn{x}, that is,
#' \eqn{|x| \ll 1}{|x| << 1}.  Similarly, compute the first and second
#' derivatives \eqn{h'(x) = x^{-1}(1+x)^{-1} - \log(1+x)/x^2} and
#' \eqn{h''(x) = 2 \log(1+x)/x^3 - 2 x^{-2}(1+x)^{-1} - x^{-1}(1+x)^{-2}}.
#'
#' @param x A numeric vector with values \eqn{x \geq -1}.
#' @param eps A numeric scalar. For values of \eqn{x} that lie in
#'   `(-eps, eps)` an approximation is used instead of a direct calculation.
#'   See **Details**.
#' @details For `dlog1pdx` and `d2log1pdx` a quadratic approximation over
#'   `(-eps, eps)`, from a Lagrangian interpolation of the values from the
#'   direct calculation for \eqn{x = } `-eps`, \eqn{0} and `eps`.
#' @return A numeric vector.
#' @examples
#' # In the limit as x tends to 0, log(1+x)/x = 1
#' log1pdx(0)
#'
#' # log1p(x) / x is fine unless x is 0 or extremely close to 0, e.g. 1e-324,
#' # when it returns NaN
#' x <- seq(from = -1, to = 2, by = 0.01)
#' y1 <- log1pdx(x)
#' y2 <- log1p(x) / x
#' y <- cbind(y1, y2)
#' matplot(x, y, type = "l")
#'
#' # In the limit as x tends to 0, 1/(x*(1+x)) - log(1+x)/x^2 tends to -1/2
#' dlog1pdx(0)
#'
#' # 1 / (x * (1 + x)) - log1p(x) / x ^ 2 is fine unless x is 0 (NaN returned)
#' # or extremely close to 0, e.g. 1e-15 returns 0.625, 1e-16 returns 0
#' # when it returns NaN
#' x <- seq(from = -1, to = 2, by = 0.01)
#' y1 <- dlog1pdx(x)
#' y2 <- 1 / (x * (1 + x)) - log1p(x) / x ^ 2
#' y <- cbind(y1, y2)
#' matplot(x, y, type = "l")
#'
#' # In the limit as x tends to 0, 2log(1+x)/x^3-2/(x^2*(1+x))-1/(x(1+x)^2)
#' # tends to 2/3
#' d2log1pdx(0)
#'
#' # 2log(1+x)/x^3-2/(x^2*(1+x))1/(x(1+x)^2) becomes unstable for
#' # abs(x) < 1e-6
#' x <- seq(from = -1, to = 2, by = 0.01)
#' y1 <- d2log1pdx(x)
#' y2 <- 2 * log1p(x) / x ^ 3 - 2 / (x ^ 2 * (1 + x)) - 1 / (x * (1 + x) ^ 2)
#' y <- cbind(y1, y2)
#' matplot(x, y, type = "l")
#' @name log1pdx
NULL
## NULL

#' @rdname log1pdx
#' @export
log1pdx <- function(x) {
  val <- log1p(x) / x
  val[x == 0] <- 1
  val[is.infinite(x) & x > 0] <- 0
  return(val)
}

#' @rdname log1pdx
#' @export
dlog1pdx <- function(x, eps = 1e-6) {
  fun <- function(x) {
    return(1 / (x * (1 + x)) - log1p(x) / x ^ 2)
  }
  val <- x
  # Deal with x = Inf case
  if (any(xIsInf <- is.infinite(x) & x > 0)) {
    val[xIsInf] <- 0
  }
  if (any(xNearZero <- abs(x) < eps & !xIsInf & !is.na(x))) {
    a <- -1 / 2
    yp <- fun(eps)
    ym <- fun(-eps)
    z <- x[xNearZero] / eps
    val[xNearZero] <- a + (yp - ym) * z / 2 + (yp + ym - 2 * a) * z ^ 2 / 2
  }
  val[!xNearZero & !xIsInf] <- fun(x[!xNearZero & !xIsInf])
  return(val)
}

#' @rdname log1pdx
#' @export
d2log1pdx <- function(x, eps = 1e-6) {
  fun <- function(x) {
    return(2 * log1p(x) / x ^ 3  - 2 / (x ^ 2 * (1 + x)) -
             1 / (x * (1 + x) ^ 2))
  }
  val <- x
  # Deal with x = Inf case
  if (any(xIsInf <- is.infinite(x) & x > 0)) {
    val[xIsInf] <- 0
  }
  if (any(xNearZero <- abs(x) < eps & !xIsInf & !is.na(x))) {
    a <- 2 / 3
    yp <- fun(eps)
    ym <- fun(-eps)
    z <- x[xNearZero] / eps
    val[xNearZero] <- a + (yp - ym) * z / 2 + (yp + ym - 2 * a) * z ^ 2 / 2
  }
  val[!xNearZero & !xIsInf] <- fun(x[!xNearZero & !xIsInf])
  return(val)
}
