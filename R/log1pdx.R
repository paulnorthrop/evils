#' Accurate `log(1+x)/x` computation
#'
#' Compute \eqn{\log(1+x)/x} accurately also for small \eqn{x}, that is,
#' \eqn{|x| \ll 1}{|x| << 1}.
#'
#' @param x A numeric vector with values \eqn{x > -1}.
#' @param minL1 A negative cutoff \eqn{m_l}. For \eqn{x \in (m_l, 1)} the
#'   computations are based on the series expansion 4.1.29 on page 68 of
#'   Abramowitz and Stegun (1972).
#' @param eps2 A non-negative cutoff \eqn{\epsilon_2}. For
#'   \eqn{|x| > \epsilon_2}, [`logcf`][`DPQ::logcf`] is used when computing
#'   \eqn{\log(1+x)/x} using the series expansion. Otherwise, that is,
#'   \eqn{|x| \leq \epsilon_2}, only a few terms of the expansion are used.
#' @param tol_logcf A non-negative number indicating the tolerance (maximal
#'   relative error) for the [`logcf`][`DPQ::logcf`] function.
#' @param trace.lcf A logical scalar. Used in
#'   [`logcf`][`DPQ::logcf`]`(.., trace = trace.lcf)`.
#' @param logCF The function to be used as [`logcf`][`DPQ::logcf`]. The default
#'   chooses the pure \R `logcfR()` when `x` is not numeric, and chooses the
#'   C-based `logcf()` when `is.numeric(x)` is `TRUE`.
#' @details For \eqn{x} in \eqn{(m_l, 1)} the computations are based on the
#'   series expansion
#' \deqn{\log(1+x) = 2 \left( t+\frac{t^3}{3}+\frac{t^5}{5}+\cdots \right),}
#'   where \eqn{t = x/(2+x)} which is formula 4.1.29 on page 68 of Abramowitz and Stegun (1972).
#'
#' Different computations are used in 3 different ranges of \eqn{x}, that is,
#'
#'  * \eqn{x < m_l} or \eqn{x > 1}: `log1pdx(x):= ` [`log1p`] `/ x`.
#'  * \eqn{|x| \leq \epsilon_2}:
#'    \eqn{\frac{2}{2+x}((((\frac19y+\frac17)y+\frac15)y+\frac13)y+1)}
#'  * \eqn{x \in (m_l, 1)} and \eqn{|x| > \epsilon_2}:
#'    \eqn{\frac{2}{2+x}} [`logcf`]\eqn{(y, 1, 2)}, where \eqn{y = t^2}.
#'
#' @return A numeric vector (with the same attributes as `x`).
#' @author Paul Northrop created this function by modifying the code in
#'   [`log1pmx`][`DPQ::log1pmx`].
#' @references Abramowitz, M. and Stegun, I. A. (1972) Handbook of Mathematical
#' Functions. New York: Dover.
#' [https://en.wikipedia.org/wiki/Abramowitz_and_Stegun](https://en.wikipedia.org/wiki/Abramowitz_and_Stegun)
#' provides links to the full text, which is in public domain.

#' @examples
#' # In the limit as x tends to 0 log(1+x)/x = 1
#' log1pdx(0)
#'
#' #
#' x <- seq(from = -1, to = 4, by = 0.01)
#' y1 <- log1pdx(x)
#' y2 <- log1p(x) / x
#' y3 <- log(1 + x) / x
#' y <- cbind(y1, y2, y3)
#' matplot(x, y, type = "l")
#'
#' #
#' ep <- 1e-8
#' x <- seq(from = -ep, to = ep, len = 10001)
#' y1 <- log1pdx(x)
#' y2 <- log1p(x) / x
#' y3 <- log(1 + x) / x # y3 is known to behave poorly near 0
#' y <- cbind(y1, y2)
#' matplot(x, y, type = "l", lty = 1, col = c("black", "red", "green"))
#' @export
log1pdx <- function(x, minL1 = -0.79149064, eps2 = 0.01, tol_logcf = 1e-14,
                    trace.lcf = FALSE,
                    logCF = if (is.numeric(x)) DPQ::logcf else DPQ::logcfR.) {

  # Check input values of eps2 and minL1
  stopifnot(is.numeric(eps2), eps2 >= 0, is.numeric(minL1),
            -1 <= minL1, minL1 < 0)

  # Create return vector of the correct length
  # There are 3 cases:
  #
  #
  #
  r <- x
  if (any(c1 <- (x > 1 | x < minL1)))
    r[c1] <- log1p(x[c1]) / x[c1]
  ## else { ## ##/* expand in [x/(2+x)]^2 */
  if (any(c2 <- !c1)) {
    x <- x[c2]
    term <- x / (2 + x)
    term2 <- 2 / (2 + x)
    y <- term * term
    r[c2] <- term2 * { ## not using ifelse(), rather what works with "mpfr"
      ## ifelse(abs(x) < eps2,
      ##        (((2 / 9 * y + 2 / 7) * y + 2 / 5) * y + 2 / 3) * y - x,
      ##        2 * y * logcf(y, 3, 2, tol_logcf) - x)
      A <- x
      if (any(isSml <- abs(x) <= eps2)) {
        i <- which(isSml)
        y. <- y[i]
        # Becomes y-like (e.g. mpfr)
        z <- if (is.numeric(x)) 1 else 1 + 0 * y.
        A[i] <- ((((z / 9 * y. + z / 7) * y. + z / 5) * y. + z / 3) * y. + 1)
      }
      if (length(iLrg <- which(!isSml))) {
        y. <- y[iLrg]
        A[iLrg] <- logCF(y., 1, 2, tol_logcf, trace = trace.lcf)
      }
      A
    }
  }
  return(r)
}