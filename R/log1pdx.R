#' Accurate `log(1+x)/x` computation
#'
#' Compute \eqn{\log(1+x)/x} accurately also for small \eqn{x}, that is,
#' \eqn{|x| \ll 1}{|x| << 1}
#' @param epsilon The desired error margin when an approximation is used.
#' @details Add details
#'
#' * Do not need to treat \eqn{\xi = 0} as a special case.
#' * The expansion should not be in \eqn{\xi} but in something else.
#' * I can generalise the existing R function for a series approximation with
#' an arbitrary numebr of terms.
#' * Create a C++ version.
#'
#' @return Returns
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
## Accurate calculation of log(1+x)-x, particularly for small x.
## See also R-interface to R's C API [src/nmath/pgamma.c] log1pmx()
##     --> log1pmxC() in >> ./utils.R >> ../src/DPQ-misc.c
log1pdx <- function(x, tol_logcf = 1e-14,
                    eps2 = 0.01,
                    minL1 = -0.79149064, ## << was hard-wired 'minLog1Value' in R's source of log1pmx()
                    trace.lcf = FALSE,
                    logCF = if(is.numeric(x)) DPQ::logcf else DPQ::logcfR.)
{
  stopifnot(is.numeric(eps2), eps2 >= 0, is.numeric(minL1), -1 <= minL1, minL1 < 0)# < -1/4 ?
  r <- x
  if(any(c1 <- (x > 1 | x < minL1)))
    r[c1] <- log1p(x[c1]) / x[c1]
  ## else { ## ##/* expand in [x/(2+x)]^2 */
  if(any(c2 <- !c1)) {
    x <- x[c2]
    term <- x / (2 + x)
    term2 <- 2 / (2 + x)
    y <- term * term
    r[c2] <- term2 * { ## not using ifelse(), rather what works with "mpfr"
      ## ifelse(abs(x) < eps2,
      ##        (((2 / 9 * y + 2 / 7) * y + 2 / 5) * y + 2 / 3) * y - x,
      ##        2 * y * logcf(y, 3, 2, tol_logcf) - x)
      A <- x
      if(any(isSml <- abs(x) <= eps2)) {
        i <- which(isSml)
        y. <- y[i]
        z <- if(is.numeric(x)) 1 else 1 + 0*y. # becomes y-like (e.g. mpfr)
        A[i] <- ((((z / 9 * y. + z / 7) * y. + z / 5) * y. + z / 3) * y. + 1)
      }
      if(length(iLrg <- which(!isSml))) {
        y. <- y[iLrg]
        A[iLrg] <- logCF(y., 1, 2, tol_logcf, trace = trace.lcf)
      }
      A
    }
  }
  return(r)
}
