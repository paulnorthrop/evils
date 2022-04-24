#' Internal evils functions
#'
#' Internal evils functions
#' @details
#' These functions are not intended to be called by the user.
#' @name evils-internal
#' @keywords internal
NULL

# ================= Functions to calculate log(1 + x) / x =================== #

#' @keywords internal
#' @rdname evils-internal
log1pxOverx <- function(x, tol = 1e-4, epsilon = 1e-15) {
  if (tol > 1) {
    stop("tol must be no larger than 1 an should be close to 0")
  }
  if (abs(x) < tol) {
    val <- log1pxOverxApprox(x, epsilon = epsilon)
  } else {
    val <- log1pxOverxDirect(x)
  }
  return(val)
}

# ------- Approximation of log(1 + x) / x  using sumR::infiniteSum() -------- #

#' @keywords internal
#' @rdname evils-internal
log1pxOverxApprox <- function(x, epsilon = 1e-15) {
  if (x <= -1 || x >= 1) {
    stop("x must be in (-1, 1)")
  }
  if (x == 0) {
    return(1)
  }
  alternate <- as.numeric(x > 0)
  if (x < 0) {
    logL <- log(abs(x))
  } else {
    logL <- NULL
  }
  temp <- sumR::infiniteSum(log1pxOverxLogFunction, parameters = x, logL = logL,
                           alternate = alternate, epsilon = epsilon)
  return(exp(temp$sum))
}

#' @keywords internal
#' @rdname evils-internal
log1pxOverxLogFunction <- function(n, x) {
  # x must be in (-1, 1) for the series to converge
  # Return the log of the absolute value of the nth term in the series in which
  # contributions from all excesses are accumulated
  return(n * log(abs(x)) - log(n + 1))
}

# -------------------- Direct evaluation of log(1 + x) / x ------------------ #

#' @keywords internal
#' @rdname evils-internal
log1pxOverxDirect <- function(x) {
  if (x < -1) {
    stop("x must be greater than or equal to -1")
  }
  # Calculate ln(1+x) / x
  if (x == 0) {
    val <- 1
  } else {
    val <- exp(log(abs(log1p(x))) - log(abs(x)))
  }
  return(val)
}
