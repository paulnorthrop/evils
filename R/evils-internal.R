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
    stop("tol must be no larger than 1 and should be close to 0")
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

# ======================= Direct functions for checking ===================== #

# ----------------------------------- GP ------------------------------------ #

#' @keywords internal
#' @rdname evils-internal
#' @export
gpLogLikDirect <- function(pars, excesses, individual = FALSE) {
  y <- excesses
  # sigma_u
  s <- pars[1]
  if (s <= 0) {
    stop("The GP scale parameter must be positive")
  }
  # xi
  x <- pars[2]
  #
  z <- x / s
  zy <- z * y
  t0 <- 1 + zy
  if (any(t0 <= 0)) {
    stop("The likelihood is 0 for this combination of data and parameters")
  }
  if (x == 0) {
    gploglik <- -log(s) - y / s
  } else {
    gploglik <- -log(s) - (1 + 1 / x) * log1p(zy)
  }
  if (!individual) {
    gploglik <- sum(gploglik)
  }
  return(gploglik)
}

# ==================== Simulation from a GP distribution ==================== #

#' @keywords internal
#' @rdname evils-internal
#' @export
rGenPareto <- function (n, loc = 0, scale = 1, shape = 0) {
  max_len <- ifelse(length(n) > 1, length(n), n)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  return(qGenPareto(stats::runif(n), loc = loc, scale = scale, shape = shape))
}

#' @keywords internal
#' @rdname evils-internal
qGenPareto <- function (p, loc = 0, scale = 1, shape = 0,
                        lower.tail = TRUE, log.p = FALSE) {
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  if (length(p) == 0) {
    return(numeric(0))
  }
  if (!log.p & any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("invalid p: p must be in [0,1].")
  }
  max_len <- max(length(p), length(loc), length(scale), length(shape))
  p <- rep_len(p, max_len)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  mult <- BoxCoxVec(x = 1 - p, lambda = -shape)
  return(loc - scale * mult)
}

#' @keywords internal
#' @rdname evils-internal
BoxCoxVec <- function(x, lambda = 1, lambda_tol = 1e-6) {
  #
  # Computes the Box-Cox transformation of a vector.  If lambda is very close
  # to zero then a first order Taylor series approximation is used.
  #
  # Args:
  #   x          : A numeric vector. (Non-negative) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  # Returns:
  #   A numeric vector.  The transformed value
  #     (x^lambda - 1) / lambda
  #
  if (any(x < 0, na.rm = TRUE)) {
    stop("Invalid x: x must be non-negative")
  }
  max_len <- max(length(x), length(lambda))
  x <- rep_len(x, max_len)
  lambda <- rep_len(lambda, max_len)
  retval <- ifelse(abs(lambda) > lambda_tol, (x ^ lambda - 1) / lambda,
                   ifelse(lambda == 0, log(x),
                          ifelse(is.infinite(x),
                                 ifelse(lambda < 0, -1 / lambda, Inf),
                                 ifelse(x == 0, ifelse(lambda > 0, -1 / lambda, -Inf),
                                        log(x) * (1 + log(x) * lambda / 2)))))
  return(retval)
}
