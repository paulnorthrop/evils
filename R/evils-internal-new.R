#' Internal evils functions
#'
#' Internal evils functions
#' @details
#' These functions are not intended to be called by the user.
#' @name evils-internal
#' @keywords internal
NULL

# =============== For calculating the GEV expected information ============== #

# --------------------- Constants and helper functions ---------------------- #

#' @keywords internal
#' @rdname evils-internal
EulersConstant <- 0.57721566490153286060651209008240243104215933593992

#' @keywords internal
#' @rdname evils-internal
AperysConstant <- 1.202056903159594285399738161511449990764986292

#' @keywords internal
#' @rdname evils-internal
pxi <- function(xi) {
  return((1 + xi) ^ 2 * gamma(1 + 2 * xi))
}

#' @keywords internal
#' @rdname evils-internal
qxi <- function(xi) {
  return(gamma(2 + xi) * (digamma(1 + xi) + (1 + xi) / xi))
}

# ------------------------------ (sigma, sigma) ----------------------------- #

#' @keywords internal
#' @rdname evils-internal
iss0Fn <- function() {
  return(pi ^ 2 / 6 + (1 - EulersConstant) ^ 2)
}

#' @keywords internal
#' @rdname evils-internal
iss0Constant <- iss0Fn()

#' @keywords internal
#' @rdname evils-internal
issFn <- function(xi) {
  return((1 - 2 * gamma(2 + xi) + pxi(xi)) / (xi ^ 2))
}

# --------------------------------- (xi, xi) -------------------------------- #

#' @keywords internal
#' @rdname evils-internal
ixx0Fn <- function() {
  val <- pi ^ 2 / 6 - pi ^ 2 * EulersConstant / 2 + EulersConstant ^ 2 -
    EulersConstant ^ 3 - 2 * AperysConstant +
    2 * EulersConstant * AperysConstant + pi ^ 2 * EulersConstant ^ 2 / 4 +
    EulersConstant ^ 4 / 4 + 3 * pi ^ 4 / 80
  return(val)
}

#' @keywords internal
#' @rdname evils-internal
ixx0Constant <- ixx0Fn()

#' @keywords internal
#' @rdname evils-internal
ixxFn <- function(xi) {
  val <- (pi ^ 2 / 6 + (1 - EulersConstant + 1 / xi) ^ 2 - 2 * qxi(xi) / xi +
            pxi(xi) / xi ^ 2) / (xi ^ 2)
  return(val)
}

# -------------------------------- (mu, sigma) ------------------------------ #

#' @keywords internal
#' @rdname evils-internal
ims0Fn <- function() {
  return(EulersConstant - 1)
}

#' @keywords internal
#' @rdname evils-internal
ims0Constant <- ims0Fn()

#' @keywords internal
#' @rdname evils-internal
imsFn <- function(xi) {
  val <- (gamma(2 + xi) - pxi(xi)) / xi
  return(val)
}

# --------------------------------- (mu, xi) -------------------------------- #

#' @keywords internal
#' @rdname evils-internal
imx0Fn <- function() {
  return((pi ^ 2 / 6 + EulersConstant ^ 2 - 2 * EulersConstant) / 2)
}

#' @keywords internal
#' @rdname evils-internal
imx0Constant <- imx0Fn()

#' @keywords internal
#' @rdname evils-internal
imxFn <- function(xi) {
  val <- (pxi(xi) / xi - qxi(xi)) / xi
  return(val)
}

# ------------------------------- (sigma, xi) ------------------------------- #

#' @keywords internal
#' @rdname evils-internal
isx0Fn <- function() {
  val <- (4 * EulersConstant + 4 * AperysConstant + pi ^ 2 * EulersConstant +
     2 * EulersConstant ^ 3 - pi ^ 2 - 6 * EulersConstant ^ 2) / 4
  return(val)
}

#' @keywords internal
#' @rdname evils-internal
isx0Constant <- isx0Fn()

#' @keywords internal
#' @rdname evils-internal
isxFn <- function(xi) {
  val <- -(1 - EulersConstant + (1 - gamma(2 + xi)) / xi - qxi(xi) +
             pxi(xi) / xi) / (xi ^ 2)
  return(val)
}

# -- Calculate a component using a quadratic approximation if xi is near 0 -- #

#' @keywords internal
#' @export
gevExpInfoComp <- function(fun, fun0, xi, eps = 3e-3) {
  eps <- abs(eps)
  val <- xi
  if (any(xiNearZero <- abs(xi) < eps)) {
    aa <- fun0
    yp <- fun(eps)
    ym <- fun(-eps)
    xiOverEps <- xi[xiNearZero] / eps
    cc <- (yp + ym - 2 * aa) * xiOverEps ^ 2 / 2
    bb <- (yp - aa) * xiOverEps
    val[xiNearZero] <- aa + bb * xi[xiNearZero] + cc * xi[xiNearZero] ^ 2
  }
  val[!xiNearZero] <- fun(xi[!xiNearZero])
  return(val)
}
# ================= Functions to calculate log(1 + x) / x =================== #

#' @keywords internal
#' @rdname evils-internal
log1pxOverxVec <- function(x, tol = 1e-4, epsilon = 1e-15) {
  if (tol > 1) {
    stop("tol must be no larger than 1 and should be close to 0")
  }
  val <- rep(NA, length(x))
  cond <- abs(x) < tol
  val[cond] <- vapply(x[cond], log1pxOverxApprox, 0.0, epsilon = epsilon)
  val[!cond] <- log1pxOverxDirect(x[!cond])
  return(val)
}

#' @keywords internal
#' @rdname evils-internal
log1pxOverxVec2 <- function(x, tol = 1e-4, epsilon = 1e-15) {
  if (tol > 1) {
    stop("tol must be no larger than 1 and should be close to 0")
  }
  cond <- abs(x) < tol
  x[cond] <- vapply(x[cond], log1pxOverxApprox, 0.0, epsilon = epsilon)
  x[!cond] <- log1pxOverxDirect(x[!cond])
  return(x)
}

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
  temp <- sumR::infiniteSum(log1pxOverxLogFunction, parameters = x,
                            logL = logL, alternate = alternate,
                            epsilon = epsilon)
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
#  if any(x < -1) {
#    stop("x must be greater than or equal to -1")
#  }
  # Calculate ln(1+x) / x
  # We will never call this function for x == 0
  return(exp(log(abs(log1p(x))) - log(abs(x))))
}

# ======================= Direct functions for checking ===================== #

# ----------------------------------- GP ------------------------------------ #

#' @keywords internal
#' @rdname evils-internal
#' @export
gpLoglikDirect <- function(pars, excesses, individual = FALSE) {
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

# ---------------------------------- GEV ------------------------------------ #

#' @keywords internal
#' @rdname evils-internal
#' @export
gevLoglikDirect <- function(pars, maxima, individual = FALSE) {
  y <- maxima
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
  z <- x / s
  w <- y - m
  v <- w / s
  zw <- z * w
  t0 <- 1 + zw
  if (any(t0 <= 0)) {
    stop("The likelihood is 0 for this combination of data and parameters")
  }
  if (x == 0) {
    gevloglik <- -log(s) - v - exp(-v)
  } else {
    log1pzwx <- log1p(zw) / x
    gevloglik <- -log(s) - (x + 1) * log1pzwx - exp(-log1pzwx)
  }
  if (!individual) {
    gevloglik <- sum(gevloglik)
  }
  return(gevloglik)
}

# ---------------------------------- NHPP ----------------------------------- #

#' @keywords internal
#' @rdname evils-internal
#' @export
nhppLoglikDirect <- function(pars, data, u, b = 1, nb = 365,
                             individual = FALSE) {
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
  zw <- x * (y - m) / s
  t0 <- 1 + zw
  zu <- x * (u - m) / s
  u0 <- 1 + zu
  if (any(t0 <= 0) || any(u0 <= 0)) {
    stop("The likelihood is 0 for this combination of data and parameters")
  }
  # Which observations is y exceed the threshold u?
  isExc <- y > u
  if (x == 0) {
    print("HERE")
    nhpploglik <- -isExc * (log(b) + log(s) + (y - m) / s) -
      exp(-(u - m) / s) / nb
  } else {
    log1pzwx <- log1p(zw) / x
    log1pzux <- log1p(zu) / x
    print(log1pzux)
    print(log1pzwx)
    nhpploglik <- -isExc * (log(b) + log(s) + (x + 1) * log1pzwx) -
      exp(-log1pzux) / nb
  }
  if (!individual) {
    nhpploglik <- sum(nhpploglik)
  }
  return(nhpploglik)
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

# ==================== Simulation from a GEV distribution =================== #

#' @keywords internal
#' @rdname evils-internal
#' @export
rGenExtremeValue <- function (n, loc = 0, scale = 1, shape = 0) {
  max_len <- ifelse(length(n) > 1, length(n), n)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  return(qGenExtremeValue(stats::runif(n), loc = loc, scale = scale,
                          shape = shape))
}

#' @keywords internal
#' @rdname evils-internal
qGenExtremeValue <- function (p, loc = 0, scale = 1, shape = 0,
                              lower.tail = TRUE, log.p = FALSE) {
  if (any(scale <= 0)) {
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
  xp <- -log(p)
  mult <- BoxCoxVec(x = xp, lambda = -shape)
  return(loc - scale * mult)
}
