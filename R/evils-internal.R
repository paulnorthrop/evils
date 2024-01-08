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
gev22e0Fn <- function() {
  return(pi ^ 2 / 6 + (1 - EulersConstant) ^ 2)
}

#' @keywords internal
#' @rdname evils-internal
gev22e0Constant <- gev22e0Fn()

#' @keywords internal
#' @rdname evils-internal
gev22eFn <- function(xi) {
  return((1 - 2 * gamma(2 + xi) + pxi(xi)) / (xi ^ 2))
}

# --------------------------------- (xi, xi) -------------------------------- #

#' @keywords internal
#' @rdname evils-internal
gev33e0Fn <- function() {
  val <- pi ^ 2 / 6 - pi ^ 2 * EulersConstant / 2 + EulersConstant ^ 2 -
    EulersConstant ^ 3 - 2 * AperysConstant +
    2 * EulersConstant * AperysConstant + pi ^ 2 * EulersConstant ^ 2 / 4 +
    EulersConstant ^ 4 / 4 + 3 * pi ^ 4 / 80
  return(val)
}

#' @keywords internal
#' @rdname evils-internal
gev33e0Constant <- gev33e0Fn()

#' @keywords internal
#' @rdname evils-internal
gev33eFn <- function(xi) {
  val <- (pi ^ 2 / 6 + (1 - EulersConstant + 1 / xi) ^ 2 - 2 * qxi(xi) / xi +
            pxi(xi) / xi ^ 2) / (xi ^ 2)
  return(val)
}

# -------------------------------- (mu, sigma) ------------------------------ #

#' @keywords internal
#' @rdname evils-internal
gev12e0Fn <- function() {
  return(EulersConstant - 1)
}

#' @keywords internal
#' @rdname evils-internal
gev12e0Constant <- gev12e0Fn()

#' @keywords internal
#' @rdname evils-internal
gev12eFn <- function(xi) {
  val <- (gamma(2 + xi) - pxi(xi)) / xi
  return(val)
}

# --------------------------------- (mu, xi) -------------------------------- #

#' @keywords internal
#' @rdname evils-internal
gev13e0Fn <- function() {
  return((pi ^ 2 / 6 + EulersConstant ^ 2 - 2 * EulersConstant) / 2)
}

#' @keywords internal
#' @rdname evils-internal
gev13e0Constant <- gev13e0Fn()

#' @keywords internal
#' @rdname evils-internal
gev13eFn <- function(xi) {
  val <- (pxi(xi) / xi - qxi(xi)) / xi
  return(val)
}

# ------------------------------- (sigma, xi) ------------------------------- #

#' @keywords internal
#' @rdname evils-internal
gev23e0Fn <- function() {
  val <- (4 * EulersConstant + 4 * AperysConstant + pi ^ 2 * EulersConstant +
     2 * EulersConstant ^ 3 - pi ^ 2 - 6 * EulersConstant ^ 2) / 4
  return(val)
}

#' @keywords internal
#' @rdname evils-internal
gev23e0Constant <- gev23e0Fn()

#' @keywords internal
#' @rdname evils-internal
gev23eFn <- function(xi) {
  val <- -(1 - EulersConstant + (1 - gamma(2 + xi)) / xi - qxi(xi) +
             pxi(xi) / xi) / (xi ^ 2)
  return(val)
}

# -- Calculate a component using a quadratic approximation if xi is near 0 -- #

#' @keywords internal
#' @rdname evils-internal
gevExpInfoComp <- function(fun, fun0, xi, eps = 3e-3) {
  eps <- abs(eps)
  val <- xi
  if (any(xiNearZero <- abs(xi) < eps)) {
    aa <- fun0
    yp <- fun(eps)
    ym <- fun(-eps)
    ff <- lagrangianInterpolation(c(-eps, 0, eps), c(ym, aa, yp))
    val[xiNearZero] <- ff(xi[xiNearZero])
  }
  val[!xiNearZero] <- fun(xi[!xiNearZero])
  return(val)
}

#' @keywords internal
#' @rdname evils-internal
lagrangianInterpolation <- function(x0, y0) {
  f <- function(x) {
    sum(y0 * sapply(seq_along(x0), \(j) {
      prod(x - x0[-j])/prod(x0[j] - x0[-j])
    }))
  }
  return(Vectorize(f, "x"))
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

#' @keywords internal
#' @rdname evils-internal
BC <- function(x, lambda, eps = 1e-6) {
  if (any(x < 0, na.rm = TRUE)) {
    stop("Invalid x: x must be non-negative")
  }
  eps <- abs(eps)
  # Recycle the vector input q, loc, scale and shape, if necessary
  maxLen <- max(length(x), length(lambda))
  x <- rep_len(x, maxLen)
  lambda <- rep_len(lambda, maxLen)
  # If either x or lambda is NA then return NA
  if (any(nas <- is.na(x) | is.na(lambda))) {
    x[nas] <- NA
  }
  # If lambda = 0 then use log(x)
  if (any(lambdaEq0 <- lambda == 0 & !nas)) {
    x[lambdaEq0] <- log(x[lambdaEq0])
  }
  # If abs(lambda) > eps or lambda = NA then use the usual formula
  if (any(large <- !nas & abs(lambda) >= eps & !lambdaEq0)) {
    x[large] <- (x[large] ^ lambda[large] - 1) / lambda[large]
  }
  # Indicators of lambda < 0 and lambda > 0
  neg <- !lambdaEq0 & !large & !nas & lambda < 0
  pos <- !lambdaEq0 & !large & !nas & lambda > 0
  # Indicators of being Inf or 0
  xInf <- is.infinite(x)
  xZero <- x == 0
  # Calculations for combinations of these indicators
  if (any(xInfNeg <- xInf & neg)) {
    x[xInfNeg] <- -1 / lambda[xInfNeg]
  }
  if (any(xInfPos <- xInf & pos)) {
    x[xInfPos] <- Inf
  }
  if (any(xZeroNeg <- xZero & neg)) {
    x[xZeroNeg] <- -Inf
  }
  if (any(xZeroPos <- xZero & pos)) {
    x[xZeroPos] <- -1 / lambda[xZeroPos]
  }
  # Use Taylor series expansion for other cases
  if (any(rest <-  !lambdaEq0 & !large & !xInf & !xZero & !nas)) {
    logxlam <- log(x[rest]) * lambda[rest]
    x[rest] <- log(x[rest]) * (1 + logxlam / 2 + logxlam ^ 2 / 6)
  }
  return(x)
}
