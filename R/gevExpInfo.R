#' GEV Distribution Expected Information
#'
#' Calculates the expected information matrix for the GEV distribution.
#'
#' @param sigma,xi Numeric vectors.
#'
#' * `gevExpInfo`: `sigma` and `xi` must have length 1.
#' * `imm`, `iss`, `ixx`, `ims`, `imx`, `isx`:
#' @param eps tolerance
#' @name gevExpInfo
NULL
## NULL

#' @rdname gevExpInfo
#' @export
imm <- function(sigma, xi) {
  return(pxi(xi) / sigma ^ 2)
}

#' @rdname gevExpInfo
#' @export
iss <- function(sigma, xi, eps = 3e-3) {
  val <- gevExpInfoComp(fun = issFn, fun0 = iss0Constant, xi = xi, eps = eps)
  return(val / sigma ^ 2)
}

#' @rdname gevExpInfo
#' @export
ixx <- function(sigma, xi, eps = 3e-3) {
  val <- gevExpInfoComp(fun = ixxFn, fun0 = ixx0Constant, xi = xi, eps = eps)
  return(val)
}

#' @rdname gevExpInfo
#' @export
ims <- function(sigma, xi, eps = 3e-3) {
  val <- gevExpInfoComp(fun = imsFn, fun0 = ims0Constant, xi = xi, eps = eps)
  return(val / sigma ^ 2)
}

#' @rdname gevExpInfo
#' @export
imx <- function(sigma, xi, eps = 3e-3) {
  val <- gevExpInfoComp(fun = imxFn, fun0 = imx0Constant, xi = xi, eps = eps)
  return(val / sigma)
}

#' @rdname gevExpInfo
#' @export
isx <- function(sigma, xi, eps = 3e-3) {
  val <- gevExpInfoComp(fun = isxFn, fun0 = isx0Constant, xi = xi, eps = eps)
  return(val / sigma)
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
  val[3, 3] <- ixx(sigma, xi, eps)
  val[2, 1] <- val[1, 2] <- ims(sigma, xi, eps)
  val[3, 1] <- val[1, 3] <- imx(sigma, xi, eps)
  val[3, 2] <- val[2, 3] <- isx(sigma, xi, eps)
  dimnames(val) <- list(c("mu", "sigma", "xi"), c("mu", "sigma", "xi"))
  return(val)
}
