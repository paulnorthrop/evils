#' Generalised Extreme Value Expected Information
#'
#' Describe
#' @name gevExpectedInformation
NULL
## NULL

#' @rdname gevExpectedInformation
#' @export
imm <- function(sigma, xi) {
  return(pxi(xi) / sigma ^ 2)
}

#' @rdname gevExpectedInformation
#' @export
iss <- function(sigma, xi, eps = 1e-3) {
  val <- expInfoComponent(fun = issFn, fun0 = iss0Constant, xi = xi, eps = eps)
  return(val / sigma ^ 2)
}

#' @rdname gevExpectedInformation
#' @export
ixx <- function(sigma, xi, eps = 1e-3) {
  val <- expInfoComponent(fun = ixxFn, fun0 = ixx0Constant, xi = xi, eps = eps)
  return(val)
}

#' @rdname gevExpectedInformation
#' @export
ims <- function(sigma, xi, eps = 1e-3) {
  val <- expInfoComponent(fun = imsFn, fun0 = ims0Constant, xi = xi, eps = eps)
  return(val / sigma ^ 2)
}

#' @rdname gevExpectedInformation
#' @export
imx <- function(sigma, xi, eps = 1e-3) {
  val <- expInfoComponent(fun = imxFn, fun0 = imx0Constant, xi = xi, eps = eps)
  return(val / sigma)
}

#' @rdname gevExpectedInformation
#' @export
isx <- function(sigma, xi, eps = 1e-3) {
  val <- expInfoComponent(fun = isxFn, fun0 = isx0Constant, xi = xi, eps = eps)
  return(val / sigma)
}
