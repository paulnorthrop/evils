% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evils-internal.R
\docType{data}
\name{evils-internal}
\alias{evils-internal}
\alias{EulersConstant}
\alias{AperysConstant}
\alias{pxi}
\alias{qxi}
\alias{gev22e0Fn}
\alias{gev22e0Constant}
\alias{gev22eFn}
\alias{gev33e0Fn}
\alias{gev33e0Constant}
\alias{gev33eFn}
\alias{gev12e0Fn}
\alias{gev12e0Constant}
\alias{gev12eFn}
\alias{gev13e0Fn}
\alias{gev13e0Constant}
\alias{gev13eFn}
\alias{gev23e0Fn}
\alias{gev23e0Constant}
\alias{gev23eFn}
\alias{gevExpInfoComp}
\alias{lagrangianInterpolation}
\alias{gpLoglikDirect}
\alias{gevLoglikDirect}
\alias{BC}
\title{Internal evils functions}
\format{
An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.
}
\usage{
EulersConstant

AperysConstant

pxi(xi)

qxi(xi)

gev22e0Fn()

gev22e0Constant

gev22eFn(xi)

gev33e0Fn()

gev33e0Constant

gev33eFn(xi)

gev12e0Fn()

gev12e0Constant

gev12eFn(xi)

gev13e0Fn()

gev13e0Constant

gev13eFn(xi)

gev23e0Fn()

gev23e0Constant

gev23eFn(xi)

gevExpInfoComp(fun, fun0, xi, eps = 0.003)

lagrangianInterpolation(x0, y0)

gpLoglikDirect(pars, excesses, individual = FALSE)

gevLoglikDirect(pars, maxima, individual = FALSE)

BC(x, lambda, eps = 1e-06)
}
\description{
Internal evils functions
}
\details{
These functions are not intended to be called by the user.
}
\keyword{internal}
