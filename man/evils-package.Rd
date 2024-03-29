% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evils-package.R
\docType{package}
\name{evils-package}
\alias{evils}
\alias{evils-package}
\title{evils: Evaluate Extreme Value Likelihoods Safely}
\description{
Provides functions to calculate contributions to the
log-likelihood function, score function and observed information matrix from
individuals observations for the GEV and GP extreme value models. If the
shape parameter is close enough to zero, direct evaluation of some of the
quantities involved can be unreliable. In these cases, the problematic
terms involve log(1+x)/x and/or the first and second derivatives of this
function. Functions are provided to evaluate these quantities reliably when
x is close to 0.
}
\details{
Add details
}
\seealso{
Useful links:
\itemize{
  \item \url{https://paulnorthrop.github.io/evils/}
  \item \url{https://github.com/paulnorthrop/evils}
  \item Report bugs at \url{https://github.com/paulnorthrop/evils/issues}
}

}
\author{
\strong{Maintainer}: Paul J. Northrop \email{p.northrop@ucl.ac.uk} [copyright holder]

}
