% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gpExpInfo.R
\name{gpExpInfo}
\alias{gpExpInfo}
\alias{gp11e}
\alias{gp22e}
\alias{gp12e}
\title{GP Distribution Expected Information}
\usage{
gp11e(scale, shape)

gp22e(shape)

gp12e(scale, shape)

gpExpInfo(scale, shape)
}
\arguments{
\item{scale, shape}{Numeric vectors. Respective values of the GP parameters
\eqn{\sigma} and \eqn{\xi}. For \code{gpExpInfo}, \code{scale} and \code{shape} must have
length 1.}
}
\value{
\code{gpExpInfo} returns a 3 by 3 numeric matrix with row and column
names \verb{loc, scale, shape}. The other functions return a numeric vector of
length equal to the maximum of the lengths of the arguments, excluding
\code{eps}.
}
\description{
Calculates the expected information matrix for the GP distribution.
}
\details{
\code{gpExpInfo} calculates, for single pair of values
\eqn{(\sigma, \xi) = } \verb{(scale, shape)}, the expected information matrix for a
single observation from a GP distribution with distribution function
\deqn{F(x) = P(X \leq x) = 1 - \left[ 1+\xi x / \sigma \right]_+^{-1/\xi},}
where \eqn{x_+ = \max(x, 0)}.
The GP expected information is defined only for \eqn{\xi > -0.5}.

The other functions are vectorised and calculate the individual
contributions to the expected information matrix. For example, \code{gp11e}
calculates the expectation \eqn{i_{\sigma\sigma}} of the negated second
derivative of the GP log-density with respect to \eqn{\sigma}, that is,
each \code{1} indicates one derivative with respect to \eqn{\sigma}. Similarly,
\code{2} denotes one derivative with respect to \eqn{\xi} so that, for example,
\code{gp12e} calculates the expectation \eqn{i_{\sigma\xi}} of the negated GP
log-density after one taking one derivative with respect to \eqn{\sigma}
and one derivative with respect to \eqn{\xi}. Note that \eqn{i_{\xi\xi}},
calculated using \code{gp22e}, depends only on \eqn{\xi}.
}
\examples{
# Expected information matrices for ...
# ... scale = 2 and shape = -0.4
gpExpInfo(2, -0.4)
# ... scale = 3 and shape = 0.001
gpExpInfo(3, 0.001)
# ... scale = 3 and shape = 0
gpExpInfo(3, 0)
# ... scale = 1 and shape = 0.1
gpExpInfo(1, 0.1)

# The individual components of the latter matrix
gp11e(1, 0.1)
gp12e(1, 0.1)
gp22e(0.1)
}
