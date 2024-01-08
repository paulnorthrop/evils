#' evils: Evaluate Extreme Value Likelihoods Safely
#'
#' Provides functions to calculate contributions to the
#' log-likelihood function, score function and observed information matrix from
#' individuals observations for the GEV and GP extreme value models. If the
#' shape parameter is close enough to zero, direct evaluation of some of the
#' quantities involved can be unreliable. In these cases, the problematic
#' terms involve log(1+x)/x and/or the first and second derivatives of this
#' function. Functions are provided to evaluate these quantities reliably when
#' x is close to 0.
#'
#' @details Add details
#' @docType package
"_PACKAGE"
