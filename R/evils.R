#' evils: Evaluate Extreme Value Likelihoods Safely
#'
#' Provides functions to calculate the log-likelihood, score and observed
#' information for extreme value models in cases where the shape parameter is
#' very close to zero.  In these cases, the quantities of interest are
#' expressed as a series expansion, whose value is approximated using the
#' \code{\link[sumR]{sumR}} package.
#'
#' @details The main functions are
#' See \code{vignette("introduction-to-evils", package = "evils")} for an
#' overview of the package.
#' @docType package
#' @name evils
#' @importFrom graphics plot
NULL
