
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evils

[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/paulnorthrop/evils?branch=main&svg=true)](https://ci.appveyor.com/project/paulnorthrop/evils)
[![R-CMD-check](https://github.com/paulnorthrop/evils/workflows/R-CMD-check/badge.svg)](https://github.com/paulnorthrop/evils/actions)
[![Coverage
Status](https://codecov.io/github/paulnorthrop/evils/coverage.svg?branch=main)](https://app.codecov.io/github/paulnorthrop/evils?branch=main)

## Evaluate Extreme Value Likelihoods Safely

Provides functions to calculate contributions to the log-likelihood
function, score function and observed information matrix from
individuals observations for the GEV and GP extreme value models. If the
shape parameter is close enough to zero, direct evaluation of some of
the quantities involved can be unreliable. In these cases, the
problematic terms involve log(1+x)/x and/or the first and second
derivatives of thisfunction. Functions are provided to evaluate these
quantities reliably when x is close to 0.

### An example

Add example

### Installation

### Vignette
