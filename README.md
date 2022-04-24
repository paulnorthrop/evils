
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evils

## Evaluate Extreme Value Likelihoods Safely

The **evils** package provides functions to calculate the
log-likelihood, score and observed information for extreme value models
in cases where the shape parameter is very close to zero. In these
cases, the quantities of interest are expressed as a series expansion,
whose value is approximated using the
[sumR](https://cran.r-project.org/package=sumR) package.

### An example

Consider a single observation *y* sampled from a generalised Pareto (GP)
distribution with scale parameter *σ* and shape parameter *ξ*. Let
*i*<sub>22</sub> be the contribution of this observation to the (2, 2)
element of the observed information matrix, that is, the negated value
of the second derivative of the log-likelihood with respect to *ξ*.

*i*<sub>22</sub> is continuous at *ξ* = 0 and can, in principle, be
evaluated for values of *ξ* that are arbitrarily close to 0. However,
direct evaluation of *i*<sub>22</sub> is unreliable if *ξ* is very close
to zero. To avoid this problem we express *i*<sub>22</sub> as a series
expansion in *ξ*, whose value is approximated using the [sumR
package](https://cran.r-project.org/package=sumR). This is illustrated
using the following code, which simulates a random sample from a GP
distribution.

``` r
library(evils)

### Simulate some data
set.seed(15042022)
y <- rGenPareto(100, 0, 1, 0)

# Approximation using sumR::infinitesum()
gpInfo(pars = c(1, 1e-7), excesses = y)[2, 2]
#> [1] 617.6421
gpInfo(pars = c(1, -1e-7), excesses = y)[2, 2]
#> [1] 617.6443
gpInfo(pars = c(1, 0), excesses = y)[2, 2]
#> [1] 617.6432

# Direct calculation, starting the break down
gpInfoDirect(pars = c(1, 1e-7), excesses = y)[2, 2]
#> [1] 615.4318
gpInfoDirect(pars = c(1, -1e-7), excesses = y)[2, 2]
#> [1] 614.5825

# Direct calculation, getting worse
gpInfoDirect(pars = c(1, 1e-10), excesses = y)[2, 2]
#> [1] 775407.3
gpInfoDirect(pars = c(1, -1e-10), excesses = y)[2, 2]
#> [1] 883236.4
```

Similar problems occur when evaluating the log-likelihood and score
function in the generalised Pareto case and for other extreme value
models.

### Installation

To get the current released version from CRAN:

``` r
install.packages("evils")
```

### Vignette

See `vignette("evils-vignette", package = "evils")` for an overview of
the package.
