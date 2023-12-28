
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evils

## Evaluate Extreme Value Likelihoods Safely

The **evils** package provides functions to calculate the
log-likelihood, score and observed information for extreme value models
in cases where the shape parameter is very close to zero. In these
cases, …

### An example

Consider a single observation $y$ sampled from a generalised Pareto (GP)
distribution with scale parameter $\sigma$ and shape parameter $\xi$.
Let $i_{22}$ be the contribution of this observation to the (2, 2)
element of the observed information matrix, that is, the negated value
of the second derivative of the log-likelihood with respect to $\xi$.

$i_{22}$ is continuous at $\xi = 0$ and can, in principle, be evaluated
for values of $\xi$ that are arbitrarily close to 0. However, direct
evaluation of $i_{22}$ is unreliable if $\xi$ is very close to zero. To
avoid this problem we express $i_{22}$ as a series expansion in $\xi$,
whose value is approximated using the [sumR
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
gpInfo(pars = c(1, -1e-7), excesses = y)[2, 2]
gpInfo(pars = c(1, 0), excesses = y)[2, 2]

# Direct calculation, starting the break down
gpInfoDirect(pars = c(1, 1e-7), excesses = y)[2, 2]
gpInfoDirect(pars = c(1, -1e-7), excesses = y)[2, 2]

# Direct calculation, getting worse
gpInfoDirect(pars = c(1, 1e-10), excesses = y)[2, 2]
gpInfoDirect(pars = c(1, -1e-10), excesses = y)[2, 2]
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
