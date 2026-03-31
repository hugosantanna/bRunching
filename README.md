
# bRunching <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

An R package implementing the Caetano, Caetano, and Nielsen (CCN) bunching correction estimator for regression models where the treatment variable is censored at zero.

When many observations "bunch" at zero, naive OLS produces biased estimates. **bRunching** corrects this by estimating the censored expectation of the latent treatment and including it as a correction term in the regression.

## Installation

```r
# install.packages("remotes")
remotes::install_github("hugosantanna/bRunching")
```

## Usage

```r
library(bRunching)

# Simulate data with censoring at zero
set.seed(42)
n <- 5000
z <- runif(n, 0, 10)
x_star <- -0.5 + 0.3 * z + rnorm(n)
x <- pmax(x_star, 0)                    # censored at zero
y <- 1 + 0.8 * x_star + rnorm(n)

df <- data.frame(y = y, x = x, z = z)

# Naive OLS (biased)
naive <- ccn("y", "x", "z", df, method = "naive")

# Corrected estimates
corrected <- ccn("y", "x", "z", df, method = "symmetric", n_bins = 10)
print(corrected)

# Discontinuity plot at the bunching point
discontinuity_plot("y", "x", df)

# Or directly from a fitted object
plot(corrected)
```

## Methods

| Method | Description |
|--------|-------------|
| `"naive"` | Uncorrected OLS with instrument-bin fixed effects |
| `"uniform"` | Assumes uniform distribution below zero |
| `"tobit"` | Pooled Tobit with instrument-bin dummies |
| `"het_tobit"` | Separate Tobit per instrument bin |
| `"symmetric"` | Symmetry-based estimator with Tobit fallback for high-censoring bins |

All methods support heterogeneous treatment effects (`het_beta`) and heterogeneous correction terms (`het_delta`).

## References

Caetano, C., Caetano, G., and Nielsen, E. (2024). "Is Video Watching Bad for Kids? The Effect of Video Watching on Children's Skills."

## Authors

Hugo Sant'Anna and Debora Mazetto
