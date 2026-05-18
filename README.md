# CxtEBTD

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An R package for **contextual empirical Bayes tensor decomposition** of sparse count tensors. CxtEBTD fits a Poisson CP decomposition with factor-specific spike-and-slab priors and incorporates observation-level covariates directly into the generative model, enabling simultaneous factor recovery, covariate effect estimation, and posterior uncertainty quantification.

## Installation

Two dependencies are hosted on GitHub and must be installed first:

```r
devtools::install_github("DongyueXie/ebpm")
devtools::install_github("DongyueXie/smashrgen")
devtools::install_github("xzhang0407/CxtEBTD")
```

## Quick Start

### Simulate a sparse count tensor

```r
library(CxtEBTD)

set.seed(1)
n <- 500; p <- 100; w <- 50; K <- 4
Xcov <- cbind(rnorm(n), rnorm(n))

gamma_true <- list(
  c( 0.6, -0.3),
  c(-0.5,  0.4),
  c( 0.3,  0.5),
  c(-0.4, -0.2)
)

sim  <- simulate_tensor(n=n, p=p, w=w, K=K,
                        Xcov=Xcov, gamma_true=gamma_true,
                        sparsity=50, seed=42)
X    <- sim$X
Xcov <- sim$Xcov
```

### Fit the model

```r
fit <- CxtEBTD(X, K=K, Xcov=Xcov,
               maxiter=50,
               convergence_criteria="factor_change")
```

### Normalize and visualize factors

```r
nf <- normalize_factors(fit)

plot_time_factors(Ef=nf$Ef)
plot_channel_factors(Ew=nf$Ew)
```

### Covariate coefficient inference

```r
# Delta-method confidence intervals for gamma
gamma_ci <- get_gamma_ci(fit, method="delta", level=0.95)
gamma_ci[[1]]  # estimates, SEs, CIs, p-values for factor 1
```

### Posterior uncertainty quantification

```r
# 95% credible intervals for observation-mode loadings
ci_L <- get_credible_interval(fit, mode="L", level=0.95)
```

### Missing data

```r
mask     <- generate_missing_mask(X, missing_rate=0.10, seed=42)
fit_miss <- CxtEBTD_missing(X=mask$X_obs, K=K, Xcov=Xcov,
                             obs_mask=mask$obs_mask,
                             maxiter=50,
                             convergence_criteria="factor_change")
```

## Questions

[Contact redacted for review]
