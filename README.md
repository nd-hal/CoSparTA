# CoSparTA

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An R package for **Covariate-Aware Sparsity-Adaptive Tensor Analysis** of sparse count tensors with an empirical Bayes approach. CoSparTA fits a Poisson CP decomposition with factor-specific spike-and-slab priors and incorporates observation-level covariates directly into the generative model, enabling simultaneous factor recovery, covariate effect estimation, and posterior uncertainty quantification.

## Installation

Two dependencies are hosted on GitHub and must be installed first:

```r
devtools::install_github("DongyueXie/ebpm")
devtools::install_github("DongyueXie/smashrgen")
devtools::install_github("nd-hal/CoSparTA")
```

## Repository Structure

```
CoSparTA/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── R/
│   ├── CoSparTA.R              # Main fitting function: CoSparTA
│   ├── missing.R               # Missing data: CoSparTA_missing,
│   │                           #   generate_missing_mask, evaluate_missing_prediction
│   ├── ebpm_covariates.R       # Covariate-dependent prior: ebpm_point_gamma_multiplier_covariates
│   ├── ebpm_wrappers.R         # Prior solver wrappers: ebpm_point_gamma_with_uq, ebps_with_uq
│   ├── inference.R             # Posterior inference: get_pip, get_posterior_quantile,
│   │                           #   get_significant_patterns, get_gamma_ci
│   ├── postprocessing.R        # Factor utilities: get_loadings, normalize_factors,
│   │                           #   match_factors, project_tensor, reconstruct_tensor,
│   │                           #   select_covariates, simulate_tensor, init_cpapr
│   ├── preprocessing.R         # Data utilities: build_tensor
│   ├── visualization.R         # Plotting: plot_time_factors, plot_channel_factors
│   ├── internals.R             # Internal CAVI update logic
│   ├── utils.R                 # Internal helpers and utility functions
│   └── RcppExports.R           # Auto-generated Rcpp bindings
├── src/
│   ├── calc_EZ_3d_cpp.cpp      # Sparse weighted aggregation for E[Z] (C++)
│   └── calc_qz_sparse_cpp.cpp  # Sparse softmax for factor responsibilities (C++)
├── data/
│   ├── demo_covariates.rda            # Demo covariate matrix
│   └── clickstream_synth_cov.rda      # Synthetic clickstream covariates
├── inst/extdata/
│   ├── demo_df.rds                    # Demo event log (long format)
│   ├── demo_channel_groups.rds        # Demo channel groupings for visualization
│   ├── clickstream_synth_tensor.rds   # Synthetic clickstream tensor (real-data demo)
│   ├── clickstream_channel_names.rds  # Channel names for clickstream demo
│   ├── demo_covariates.csv            # Demo covariate matrix (raw CSV)
│   └── clickstream_synth_cov.csv      # Synthetic clickstream covariates (raw CSV)
├── examples/
│   ├── demo.R                  # Full pipeline walkthrough (synthetic data)
│   └── demo_clickstream.R      # Real clickstream data example
└── tests/
    └── testthat/
        ├── test_core.R         # Tests: decomposition, missing data
        └── test_inference.R    # Tests: PIPs, credible intervals, gamma CI
```

## Quick Start

The following example walks through the complete analysis pipeline using a
synthetic clickstream dataset bundled with the package. The full runnable
script is at `examples/demo.R`.

### Stage 1: Build the tensor

Raw data arrive as a long-format data frame with one row per observed event.
`build_tensor` converts this to the required n x p x w integer count tensor.
Dimension labels can be supplied explicitly to ensure all time bins are
retained even if some are unobserved.

```r
library(CoSparTA)

# Load bundled demo data
df             <- readRDS(system.file("extdata", "demo_df.rds", package = "CoSparTA"))
data("demo_covariates", package = "CoSparTA"); cov_df <- demo_covariates
channel_groups <- readRDS(system.file("extdata", "demo_channel_groups.rds", package = "CoSparTA"))

Xcov_mat <- as.matrix(cov_df[, c("cov1", "cov2")])

# Reconstruct universe of labels
session_ids   <- sort(unique(df$session_id))
hour_labels   <- paste0("H", sprintf("%02d", 0:99))
channel_names <- c("Retail", "Social", "Review", "Deal", "Search")
website_names <- as.vector(
  outer(channel_names, paste0("_", sprintf("%02d", 1:10)), paste0))

# Build 1000 x 100 x 50 tensor
tensor_out <- build_tensor(
  data         = df,
  row          = "session_id",
  col          = "hour",
  slice        = "website",
  value        = "count",
  row_levels   = session_ids,
  col_levels   = hour_labels,
  slice_levels = website_names
)
X <- tensor_out$X
```

### Stage 2: Fit the model

`CoSparTA` is the main entry point. Supplying `Xcov` activates
covariate-dependent priors on the observation mode. An optional CP-APR
warm-start is available via `init_cpapr` when Python is accessible.

```r
# Optional CP-APR warm-start (requires Python 3.9+ with pyCP_APR)
# use_virtualenv("cosparta_env", required = TRUE)
init_vals <- tryCatch(
  init_cpapr(X, K = 4, virtualenv = "cosparta_env"),
  error = function(e) "random_gamma"
)

fit <- CoSparTA(
  X                    = X,
  K                    = 4,
  Xcov                 = Xcov_mat,
  init                 = init_vals,
  maxiter              = 20,
  convergence_criteria = "factor_change",
  tol                  = 1e-6,
  verbose              = TRUE
)
# Retrieve outputs via the unified getter:
U1     <- get_loadings(fit, "U1")      # n x K normalized loadings
U2     <- get_loadings(fit, "U2")      # p x K time factors
U3     <- get_loadings(fit, "U3")      # w x K channel weights
lambda <- get_loadings(fit, "weight")  # length-K component weights
gamma  <- get_loadings(fit, "gamma")   # list of K covariate coefficient vectors
```

Rank-specific covariate sets and missing data are also supported:

```r
# Rank-specific covariates
fit_rankspec <- CoSparTA(X, K = 4,
  Xcov = list(Xcov_mat, NULL, Xcov_mat[,1,drop=FALSE], Xcov_mat[,2,drop=FALSE]),
  init = init_vals, convergence_criteria = "factor_change")

# Missing data
mask     <- generate_missing_mask(X, missing_rate = 0.10)
fit_miss <- CoSparTA_missing(X, K = 4, Xcov = Xcov_mat,
                             obs_mask = mask$obs_mask,
                             convergence_criteria = "factor_change")
```

### Stage 3: Uncertainty quantification

```r
# 95% posterior credible intervals for channel mode (U3)
ci_w <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = "W")

# 95% posterior credible intervals for observation mode (U1)
ci_l <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = "L")

# Delta-method CIs for covariate coefficients
ci_gamma <- get_gamma_ci(fit, level = 0.95)
print(ci_gamma[[1]])   # factor 1 results
```

### Stage 4: Covariate screening

```r
# Fit covariate-free model, screen covariates against log-transformed active loadings
fit_unsup <- CoSparTA(X, K = 4, Xcov = NULL, init = init_vals,
                      convergence_criteria = "factor_change")
sel <- select_covariates(K = 4,
                         covariate_data = as.data.frame(Xcov_mat),
                         fit = fit_unsup)
sel$selected   # list of selected covariate names per factor
```

### Stage 5: Post-processing

```r
# Normalized factors are available directly from the fit object
U1 <- get_loadings(fit, "U1")   # already normalized and ordered

# Project new observations onto the factor space (no refit)
X_new <- X[991:1000, , ]
F_new <- project_tensor(X_new, fit)

# Reconstruct the denoised Poisson mean tensor
X_hat <- reconstruct_tensor(fit)
```

### Stage 6: Visualization

```r
# Faceted line plot of K time factors
plot_time_factors(fit, time_labels = 1:100)

# Faceted bar plot of K channel factors with category grouping
plot_channel_factors(fit,
                     channel_names  = website_names,
                     channel_groups = channel_groups)
```

## Questions?

Let me know if you have any requests, bugs, etc.
Email: xzhang38@nd.edu
