# CoSparTA

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An R package for **contextual empirical Bayes tensor decomposition** of sparse count tensors. CoSparTA fits a Poisson CP decomposition with factor-specific spike-and-slab priors and incorporates observation-level covariates directly into the generative model, enabling simultaneous factor recovery, covariate effect estimation, and posterior uncertainty quantification.

## Installation

Two dependencies are hosted on GitHub and must be installed first:

```r
devtools::install_github("DongyueXie/ebpm")
devtools::install_github("DongyueXie/smashrgen")
devtools::install_github("xzhang0407/CoSparTA")
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
df             <- readRDS("data/demo_df.rds")
cov_df         <- read.csv("data/demo_covariates.csv")
channel_groups <- readRDS("data/demo_channel_groups.rds")

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
# use_virtualenv("cxtebtd_env", required = TRUE)
init_vals <- tryCatch(
  init_cpapr(X, K = 4, virtualenv = "cxtebtd_env"),
  error = function(e) "random_gamma"
)

fit <- CoSparTA(
  X                    = X,
  K                    = 4,
  Xcov                 = Xcov_mat,
  init                 = init_vals,
  maxiter              = 20,
  convergence_criteria = "ELBO",
  tol                  = 1e-6,
  verbose              = TRUE
)
# Key outputs:
# fit$res$ql$El -- n x K posterior mean loadings (U1)
# fit$res$qf$Ef -- p x K time factors            (U2)
# fit$res$qw$Ew -- w x K channel weights         (U3)
# fit$res$gl    -- list of K covariate coefficient vectors
```

Rank-specific covariate sets and missing data are also supported:

```r
# Rank-specific covariates
fit_rankspec <- CoSparTA(X, K = 4,
  Xcov = list(Xcov_mat, NULL, Xcov_mat[,1,drop=FALSE], Xcov_mat[,2,drop=FALSE]),
  init = init_vals, convergence_criteria = "ELBO")

# Missing data
mask     <- generate_missing_mask(X, missing_rate = 0.10)
fit_miss <- CoSparTA_missing(X, K = 4, Xcov = Xcov_mat,
                             obs_mask = mask$obs_mask,
                             convergence_criteria = "ELBO")
```

### Stage 3: Uncertainty quantification

```r
# 95% posterior credible intervals for channel mode (U3)
ci_w <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = "W")

# 95% posterior credible intervals for observation mode (U1)
ci_l <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = "L")

# Delta-method CIs for covariate coefficients
ci_gamma <- get_gamma_ci(fit, method = "delta", level = 0.95)
print(ci_gamma[[1]])   # factor 1 results
```

### Stage 4: Covariate screening

```r
# Fit covariate-free model, then screen covariates
fit_unsup  <- CoSparTA(X, K = 4, Xcov = NULL, init = init_vals,
                       convergence_criteria = "ELBO")
nf_unsup   <- normalize_factors(fit_unsup)
sel        <- select_covariates(K = 4,
                                covariate_data = as.data.frame(Xcov_mat),
                                El             = nf_unsup$El)
sel$selected   # selected covariate indices per factor
```

### Stage 5: Post-processing

```r
# Normalize to unit norm, sort by weight lambda
nf <- normalize_factors(fit)

# Project new observations onto the factor space (no refit)
X_new <- X[991:1000, , ]   # held-out sessions for illustration
F_new <- project_tensor(X_new, fit)

# Reconstruct the denoised Poisson mean tensor
X_hat <- reconstruct_tensor(fit)
```

### Stage 6: Visualization

```r
# Faceted line plot of K time factors
plot_time_factors(Ef = nf$Ef, time_labels = 1:100)

# Faceted bar plot of K channel factors with category grouping
plot_channel_factors(Ew             = nf$Ew,
                     channel_names  = website_names,
                     channel_groups = channel_groups)
```

## Questions

[Contact redacted for review]
