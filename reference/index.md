# Package index

## Decomposition

Main model fitting functions

- [`CoSparTA()`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  : Covariate-aware Empirical Bayes Tensor Decomposition
- [`CoSparTA_missing()`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md)
  : Covariate-aware Empirical Bayes Tensor Decomposition with Missing
  Data

## Core Prior

Empirical Bayes prior solvers

- [`ebpm_point_gamma_multiplier_covariates()`](https://nd-hal.github.io/CoSparTA/reference/ebpm_point_gamma_multiplier_covariates.md)
  : Empirical Bayes Poisson Mean Estimation with Covariate-Dependent
  Rates
- [`ebpm_point_gamma_with_uq()`](https://nd-hal.github.io/CoSparTA/reference/ebpm_point_gamma_with_uq.md)
  : Empirical Bayes Point-Gamma with Posterior Variance and PIP
- [`ebps_with_uq()`](https://nd-hal.github.io/CoSparTA/reference/ebps_with_uq.md)
  : Empirical Bayes Poisson Smoothing with Posterior Variance

## Inference

Posterior uncertainty quantification

- [`get_pip()`](https://nd-hal.github.io/CoSparTA/reference/get_pip.md)
  : Extract Posterior Inclusion Probabilities from a CoSparTA fit
- [`get_significant_patterns()`](https://nd-hal.github.io/CoSparTA/reference/get_significant_patterns.md)
  : Identify Significant Channels and Time Points per Factor via lFDR
  Control
- [`get_posterior_quantile()`](https://nd-hal.github.io/CoSparTA/reference/get_posterior_quantile.md)
  : Compute Exact Posterior Quantiles from a CoSparTA fit
- [`get_gamma_ci()`](https://nd-hal.github.io/CoSparTA/reference/get_gamma_ci.md)
  : Confidence Intervals for Covariate Coefficients

## Post-processing

Factor extraction and manipulation

- [`get_loadings()`](https://nd-hal.github.io/CoSparTA/reference/get_loadings.md)
  : Retrieve Factor Matrices or Weights from a CoSparTA Fit
- [`normalize_factors()`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md)
  : Normalize and Sort Components of a CoSparTA fit
- [`project_tensor()`](https://nd-hal.github.io/CoSparTA/reference/project_tensor.md)
  : Project a New Tensor onto the Learned Factor Space
- [`reconstruct_tensor()`](https://nd-hal.github.io/CoSparTA/reference/reconstruct_tensor.md)
  : Reconstruct the Denoised Mean Tensor from a CoSparTA fit
- [`match_factors()`](https://nd-hal.github.io/CoSparTA/reference/match_factors.md)
  : Match and Align Estimated Factors to Reference Factors
- [`select_covariates()`](https://nd-hal.github.io/CoSparTA/reference/select_covariates.md)
  : Two-Step Covariate Selection via Unsupervised Decomposition

## Data Utilities

Tensor construction and simulation

- [`build_tensor()`](https://nd-hal.github.io/CoSparTA/reference/build_tensor.md)
  : Build a 3-Way Tensor from Long-Format Data
- [`simulate_tensor()`](https://nd-hal.github.io/CoSparTA/reference/simulate_tensor.md)
  : Simulate a Count Tensor with Known Ground Truth
- [`init_cpapr()`](https://nd-hal.github.io/CoSparTA/reference/init_cpapr.md)
  : Initialize Factor Matrices via CP-APR (Poisson Tensor Factorization)

## Missing Data

Missing data utilities

- [`generate_missing_mask()`](https://nd-hal.github.io/CoSparTA/reference/generate_missing_mask.md)
  : Generate Missing Mask for Tensor Simulation
- [`evaluate_missing_prediction()`](https://nd-hal.github.io/CoSparTA/reference/evaluate_missing_prediction.md)
  : Evaluate Prediction Quality on Held-out Missing Entries

## Visualization

Factor plotting functions

- [`plot_time_factors()`](https://nd-hal.github.io/CoSparTA/reference/plot_time_factors.md)
  : Plot Time-Mode Factors
- [`plot_channel_factors()`](https://nd-hal.github.io/CoSparTA/reference/plot_channel_factors.md)
  : Plot Channel-Mode Factors

## Data

Built-in datasets

- [`clickstream_synth_cov`](https://nd-hal.github.io/CoSparTA/reference/clickstream_synth_cov.md)
  : Synthetic clickstream covariate matrix
- [`demo_covariates`](https://nd-hal.github.io/CoSparTA/reference/demo_covariates.md)
  : Demo covariate matrix
