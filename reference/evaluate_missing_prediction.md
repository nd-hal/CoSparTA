# Evaluate Prediction Quality on Held-out Missing Entries

Computes prediction metrics for held-out entries after fitting
[`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).
Reconstructs the predicted Poisson rate at masked positions using the
fitted factor matrices and compares against the true values.

## Usage

``` r
evaluate_missing_prediction(fit, missing_info)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).

- missing_info:

  Output from
  [`generate_missing_mask`](https://nd-hal.github.io/CoSparTA/reference/generate_missing_mask.md).

## Value

A named list with:

- rmse:

  Root mean squared error.

- mae:

  Mean absolute error.

- correlation:

  Pearson correlation between predicted and true values.

- deviance:

  Poisson deviance.

- predicted:

  Numeric vector of predicted values at masked positions.

- true_values:

  Numeric vector of true values at masked positions.

## Examples

``` r
if (FALSE) { # \dontrun{
X <- array(rpois(20 * 12 * 4, lambda = 1.5), dim = c(20, 12, 4))
mask_info <- generate_missing_mask(X, missing_rate = 0.1, seed = 42)
fit <- CoSparTA_missing(mask_info$X_obs, K = 3,
                        obs_mask = mask_info$obs_mask)
metrics <- evaluate_missing_prediction(fit, mask_info)
} # }
```
