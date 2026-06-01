# Generate Missing Mask for Tensor Simulation

Randomly masks a proportion of non-zero entries in a tensor to simulate
missing data. Returns the observed tensor, observation mask, and
metadata needed for evaluation of imputation quality.

## Usage

``` r
generate_missing_mask(X, missing_rate = 0.1, seed = NULL)
```

## Arguments

- X:

  Non-negative integer array of dimensions `n x p x w`.

- missing_rate:

  Numeric proportion of non-zero entries to mask. Default `0.1` (10%).

- seed:

  Optional integer random seed for reproducibility. Default `NULL`.

## Value

A named list with:

- X_obs:

  Observed tensor with masked entries set to `NA`.

- obs_mask:

  Logical array of same dimensions as `X`, where `TRUE` indicates
  observed entries.

- obs_indices:

  data.table of observed entry indices (V1, V2, V3).

- missing_nonzero_indices:

  data.table of masked non-zero entry indices (V1, V2, V3).

- true_values:

  Numeric vector of true values at masked positions.

- n_missing:

  Number of masked entries.

## Examples

``` r
if (FALSE) { # \dontrun{
X <- array(rpois(20 * 12 * 4, lambda = 1.5), dim = c(20, 12, 4))
mask_info <- generate_missing_mask(X, missing_rate = 0.1, seed = 42)
} # }
```
