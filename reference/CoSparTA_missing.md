# Covariate-aware Empirical Bayes Tensor Decomposition with Missing Data

Extends
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md) to
handle tensors with missing entries. Missing entries are excluded from
the likelihood and the ELBO, and per-observation scales are computed
using only observed entries for each mode. The generative model and
priors are identical to
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).

## Usage

``` r
CoSparTA_missing(
  X,
  K,
  Xcov = NULL,
  obs_mask = NULL,
  lib_size = NULL,
  init = "random_gamma",
  maxiter = 100,
  maxiter_init = 100,
  tol = 1e-06,
  compute_elbo_final = FALSE,
  ebpm.fn = c(ebpm_point_gamma_multiplier_covariates, ebps_with_uq,
    ebpm_point_gamma_with_uq),
  fix_L = FALSE,
  fix_F = FALSE,
  fix_W = FALSE,
  smooth_F = TRUE,
  printevery = 10,
  verbose = TRUE,
  convergence_criteria = "ELBO",
  n_stable = 3L,
  U1_true = NULL,
  U2_true = NULL,
  U3_true = NULL
)
```

## Arguments

- X:

  A 3-dimensional non-negative integer array of dimensions `n x p x w`.
  May contain `NA` for missing entries.

- K:

  Integer. Number of components (CP rank).

- Xcov:

  Covariate input for the observation mode. Can be: (1) `NULL` for fully
  unsupervised decomposition; (2) a numeric matrix of dimension `n x q`,
  applied identically to all K components; or (3) a list of length `K`,
  where each element is either a numeric covariate matrix (dimensions
  `n x q_k`, potentially different numbers of covariates per component)
  or `NULL` for unsupervised components. Default `NULL`.

- obs_mask:

  Optional logical array of same dimensions as `X`, where `TRUE`
  indicates observed entries. If `NULL`, inferred from `is.na(X)`.
  Default `NULL`.

- lib_size:

  Numeric vector of length `n`. Default `NULL`.

- init:

  Initialization method. Default `'random_gamma'`.

- maxiter:

  Maximum EM iterations. Default `100`.

- maxiter_init:

  Maximum initialization iterations. Default `100`.

- tol:

  Convergence tolerance. Default `1e-6`.

- compute_elbo_final:

  Logical. If `TRUE`, computes the ELBO after the final iteration even
  when `convergence_criteria != 'ELBO'`. Default `FALSE`.

- ebpm.fn:

  List of three EBPM functions for L, F, W modes respectively. Default
  uses `ebpm_point_gamma_multiplier_covariates` for L, `ebps_with_uq`
  for F, and `ebpm_point_gamma_with_uq` for W.

- fix_L:

  Logical. Fix L at initialization. Default `FALSE`.

- fix_F:

  Logical. Fix F at initialization. Default `FALSE`.

- fix_W:

  Logical. Fix W at initialization. Default `FALSE`.

- smooth_F:

  Logical. Apply smoothing to F. Default `TRUE`.

- printevery:

  Print progress every this many iterations. Default `10`.

- verbose:

  Logical. Print progress. Default `TRUE`.

- convergence_criteria:

  Convergence criterion: `'ELBO'` (default) or `'mKLabs'`.

- n_stable:

  Integer. Number of consecutive iterations with factor change below
  `tol` required before declaring convergence. Only used when
  `convergence_criteria = 'factor_change'`. Default `3`.

- U1_true:

  Optional true L matrix for simulation evaluation.

- U2_true:

  Optional true F matrix for simulation evaluation.

- U3_true:

  Optional true W matrix for simulation evaluation.

## Value

Same structure as
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
(including all `_normed` fields) with one additional field:

- obs_structure:

  Internal observation structure used during fitting, required by
  [`evaluate_missing_prediction`](https://nd-hal.github.io/CoSparTA/reference/evaluate_missing_prediction.md).

## See also

[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md),
[`generate_missing_mask`](https://nd-hal.github.io/CoSparTA/reference/generate_missing_mask.md),
[`evaluate_missing_prediction`](https://nd-hal.github.io/CoSparTA/reference/evaluate_missing_prediction.md)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
X <- array(rpois(20 * 12 * 4, lambda = 1.5), dim = c(20, 12, 4))
mask_info <- generate_missing_mask(X, missing_rate = 0.1, seed = 42)

# With covariates
Xcov <- matrix(rnorm(20 * 2), nrow = 20)
fit <- CoSparTA_missing(mask_info$X_obs, K = 3, Xcov = Xcov,
                        obs_mask = mask_info$obs_mask)

# Without covariates
fit0 <- CoSparTA_missing(mask_info$X_obs, K = 3,
                         obs_mask = mask_info$obs_mask)
} # }
```
