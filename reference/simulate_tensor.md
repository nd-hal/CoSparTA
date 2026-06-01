# Simulate a Count Tensor with Known Ground Truth

Generates a synthetic Poisson count tensor from a CP structure with
known factor matrices, optional covariate effects, and spike-and-slab
sparsity. Useful for simulation studies evaluating decomposition
accuracy.

## Usage

``` r
simulate_tensor(
  n = 100,
  p = 20,
  w = 10,
  K = 3,
  Xcov = NULL,
  gamma_true = NULL,
  pi0 = 0.2,
  alpha_true = 3,
  beta_true = 2,
  weights = NULL,
  sparsity = 20,
  seed = 42
)
```

## Arguments

- n:

  Integer. Number of observations (mode 1). Default `100`.

- p:

  Integer. Number of time points (mode 2). Default `20`.

- w:

  Integer. Number of channels (mode 3). Default `10`.

- K:

  Integer. Number of components. Default `3`.

- Xcov:

  Optional numeric covariate matrix of dimension `n x q`. If `NULL`,
  covariates are generated automatically: an intercept, a binary
  covariate, and a continuous covariate. Default `NULL`.

- gamma_true:

  Optional list of length K, where each element is a numeric vector of
  covariate coefficients for that component. If `NULL` and `Xcov` is
  also `NULL`, random coefficients are generated. If `Xcov` is provided,
  must be supplied. Set to `FALSE` to generate an unsupervised tensor
  (no covariate effects). Default `NULL`.

- pi0:

  Numeric. Spike (structural zero) probability for mode-1 factors.
  Default `0.2`.

- alpha_true:

  Numeric. Gamma shape parameter for mode-1 slab. Default `3`.

- beta_true:

  Numeric. Gamma rate parameter for mode-1 slab. Default `2`.

- weights:

  Numeric vector of length K giving component weights. Default `NULL`
  generates decreasing weights.

- sparsity:

  Numeric scaling factor controlling overall tensor sparsity (higher =
  sparser). Default `20`.

- seed:

  Integer random seed. Default `42`.

## Value

A named list with:

- X:

  The observed Poisson count tensor, `n x p x w`.

- lambda_true:

  The true Poisson rate tensor, `n x p x w`.

- U1_true:

  True mode-1 factor matrix (normalized to unit norm), `n x K`.

- U2_true:

  True mode-2 factor matrix (normalized), `p x K`.

- U3_true:

  True mode-3 factor matrix (normalized), `w x K`.

- weights:

  Component weights, length K.

- Xcov:

  Covariate matrix used, `n x q`. `NULL` if unsupervised.

- gamma_true:

  List of true covariate coefficient vectors. `NULL` if unsupervised.

- sparsity_pct:

  Percentage of zeros in `X`.

## See also

[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md),
[`match_factors`](https://nd-hal.github.io/CoSparTA/reference/match_factors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Supervised simulation
sim <- simulate_tensor(n = 100, p = 20, w = 10, K = 3)
fit <- CoSparTA(sim$X, K = 3, Xcov = sim$Xcov)

# Unsupervised simulation
sim0 <- simulate_tensor(n = 100, p = 20, w = 10, K = 2, gamma_true = FALSE)
fit0 <- CoSparTA(sim0$X, K = 2)
} # }
```
