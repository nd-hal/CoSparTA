# Initialize Factor Matrices via CP-APR (Poisson Tensor Factorization)

Computes non-negative CP-APR (Alternating Poisson Regression) factor
matrices to use as warm-start initialization for
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).
Calls `pyCP_APR` via `reticulate`, which requires Python 3.9+ with
`pyCP_APR` and `numpy` installed.

## Usage

``` r
init_cpapr(
  X,
  K,
  n_iters = 150,
  method = "torch",
  random_state = 42,
  virtualenv = "ebtd1"
)
```

## Arguments

- X:

  A 3-dimensional non-negative integer array of dimensions `n x p x w`.

- K:

  Integer. Number of components (CP rank).

- n_iters:

  Integer. Maximum number of CP-APR iterations. Default `150`.

- method:

  Character string. Optimization backend for pyCP_APR: `'torch'`
  (default) or `'numpy'`.

- random_state:

  Integer. Random seed for reproducibility. Default `42`.

- virtualenv:

  Character string. Name of the Python virtual environment containing
  pyCP_APR. Default `'ebtd1'`.

## Value

A list of three matrices `list(L, F, W)` with dimensions `n x K`,
`p x K`, `w x K`, suitable for passing as the `init` argument to
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md) or
[`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).

## Details

Requires the `reticulate` R package and a Python environment with
`pyCP_APR` installed. To set up:


    pip install pyCP_APR numpy

The function adds a tiny constant (1e-7) to the last tensor slice
`X[, p, w]` to prevent pyCP_APR from failing on all-zero boundary
slices. pyCP_APR returns factor matrices with a zero-padded first row
(0-indexed), which is automatically removed before returning.

## See also

[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
X <- array(rpois(100 * 20 * 10, lambda = 1.5), dim = c(100, 20, 10))
init <- init_cpapr(X, K = 3)
fit <- CoSparTA(X, K = 3, init = init)
} # }
```
