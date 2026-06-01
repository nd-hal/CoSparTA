# Reconstruct the Denoised Mean Tensor from a CoSparTA fit

Reconstructs the denoised Poisson mean tensor \\\hat{X}\[i,j,m\] =
\sum\_{k=1}^{K} L\_{ik} F\_{jk} W\_{mk}\\ from the raw (unnormalized)
posterior mean factor matrices stored in a fitted
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
object. The result is the best rank-`K` approximation to the observed
tensor under the fitted model.

Implementation uses a loop over `K` components. For each \\k\\ the outer
product \\\mathbf{l}\_k \otimes \mathbf{f}\_k \otimes \mathbf{w}\_k\\ is
added to the accumulator via `tcrossprod` and a sweep, avoiding
allocation of `K` separate `n x p x w` arrays.

## Usage

``` r
reconstruct_tensor(fit, normalized = TRUE)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  or
  [`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).

- normalized:

  Logical. If `TRUE` (default), reads from the `_normed` fields
  (`El_normed`, `Ef_normed`, `Ew_normed`) and scales each component by
  `lambda_normed[k]`. If `FALSE`, uses raw posterior means directly
  (current behavior).

## Value

A numeric array of dimensions `n x p x w` containing the reconstructed
denoised mean tensor \\\hat{X}\\.

## See also

[`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md),
[`project_tensor`](https://nd-hal.github.io/CoSparTA/reference/project_tensor.md),
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3)

X_hat <- reconstruct_tensor(fit)
dim(X_hat)          # same as dim(X)

# Mean squared reconstruction error
mean((X - X_hat)^2)
} # }
```
