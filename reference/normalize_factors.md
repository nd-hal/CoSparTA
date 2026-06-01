# Normalize and Sort Components of a CoSparTA fit

Extracts the factor matrices from a fitted
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
object, normalizes each column to unit Frobenius norm, computes a
per-component weight \\\lambda_k = \\\mathbf{l}\_k\\ \\\mathbf{f}\_k\\
\\\mathbf{w}\_k\\\\, and returns components sorted by \\\lambda\\
descending. This is the canonical form for comparing decompositions and
for downstream use with `project_tensor` and `reconstruct_tensor`.

Since
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
now produces normalized outputs directly in `_normed` fields (e.g.
`fit$res$ql$El_normed`), this function is primarily useful for legacy
code or custom normalization workflows applied to fit objects that lack
those fields.

## Usage

``` r
normalize_factors(fit)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  or
  [`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).

## Value

A named list with:

- El:

  Normalized and reordered observation loading matrix (`n x K`). Each
  column has unit Frobenius norm.

- Ef:

  Normalized and reordered time factor matrix (`p x K`). Each column has
  unit Frobenius norm.

- Ew:

  Normalized and reordered channel weight matrix (`w x K`). Each column
  has unit Frobenius norm.

- lambda:

  Numeric vector of length `K` giving the component weights \\\lambda_k
  = \\\mathbf{l}\_k\\ \\\mathbf{f}\_k\\ \\\mathbf{w}\_k\\\\, sorted
  descending.

- order:

  Integer vector of length `K` giving the permutation of original
  component indices in new (descending \\\lambda\\) order.

- gamma_list:

  A list of length `K`. Each element is the gamma coefficient vector for
  that component after reordering by descending lambda. `NULL` if no
  components have gamma estimates (i.e. the model was fit without
  covariates).

- lambda_order:

  Integer vector of length `K` giving the permutation used to reorder
  components by descending lambda. Useful for reordering other
  component-indexed objects (e.g. `Xcov_list`) to match the normalized
  factor ordering.

## See also

[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md),
[`project_tensor`](https://nd-hal.github.io/CoSparTA/reference/project_tensor.md),
[`reconstruct_tensor`](https://nd-hal.github.io/CoSparTA/reference/reconstruct_tensor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3, Xcov = Xcov)

norm_fit <- normalize_factors(fit)
norm_fit$lambda   # component weights, descending
norm_fit$order    # e.g. c(2, 1, 3) — original rank 2 is now rank 1
norm_fit$El       # normalized loading matrix (n x K)
} # }
```
