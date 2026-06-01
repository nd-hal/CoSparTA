# Extract Posterior Inclusion Probabilities from a CoSparTA fit

Extracts the posterior inclusion probability (PIP) matrices from a
fitted
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
object. PIPs give the probability that each factor element is truly
non-zero vs. noise, derived from the spike-and-slab prior. A PIP close
to 1 indicates strong evidence for a non-zero loading; a PIP close to 0
indicates the loading is likely noise.

## Usage

``` r
get_pip(fit, mode = "L", threshold = NULL, normalized = TRUE)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).

- mode:

  Character string specifying which mode to extract: `'L'` for
  observation loadings, `'F'` for time factors, `'W'` for channel
  weights. Default `'L'`.

- threshold:

  Numeric value in `[0, 1]`. If supplied, returns a logical matrix
  indicating which elements exceed the threshold. If `NULL`, returns the
  raw PIP matrix. Default `NULL`.

- normalized:

  Logical. If `TRUE` (default), reads from the `_normed` PIP fields
  (unit-norm columns, descending \\\lambda\\ order). If `FALSE`, reads
  from the raw PIP fields.

## Value

If `threshold = NULL`, a numeric matrix of PIPs with the same dimensions
as the corresponding factor matrix (`n x K` for L, `p x K` for F,
`w x K` for W). If `threshold` is supplied, a logical matrix of the same
dimensions where `TRUE` indicates PIP \> threshold. Returns `NULL` with
a warning if PIPs are not available for the requested mode (e.g. when
`ebps_with_uq` is used for F, which has no spike component).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3, Xcov = Xcov)

# Raw PIP matrix for observation loadings
pip_L <- get_pip(fit, mode = 'L')

# Which loadings exceed 0.9 PIP threshold
sig_L <- get_pip(fit, mode = 'L', threshold = 0.9)
} # }
```
