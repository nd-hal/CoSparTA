# Retrieve Factor Matrices or Weights from a CoSparTA Fit

Unified getter that retrieves a named output from a fitted
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
object by logical name, so users do not need to navigate the internal
slot structure (e.g.\\ `fit$res$ql$El_normed`).

## Usage

``` r
get_loadings(fit, what = "U1", normalized = TRUE)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  or
  [`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).

- what:

  Character string specifying what to retrieve. One of:

  `"U1"`

  :   Observation loading matrix (`n x K`).

  `"U2"`

  :   Time factor matrix (`p x K`).

  `"U3"`

  :   Channel weight matrix (`w x K`).

  `"weight"`

  :   Component weight vector \\\lambda\\ of length `K` (product of mode
      norms; sorted descending when `normalized = TRUE`).

  `"gamma"`

  :   List of length `K` of covariate coefficient vectors (one per
      component). `NULL` entries indicate unsupervised components.

- normalized:

  Logical. If `TRUE` (default), returns the normalized and reordered
  version stored in the `_normed` fields (unit-norm columns, descending
  \\\lambda\\ order). If `FALSE`, returns the raw posterior estimates in
  original component order.

## Value

The requested object:

- `"U1"`, `"U2"`, `"U3"`:

  A numeric matrix.

- `"weight"`:

  A numeric vector of length `K`.

- `"gamma"`:

  A list of length `K`.

## See also

[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md),
[`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3, Xcov = Xcov)

# Normalized observation loadings (n x 3, unit-norm columns)
U1 <- get_loadings(fit, "U1")

# Raw observation loadings
U1_raw <- get_loadings(fit, "U1", normalized = FALSE)

# Component weights (descending)
w <- get_loadings(fit, "weight")

# Covariate coefficient lists
gamma <- get_loadings(fit, "gamma")
gamma[[1]]  # gamma for component 1
} # }
```
