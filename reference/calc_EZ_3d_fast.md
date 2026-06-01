# Fast Sparse EZ Aggregation Wrapper

Thin R wrapper around
[`calc_EZ_3d_cpp`](https://nd-hal.github.io/CoSparTA/reference/calc_EZ_3d_cpp.md)
that accepts the same sparse coordinate-format data.table `x` used by
`calc_EZ_3d()`, making the Rcpp replacement a drop-in at the call site
in
[`CoSparTA()`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
and
[`CoSparTA_missing()`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).

## Usage

``` r
calc_EZ_3d_fast(x, alpha_k, n, p, w)
```

## Arguments

- x:

  A data.table with integer columns `V1`, `V2`, `V3` (1-based indices)
  and a numeric column `x` (observed counts).

- alpha_k:

  Numeric vector of length `nrow(x)`: rank-k responsibilities.

- n:

  Integer. L-mode (user) dimension.

- p:

  Integer. F-mode (time) dimension.

- w:

  Integer. W-mode (channel) dimension.

## Value

Named list with `rs` (length n), `cs` (length p), `zs` (length w).
