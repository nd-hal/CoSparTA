# Sparse Weighted Aggregation for Tensor EM

For a sparse 3-way tensor stored in coordinate format, computes the
rank-k sufficient statistics (row sums, column sums, slice sums) needed
for the EM update of factor matrices L, F, W. Replaces the dplyr
group-by approach in `calc_EZ_3d()`.

## Usage

``` r
calc_EZ_3d_cpp(V1, V2, V3, x_vals, alpha_k, n, p, w)
```

## Arguments

- V1:

  Integer vector of user (L-mode) indices (1-based).

- V2:

  Integer vector of time (F-mode) indices (1-based).

- V3:

  Integer vector of channel (W-mode) indices (1-based).

- x_vals:

  Numeric vector of observed counts at each nonzero position.

- alpha_k:

  Numeric vector of length `nnz`: the rank-k responsibility weights
  (alpha\[, k\]).

- n:

  Integer. Number of users (L-mode dimension).

- p:

  Integer. Number of time points (F-mode dimension).

- w:

  Integer. Number of channels (W-mode dimension).

## Value

A named list with three numeric vectors:

- rs:

  Length-`n` vector: row sums (L-mode aggregation).

- cs:

  Length-`p` vector: column sums (F-mode aggregation).

- zs:

  Length-`w` vector: slice sums (W-mode aggregation).
