# Sparse Softmax (Responsibility) Computation

Computes the softmax responsibility matrix alpha (nnz x K) at nonzero
positions of a sparse 3-way tensor, using the log-sum-exp trick for
numerical stability. Replaces the five-line R block:


      alpha = Elogl[V1,] + Elogf[V2,] + Elogw[V3,]
      exp_offset = rowMaxs(alpha)
      alpha = alpha - outer(exp_offset, rep(1,K))
      alpha = exp(alpha)
      alpha = alpha / rowsums(alpha)

## Usage

``` r
calc_qz_sparse_cpp(V1, V2, V3, Elogl, Elogf, Elogw)
```

## Arguments

- V1:

  Integer vector of user (L-mode) indices (1-based), length nnz.

- V2:

  Integer vector of time (F-mode) indices (1-based), length nnz.

- V3:

  Integer vector of channel (W-mode) indices (1-based), length nnz.

- Elogl:

  Numeric matrix (n x K): expected log loadings for L-mode.

- Elogf:

  Numeric matrix (p x K): expected log loadings for F-mode.

- Elogw:

  Numeric matrix (w x K): expected log loadings for W-mode.

## Value

A numeric matrix of size (nnz x K): the row-normalized softmax
responsibilities at each nonzero position.
