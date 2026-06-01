# Rescale loadings and factors to similar column norms

Adjusts the column scales of L and F so that their geometric mean is
preserved. Used for numerical stability when K \> 1.

## Usage

``` r
adjLF(L, FF)
```

## Arguments

- L:

  Numeric matrix of dimensions `n x K` (observation loadings).

- FF:

  Numeric matrix of dimensions `p x K` (time factors).

## Value

A named list with:

- L_init:

  Rescaled loading matrix, same dimensions as `L`.

- F_init:

  Rescaled factor matrix, same dimensions as `FF`.
