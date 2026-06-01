# Standardize Poisson matrix factorization to multinomial parameterization

Converts a Poisson factorization \\X \approx L F^\top\\ into a
standardized form where columns of F sum to 1 and L absorbs the scale.

## Usage

``` r
poisson_to_multinom(FF, L)
```

## Arguments

- FF:

  Numeric matrix of dimensions `p x K` (factors).

- L:

  Numeric matrix of dimensions `n x K` (loadings).

## Value

A named list with:

- FF:

  Column-normalized factor matrix (`p x K`).

- L:

  Rescaled loading matrix (`n x K`).

- s:

  Numeric vector of length `n` giving row scales.
