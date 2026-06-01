# Log-transformation of count matrix for empirical Bayes matrix factorization

Applies a library-size normalized log transformation suitable as
preprocessing input for EBMF methods.

## Usage

``` r
log_for_ebmf(Y)
```

## Arguments

- Y:

  Non-negative integer count matrix (`n x p`).

## Value

A numeric matrix of the same dimensions as `Y` with log-transformed,
library-size normalized values.
