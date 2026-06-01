# Mean KL divergence between two non-negative matrices

Computes the mean generalized KL divergence \\D(A \\ B) = A \log(A/B) -
A + B\\, averaged over all elements. Used as a convergence criterion in
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).

## Usage

``` r
mKL(A, B)
```

## Arguments

- A:

  Non-negative numeric matrix or vector (observed).

- B:

  Non-negative numeric matrix or vector (fitted), same dimensions as
  `A`.

## Value

A single numeric value giving the mean KL divergence.
