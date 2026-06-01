# Match and Align Estimated Factors to Reference Factors

Finds the optimal permutation that aligns the columns of estimated
factor matrices to reference (true) factor matrices using the Tucker
congruence coefficient (cosine similarity) and the Hungarian algorithm.
This is the standard method for comparing tensor decompositions in
simulation studies, where the component ordering is arbitrary.

For each pair of columns (one from reference, one from estimated), the
congruence coefficient is \\a^\top b / (\\a\\ \\b\\)\\. When multiple
modes are provided, the joint congruence is the product of per-mode
congruences, so a good match must be good on all modes simultaneously.

## Usage

``` r
match_factors(ref, est, absolute_value = TRUE)
```

## Arguments

- ref:

  A list of reference factor matrices, one per mode. Each matrix has
  dimensions `d_m x K`. Typically the true factors from a simulation.

- est:

  A list of estimated factor matrices, same structure as `ref`.
  Dimensions must match: same number of modes, same `d_m` per mode, same
  K.

- absolute_value:

  Logical. If `TRUE` (default), uses absolute cosine similarity,
  ignoring sign differences. Appropriate for non-negative
  decompositions.

## Value

A named list with:

- permutation:

  Integer vector of length K giving the optimal column permutation.
  `permutation[i]` is the column index in `est` that best matches column
  `i` in `ref`.

- mean_congruence:

  Numeric scalar: mean congruence coefficient across all K matched
  pairs. 1.0 = perfect recovery.

- per_component:

  Numeric vector of length K giving the congruence for each matched
  pair.

## Details

Requires the `clue` package for the Hungarian algorithm
([`clue::solve_LSAP`](https://rdrr.io/pkg/clue/man/solve_LSAP.html)). If
not installed, falls back to a greedy matching algorithm with a warning.

## See also

[`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md),
[`simulate_tensor`](https://nd-hal.github.io/CoSparTA/reference/simulate_tensor.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulation: compare estimated to true factors
sim <- simulate_tensor(n = 100, p = 20, w = 10, K = 3)
fit <- CoSparTA(sim$X, K = 3, Xcov = sim$Xcov)
nf <- normalize_factors(fit)

result <- match_factors(
  ref = list(sim$U1_true, sim$U2_true, sim$U3_true),
  est = list(nf$El, nf$Ef, nf$Ew)
)
result$mean_congruence   # overall recovery quality
result$permutation       # which estimated component matches which true one
} # }
```
