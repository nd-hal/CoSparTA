# Identify Significant Channels and Time Points per Factor via lFDR Control

Implements Bayesian local-FDR control (Algorithms 1 and 2 from the EBTD
framework) to identify which channels and time points are truly active
for each factor. For each factor `l`, the posterior spike probabilities
(local-fdr values) from the F and W modes are sorted and thresholded to
produce a discovery set at a controlled FDR level.

The local-fdr value for entry \\(j, l)\\ is \\\rho\_{jl} = 1 -
\text{PIP}\_{jl}\\, i.e., the posterior probability that the entry is
truly zero. Algorithm 1 finds the largest discovery set such that the
mean local-fdr does not exceed `alpha`.

## Usage

``` r
get_significant_patterns(fit, alpha = 0.05, mode = "both", normalized = TRUE)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).

- alpha:

  Numeric FDR level in `(0, 1)`. Default `0.05`.

- mode:

  Character string specifying which mode(s) to run discovery on: `'F'`
  for time factors, `'W'` for channel weights, or `'both'`. Default
  `'both'`.

- normalized:

  Logical. If `TRUE` (default), reads from `_normed` PIP fields. If
  `FALSE`, reads from raw PIP fields.

## Value

A list of length `K` (one element per factor). Each element is a named
list with:

- factor:

  Integer factor index.

- active_times:

  Integer vector of active time point indices for this factor (from F
  mode). `NULL` if mode is `'W'` or PIPs unavailable.

- active_channels:

  Integer vector of active channel indices for this factor (from W
  mode). `NULL` if mode is `'F'` or PIPs unavailable.

- n_active_times:

  Number of active time points discovered.

- n_active_channels:

  Number of active channels discovered.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3, Xcov = Xcov)

# Discover active channels and time points at 5% FDR
patterns <- get_significant_patterns(fit, alpha = 0.05)

# Factor 1 active time points
patterns[[1]]$active_times

# Factor 2 active channels
patterns[[2]]$active_channels
} # }
```
