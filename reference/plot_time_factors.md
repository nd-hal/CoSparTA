# Plot Time-Mode Factors

Produces a faceted line plot of the time-mode factor matrix (Ef), with
one panel per component. Useful for visualizing temporal patterns
learned by
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).

## Usage

``` r
plot_time_factors(
  fit = NULL,
  ranks = NULL,
  time_labels = NULL,
  normalize = TRUE,
  ncol = 1,
  Ef = NULL,
  lambda = NULL,
  xlim = NULL
)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).
  Either `fit` or `Ef` must be supplied.

- ranks:

  Integer vector specifying which components to plot. Default `NULL`
  plots all K components.

- time_labels:

  Numeric vector of length `p` giving x-axis values (e.g., seconds, time
  indices). Default `NULL` uses `1:p`.

- normalize:

  Logical. If `TRUE` (default) and `fit` is supplied, reads
  `fit$res$qf$Ef_normed` (unit-norm, descending \\\lambda\\ order).
  Falls back to
  [`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md)
  for legacy fits that lack `_normed` fields. Ignored when `Ef` is
  provided directly.

- ncol:

  Integer. Number of columns in the facet layout. Default `1`.

- Ef:

  Numeric matrix of dimensions `p x K` (time factor matrix). If
  provided, used directly and `fit` is ignored.

- lambda:

  Numeric vector of length `K` (component weights). Accepted for API
  symmetry when `Ef` is supplied directly; currently unused by this
  function.

- xlim:

  Numeric vector of length 2 passed to
  `ggplot2::coord_cartesian(xlim = xlim)` to restrict the x-axis range.
  Default `NULL` applies no restriction.

## Value

A `ggplot` object.

## See also

[`plot_channel_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_channel_factors.md),
[`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3, Xcov = Xcov)
plot_time_factors(fit)
plot_time_factors(fit, ranks = c(1, 3), time_labels = seq(0, 600, length.out = 20))

# Using raw factor matrices directly
nf <- normalize_factors(fit)
plot_time_factors(Ef = nf$Ef, lambda = nf$lambda)
} # }
```
