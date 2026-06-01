# Plot Channel-Mode Factors

Produces a faceted bar plot of the channel-mode factor matrix (Ew), with
one row of panels per component. Supports two display modes: individual
channel names on the x-axis, or channels grouped into higher-level
categories with group labels as facet columns.

## Usage

``` r
plot_channel_factors(
  fit = NULL,
  ranks = NULL,
  channel_names = NULL,
  channel_groups = NULL,
  normalize = TRUE,
  show_names = FALSE,
  Ew = NULL,
  lambda = NULL
)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).
  Either `fit` or `Ew` must be supplied.

- ranks:

  Integer vector specifying which components to plot. Default `NULL`
  plots all K components.

- channel_names:

  Character vector of length `w` giving individual channel names.
  Default `NULL`.

- channel_groups:

  Character vector of length `w` giving group labels for each channel
  (e.g., `c(rep("TextEmo", 5), rep("Gaze", 25))`). When supplied,
  channels are grouped into facet columns. Default `NULL`.

- normalize:

  Logical. If `TRUE` (default) and `fit` is supplied, reads
  `fit$res$qw$Ew_normed` (unit-norm, descending \\\lambda\\ order).
  Falls back to
  [`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md)
  for legacy fits that lack `_normed` fields. Ignored when `Ew` is
  provided directly.

- show_names:

  Logical. If `TRUE` and `channel_names` is supplied, shows individual
  channel names on the x-axis. Default `FALSE`.

- Ew:

  Numeric matrix of dimensions `w x K` (channel factor matrix). If
  provided, used directly and `fit` is ignored.

- lambda:

  Numeric vector of length `K` (component weights). Accepted for API
  symmetry when `Ew` is supplied directly; currently unused by this
  function.

## Value

A `ggplot` object.

## See also

[`plot_time_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_time_factors.md),
[`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3, Xcov = Xcov)

# Simple: channels labeled 1 to w
plot_channel_factors(fit)

# With channel names
plot_channel_factors(fit, channel_names = c("page_view", "cart", "purchase"))

# With channel groups (faceted columns)
groups <- c(rep("TextEmo", 5), rep("Gaze", 25), rep("AU_c", 6), rep("AU_r", 14))
plot_channel_factors(fit, channel_groups = groups)

# Using raw factor matrices directly
nf <- normalize_factors(fit)
plot_channel_factors(Ew = nf$Ew, lambda = nf$lambda, channel_names = c("a", "b"))
} # }
```
