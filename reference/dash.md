# Launch an Interactive Fit-Explorer Dashboard

Launches a local, interactive Quarto dashboard for exploring a fitted
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md) /
[`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md)
object, or a set of raw (already normalized) factor matrices. The
dashboard shows the temporal factor pattern
([`plot_time_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_time_factors.md)),
the channel factor pattern
([`plot_channel_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_channel_factors.md)),
a component-weight chart, and a table of covariate (gamma) coefficients,
all restricted to a set of components selected via checkboxes.

## Usage

``` r
dash(
  fit = NULL,
  channel_names = NULL,
  time_labels = NULL,
  channel_groups = NULL,
  covariate_names = NULL,
  Ef = NULL,
  Ew = NULL,
  lambda = NULL,
  gamma_list = NULL,
  intercept = TRUE
)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  or
  [`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).
  Either `fit` or all of `Ef`, `Ew`, and `lambda` must be supplied.
  `fit$res$ql$El_normed` (per-observation loadings) is never read — the
  dashboard has no observation-mode panel.

- channel_names:

  Character vector of channel labels, forwarded to
  [`plot_channel_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_channel_factors.md).
  Fixed at launch; not adjustable from the dashboard. Default `NULL`.

- time_labels:

  Numeric vector of time-axis labels, forwarded to
  [`plot_time_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_time_factors.md).
  Fixed at launch; not adjustable from the dashboard. Default `NULL`.

- channel_groups:

  Character vector of group labels for each channel, forwarded to
  [`plot_channel_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_channel_factors.md)'s
  `channel_groups` argument (e.g.
  `c(rep("TextEmo", 5), rep("Gaze", 25))`). Fixed at launch; not
  adjustable from the dashboard. Default `NULL`.

- covariate_names:

  Names for the gamma table's covariate columns. When `intercept = TRUE`
  (the default), these describe every column except the always-present,
  always-first `"Intercept"` column. When `intercept = FALSE`, there is
  no `"Intercept"` column and these describe every column. Two forms are
  accepted, mirroring the `Xcov` convention used elsewhere in the
  package:

  a character vector

  :   applied to every supervised component. Requires all supervised
      components to have the same number of covariates; errors with a
      clear message otherwise (use the list form instead when they
      differ).

  a list of length `K`

  :   per-component names. `covariate_names[[k]]` must have length
      `length(gamma_list[[k]]) - 1` when `intercept = TRUE`, or exactly
      `length(gamma_list[[k]])` when `intercept = FALSE`. Entries for
      unsupervised components are ignored.

  Default `NULL` uses generic `"V1"`, `"V2"`, ... labels. Works with
  either input path, including the `fit` path.

- Ef:

  Numeric matrix of dimensions `p x K` (time factor matrix, already
  normalized — e.g. `fit$res$qf$Ef_normed` or
  `normalize_factors(fit)$Ef`). Supply together with `Ew` and `lambda`
  to bypass `fit`, mirroring the `Ef`/`Ew` raw-matrix path already
  supported by
  [`plot_time_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_time_factors.md)
  and
  [`plot_channel_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_channel_factors.md).

- Ew:

  Numeric matrix of dimensions `w x K` (channel factor matrix, already
  normalized). Supply together with `Ef` and `lambda` to bypass `fit`.

- lambda:

  Numeric vector of length `K` giving component weights (used by the
  component-weight panel). Supply together with `Ef` and `Ew` to bypass
  `fit`.

- gamma_list:

  List of length `K` of covariate coefficient vectors, one per
  component, in the same format returned by
  [`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md)'s
  `gamma_list`: each element is either a numeric vector (whose first
  entry is the intercept when `intercept = TRUE`, or is itself a
  covariate coefficient when `intercept = FALSE`) or a scalar `NA` for a
  component with no covariates. `NULL` (the default) means every
  component is unsupervised. Only used when bypassing `fit`; ignored
  (and derived from `fit$res$gl_normed` instead) when `fit` is supplied.

- intercept:

  Logical. Whether each element of `gamma_list` begins with an intercept
  term. Default `TRUE`: `gamma_list[[k]][1]` is the intercept,
  `covariate_names[[k]]` describes the remaining entries, and the gamma
  table's first column is always labeled `"Intercept"`. Set to `FALSE`
  for fits where `Xcov` was supplied without an intercept column, so
  `gamma_list[[k]]` is entirely covariate coefficients:
  `covariate_names[[k]]` must then match its length exactly, and the
  gamma table has no `"Intercept"` column. If left at its default and
  `covariate_names` lengths match `gamma_list` lengths exactly rather
  than length `- 1`, this is unambiguous evidence of a no-intercept fit:
  `dash()` switches to `intercept = FALSE` automatically and warns,
  rather than erroring. Pass `intercept = FALSE` explicitly to silence
  that warning.

## Value

Invisibly, the path to the temporary `.qmd` file used to launch the
dashboard. Called primarily for its side effect of starting the
dashboard server.

## Details

The dashboard is a Quarto document (`format: dashboard`,
`server: shiny`) bundled with the package at
`inst/dashboard/cosparta_dashboard.qmd`. `dash()` resolves `Ef`, `Ew`,
`lambda`, and `gamma_list` — either directly from the arguments of the
same name, or extracted from `fit` — saves them (together with
`channel_names`, `time_labels`, `channel_groups`, and the resolved
`covariate_names`) to a temporary `.rds` file next to a copy of the
template, then calls
[`quarto_preview`](https://quarto-dev.github.io/quarto-r/reference/quarto_preview.html),
which starts a local Shiny-backed server and opens a browser tab.

This call blocks the R console for as long as the dashboard is open —
the same behavior as
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html). That is
expected, not a bug. Stop the server (e.g. Escape in RStudio, or Ctrl+C
in a terminal) to return control to the console.

The only interactive control is a checkbox group for selecting which
components (ranks) to display; all four panels re-render from that
selection. `channel_names`, `time_labels`, `channel_groups`, and
`covariate_names` are read once at launch and are not exposed as
dashboard inputs.

Requires a local Quarto CLI installation. If it is not found, install it
from <https://quarto.org/docs/get-started/>.

## See also

[`plot_time_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_time_factors.md),
[`plot_channel_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_channel_factors.md),
[`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md),
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md),
[`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3, Xcov = Xcov)
dash(fit, channel_names = c("Retail", "Social", "Search"),
     time_labels = seq(0, 600, length.out = 20))

# Raw-matrix path, bypassing fit (e.g. after normalize_factors())
nf <- normalize_factors(fit)
dash(Ef = nf$Ef, Ew = nf$Ew, lambda = nf$lambda, gamma_list = nf$gamma_list,
     covariate_names = c("age", "income"))

# Per-component covariate names, differing counts across components
dash(fit, covariate_names = list(c("age", "income"), NULL, c("region")))

# Fit without an intercept in Xcov: gamma_list[[k]] is all covariates
dash(Ef = nf$Ef, Ew = nf$Ew, lambda = nf$lambda, gamma_list = nf$gamma_list,
     covariate_names = c("age", "income"), intercept = FALSE)
} # }
```
