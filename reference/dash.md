# Launch an Interactive Fit-Explorer Dashboard

Launches a local, interactive Quarto dashboard for exploring a fitted
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md) /
[`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md)
object. The dashboard shows the temporal factor pattern
([`plot_time_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_time_factors.md)),
the channel factor pattern
([`plot_channel_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_channel_factors.md)),
a component-weight chart, and a table of covariate (gamma) coefficients,
all restricted to a set of components selected via checkboxes.

## Usage

``` r
dash(fit, channel_names = NULL, time_labels = NULL)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  or
  [`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).

- channel_names:

  Character vector of channel labels, forwarded to
  [`plot_channel_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_channel_factors.md).
  Fixed at launch; not adjustable from the dashboard. Default `NULL`.

- time_labels:

  Numeric vector of time-axis labels, forwarded to
  [`plot_time_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_time_factors.md).
  Fixed at launch; not adjustable from the dashboard. Default `NULL`.

## Value

Invisibly, the path to the temporary `.qmd` file used to launch the
dashboard. Called primarily for its side effect of starting the
dashboard server.

## Details

The dashboard is a Quarto document (`format: dashboard`,
`server: shiny`) bundled with the package at
`inst/dashboard/cosparta_dashboard.qmd`. `dash()` saves `fit` (together
with `channel_names` and `time_labels`) to a temporary `.rds` file next
to a copy of the template, then calls
[`quarto_preview`](https://quarto-dev.github.io/quarto-r/reference/quarto_preview.html),
which starts a local Shiny-backed server and opens a browser tab.

This call blocks the R console for as long as the dashboard is open —
the same behavior as
[`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html). That is
expected, not a bug. Stop the server (e.g. Escape in RStudio, or Ctrl+C
in a terminal) to return control to the console.

The only interactive control is a checkbox group for selecting which
components (ranks) to display; all four panels re-render from that
selection. `channel_names` and `time_labels` are read once at launch and
are not exposed as dashboard inputs.

Requires a local Quarto CLI installation. If it is not found, install it
from <https://quarto.org/docs/get-started/>.

## See also

[`plot_time_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_time_factors.md),
[`plot_channel_factors`](https://nd-hal.github.io/CoSparTA/reference/plot_channel_factors.md),
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md),
[`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3, Xcov = Xcov)
dash(fit, channel_names = c("Retail", "Social", "Search"),
     time_labels = seq(0, 600, length.out = 20))
} # }
```
