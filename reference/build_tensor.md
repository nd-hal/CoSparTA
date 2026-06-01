# Build a 3-Way Tensor from Long-Format Data

Converts a long-format data frame into a 3-dimensional array suitable
for
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).
Supports two input styles: (1) raw event logs where each row is a single
event and counts are computed automatically, or (2) pre-aggregated data
where a value column is already provided. Optionally bins a numeric or
datetime time column into intervals.

## Usage

``` r
build_tensor(
  data,
  row,
  col,
  slice,
  value = NULL,
  time_bins = NULL,
  fill = 0L,
  row_levels = NULL,
  col_levels = NULL,
  slice_levels = NULL
)
```

## Arguments

- data:

  A data frame in long format.

- row:

  Character string. Column name for mode 1 (e.g., users/sessions).

- col:

  Character string. Column name for mode 2 (e.g., time).

- slice:

  Character string. Column name for mode 3 (e.g., channels).

- value:

  Character string or `NULL`. Column name containing counts or values.
  If `NULL`, each row is treated as a single event and occurrences are
  counted per `(row, col, slice)` combination. Default `NULL`.

- time_bins:

  Controls time binning for the `col` dimension:

  `NULL`

  :   No binning. Treats `col` as categorical; each unique value becomes
      a time index. Default.

  Integer

  :   Number of equal-width bins. Works with numeric and
      `POSIXct`/`POSIXlt` datetime columns.

  Numeric vector

  :   Explicit break points for
      [`cut()`](https://rdrr.io/r/base/cut.html). Must cover the full
      range of the `col` column.

- fill:

  Numeric. Value to fill for unobserved `(row, col, slice)`
  combinations. Default `0L`.

- row_levels:

  Optional character or factor vector specifying the complete set of
  levels for the row dimension (e.g. all session IDs). If `NULL`, levels
  are inferred from the data.

- col_levels:

  Optional character or factor vector specifying the complete set of
  levels for the column dimension (e.g. all time slots). If `NULL`,
  levels are inferred from the data. Recommended for time dimensions
  where empty slots should be preserved.

- slice_levels:

  Optional character or factor vector specifying the complete set of
  levels for the slice dimension (e.g. all channels/websites). If
  `NULL`, levels are inferred from the data.

## Value

A named list with:

- X:

  Integer array of dimensions `n x p x w`.

- dim1_labels:

  Character or factor vector of length `n` giving the unique values of
  `row`, in the order they appear as rows of `X`.

- dim2_labels:

  Character or factor vector of length `p` giving the unique values (or
  bin labels) of `col`.

- dim3_labels:

  Character or factor vector of length `w` giving the unique values of
  `slice`.

## Details

The function handles the following time column types:

- Character or factor:

  Treated as categorical. Sorted alphabetically unless already a factor
  with levels.

- Numeric or integer:

  When `time_bins` is an integer, binned into equal-width intervals via
  [`cut()`](https://rdrr.io/r/base/cut.html).

- POSIXct or POSIXlt:

  When `time_bins` is an integer, the time range is divided into
  equal-width intervals. Bin labels show the interval start times. When
  `time_bins` is a vector of `POSIXct` break points, those are used
  directly.

When `value = NULL` (event counting mode), duplicate `(row, col, slice)`
combinations are counted. When `value` is supplied, values for duplicate
combinations are summed.

## See also

[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Raw event log: count events per (user, time_bin, channel)
events <- data.frame(
  user = rep(paste0("u", 1:10), each = 50),
  timestamp = runif(500, 0, 600),
  channel = sample(c("view", "click", "purchase"), 500, replace = TRUE)
)
tensor_data <- build_tensor(events, row = "user", col = "timestamp",
                             slice = "channel", time_bins = 12)
dim(tensor_data$X)       # 10 x 12 x 3
tensor_data$dim2_labels  # bin labels

# Pre-aggregated monthly data
monthly <- data.frame(
  user = rep(paste0("u", 1:5), each = 12),
  month = rep(month.abb, 5),
  channel = sample(c("web", "app"), 60, replace = TRUE),
  count = rpois(60, 3)
)
tensor_data <- build_tensor(monthly, row = "user", col = "month",
                             slice = "channel", value = "count")

# Datetime timestamps with explicit bin breaks
events_dt <- data.frame(
  session = rep(1:20, each = 30),
  time = as.POSIXct("2024-01-01") + runif(600, 0, 86400),
  action = sample(c("browse", "cart", "buy"), 600, replace = TRUE)
)
tensor_data <- build_tensor(events_dt, row = "session", col = "time",
                             slice = "action", time_bins = 24)

# Feed directly into CoSparTA
fit <- CoSparTA(tensor_data$X, K = 3)

# Ensure all 48 time slots are present even if some are empty
tensor_out <- build_tensor(
  data       = df,
  row        = "session_id",
  col        = "hour",
  slice      = "website",
  value      = "count",
  col_levels = paste0("H", sprintf("%02d", 0:47))
)
} # }
```
