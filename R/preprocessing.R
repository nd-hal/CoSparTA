#' Build a 3-Way Tensor from Long-Format Data
#'
#' @description
#' Converts a long-format data frame into a 3-dimensional array suitable for
#' \code{\link{CoSparTA}}. Supports two input styles: (1) raw event logs where
#' each row is a single event and counts are computed automatically, or
#' (2) pre-aggregated data where a value column is already provided.
#' Optionally bins a numeric or datetime time column into intervals.
#'
#' @param data A data frame in long format.
#' @param row Character string. Column name for mode 1 (e.g., users/sessions).
#' @param col Character string. Column name for mode 2 (e.g., time).
#' @param slice Character string. Column name for mode 3 (e.g., channels).
#' @param value Character string or \code{NULL}. Column name containing counts
#'   or values. If \code{NULL}, each row is treated as a single event and
#'   occurrences are counted per \code{(row, col, slice)} combination.
#'   Default \code{NULL}.
#' @param time_bins Controls time binning for the \code{col} dimension:
#'   \describe{
#'     \item{\code{NULL}}{No binning. Treats \code{col} as categorical; each
#'       unique value becomes a time index. Default.}
#'     \item{Integer}{Number of equal-width bins. Works with numeric and
#'       \code{POSIXct}/\code{POSIXlt} datetime columns.}
#'     \item{Numeric vector}{Explicit break points for \code{cut()}.
#'       Must cover the full range of the \code{col} column.}
#'   }
#' @param row_levels Optional character or factor vector specifying the complete
#'   set of levels for the row dimension (e.g. all session IDs). If \code{NULL},
#'   levels are inferred from the data.
#' @param col_levels Optional character or factor vector specifying the complete
#'   set of levels for the column dimension (e.g. all time slots). If
#'   \code{NULL}, levels are inferred from the data. Recommended for time
#'   dimensions where empty slots should be preserved.
#' @param slice_levels Optional character or factor vector specifying the
#'   complete set of levels for the slice dimension (e.g. all
#'   channels/websites). If \code{NULL}, levels are inferred from the data.
#' @param fill Numeric. Value to fill for unobserved \code{(row, col, slice)}
#'   combinations. Default \code{0L}.
#'
#' @return A named list with:
#' \describe{
#'   \item{X}{Integer array of dimensions \code{n x p x w}.}
#'   \item{dim1_labels}{Character or factor vector of length \code{n} giving
#'     the unique values of \code{row}, in the order they appear as rows of
#'     \code{X}.}
#'   \item{dim2_labels}{Character or factor vector of length \code{p} giving
#'     the unique values (or bin labels) of \code{col}.}
#'   \item{dim3_labels}{Character or factor vector of length \code{w} giving
#'     the unique values of \code{slice}.}
#' }
#'
#' @details
#' The function handles the following time column types:
#' \describe{
#'   \item{Character or factor}{Treated as categorical. Sorted
#'     alphabetically unless already a factor with levels.}
#'   \item{Numeric or integer}{When \code{time_bins} is an integer, binned
#'     into equal-width intervals via \code{cut()}.}
#'   \item{POSIXct or POSIXlt}{When \code{time_bins} is an integer, the
#'     time range is divided into equal-width intervals. Bin labels show
#'     the interval start times. When \code{time_bins} is a vector of
#'     \code{POSIXct} break points, those are used directly.}
#' }
#'
#' When \code{value = NULL} (event counting mode), duplicate
#' \code{(row, col, slice)} combinations are counted. When \code{value} is
#' supplied, values for duplicate combinations are summed.
#'
#' @examples
#' \dontrun{
#' # Raw event log: count events per (user, time_bin, channel)
#' events <- data.frame(
#'   user = rep(paste0("u", 1:10), each = 50),
#'   timestamp = runif(500, 0, 600),
#'   channel = sample(c("view", "click", "purchase"), 500, replace = TRUE)
#' )
#' tensor_data <- build_tensor(events, row = "user", col = "timestamp",
#'                              slice = "channel", time_bins = 12)
#' dim(tensor_data$X)       # 10 x 12 x 3
#' tensor_data$dim2_labels  # bin labels
#'
#' # Pre-aggregated monthly data
#' monthly <- data.frame(
#'   user = rep(paste0("u", 1:5), each = 12),
#'   month = rep(month.abb, 5),
#'   channel = sample(c("web", "app"), 60, replace = TRUE),
#'   count = rpois(60, 3)
#' )
#' tensor_data <- build_tensor(monthly, row = "user", col = "month",
#'                              slice = "channel", value = "count")
#'
#' # Datetime timestamps with explicit bin breaks
#' events_dt <- data.frame(
#'   session = rep(1:20, each = 30),
#'   time = as.POSIXct("2024-01-01") + runif(600, 0, 86400),
#'   action = sample(c("browse", "cart", "buy"), 600, replace = TRUE)
#' )
#' tensor_data <- build_tensor(events_dt, row = "session", col = "time",
#'                              slice = "action", time_bins = 24)
#'
#' # Feed directly into CoSparTA
#' fit <- CoSparTA(tensor_data$X, K = 3)
#'
#' # Ensure all 48 time slots are present even if some are empty
#' tensor_out <- build_tensor(
#'   data       = df,
#'   user_col   = "session_id",
#'   time_col   = "hour",
#'   chan_col   = "website",
#'   count_col  = "count",
#'   col_levels = paste0("H", sprintf("%02d", 0:47))
#' )
#' }
#'
#' @seealso \code{\link{CoSparTA}}
#' @export
build_tensor <- function(data, row, col, slice, value = NULL,
                         time_bins = NULL, fill = 0L,
                         row_levels = NULL, col_levels = NULL,
                         slice_levels = NULL) {

  # --- 1. Input validation ---
  if (!is.data.frame(data)) stop("'data' must be a data frame.")
  needed <- c(row, col, slice)
  if (!is.null(value)) needed <- c(needed, value)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in 'data': ",
         paste(missing_cols, collapse = ", "))
  }

  # Working copy of relevant columns only
  d <- data.frame(
    r = data[[row]],
    c = data[[col]],
    s = data[[slice]],
    stringsAsFactors = FALSE
  )
  if (!is.null(value)) d$v <- data[[value]]

  # --- 2. Time binning ---
  if (!is.null(time_bins)) {
    col_vals <- d$c

    if (length(time_bins) == 1L) {
      # Single integer: equal-width bins
      if (inherits(col_vals, c("POSIXct", "POSIXlt"))) {
        rng    <- range(col_vals, na.rm = TRUE)
        breaks <- seq(rng[1], rng[2], length.out = time_bins + 1L)
        labels <- format(breaks[-length(breaks)])
        d$c    <- cut(col_vals, breaks = breaks, labels = labels,
                      include.lowest = TRUE)
      } else if (is.numeric(col_vals) || is.integer(col_vals)) {
        d$c <- cut(col_vals, breaks = time_bins, include.lowest = TRUE)
      } else {
        stop("'time_bins' as a single integer requires a numeric or ",
             "POSIXct/POSIXlt time column.")
      }
    } else {
      # Numeric vector: explicit break points
      d$c <- cut(col_vals, breaks = time_bins, include.lowest = TRUE)
    }
  }

  # --- 3. Factor levels to define the full universe ---
  d$r <- if (!is.null(row_levels))   factor(d$r, levels = row_levels)   else factor(d$r, levels = unique(d$r))
  d$s <- if (!is.null(slice_levels)) factor(d$s, levels = slice_levels) else factor(d$s, levels = unique(d$s))

  # col: user-supplied levels take priority; then cut() levels; otherwise sort
  if (!is.null(col_levels)) {
    d$c <- factor(d$c, levels = col_levels)
  } else if (is.factor(d$c)) {
    # already has ordered levels (from cut() or original factor)
  } else {
    d$c <- factor(d$c, levels = sort(unique(d$c)))
  }

  # --- 4. Aggregate via data.table ---
  dt <- data.table::data.table(r = d$r, c = d$c, s = d$s)

  if (is.null(value)) {
    agg <- dt[, .(count = .N), by = .(r, c, s)]
  } else {
    dt[, v := d$v]
    agg <- dt[, .(count = sum(v)), by = .(r, c, s)]
  }

  # --- 5. Build the 3D array ---
  dim1_labels <- levels(agg$r)
  dim2_labels <- levels(agg$c)
  dim3_labels <- levels(agg$s)

  n <- length(dim1_labels)
  p <- length(dim2_labels)
  w <- length(dim3_labels)

  X <- array(as.integer(fill), dim = c(n, p, w))

  i_idx <- match(as.character(agg$r), dim1_labels)
  j_idx <- match(as.character(agg$c), dim2_labels)
  m_idx <- match(as.character(agg$s), dim3_labels)
  X[cbind(i_idx, j_idx, m_idx)] <- as.integer(agg$count)

  cat(sprintf("Tensor built: %d x %d x %d (%.1f%% zeros)\n",
              n, p, w, 100 * mean(X == 0L)))

  list(
    tensor      = X,
    X           = X,
    dim1_labels = dim1_labels,
    dim2_labels = dim2_labels,
    dim3_labels = dim3_labels
  )
}
