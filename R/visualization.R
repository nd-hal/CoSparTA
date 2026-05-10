#' Plot Time-Mode Factors
#'
#' @description
#' Produces a faceted line plot of the time-mode factor matrix (Ef), with one
#' panel per component. Useful for visualizing temporal patterns learned by
#' \code{\link{CxtEBTD}}.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD}}. Either
#'   \code{fit} or \code{Ef} must be supplied.
#' @param Ef Numeric matrix of dimensions \code{p x K} (time factor matrix).
#'   If provided, used directly and \code{fit} is ignored.
#' @param lambda Numeric vector of length \code{K} (component weights). Accepted
#'   for API symmetry when \code{Ef} is supplied directly; currently unused by
#'   this function.
#' @param ranks Integer vector specifying which components to plot. Default
#'   \code{NULL} plots all K components.
#' @param time_labels Numeric vector of length \code{p} giving x-axis values
#'   (e.g., seconds, time indices). Default \code{NULL} uses \code{1:p}.
#' @param normalize Logical. If \code{TRUE} and \code{fit} is supplied, uses
#'   \code{\link{normalize_factors}} to normalize columns to unit norm. Ignored
#'   when \code{Ef} is provided directly. Default \code{TRUE}.
#' @param ncol Integer. Number of columns in the facet layout. Default \code{1}.
#' @param xlim Numeric vector of length 2 passed to
#'   \code{ggplot2::coord_cartesian(xlim = xlim)} to restrict the x-axis range.
#'   Default \code{NULL} applies no restriction.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X, K = 3, Xcov = Xcov)
#' plot_time_factors(fit)
#' plot_time_factors(fit, ranks = c(1, 3), time_labels = seq(0, 600, length.out = 20))
#'
#' # Using raw factor matrices directly
#' nf <- normalize_factors(fit)
#' plot_time_factors(Ef = nf$Ef, lambda = nf$lambda)
#' }
#'
#' @seealso \code{\link{plot_channel_factors}}, \code{\link{normalize_factors}}
#' @import ggplot2
#' @export
plot_time_factors <- function(fit = NULL, ranks = NULL, time_labels = NULL,
                               normalize = TRUE, ncol = 1,
                               Ef = NULL, lambda = NULL, xlim = NULL) {

  if (!is.null(Ef)) {
    # use directly supplied matrix
  } else if (!is.null(fit)) {
    if (normalize) {
      nf <- normalize_factors(fit)
      Ef <- nf$Ef
    } else {
      Ef <- fit$res$qf$Ef
    }
  } else {
    stop("Either 'fit' or 'Ef' must be provided.")
  }

  p <- nrow(Ef)
  K <- ncol(Ef)

  if (is.null(ranks)) ranks <- seq_len(K)
  Ef <- Ef[, ranks, drop = FALSE]

  if (is.null(time_labels)) time_labels <- seq_len(p)

  # Build long-format data frame
  df <- do.call(rbind, lapply(seq_along(ranks), function(j) {
    data.frame(
      time  = time_labels,
      value = Ef[, j],
      rank  = factor(paste0("R", ranks[j]), levels = paste0("R", ranks))
    )
  }))

  p <- ggplot(df, aes(x = time, y = value, colour = rank)) +
    geom_line(size = 0.8) +
    facet_wrap(~ rank, scales = "free_y", ncol = ncol) +
    xlab("Time") +
    ylab("Loading") +
    theme_minimal() +
    theme(legend.position = "none")

  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  p
}


#' Plot Channel-Mode Factors
#'
#' @description
#' Produces a faceted bar plot of the channel-mode factor matrix (Ew), with one
#' row of panels per component. Supports two display modes: individual channel
#' names on the x-axis, or channels grouped into higher-level categories with
#' group labels as facet columns.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD}}. Either
#'   \code{fit} or \code{Ew} must be supplied.
#' @param Ew Numeric matrix of dimensions \code{w x K} (channel factor matrix).
#'   If provided, used directly and \code{fit} is ignored.
#' @param lambda Numeric vector of length \code{K} (component weights). Accepted
#'   for API symmetry when \code{Ew} is supplied directly; currently unused by
#'   this function.
#' @param ranks Integer vector specifying which components to plot. Default
#'   \code{NULL} plots all K components.
#' @param channel_names Character vector of length \code{w} giving individual
#'   channel names. Default \code{NULL}.
#' @param channel_groups Character vector of length \code{w} giving group labels
#'   for each channel (e.g., \code{c(rep("TextEmo", 5), rep("Gaze", 25))}). When
#'   supplied, channels are grouped into facet columns. Default \code{NULL}.
#' @param normalize Logical. If \code{TRUE} and \code{fit} is supplied, uses
#'   \code{\link{normalize_factors}} to normalize columns to unit norm. Ignored
#'   when \code{Ew} is provided directly. Default \code{TRUE}.
#' @param show_names Logical. If \code{TRUE} and \code{channel_names} is supplied,
#'   shows individual channel names on the x-axis. Default \code{FALSE}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X, K = 3, Xcov = Xcov)
#'
#' # Simple: channels labeled 1 to w
#' plot_channel_factors(fit)
#'
#' # With channel names
#' plot_channel_factors(fit, channel_names = c("page_view", "cart", "purchase"))
#'
#' # With channel groups (faceted columns)
#' groups <- c(rep("TextEmo", 5), rep("Gaze", 25), rep("AU_c", 6), rep("AU_r", 14))
#' plot_channel_factors(fit, channel_groups = groups)
#'
#' # Using raw factor matrices directly
#' nf <- normalize_factors(fit)
#' plot_channel_factors(Ew = nf$Ew, lambda = nf$lambda, channel_names = c("a", "b"))
#' }
#'
#' @seealso \code{\link{plot_time_factors}}, \code{\link{normalize_factors}}
#' @export
plot_channel_factors <- function(fit = NULL, ranks = NULL, channel_names = NULL,
                                  channel_groups = NULL, normalize = TRUE,
                                  show_names = FALSE, Ew = NULL, lambda = NULL) {

  if (!is.null(Ew)) {
    # use directly supplied matrix
  } else if (!is.null(fit)) {
    if (normalize) {
      nf <- normalize_factors(fit)
      Ew <- nf$Ew
    } else {
      Ew <- fit$res$qw$Ew
    }
  } else {
    stop("Either 'fit' or 'Ew' must be provided.")
  }

  w <- nrow(Ew)
  K <- ncol(Ew)

  if (is.null(ranks)) ranks <- seq_len(K)
  Ew <- Ew[, ranks, drop = FALSE]

  # Build long-format data frame
  df <- do.call(rbind, lapply(seq_along(ranks), function(j) {
    d <- data.frame(
      channel_idx = seq_len(w),
      value       = Ew[, j],
      rank        = factor(paste0("R", ranks[j]), levels = paste0("R", ranks))
    )
    if (!is.null(channel_names)) d$channel_name  <- channel_names
    if (!is.null(channel_groups)) {
      # Preserve first-appearance order for group levels
      grp_levels <- unique(channel_groups)
      d$channel_group <- factor(channel_groups, levels = grp_levels)
    }
    d
  }))

  base_theme <- theme_minimal() +
    theme(
      legend.position  = "none",
      strip.text.y     = element_text(size = 10),
      axis.text.y      = element_text(size = 5)
    )

  if (!is.null(channel_groups)) {
    # Compute within-group sequential index for x-axis
    df$within_idx <- ave(df$channel_idx, df$rank, df$channel_group,
                         FUN = function(i) seq_along(i))

    p <- ggplot(df, aes(x = within_idx, y = value, fill = channel_group)) +
      geom_bar(stat = "identity") +
      facet_grid(rank ~ channel_group, scales = "free") +
      xlab("Channels") +
      ylab("Loading") +
      base_theme +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
      )

  } else {
    if (!is.null(channel_names) && show_names) {
      x_var  <- "channel_name"
      x_text <- element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
    } else {
      x_var  <- "channel_idx"
      x_text <- element_text()
    }

    p <- ggplot(df, aes_string(x = x_var, y = "value", fill = "rank")) +
      geom_bar(stat = "identity") +
      facet_wrap(~ rank, scales = "free_y", ncol = 1) +
      xlab("Channels") +
      ylab("Loading") +
      base_theme +
      theme(axis.text.x = x_text)
  }

  p
}
