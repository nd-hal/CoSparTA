#' Launch an Interactive Fit-Explorer Dashboard
#'
#' @description
#' Launches a local, interactive Quarto dashboard for exploring a fitted
#' \code{\link{CoSparTA}} / \code{\link{CoSparTA_missing}} object. The
#' dashboard shows the temporal factor pattern
#' (\code{\link{plot_time_factors}}), the channel factor pattern
#' (\code{\link{plot_channel_factors}}), a component-weight chart, and a
#' table of covariate (gamma) coefficients, all restricted to a set of
#' components selected via checkboxes.
#'
#' @param fit A fitted object returned by \code{\link{CoSparTA}} or
#'   \code{\link{CoSparTA_missing}}.
#' @param channel_names Character vector of channel labels, forwarded to
#'   \code{\link{plot_channel_factors}}. Fixed at launch; not adjustable from
#'   the dashboard. Default \code{NULL}.
#' @param time_labels Numeric vector of time-axis labels, forwarded to
#'   \code{\link{plot_time_factors}}. Fixed at launch; not adjustable from
#'   the dashboard. Default \code{NULL}.
#'
#' @details
#' The dashboard is a Quarto document (\code{format: dashboard},
#' \code{server: shiny}) bundled with the package at
#' \code{inst/dashboard/cosparta_dashboard.qmd}. \code{dash()} saves
#' \code{fit} (together with \code{channel_names} and \code{time_labels}) to
#' a temporary \code{.rds} file next to a copy of the template, then calls
#' \code{\link[quarto]{quarto_preview}}, which starts a local Shiny-backed
#' server and opens a browser tab.
#'
#' This call blocks the R console for as long as the dashboard is open —
#' the same behavior as \code{shiny::runApp()}. That is expected, not a bug.
#' Stop the server (e.g. Escape in RStudio, or Ctrl+C in a terminal) to
#' return control to the console.
#'
#' The only interactive control is a checkbox group for selecting which
#' components (ranks) to display; all four panels re-render from that
#' selection. \code{channel_names} and \code{time_labels} are read once at
#' launch and are not exposed as dashboard inputs.
#'
#' Requires a local Quarto CLI installation. If it is not found, install it
#' from \url{https://quarto.org/docs/get-started/}.
#'
#' @return Invisibly, the path to the temporary \code{.qmd} file used to
#'   launch the dashboard. Called primarily for its side effect of starting
#'   the dashboard server.
#'
#' @examples
#' \dontrun{
#' fit <- CoSparTA(X, K = 3, Xcov = Xcov)
#' dash(fit, channel_names = c("Retail", "Social", "Search"),
#'      time_labels = seq(0, 600, length.out = 20))
#' }
#'
#' @seealso \code{\link{plot_time_factors}}, \code{\link{plot_channel_factors}},
#'   \code{\link{CoSparTA}}, \code{\link{CoSparTA_missing}}
#' @export
dash <- function(fit, channel_names = NULL, time_labels = NULL) {

  if (!requireNamespace("quarto", quietly = TRUE)) {
    stop("Package 'quarto' is required to launch the dashboard. ",
         "Install it with install.packages('quarto').")
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to launch the dashboard. ",
         "Install it with install.packages('shiny').")
  }
  if (is.null(quarto::quarto_path())) {
    stop(
      "Quarto CLI was not found on this system. The CoSparTA dashboard ",
      "requires a local Quarto installation.\n",
      "Install it from https://quarto.org/docs/get-started/, restart R, ",
      "and try dash() again."
    )
  }

  # --- validate fit structure -------------------------------------------
  El <- fit$res$ql$El_normed
  if (is.null(El) || !is.matrix(El)) {
    stop(
      "`fit` does not look like a fitted CoSparTA / CoSparTA_missing ",
      "object: expected a matrix at fit$res$ql$El_normed. Pass the object ",
      "returned by CoSparTA() or CoSparTA_missing()."
    )
  }
  K  <- ncol(El)
  Ef <- fit$res$qf$Ef_normed
  Ew <- fit$res$qw$Ew_normed
  if (is.null(Ef) || !is.matrix(Ef) || ncol(Ef) != K) {
    stop("`fit` is missing a valid fit$res$qf$Ef_normed (p x K) matrix.")
  }
  if (is.null(Ew) || !is.matrix(Ew) || ncol(Ew) != K) {
    stop("`fit` is missing a valid fit$res$qw$Ew_normed (w x K) matrix.")
  }
  if (is.null(fit$res$lambda_normed) || length(fit$res$lambda_normed) != K) {
    stop("`fit` is missing fit$res$lambda_normed (component weights) of length K.")
  }
  if (is.null(fit$res$gl_normed) || length(fit$res$gl_normed) != K) {
    stop("`fit` is missing fit$res$gl_normed (gamma estimates) of length K.")
  }

  if (!is.null(channel_names) && length(channel_names) != nrow(Ew)) {
    stop(sprintf("`channel_names` has length %d but fit has %d channels.",
                  length(channel_names), nrow(Ew)))
  }
  if (!is.null(time_labels) && length(time_labels) != nrow(Ef)) {
    stop(sprintf("`time_labels` has length %d but fit has %d time points.",
                  length(time_labels), nrow(Ef)))
  }

  # --- stage a private copy of the template + data -----------------------
  template <- system.file("dashboard", "cosparta_dashboard.qmd", package = "CoSparTA")
  if (!nzchar(template)) {
    stop("Could not locate the bundled dashboard template ",
         "(inst/dashboard/cosparta_dashboard.qmd). Is CoSparTA installed correctly?")
  }

  dash_dir <- tempfile("cosparta_dash_")
  dir.create(dash_dir)
  qmd_path <- file.path(dash_dir, "cosparta_dashboard.qmd")
  file.copy(template, qmd_path, overwrite = TRUE)

  saveRDS(
    list(fit = fit, channel_names = channel_names, time_labels = time_labels),
    file.path(dash_dir, "dash_data.rds")
  )

  message("Launching CoSparTA dashboard (this blocks the console until the ",
          "dashboard is stopped, like shiny::runApp())...")
  quarto::quarto_preview(qmd_path)

  invisible(qmd_path)
}
