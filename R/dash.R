#' Launch an Interactive Fit-Explorer Dashboard
#'
#' @description
#' Launches a local, interactive Quarto dashboard for exploring a fitted
#' \code{\link{CoSparTA}} / \code{\link{CoSparTA_missing}} object, or a set
#' of raw (already normalized) factor matrices. The dashboard shows the
#' temporal factor pattern (\code{\link{plot_time_factors}}), the channel
#' factor pattern (\code{\link{plot_channel_factors}}), a component-weight
#' chart, and a table of covariate (gamma) coefficients, all restricted to
#' a set of components selected via checkboxes.
#'
#' @param fit A fitted object returned by \code{\link{CoSparTA}} or
#'   \code{\link{CoSparTA_missing}}. Either \code{fit} or all of \code{Ef},
#'   \code{Ew}, and \code{lambda} must be supplied. \code{fit$res$ql$El_normed}
#'   (per-observation loadings) is never read — the dashboard has no
#'   observation-mode panel.
#' @param channel_names Character vector of channel labels, forwarded to
#'   \code{\link{plot_channel_factors}}. Fixed at launch; not adjustable from
#'   the dashboard. Default \code{NULL}.
#' @param time_labels Numeric vector of time-axis labels, forwarded to
#'   \code{\link{plot_time_factors}}. Fixed at launch; not adjustable from
#'   the dashboard. Default \code{NULL}.
#' @param channel_groups Character vector of group labels for each channel,
#'   forwarded to \code{\link{plot_channel_factors}}'s \code{channel_groups}
#'   argument (e.g. \code{c(rep("TextEmo", 5), rep("Gaze", 25))}). Fixed at
#'   launch; not adjustable from the dashboard. Default \code{NULL}.
#' @param covariate_names Names for the gamma table's covariate columns.
#'   When \code{intercept = TRUE} (the default), these describe every
#'   column except the always-present, always-first \code{"Intercept"}
#'   column. When \code{intercept = FALSE}, there is no \code{"Intercept"}
#'   column and these describe every column. Two forms are accepted,
#'   mirroring the \code{Xcov} convention used elsewhere in the package:
#'   \describe{
#'     \item{a character vector}{applied to every supervised component.
#'       Requires all supervised components to have the same number of
#'       covariates; errors with a clear message otherwise (use the list
#'       form instead when they differ).}
#'     \item{a list of length \code{K}}{per-component names.
#'       \code{covariate_names[[k]]} must have length
#'       \code{length(gamma_list[[k]]) - 1} when \code{intercept = TRUE},
#'       or exactly \code{length(gamma_list[[k]])} when
#'       \code{intercept = FALSE}. Entries for unsupervised components are
#'       ignored.}
#'   }
#'   Default \code{NULL} uses generic \code{"V1"}, \code{"V2"}, ... labels.
#'   Works with either input path, including the \code{fit} path.
#' @param Ef Numeric matrix of dimensions \code{p x K} (time factor matrix,
#'   already normalized — e.g. \code{fit$res$qf$Ef_normed} or
#'   \code{normalize_factors(fit)$Ef}). Supply together with \code{Ew} and
#'   \code{lambda} to bypass \code{fit}, mirroring the \code{Ef}/\code{Ew}
#'   raw-matrix path already supported by \code{\link{plot_time_factors}}
#'   and \code{\link{plot_channel_factors}}.
#' @param Ew Numeric matrix of dimensions \code{w x K} (channel factor
#'   matrix, already normalized). Supply together with \code{Ef} and
#'   \code{lambda} to bypass \code{fit}.
#' @param lambda Numeric vector of length \code{K} giving component weights
#'   (used by the component-weight panel). Supply together with \code{Ef}
#'   and \code{Ew} to bypass \code{fit}.
#' @param gamma_list List of length \code{K} of covariate coefficient
#'   vectors, one per component, in the same format returned by
#'   \code{\link{normalize_factors}}'s \code{gamma_list}: each element is
#'   either a numeric vector (whose first entry is the intercept when
#'   \code{intercept = TRUE}, or is itself a covariate coefficient when
#'   \code{intercept = FALSE}) or a scalar \code{NA} for a component with no
#'   covariates. \code{NULL} (the default) means every component is
#'   unsupervised. Only used when bypassing \code{fit}; ignored (and derived
#'   from \code{fit$res$gl_normed} instead) when \code{fit} is supplied.
#' @param intercept Logical. Whether each element of \code{gamma_list}
#'   begins with an intercept term. Default \code{TRUE}: \code{gamma_list[[k]][1]}
#'   is the intercept, \code{covariate_names[[k]]} describes the remaining
#'   entries, and the gamma table's first column is always labeled
#'   \code{"Intercept"}. Set to \code{FALSE} for fits where \code{Xcov} was
#'   supplied without an intercept column, so \code{gamma_list[[k]]} is
#'   entirely covariate coefficients: \code{covariate_names[[k]]} must then
#'   match its length exactly, and the gamma table has no \code{"Intercept"}
#'   column. If left at its default and \code{covariate_names} lengths match
#'   \code{gamma_list} lengths exactly rather than length \code{- 1}, this is
#'   unambiguous evidence of a no-intercept fit: \code{dash()} switches to
#'   \code{intercept = FALSE} automatically and warns, rather than erroring.
#'   Pass \code{intercept = FALSE} explicitly to silence that warning.
#'
#' @details
#' The dashboard is a Quarto document (\code{format: dashboard},
#' \code{server: shiny}) bundled with the package at
#' \code{inst/dashboard/cosparta_dashboard.qmd}. \code{dash()} resolves
#' \code{Ef}, \code{Ew}, \code{lambda}, and \code{gamma_list} — either
#' directly from the arguments of the same name, or extracted from
#' \code{fit} — saves them (together with \code{channel_names},
#' \code{time_labels}, \code{channel_groups}, and the resolved
#' \code{covariate_names}) to a temporary \code{.rds} file next to a copy
#' of the template, then calls \code{\link[quarto]{quarto_preview}}, which
#' starts a local Shiny-backed server and opens a browser tab.
#'
#' This call blocks the R console for as long as the dashboard is open —
#' the same behavior as \code{shiny::runApp()}. That is expected, not a
#' bug. Stop the server (e.g. Escape in RStudio, or Ctrl+C in a terminal)
#' to return control to the console.
#'
#' The only interactive control is a checkbox group for selecting which
#' components (ranks) to display; all four panels re-render from that
#' selection. \code{channel_names}, \code{time_labels}, \code{channel_groups},
#' and \code{covariate_names} are read once at launch and are not exposed
#' as dashboard inputs.
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
#'
#' # Raw-matrix path, bypassing fit (e.g. after normalize_factors())
#' nf <- normalize_factors(fit)
#' dash(Ef = nf$Ef, Ew = nf$Ew, lambda = nf$lambda, gamma_list = nf$gamma_list,
#'      covariate_names = c("age", "income"))
#'
#' # Per-component covariate names, differing counts across components
#' dash(fit, covariate_names = list(c("age", "income"), NULL, c("region")))
#'
#' # Fit without an intercept in Xcov: gamma_list[[k]] is all covariates
#' dash(Ef = nf$Ef, Ew = nf$Ew, lambda = nf$lambda, gamma_list = nf$gamma_list,
#'      covariate_names = c("age", "income"), intercept = FALSE)
#' }
#'
#' @seealso \code{\link{plot_time_factors}}, \code{\link{plot_channel_factors}},
#'   \code{\link{normalize_factors}}, \code{\link{CoSparTA}},
#'   \code{\link{CoSparTA_missing}}
#' @export
dash <- function(fit = NULL, channel_names = NULL, time_labels = NULL,
                  channel_groups = NULL, covariate_names = NULL,
                  Ef = NULL, Ew = NULL, lambda = NULL, gamma_list = NULL,
                  intercept = TRUE) {

  intercept_specified <- !missing(intercept)

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

  # --- resolve Ef / Ew / lambda / gamma_list, from raw matrices or fit ---
  raw_supplied <- !is.null(Ef) || !is.null(Ew) || !is.null(lambda)
  if (raw_supplied) {
    if (is.null(Ef) || is.null(Ew) || is.null(lambda)) {
      stop("When bypassing `fit`, all of 'Ef', 'Ew', and 'lambda' must be ",
           "supplied together.")
    }
    if (!is.matrix(Ef)) stop("'Ef' must be a matrix (p x K).")
    if (!is.matrix(Ew)) stop("'Ew' must be a matrix (w x K).")
    K <- ncol(Ef)
    if (ncol(Ew) != K) {
      stop(sprintf("'Ef' has %d components but 'Ew' has %d.", K, ncol(Ew)))
    }
    if (length(lambda) != K) {
      stop(sprintf("'lambda' has length %d but Ef/Ew have %d components.",
                    length(lambda), K))
    }
    if (!is.null(gamma_list) && length(gamma_list) != K) {
      stop(sprintf("'gamma_list' has length %d but Ef/Ew have %d components.",
                    length(gamma_list), K))
    }
  } else if (!is.null(fit)) {
    Ef <- fit$res$qf$Ef_normed
    if (is.null(Ef) || !is.matrix(Ef)) {
      stop(
        "`fit` does not look like a fitted CoSparTA / CoSparTA_missing ",
        "object: expected a matrix at fit$res$qf$Ef_normed. Pass the object ",
        "returned by CoSparTA() or CoSparTA_missing()."
      )
    }
    K  <- ncol(Ef)
    Ew <- fit$res$qw$Ew_normed
    if (is.null(Ew) || !is.matrix(Ew) || ncol(Ew) != K) {
      stop("`fit` is missing a valid fit$res$qw$Ew_normed (w x K) matrix.")
    }
    lambda <- fit$res$lambda_normed
    if (is.null(lambda) || length(lambda) != K) {
      stop("`fit` is missing fit$res$lambda_normed (component weights) of length K.")
    }
    gl <- fit$res$gl_normed
    if (is.null(gl) || length(gl) != K) {
      stop("`fit` is missing fit$res$gl_normed (gamma estimates) of length K.")
    }
    # Translate fit's internal $gamma/$type structure into the same
    # (numeric vector | scalar NA) convention normalize_factors()$gamma_list
    # uses, so the dashboard template only ever has to understand one
    # representation regardless of which input path was used.
    gamma_list <- lapply(gl, function(g) {
      if (!is.null(g) && !is.null(g$type) && identical(g$type, "covariate_dependent")) {
        g$gamma
      } else {
        NA
      }
    })
  } else {
    stop("Either 'fit' or all of 'Ef', 'Ew', and 'lambda' must be provided.")
  }

  if (!is.null(channel_names) && length(channel_names) != nrow(Ew)) {
    stop(sprintf("`channel_names` has length %d but there are %d channels.",
                  length(channel_names), nrow(Ew)))
  }
  if (!is.null(time_labels) && length(time_labels) != nrow(Ef)) {
    stop(sprintf("`time_labels` has length %d but there are %d time points.",
                  length(time_labels), nrow(Ef)))
  }
  if (!is.null(channel_groups) && length(channel_groups) != nrow(Ew)) {
    stop(sprintf("`channel_groups` has length %d but there are %d channels.",
                  length(channel_groups), nrow(Ew)))
  }

  resolved <- .dash_resolve_covariate_names(covariate_names, gamma_list, K,
                                              intercept, intercept_specified)
  covariate_names_resolved <- resolved$names
  intercept                <- resolved$intercept

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
    list(
      Ef = Ef, Ew = Ew, lambda = lambda, gamma_list = gamma_list,
      channel_names = channel_names, time_labels = time_labels,
      channel_groups = channel_groups, covariate_names = covariate_names_resolved,
      intercept = intercept
    ),
    file.path(dash_dir, "dash_data.rds")
  )

  message("Launching CoSparTA dashboard (this blocks the console until the ",
          "dashboard is stopped, like shiny::runApp())...")
  quarto::quarto_preview(qmd_path)

  invisible(qmd_path)
}

#' @keywords internal
.dash_gamma_len <- function(g) {
  if (is.null(g) || (length(g) == 1L && is.na(g))) 0L else length(g)
}

# Validates and normalizes covariate_names (either a single character vector
# applied to every supervised component, or a list of length K giving
# per-component names) into a length-K list where element k is NULL (no
# names / unsupervised) or a character vector describing component k's
# covariates. When `intercept` is TRUE, gamma_list[[k]][1] is an intercept
# term and the expected name length is length(gamma_list[[k]]) - 1; when
# FALSE, gamma_list[[k]] is entirely covariate coefficients and the expected
# name length is length(gamma_list[[k]]) exactly.
#
# Returns list(names = <resolved per-component name list>, intercept =
# <possibly auto-detected intercept flag>). When the caller did not
# explicitly pass `intercept` (intercept_specified == FALSE) and every
# supplied name vector's length matches gamma_list lengths exactly (never
# length - 1), that unambiguously means the fit has no intercept term, so
# `intercept` is switched to FALSE (with a warning) instead of erroring.
#' @keywords internal
.dash_resolve_covariate_names <- function(covariate_names, gamma_list, K,
                                            intercept, intercept_specified) {
  if (is.null(covariate_names)) {
    return(list(names = vector("list", K), intercept = intercept))
  }

  gamma_lens <- if (is.null(gamma_list)) {
    rep(0L, K)
  } else {
    vapply(gamma_list, .dash_gamma_len, integer(1))
  }
  supervised <- gamma_lens > 0L

  no_intercept_msg <- function(observed, expected, where) {
    paste0(where, " has length ", observed, " but expected ", expected, " ",
           if (intercept) {
             "non-intercept covariate(s) (gamma length minus 1 for the intercept)."
           } else {
             "covariate(s) (gamma length, since `intercept = FALSE`)."
           })
  }

  auto_detect_warning <- function() {
    warning(
      "`covariate_names` length(s) match `gamma_list` length(s) exactly ",
      "(not length - 1), implying these components were fit without an ",
      "intercept term. Proceeding with `intercept = FALSE`. Pass ",
      "`intercept = FALSE` explicitly to silence this warning.",
      call. = FALSE
    )
  }

  if (is.list(covariate_names)) {
    if (length(covariate_names) != K) {
      stop(sprintf(
        "`covariate_names` is a list of length %d but there are %d components; ",
        length(covariate_names), K),
        "supply one element per component."
      )
    }

    if (!intercept_specified) {
      named_k <- which(supervised & !vapply(covariate_names, is.null, logical(1)))
      if (length(named_k) > 0L) {
        matches_with_icpt <- all(vapply(named_k, function(k) {
          length(covariate_names[[k]]) == gamma_lens[k] - 1L
        }, logical(1)))
        matches_without_icpt <- all(vapply(named_k, function(k) {
          length(covariate_names[[k]]) == gamma_lens[k]
        }, logical(1)))
        if (!matches_with_icpt && matches_without_icpt) {
          intercept <- FALSE
          auto_detect_warning()
        }
      }
    }

    out <- vector("list", K)
    for (k in which(supervised)) {
      nm <- covariate_names[[k]]
      if (is.null(nm)) next
      expected <- if (intercept) gamma_lens[k] - 1L else gamma_lens[k]
      if (length(nm) != expected) {
        stop(no_intercept_msg(length(nm), expected,
                               sprintf("`covariate_names[[%d]]`", k)))
      }
      out[[k]] <- nm
    }
    return(list(names = out, intercept = intercept))
  }

  if (!is.character(covariate_names)) {
    stop("`covariate_names` must be a character vector or a list of length K.")
  }

  supervised_lens <- unique(gamma_lens[supervised])
  if (length(supervised_lens) > 1L) {
    stop(
      "`covariate_names` was supplied as a single character vector, but ",
      "components have different numbers of covariates (",
      paste(sort(supervised_lens), collapse = ", "),
      "). Supply `covariate_names` as a list of length K instead, with ",
      "one name vector per component."
    )
  }
  if (length(supervised_lens) == 1L) {
    if (!intercept_specified &&
        length(covariate_names) != supervised_lens - 1L &&
        length(covariate_names) == supervised_lens) {
      intercept <- FALSE
      auto_detect_warning()
    }
    expected <- if (intercept) supervised_lens - 1L else supervised_lens
    if (length(covariate_names) != expected) {
      stop(no_intercept_msg(length(covariate_names), expected,
                             "`covariate_names`"))
    }
  }
  list(names = lapply(supervised, function(s) if (s) covariate_names else NULL),
       intercept = intercept)
}
