#' Extract Posterior Inclusion Probabilities from a CxtEBTD fit
#'
#' @description
#' Extracts the posterior inclusion probability (PIP) matrices from a fitted
#' \code{\link{CxtEBTD}} object. PIPs give the probability that each factor
#' element is truly non-zero vs. noise, derived from the spike-and-slab prior.
#' A PIP close to 1 indicates strong evidence for a non-zero loading; a PIP
#' close to 0 indicates the loading is likely noise.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD}}.
#' @param mode Character string specifying which mode to extract:
#'   \code{'L'} for observation loadings, \code{'F'} for time factors,
#'   \code{'W'} for channel weights. Default \code{'L'}.
#' @param threshold Numeric value in \code{[0, 1]}. If supplied, returns a
#'   logical matrix indicating which elements exceed the threshold. If
#'   \code{NULL}, returns the raw PIP matrix. Default \code{NULL}.
#'
#' @return If \code{threshold = NULL}, a numeric matrix of PIPs with the same
#'   dimensions as the corresponding factor matrix (\code{n x K} for L,
#'   \code{p x K} for F, \code{w x K} for W). If \code{threshold} is
#'   supplied, a logical matrix of the same dimensions where \code{TRUE}
#'   indicates PIP > threshold. Returns \code{NULL} with a warning if PIPs
#'   are not available for the requested mode (e.g. when \code{ebps_with_uq}
#'   is used for F, which has no spike component).
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X, K = 3, Xcov = Xcov)
#'
#' # Raw PIP matrix for observation loadings
#' pip_L <- get_pip(fit, mode = 'L')
#'
#' # Which loadings exceed 0.9 PIP threshold
#' sig_L <- get_pip(fit, mode = 'L', threshold = 0.9)
#' }
#'
#' @export
get_pip <- function(fit, mode = 'L', threshold = NULL) {

  pip_mat <- switch(mode,
    'L' = fit$res$ql$PIPl,
    'F' = fit$res$qf$PIPf,
    'W' = fit$res$qw$PIPw,
    stop("mode must be one of 'L', 'F', or 'W'")
  )

  if (is.null(pip_mat) || all(is.na(pip_mat))) {
    warning(paste("PIPs not available for mode", mode,
                  "-- the prior used for this mode has no spike component."))
    return(NULL)
  }

  if (!is.null(threshold)) {
    if (threshold < 0 || threshold > 1) {
      stop("threshold must be between 0 and 1")
    }
    return(pip_mat > threshold)
  }

  return(pip_mat)
}


#' Compute Credible Intervals from a CxtEBTD fit
#'
#' @description
#' Computes approximate credible intervals for factor elements using the
#' posterior mean and variance stored in a fitted \code{\link{CxtEBTD}}
#' object. Intervals are computed under a normal approximation to the
#' posterior, which is reasonable for elements with non-negligible PIP.
#' For elements with low PIP (likely zero), intervals should be interpreted
#' cautiously.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD}}.
#' @param mode Character string specifying which mode to extract:
#'   \code{'L'} for observation loadings, \code{'F'} for time factors,
#'   \code{'W'} for channel weights. Default \code{'L'}.
#' @param level Numeric credible level in \code{(0, 1)}. Default \code{0.95}
#'   for 95\% credible intervals.
#'
#' @return A list with three matrices of the same dimensions as the
#'   corresponding factor matrix:
#' \describe{
#'   \item{mean}{Posterior mean matrix.}
#'   \item{lower}{Lower bound of the credible interval.}
#'   \item{upper}{Upper bound of the credible interval.}
#' }
#' Returns \code{NULL} with a warning if posterior variance is not available
#' for the requested mode.
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X, K = 3, Xcov = Xcov)
#'
#' # 95% credible intervals for observation loadings
#' ci_L <- get_credible_interval(fit, mode = 'L', level = 0.95)
#' ci_L$mean    # posterior means
#' ci_L$lower   # lower bounds
#' ci_L$upper   # upper bounds
#'
#' # 90% credible intervals for channel weights
#' ci_W <- get_credible_interval(fit, mode = 'W', level = 0.90)
#' }
#'
#' @export
get_credible_interval <- function(fit, mode = 'L', level = 0.95) {

  if (level <= 0 || level >= 1) {
    stop("level must be between 0 and 1")
  }

  mean_mat <- switch(mode,
    'L' = fit$res$ql$El,
    'F' = fit$res$qf$Ef,
    'W' = fit$res$qw$Ew,
    stop("mode must be one of 'L', 'F', or 'W'")
  )

  var_mat <- switch(mode,
    'L' = fit$res$ql$Varl,
    'F' = fit$res$qf$Varf,
    'W' = fit$res$qw$Varw
  )

  if (is.null(var_mat) || all(is.na(var_mat))) {
    warning(paste("Posterior variance not available for mode", mode,
                  "-- credible intervals cannot be computed."))
    return(NULL)
  }

  z <- qnorm((1 + level) / 2)
  sd_mat <- sqrt(var_mat)

  return(list(
    mean  = mean_mat,
    lower = pmax(0, mean_mat - z * sd_mat),  # non-negative since loadings >= 0
    upper = mean_mat + z * sd_mat
  ))
}

#' Identify Significant Channels and Time Points per Factor via lFDR Control
#'
#' @description
#' Implements Bayesian local-FDR control (Algorithms 1 and 2 from the EBTD
#' framework) to identify which channels and time points are truly active for
#' each factor. For each factor \code{l}, the posterior spike probabilities
#' (local-fdr values) from the F and W modes are sorted and thresholded to
#' produce a discovery set at a controlled FDR level.
#'
#' The local-fdr value for entry \eqn{(j, l)} is
#' \eqn{\rho_{jl} = 1 - \text{PIP}_{jl}}, i.e., the posterior probability
#' that the entry is truly zero. Algorithm 1 finds the largest discovery set
#' such that the mean local-fdr does not exceed \code{alpha}.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD}}.
#' @param alpha Numeric FDR level in \code{(0, 1)}. Default \code{0.05}.
#' @param mode Character string specifying which mode(s) to run discovery on:
#'   \code{'F'} for time factors, \code{'W'} for channel weights, or
#'   \code{'both'}. Default \code{'both'}.
#'
#' @return A list of length \code{K} (one element per factor). Each element
#'   is a named list with:
#' \describe{
#'   \item{factor}{Integer factor index.}
#'   \item{active_times}{Integer vector of active time point indices for this
#'     factor (from F mode). \code{NULL} if mode is \code{'W'} or PIPs
#'     unavailable.}
#'   \item{active_channels}{Integer vector of active channel indices for this
#'     factor (from W mode). \code{NULL} if mode is \code{'F'} or PIPs
#'     unavailable.}
#'   \item{n_active_times}{Number of active time points discovered.}
#'   \item{n_active_channels}{Number of active channels discovered.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X, K = 3, Xcov = Xcov)
#'
#' # Discover active channels and time points at 5% FDR
#' patterns <- get_significant_patterns(fit, alpha = 0.05)
#'
#' # Factor 1 active time points
#' patterns[[1]]$active_times
#'
#' # Factor 2 active channels
#' patterns[[2]]$active_channels
#' }
#'
#' @export
get_significant_patterns <- function(fit, alpha = 0.05, mode = 'both') {

  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (!mode %in% c('F', 'W', 'both')) stop("mode must be 'F', 'W', or 'both'")

  K <- ncol(fit$res$ql$El)
  results <- vector("list", K)

  for (l in 1:K) {

    active_times    <- NULL
    active_channels <- NULL

    # -- F mode: time points --
    if (mode %in% c('F', 'both')) {
      pip_f <- fit$res$qf$PIPf
      if (!is.null(pip_f) && !all(is.na(pip_f[, l]))) {
        active_times <- lfdr_discovery(1 - pip_f[, l], alpha)
      } else {
        warning(paste("PIPs not available for F mode, factor", l))
      }
    }

    # -- W mode: channels --
    if (mode %in% c('W', 'both')) {
      pip_w <- fit$res$qw$PIPw
      if (!is.null(pip_w) && !all(is.na(pip_w[, l]))) {
        active_channels <- lfdr_discovery(1 - pip_w[, l], alpha)
      } else {
        warning(paste("PIPs not available for W mode, factor", l))
      }
    }

    results[[l]] <- list(
      factor           = l,
      active_times     = active_times,
      active_channels  = active_channels,
      n_active_times   = length(active_times),
      n_active_channels = length(active_channels)
    )
  }

  return(results)
}


#' Compute Exact Posterior Quantiles from a CxtEBTD fit
#'
#' @description
#' Computes exact quantiles of the marginal posterior distribution for each
#' factor element, using the closed-form spike-and-slab structure of the
#' point-gamma posterior. The marginal CDF is:
#' \deqn{F(\theta) = \hat\pi_i \cdot \mathbf{1}[\theta \geq 0] +
#'   (1 - \hat\pi_i) \cdot F_{\text{Gamma}}(\theta)}
#' so the quantile at probability \eqn{\tau} is 0 when
#' \eqn{\tau \leq \hat\pi_i}, and otherwise
#' \eqn{F_{\text{Gamma}}^{-1}\!\left((\tau - \hat\pi_i)\,/\,(1 - \hat\pi_i)\right)}
#' where the gamma is parameterized by the posterior shape and rate.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD}} or
#'   \code{\link{CxtEBTD_missing}}.
#' @param probs Numeric vector of probabilities in \code{[0, 1]} at which to
#'   evaluate the quantile function. Default \code{c(0.025, 0.975)}.
#' @param mode Character string specifying which mode to extract:
#'   \code{'L'} for observation loadings, \code{'F'} for time factors,
#'   \code{'W'} for channel weights. Default \code{'L'}.
#'
#' @return A named list with one matrix per entry of \code{probs}. Each
#'   matrix has the same dimensions as the corresponding factor matrix
#'   (\code{n x K} for L, \code{p x K} for F, \code{w x K} for W). Names
#'   are formatted as \code{paste0("q", round(probs * 100, 1))}, e.g.
#'   \code{"q2.5"} and \code{"q97.5"} for the default \code{probs}.
#'
#' @seealso \code{\link{get_credible_interval}}, \code{\link{get_pip}}
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X, K = 3, Xcov = Xcov)
#'
#' # 95% equal-tailed posterior intervals for observation loadings
#' q_L <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = 'L')
#' q_L$q2.5    # lower bound matrix (n x K)
#' q_L$q97.5   # upper bound matrix (n x K)
#'
#' # Posterior median for channel weights
#' q_W <- get_posterior_quantile(fit, probs = 0.5, mode = 'W')
#' q_W$q50
#' }
#'
#' @export
get_posterior_quantile <- function(fit, probs = c(0.025, 0.975), mode = 'L') {

  pip_mat <- switch(mode,
    'L' = fit$res$ql$PIPl,
    'F' = fit$res$qf$PIPf,
    'W' = fit$res$qw$PIPw,
    stop("mode must be one of 'L', 'F', or 'W'")
  )

  shape_mat <- switch(mode,
    'L' = fit$res$ql$shape_post_l,
    'F' = fit$res$qf$shape_post_f,
    'W' = fit$res$qw$shape_post_w
  )

  rate_mat <- switch(mode,
    'L' = fit$res$ql$rate_post_l,
    'F' = fit$res$qf$rate_post_f,
    'W' = fit$res$qw$rate_post_w
  )

  if (is.null(shape_mat) || is.null(rate_mat)) {
    stop("Posterior Gamma parameters not available. Refit the model with the current package version.")
  }

  pi_hat_mat <- 1 - pip_mat

  out <- lapply(probs, function(tau) {
    q_mat <- ifelse(
      tau <= pi_hat_mat,
      0,
      qgamma((tau - pi_hat_mat) / (1 - pi_hat_mat),
             shape = shape_mat, rate = rate_mat)
    )
    q_mat
  })

  names(out) <- paste0("q", round(probs * 100, 1))
  out
}


#' Algorithm 1: LFDR Discovery Set
#'
#' @description
#' Given a vector of local-fdr values, finds the largest discovery set such
#' that the mean local-fdr does not exceed alpha.
#'
#' @param lfdr_vals Numeric vector of local-fdr values in \code{[0, 1]}.
#' @param alpha Numeric FDR level.
#'
#' @return Integer vector of indices in the discovery set (in original order).
#'
#' @keywords internal
lfdr_discovery <- function(lfdr_vals, alpha) {

  M <- length(lfdr_vals)
  if (M == 0) return(integer(0))

  # sort by local-fdr ascending (most likely signal first)
  sorted_idx <- order(lfdr_vals)
  sorted_lfdr <- lfdr_vals[sorted_idx]

  # find k* = max k such that mean of top k <= alpha
  cumulative_mean <- cumsum(sorted_lfdr) / seq_len(M)
  k_star <- max(c(0, which(cumulative_mean <= alpha)))

  if (k_star == 0) return(integer(0))

  # return original indices of discovery set
  sort(sorted_idx[1:k_star])
}
