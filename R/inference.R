#' Extract Posterior Inclusion Probabilities from a CoSparTA fit
#'
#' @description
#' Extracts the posterior inclusion probability (PIP) matrices from a fitted
#' \code{\link{CoSparTA}} object. PIPs give the probability that each factor
#' element is truly non-zero vs. noise, derived from the spike-and-slab prior.
#' A PIP close to 1 indicates strong evidence for a non-zero loading; a PIP
#' close to 0 indicates the loading is likely noise.
#'
#' @param fit A fitted object returned by \code{\link{CoSparTA}}.
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
#' fit <- CoSparTA(X, K = 3, Xcov = Xcov)
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


#' Compute Credible Intervals from a CoSparTA fit
#'
#' @description
#' Computes approximate credible intervals for factor elements using the
#' posterior mean and variance stored in a fitted \code{\link{CoSparTA}}
#' object. Intervals are computed under a normal approximation to the
#' posterior, which is reasonable for elements with non-negligible PIP.
#' For elements with low PIP (likely zero), intervals should be interpreted
#' cautiously.
#'
#' @param fit A fitted object returned by \code{\link{CoSparTA}}.
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
#' fit <- CoSparTA(X, K = 3, Xcov = Xcov)
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
#' @param fit A fitted object returned by \code{\link{CoSparTA}}.
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
#' fit <- CoSparTA(X, K = 3, Xcov = Xcov)
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


#' Compute Exact Posterior Quantiles from a CoSparTA fit
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
#' @param fit A fitted object returned by \code{\link{CoSparTA}} or
#'   \code{\link{CoSparTA_missing}}.
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
#' fit <- CoSparTA(X, K = 3, Xcov = Xcov)
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

  if (is.null(shape_mat) || is.null(rate_mat) ||
      all(is.na(shape_mat)) || all(is.na(rate_mat))) {
    stop("Posterior Gamma parameters not available. Refit the model with the current package version.")
  }

  pi_hat_mat <- 1 - pip_mat

  out <- lapply(probs, function(tau) {
    q_mat <- matrix(0, nrow = nrow(pip_mat), ncol = ncol(pip_mat))
    slab_idx <- !is.na(pi_hat_mat) & tau > pi_hat_mat & pi_hat_mat < 1
    if (any(slab_idx)) {
      adjusted_prob <- (tau - pi_hat_mat[slab_idx]) / (1 - pi_hat_mat[slab_idx])
      q_mat[slab_idx] <- qgamma(adjusted_prob,
                                 shape = shape_mat[slab_idx],
                                 rate = rate_mat[slab_idx])
    }
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
#' Confidence Intervals for Covariate Coefficients (Gamma)
#'
#' @description
#' Computes confidence intervals for the covariate effect parameters
#' \eqn{\gamma} estimated by \code{\link{CoSparTA}}. Two methods are available:
#'
#' \describe{
#'   \item{\code{"delta"}}{Fast asymptotic intervals based on the Hessian of
#'     the marginal log-likelihood from the EBPM optimization. The full
#'     Hessian (covering \eqn{\gamma} together with nuisance parameters
#'     \eqn{\pi_0}, \eqn{\alpha}, \eqn{\beta}) is inverted first, and then
#'     the \eqn{\gamma} submatrix is extracted, implementing
#'     \eqn{\mathrm{Var}(\hat\gamma_k) = [H_k^{-1}]_{\gamma,\gamma}} per
#'     equation (21) of the paper. This correctly accounts for uncertainty in
#'     the nuisance parameters. These are \strong{conditional} standard errors
#'     -- they measure uncertainty in \eqn{\gamma} given the current F and W
#'     estimates, and may underestimate the true uncertainty. Recommended for
#'     exploratory screening.}
#'   \item{\code{"bootstrap"}}{Parametric bootstrap intervals. Generates
#'     \code{B} synthetic tensors from the fitted Poisson rates, refits
#'     \code{\link{CoSparTA}} on each, and collects the bootstrap distribution
#'     of \eqn{\gamma}. This correctly accounts for all sources of estimation
#'     uncertainty (including uncertainty in F and W). Recommended for
#'     publication-quality inference. Computationally expensive: requires
#'     \code{B} full model refits.}
#' }
#'
#' @param fit A fitted object returned by \code{\link{CoSparTA}}.
#' @param method Character: \code{"delta"} or \code{"bootstrap"}.
#'   Default \code{"bootstrap"}.
#' @param level Numeric confidence level in \code{(0, 1)}. Default \code{0.95}.
#' @param B Integer. Number of bootstrap replicates (only used when
#'   \code{method = "bootstrap"}). Default \code{200}.
#' @param X The original tensor used to fit the model (required for
#'   \code{method = "bootstrap"}).
#' @param K Integer. Number of components (required for bootstrap).
#' @param Xcov Covariate matrix or list used in the original fit (required
#'   for bootstrap).
#' @param init_fn Optional function with signature
#'   \code{function(X_star, K)} that returns an initialization object
#'   (e.g., a list of three matrices from \code{\link{init_cpapr}}).
#'   If provided, this function is called on each bootstrap replicate
#'   to generate a fresh initialization. If \code{NULL}, the
#'   \code{init} argument from \code{...} is reused for all
#'   replicates. Using per-replicate initialization (e.g., via
#'   \code{init_cpapr}) is recommended for publication-quality
#'   bootstrap inference.
#' @param verbose Logical. If \code{TRUE}, prints bootstrap progress.
#'   Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{CoSparTA}} during
#'   bootstrap refitting (e.g., \code{maxiter}, \code{init},
#'   \code{convergence_criteria}).
#'
#' @return A list of length K (one element per component). Each element is
#'   a named list with:
#' \describe{
#'   \item{estimate}{Numeric vector of gamma estimates from the original fit.}
#'   \item{se}{Standard errors.}
#'   \item{lower}{Lower confidence bound.}
#'   \item{upper}{Upper confidence bound.}
#'   \item{pvalue}{Two-sided p-values (Wald test for delta method;
#'     bootstrap percentile-based for bootstrap, computed as
#'     \code{2 * min(prop >= 0, prop <= 0)}).}
#'   \item{method}{Character string indicating which method was used.}
#' }
#' If a component was fitted without covariates (unsupervised rank),
#' its entry is \code{NULL}.
#'
#' @examples
#' \dontrun{
#' fit <- CoSparTA(X, K = 3, Xcov = Xcov, maxiter = 50)
#'
#' # Fast delta method (exploratory)
#' ci_delta <- get_gamma_ci(fit, method = "delta")
#'
#' # Parametric bootstrap (publication quality)
#' ci_boot <- get_gamma_ci(fit, method = "bootstrap", B = 200,
#'                          X = X, K = 3, Xcov = Xcov,
#'                          maxiter = 50, convergence_criteria = "ELBO")
#' }
#'
#' @seealso \code{\link{CoSparTA}},
#'   \code{\link{ebpm_point_gamma_multiplier_covariates}}
#' @export
get_gamma_ci <- function(fit, method = "bootstrap", level = 0.95,
                          B = 200, X = NULL, K = NULL, Xcov = NULL,
                          init_fn = NULL, verbose = TRUE, ...) {

  method <- match.arg(method, c("delta", "bootstrap"))
  gl     <- fit$res$gl
  K_fit  <- length(gl)

  # ------------------------------------------------------------------ #
  #  Delta method                                                        #
  # ------------------------------------------------------------------ #
  if (method == "delta") {
    results <- vector("list", K_fit)

    for (k in seq_len(K_fit)) {
      if (is.null(gl[[k]]) || isTRUE(gl[[k]]$type != "covariate_dependent")) {
        results[[k]] <- NULL
        next
      }

      gamma_est <- gl[[k]]$gamma
      q         <- length(gamma_est)
      hessian   <- gl[[k]]$hessian

      if (is.null(hessian)) {
        stop("Hessian not available for component ", k, ". ",
             "Refit the model with the current package version.")
      }

      # Invert the full Hessian first, then extract the gamma submatrix.
      # This implements Var(gamma_hat_k) = [H_k^{-1}]_{gamma,gamma} per
      # eq(21), correctly propagating uncertainty from nuisance parameters
      # (positions 1-3: logit_pi0, log_alpha, log_beta).
      gamma_idx  <- 4:(3 + q)
      vcov_gamma <- tryCatch({
        V_full <- solve(hessian)
        V_full[gamma_idx, gamma_idx, drop = FALSE]
      }, error = function(e) {
        warning(sprintf("Hessian inversion failed for component %d: %s. ",
                        k, conditionMessage(e)),
                "Standard errors set to NA.")
        NULL
      })

      if (is.null(vcov_gamma)) {
        se <- rep(NA_real_, q)
      } else {
        se <- sqrt(pmax(diag(vcov_gamma), 0))
      }

      z     <- qnorm((1 + level) / 2)
      lower  <- gamma_est - z * se
      upper  <- gamma_est + z * se
      pvalue <- 2 * pnorm(-abs(gamma_est / se))

      results[[k]] <- list(
        estimate = gamma_est,
        se       = se,
        lower    = lower,
        upper    = upper,
        pvalue   = pvalue,
        method   = "delta"
      )
    }

    return(results)
  }

  # ------------------------------------------------------------------ #
  #  Bootstrap method                                                    #
  # ------------------------------------------------------------------ #
  if (is.null(X) || is.null(K) || is.null(Xcov)) {
    stop("'X', 'K', and 'Xcov' must be provided for method = 'bootstrap'.")
  }

  lambda_hat <- reconstruct_tensor(fit)

  # Capture extra args and enforce adj_LF_scale = FALSE, verbose = FALSE
  extra_args <- list(...)
  if (is.null(extra_args$adj_LF_scale)) extra_args$adj_LF_scale <- FALSE
  extra_args$verbose <- FALSE

  # Determine q for each rank; initialise storage
  gamma_boot <- vector("list", K_fit)
  for (k in seq_len(K_fit)) {
    if (!is.null(gl[[k]]) && isTRUE(gl[[k]]$type == "covariate_dependent")) {
      q <- length(gl[[k]]$gamma)
      gamma_boot[[k]] <- matrix(NA_real_, nrow = B, ncol = q)
    }
  }

  nf_orig <- normalize_factors(fit)

  for (b in seq_len(B)) {
    if (verbose && b %% 10 == 0) {
      cat(sprintf("Bootstrap %d/%d\n", b, B))
    }

    X_star <- array(rpois(length(lambda_hat), lambda_hat), dim = dim(lambda_hat))

    # Per-replicate initialization if init_fn provided
    if (!is.null(init_fn)) {
      boot_init <- tryCatch(
        init_fn(X_star, K),
        error = function(e) {
          warning(sprintf("Bootstrap %d: init_fn failed (%s), using fallback.",
                          b, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(boot_init)) {
        extra_args$init <- boot_init
      }
    }

    fit_star <- tryCatch(
      do.call(CoSparTA, c(list(X = X_star, K = K, Xcov = Xcov), extra_args)),
      error = function(e) {
        warning(sprintf("Bootstrap replicate %d failed: %s -- row filled with NA.", b,
                        conditionMessage(e)))
        NULL
      }
    )

    if (is.null(fit_star)) next

    nf_star <- normalize_factors(fit_star)
    mf_b <- match_factors(
      ref = list(nf_orig$El, nf_orig$Ef, nf_orig$Ew),
      est = list(nf_star$El, nf_star$Ef, nf_star$Ew)
    )
    perm_b <- mf_b$permutation

    for (k in seq_len(K_fit)) {
      k_boot <- perm_b[k]
      if (!is.null(gamma_boot[[k]]) &&
          !is.null(fit_star$res$gl[[k_boot]]) &&
          isTRUE(fit_star$res$gl[[k_boot]]$type == "covariate_dependent")) {
        gamma_k <- fit_star$res$gl[[k_boot]]$gamma
        if (length(gamma_k) == ncol(gamma_boot[[k]])) {
          gamma_boot[[k]][b, ] <- gamma_k
        }
      }
    }
  }

  # Summarise bootstrap distributions
  results <- vector("list", K_fit)
  for (k in seq_len(K_fit)) {
    if (is.null(gamma_boot[[k]])) {
      results[[k]] <- NULL
      next
    }

    gamma_est <- gl[[k]]$gamma
    q         <- length(gamma_est)
    boot_mat  <- gamma_boot[[k]]

    se    <- apply(boot_mat, 2, sd,       na.rm = TRUE)
    lower <- apply(boot_mat, 2, quantile, probs = (1 - level) / 2, na.rm = TRUE)
    upper <- apply(boot_mat, 2, quantile, probs = (1 + level) / 2, na.rm = TRUE)

    pvalue <- vapply(seq_len(q), function(j) {
      col_j <- boot_mat[, j]
      2 * min(mean(col_j >= 0, na.rm = TRUE),
              mean(col_j <= 0, na.rm = TRUE))
    }, numeric(1))

    results[[k]] <- list(
      estimate = gamma_est,
      se       = se,
      lower    = lower,
      upper    = upper,
      pvalue   = pvalue,
      method   = "bootstrap"
    )
  }

  results
}


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
