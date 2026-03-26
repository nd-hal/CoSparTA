#' Empirical Bayes Point-Gamma with Posterior Variance and PIP
#'
#' @description
#' A wrapper around \code{\link[ebpm]{ebpm_point_gamma}} that extends the
#' posterior output with posterior variance and posterior inclusion
#' probabilities (PIPs). All optimization is delegated to the original
#' function; only the posterior summary is augmented.
#'
#' @param x Non-negative integer vector of observed counts.
#' @param s Numeric scalar or vector of length \code{length(x)} giving
#'   per-observation exposure/scaling factors. Default \code{1}.
#' @param ... Additional arguments passed to
#'   \code{\link[ebpm]{ebpm_point_gamma}} (e.g. \code{g_init},
#'   \code{fix_g}, \code{pi0}, \code{control}).
#'
#' @return Same structure as \code{\link[ebpm]{ebpm_point_gamma}} but with
#'   two additional columns in \code{posterior}:
#' \describe{
#'   \item{var}{Posterior variance \eqn{\text{Var}(\theta_i \mid x_i)},
#'     decomposed as within-component variance plus between-component
#'     variance from spike uncertainty.}
#'   \item{pip}{Posterior inclusion probability
#'     \eqn{P(\theta_i \neq 0 \mid x_i) = 1 - \hat{\pi}_i}.}
#' }
#'
#' @seealso \code{\link[ebpm]{ebpm_point_gamma}},
#'   \code{\link{ebpm_point_gamma_multiplier_covariates}}
#'
#' @export
ebpm_point_gamma_with_uq <- function(x, s = 1, ...) {

  if (length(s) == 1) s <- rep(s, length(x))

  # delegate all optimization to the original function
  fit <- ebpm::ebpm_point_gamma(x, s, ...)

  # extract fitted hyperparameters
  pi0   <- fit$fitted_g$pi0
  a     <- fit$fitted_g$shape
  b     <- 1 / fit$fitted_g$scale  # convert scale to rate

  # recompute pi_hat (mirrors logic in ebpm_point_gamma)
  nb_prob <- b / (b + s)
  nb      <- exp(dnbinom_cts_log_vec(x, a, prob = nb_prob))

  if (pi0 == 0) {
    pi_hat <- rep(0, length(x))
  } else {
    pi_hat <- pi0 * as.integer(x == 0) /
              (pi0 * as.integer(x == 0) + (1 - pi0) * nb)
  }

  # posterior variance (law of total variance)
  post_mean_nospike <- (a + x) / (b + s)
  var_within   <- (1 - pi_hat) * post_mean_nospike / (b + s)
  var_between  <- pi_hat * (1 - pi_hat) * post_mean_nospike^2
  mu_var       <- var_within + var_between

  # augment posterior
  fit$posterior$var <- mu_var
  fit$posterior$pip <- 1 - pi_hat

  return(fit)
}

#' Empirical Bayes Poisson Smoothing with Posterior Variance
#'
#' @description
#' A wrapper around \code{\link[smashrgen]{ebps}} that augments the posterior
#' output with posterior variance when available. Posterior variance is
#' computed from the log-normal formula using the variational posterior
#' parameters returned by \code{ebps}. Note that \code{ebps} uses a
#' Gaussian smoothing prior with no spike component, so posterior inclusion
#' probabilities (PIPs) are not applicable and returned as \code{NA}.
#'
#' @param x Non-negative integer vector of observed counts.
#' @param s Numeric scalar or vector of scaling factors. Default \code{NULL}.
#' @param ... Additional arguments passed to \code{\link[smashrgen]{ebps}}
#'   (e.g. \code{g_init}, \code{general_control}, \code{smooth_control}).
#'
#' @return Same structure as \code{\link[smashrgen]{ebps}} but with two
#'   additional columns in \code{posterior}:
#' \describe{
#'   \item{var}{Posterior variance \eqn{\text{Var}(\lambda_i \mid x_i)}
#'     computed via the log-normal formula
#'     \eqn{\exp(2m_i + v_i)(\exp(v_i) - 1)} when \code{var_log} is
#'     available (i.e. \code{wave_trans = 'dwt'}). \code{NA} otherwise.}
#'   \item{pip}{Not applicable for \code{ebps} (no spike-slab component).
#'     Always \code{NA}.}
#' }
#'
#' @seealso \code{\link[smashrgen]{ebps}},
#'   \code{\link{ebpm_point_gamma_with_uq}}
#'
#' @export
ebps_with_uq <- function(x, s = NULL, ...) {

  fit <- smashrgen::ebps(x, s, ...)

  # var_log only available for wave_trans='dwt'
  if (!is.null(fit$posterior$var_log)) {
    m_i <- fit$posterior$mean_log
    v_i <- fit$posterior$var_log
    fit$posterior$var <- exp(2 * m_i + v_i) * (exp(v_i) - 1)
  } else {
    fit$posterior$var <- NA
  }

  # PIP not applicable for smoothing prior
  fit$posterior$pip <- NA

  return(fit)
}
