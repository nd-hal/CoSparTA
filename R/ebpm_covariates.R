# Taylor approximation for lgamma(a+x) - lgamma(a) for small x
# Adapted from ebpm package (DongyueXie/ebpm) to avoid ::: call

#' @keywords internal
.lgamma_diff_taylor_local <- function(x, dx) {
  c <- x
  out <- digamma(x) * dx + 1/2 * psigamma(c, deriv = 1) * dx^2
  return(out)
}

#' Empirical Bayes Poisson Mean Estimation with Covariate-Dependent Rates
#'
#' @description
#' Fits a spike-and-slab point-gamma prior to a vector of counts \code{x}
#' where the Poisson rate is modulated by individual-level covariates. The
#' effective rate for observation \eqn{i} is
#' \eqn{\lambda_i = \beta \cdot \exp(X_i^\top \gamma)},
#' so that the prior mean scales multiplicatively with covariates. Prior
#' parameters \eqn{(\pi_0, \alpha, \beta, \gamma)} are estimated by maximizing
#' the marginal likelihood via \code{nlm} (with \code{optim} as fallback).
#' This function is the core building block of the covariate-aware mode update
#' in \code{\link{CoSparTA}}.
#'
#' @param x Non-negative integer vector of observed counts.
#' @param s Numeric scalar or vector of length \code{length(x)} giving
#'   per-observation exposure/scaling factors. Default \code{1}.
#' @param X Numeric covariate matrix of dimension \code{n x q}, where rows
#'   correspond to observations in \code{x} and columns to covariates.
#' @param g_init Optional numeric vector of starting values in the order
#'   \code{c(pi0, alpha, beta, gamma_1, ..., gamma_q)}. If \code{NULL},
#'   initialization uses a baseline \code{ebpm::ebpm_point_gamma} fit with
#'   \code{gamma = 0}. Default \code{NULL}.
#' @param control Optional list of control arguments passed to \code{nlm}.
#'   Defaults: \code{stepmax = 1}, \code{gradtol = 1e-6},
#'   \code{steptol = 1e-8}, \code{iterlim = 1000},
#'   \code{check.analyticals = FALSE}.
#'
#' @return A named list with:
#' \describe{
#'   \item{fitted_g}{List of estimated prior parameters:
#'     \code{pi0} — spike (point-mass at zero) probability;
#'     \code{shape} — gamma shape parameter \eqn{\alpha};
#'     \code{scale} — gamma scale parameter \eqn{1/\beta};
#'     \code{gamma} — covariate coefficient vector of length \code{q};
#'     \code{type} — always \code{"covariate_dependent"}.}
#'   \item{posterior}{Data frame with one row per observation and columns:
#'     \code{mean} — posterior mean \eqn{E[\theta_i \mid x_i]};
#'     \code{mean_log} — posterior log-mean \eqn{E[\log\theta_i \mid x_i]}
#'     (\code{-Inf} for observations with \eqn{x_i = 0});
#'     \code{var} — posterior variance
#'     \eqn{\text{Var}(\theta_i \mid x_i)}, decomposed as within-component
#'     variance (uncertainty in the gamma component) plus between-component
#'     variance (uncertainty about whether the spike fired);
#'     \code{pip} — posterior inclusion probability
#'     \eqn{P(\theta_i \neq 0 \mid x_i) = 1 - \hat{\pi}_i}, giving the
#'     probability that observation \eqn{i} has a truly non-zero loading;
#'     \code{shape_post} — posterior gamma shape parameter
#'     \eqn{\alpha + x_i}, the shape of the gamma kernel in the non-spike
#'     component; \code{rate_post} — posterior gamma rate parameter
#'     \eqn{\beta / \lambda_i + s_i}, the rate of the gamma kernel.}
#'   \item{log_likelihood}{Maximized marginal log-likelihood.}
#'   \item{convergence_code}{Optimizer convergence code. For \code{nlm}:
#'     1--2 indicate convergence; 3--5 indicate potential issues (see
#'     \code{?nlm}). For \code{optim} fallback: 0 = converged.}
#' }
#'
#' @seealso \code{\link{CoSparTA}}, \code{\link[ebpm]{ebpm_point_gamma}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 200
#' X <- matrix(rnorm(n * 2), nrow = n)
#' true_gamma <- c(0.4, -0.3)
#' x <- rpois(n, lambda = exp(0.5 + X %*% true_gamma))
#'
#' fit <- ebpm_point_gamma_multiplier_covariates(x, s = 1, X = X)
#' fit$fitted_g$gamma     # estimated covariate coefficients
#' fit$fitted_g$pi0       # estimated spike probability
#' head(fit$posterior)    # posterior means and log-means
#' }
#'
#' @export
ebpm_point_gamma_multiplier_covariates <- function(x, s = 1, X, g_init = NULL, control = NULL) {

  # Handle inputs
  if(length(s) == 1) {s <- rep(s, length(x))}
  if(is.null(control)) {
    control <- list(
      stepmax = 1,
      gradtol = 1e-6,
      steptol = 1e-8,
      iterlim = 1000,
      check.analyticals = FALSE
    )
  }

  n_gamma <- ncol(X)

  # Initialization
  if(is.null(g_init)) {
    baseline <- ebpm::ebpm_point_gamma(x, s)
    pi0_init <- sum(x == 0) / length(x)
    pi0_init <- pmax(0.01, pmin(0.99, pi0_init))
    alpha_init <- baseline$fitted_g$shape
    beta_init <- 1 / baseline$fitted_g$scale  # Convert scale to rate
    gamma_init <- rep(0, n_gamma)

    g_init <- c(pi0_init, alpha_init, beta_init, gamma_init)
  }

  # Optimization
  fn_params <- list(x = x, s = s, X = X, n_gamma = n_gamma)

  opt_result <- try({
    do.call(nlm, c(
      list(pg_multiplier_nlm_fn, transform_param_multiplier(g_init, n_gamma)),
      fn_params,
      control,
      list(hessian = TRUE)
    ))
  }, silent = TRUE)

  # Fallback to optim if nlm fails
  if(inherits(opt_result, "try-error") || opt_result$code > 3) {
    opt_result <- try({
      opt_temp <- optim(
        par = transform_param_multiplier(g_init, n_gamma),
        fn = pg_multiplier_nlm_fn,
        x = x, s = s, X = X, n_gamma = n_gamma,
        method = "L-BFGS-B",
        control = list(maxit = 1000)
      )
      list(estimate = opt_temp$par,
           minimum = opt_temp$value,
           code = ifelse(opt_temp$convergence == 0, 1, 4))
    }, silent = TRUE)
  }

  if(inherits(opt_result, "try-error")) {
    stop("Optimization failed")
  }

  # Extract parameters
  opt_params <- transform_param_back_multiplier(opt_result$estimate, n_gamma)
  pi0_est <- opt_params[1]
  alpha_est <- opt_params[2]
  beta_est <- opt_params[3]
  gamma_est <- opt_params[4:(3 + n_gamma)]

  # Compute posteriors
  linear_pred <- X %*% gamma_est
  lambda_i <- exp(linear_pred)

  # Posterior spike probabilities
  pi_hat <- rep(0, length(x))
  zero_idx <- (x == 0)

  if(any(zero_idx)) {
    nb_prob <- beta_est / (beta_est + s[zero_idx] * lambda_i[zero_idx])
    nb_zero_prob <- dnbinom(0, size = alpha_est, prob = nb_prob, log = FALSE)
    denom <- pi0_est + (1 - pi0_est) * nb_zero_prob
    pi_hat[zero_idx] <- pi0_est / denom
  }

  # Posterior means: (1 - pi_hat) * (alpha + Y) / (beta/lambda + s)
  beta_eff_i <- beta_est / lambda_i
  mu_pm <- (1 - pi_hat) * (alpha_est + x) / (beta_eff_i + s)

  # Posterior log-means
  mu_log_pm <- rep(-Inf, length(x))
  nonzero_idx <- (x > 0)
  if(any(nonzero_idx)) {
    mu_log_pm[nonzero_idx] <- (1 - pi_hat[nonzero_idx]) *
      (digamma(alpha_est + x[nonzero_idx]) - log(beta_eff_i[nonzero_idx] + s[nonzero_idx]))
  }

  var_within <- (1 - pi_hat) * (alpha_est + x) / (beta_eff_i + s)^2
  var_between <- pi_hat * (1 - pi_hat) * ((alpha_est + x) / (beta_eff_i + s))^2
  mu_var <- var_within + var_between
  
  posterior <- data.frame(
    mean       = mu_pm,
    mean_log   = mu_log_pm,
    var        = mu_var,
    pip        = 1 - pi_hat,
    shape_post = alpha_est + x,
    rate_post  = beta_eff_i + s
  )

  fitted_g = list(
  pi0 = pi0_est,
  shape = alpha_est,
  scale = 1/beta_est,  # Note: convert rate to scale
  gamma = gamma_est,
  type = "covariate_dependent")

  return(list(
    #fitted_gamma = gamma_est,
    fitted_g = fitted_g,
    posterior = posterior,
    log_likelihood = -opt_result$minimum,
    convergence_code = opt_result$code,
    hessian = if (!is.null(opt_result$hessian)) opt_result$hessian else NULL,
    family = "point_gamma"
  ))
}

#' @keywords internal
# Objective function
pg_multiplier_nlm_fn <- function(par, x, s, X, n_gamma) {

  pi0 <- plogis(par[1])
  alpha <- exp(par[2])
  beta <- exp(par[3])
  gamma <- par[4:(3 + n_gamma)]

  # Parameter bounds
  if(pi0 < 0.001 || pi0 > 0.999 || alpha < 0.001 || alpha > 1000 ||
     beta < 0.001 || beta > 1000) {
    return(1e10)
  }

  # Compute linear predictors
  linear_pred <- X %*% gamma
  linear_pred <- pmax(-100, pmin(100, linear_pred))
  lambda_i <- exp(linear_pred)

  # Compute negative binomial probabilities: p = beta / (beta + s * lambda)
  nb_prob <- beta / (beta + s * lambda_i)
  #d_log <- dnbinom(x, size = alpha, prob = nb_prob, log = TRUE)
  d_log <- dnbinom_cts_log_vec(x, alpha, prob = nb_prob)

  if(any(!is.finite(d_log))) {
    return(1e10)
  }

  zero_idx <- (x == 0)
  nonzero_idx <- !zero_idx

  # Likelihood
  ll_nonzero <- if(any(nonzero_idx)) {
    sum(log(1 - pi0) + d_log[nonzero_idx])
  } else {
    0
  }

  ll_zero <- if(any(zero_idx)) {
    sum(log(pi0 + (1 - pi0) * exp(d_log[zero_idx])))
  } else {
    0
  }

  total_ll <- ll_nonzero + ll_zero

  if(!is.finite(total_ll)) {
    return(1e10)
  }

  return(-total_ll)
}

#' @keywords internal
# Parameter transformation functions
transform_param_multiplier <- function(par, n_gamma) {
  # par = (pi0, alpha, beta, gamma)
  transformed <- par
  transformed[1] <- log(par[1] / (1 - par[1]))  # logit(pi0)
  transformed[2] <- log(par[2])                  # log(alpha)
  transformed[3] <- log(par[3])                  # log(beta)
  # gamma remains untransformed
  return(transformed)
}

#' @keywords internal
transform_param_back_multiplier <- function(par, n_gamma) {
  # par = (logit_pi0, log_alpha, log_beta, gamma)
  original <- par
  original[1] <- plogis(par[1])  # pi0
  original[2] <- exp(par[2])     # alpha
  original[3] <- exp(par[3])     # beta
  # gamma remains untransformed
  return(original)
}

#' @keywords internal
# it is equivalent to dnbinom in R wiht log = T when X is integer; I allow  it  to compute when x is not integer
dnbinom_cts_log_vec <- function(x, a, prob){
  if(length(x) > 1 && length(a) == 1){a = replicate(length(x), a)}
  tmp = x*log(1-prob)
  tmp[x == 0] = 0 ## R says 0*-Inf = NaN
  ## compute lgamma(a + x) - lgamma(a)
  lgamma_diff = lgamma(a + x)  - lgamma(a)

  subset = (lgamma_diff==0)  ## can occur when x very small compared with a
  lgamma_diff[subset] = .lgamma_diff_taylor_local(a[subset], x[subset])

  return(a*log(prob) + tmp + lgamma_diff - lgamma(x+1))
}
