# Taylor approximation for lgamma(a+x) - lgamma(a) for small x
# Adapted from ebpm package (DongyueXie/ebpm) to avoid ::: call
lgamma_diff_taylor_local <- function(x, dx) {
  c <- x
  out <- digamma(x) * dx + 1/2 * psigamma(c, deriv = 1) * dx^2
  return(out)
}

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
      control
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

  posterior <- data.frame(mean = mu_pm, mean_log = mu_log_pm)

  fitted_g = list(
  pi0 = pi0_est,
  shape = alpha_est,
  scale = 1/beta_est,  # Note: convert rate to scale
  gamma = gamma_est,
  type = "covariate_dependent")

  return(list(
    #fitted_pi0 = pi0_est,
    #fitted_alpha = alpha_est,
    #fitted_beta = beta_est,
    #fitted_gamma = gamma_est,
    fitted_g = fitted_g,
    posterior = posterior,
    log_likelihood = -opt_result$minimum,
    convergence_code = opt_result$code
  ))
}

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

transform_param_back_multiplier <- function(par, n_gamma) {
  # par = (logit_pi0, log_alpha, log_beta, gamma)
  original <- par
  original[1] <- plogis(par[1])  # pi0
  original[2] <- exp(par[2])     # alpha
  original[3] <- exp(par[3])     # beta
  # gamma remains untransformed
  return(original)
}

# it is equivalent to dnbinom in R wiht log = T when X is integer; I allow  it  to compute when x is not integer
dnbinom_cts_log_vec <- function(x, a, prob){
  #browser()
  if(length(x) > 1 && length(a) == 1){a = replicate(length(x), a)}
  tmp = x*log(1-prob)
  tmp[x == 0] = 0 ## R says 0*-Inf = NaN
  ## compute lgamma(a + x) - lgamma(a)
  lgamma_diff = lgamma(a + x)  - lgamma(a)

  subset = (lgamma_diff==0)  ## can occur when x very small compared with a
  lgamma_diff[subset] = lgamma_diff_taylor_local(a[subset], x[subset])

  return(a*log(prob) + tmp + lgamma_diff - lgamma(x+1))
}
