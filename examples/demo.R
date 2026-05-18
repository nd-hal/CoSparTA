# =============================================================================
# CxtEBTD — end-to-end demo
# Runs top-to-bottom; all objects are self-contained.
# =============================================================================

devtools::load_all("~/Desktop/CxtEBTD")

# =============================================================================
# STEP 1: Simulate a tensor
#
# 500 observations x 100 time points x 50 channels, K = 4 components.
# Two continuous covariates (no intercept column).
# sparsity = 50 targets ~99 % zeros for a tensor of this size.
# =============================================================================

n <- 500; p <- 100; w <- 50; K <- 4

# Two standardised continuous covariates
set.seed(1)
Xcov_demo <- cbind(rnorm(n), rnorm(n))

# Fixed covariate effects: one length-2 gamma vector per component
gamma_demo <- list(
  c( 0.6, -0.3),
  c(-0.5,  0.4),
  c( 0.3,  0.5),
  c(-0.4, -0.2)
)

sim <- simulate_tensor(
  n         = n,
  p         = p,
  w         = w,
  K         = K,
  Xcov      = Xcov_demo,
  gamma_true = gamma_demo,
  sparsity  = 50,
  seed      = 42
)

cat(sprintf("Observed sparsity: %.1f%% zeros\n", sim$sparsity_pct))

X_obs <- sim$X
Xcov  <- sim$Xcov      # pass this to the fitter

# =============================================================================
# STEP 2: Fit CxtEBTD
# =============================================================================

fit <- CxtEBTD(
  X                    = X_obs,
  K                    = K,
  Xcov                 = Xcov,
  init                 = "random_gamma",
  maxiter              = 50,
  convergence_criteria = "factor_change",
  verbose              = TRUE
)

# =============================================================================
# STEP 3: Normalize and visualize
# =============================================================================

nf <- normalize_factors(fit)

cat("\nComponent weights (lambda), descending:\n")
print(round(nf$lambda, 4))

# Time-factor plot (one line per component)
p_time <- plot_time_factors(Ef = nf$Ef)
print(p_time)

# Channel-factor plot (one bar per channel)
p_chan <- plot_channel_factors(Ew = nf$Ew)
print(p_chan)

# =============================================================================
# STEP 4: Covariate inference (delta method)
# =============================================================================

cat("\n--- Covariate coefficient CIs (delta method) ---\n")
gamma_ci <- get_gamma_ci(fit, method = "delta", level = 0.95)

for (k in seq_len(K)) {
  ci <- gamma_ci[[k]]
  if (is.null(ci)) {
    cat(sprintf("Factor %d: unsupervised (no gamma)\n", k))
  } else {
    cat(sprintf("\nFactor %d:\n", k))
    print(data.frame(
      estimate = round(ci$estimate, 4),
      se       = round(ci$se,       4),
      lower    = round(ci$lower,    4),
      upper    = round(ci$upper,    4),
      pvalue   = signif(ci$pvalue,  3)
    ))
  }
}

# =============================================================================
# STEP 5: Uncertainty quantification (L-mode credible intervals)
# =============================================================================

cat("\n--- L-mode 95% credible intervals ---\n")
ci_L <- get_credible_interval(fit, mode = "L", level = 0.95)
cat(sprintf("lower range: [%.4f, %.4f]\n", min(ci_L$lower, na.rm = TRUE),
                                             max(ci_L$lower, na.rm = TRUE)))
cat(sprintf("upper range: [%.4f, %.4f]\n", min(ci_L$upper, na.rm = TRUE),
                                             max(ci_L$upper, na.rm = TRUE)))

# =============================================================================
# STEP 6: Missing data
# =============================================================================

cat("\n--- Missing data ---\n")
mask <- generate_missing_mask(X_obs, missing_rate = 0.10, seed = 42)

fit_miss <- CxtEBTD_missing(
  X                    = mask$X_obs,
  K                    = K,
  Xcov                 = Xcov,
  obs_mask             = mask$obs_mask,
  init                 = "random_gamma",
  maxiter              = 50,
  convergence_criteria = "factor_change",
  verbose              = TRUE
)

cat("\nMissing-data fit factor dimensions:\n")
cat(sprintf("  El: %d x %d\n", nrow(fit_miss$res$ql$El), ncol(fit_miss$res$ql$El)))
cat(sprintf("  Ef: %d x %d\n", nrow(fit_miss$res$qf$Ef), ncol(fit_miss$res$qf$Ef)))
cat(sprintf("  Ew: %d x %d\n", nrow(fit_miss$res$qw$Ew), ncol(fit_miss$res$qw$Ew)))

cat("\nDemo complete.\n")
