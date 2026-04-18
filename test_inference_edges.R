# =============================================================================
# Inference edge-case tests for CxtEBTD
# Tests I–O: covers get_pip, get_credible_interval, get_significant_patterns,
# get_posterior_quantile, normalize_factors, project_tensor, reconstruct_tensor
# under edge conditions not exercised by test_simulation.R
# =============================================================================

remotes::install_github("xzhang0407/CxtEBTD", ref = "dev")
devtools::load_all("~/Desktop/CxtEBTD")
library(ebpm)
library(smashrgen)
library(slam)
library(data.table)
library(matrixStats)

# =============================================================================
# Helper functions
# =============================================================================

#' Generate covariate-dependent U1 (loading matrix) for one component
gen_U1_with_covariates <- function(n, X, gamma_true, pi0_true, alpha_true, beta_true, seed) {
  set.seed(seed)
  lambda_i   <- as.vector(exp(X %*% gamma_true))
  active     <- rbinom(n, 1, 1 - pi0_true)
  u          <- numeric(n)
  idx        <- which(active == 1)
  if (length(idx) > 0) {
    u[idx]   <- rgamma(length(idx), shape = alpha_true, rate = beta_true / lambda_i[idx])
  }
  u <- u / sqrt(sum(u^2))
  u
}

# =============================================================================
# Simulation setup (identical to test_simulation.R)
# =============================================================================

set.seed(42)
K <- 2
n <- 100
p <- 20
w <- 10

X_cov <- cbind(1, rbinom(n, 1, 0.5), runif(n, 20, 80))

gamma1_true <- c(0,  0.8,  0.15)
gamma2_true <- c(0, -0.5,  0.25)

U1_true <- cbind(
  gen_U1_with_covariates(n, X_cov, gamma1_true, pi0_true = 0.2,
                         alpha_true = 3,   beta_true = 2,   seed = 1),
  gen_U1_with_covariates(n, X_cov, gamma2_true, pi0_true = 0.2,
                         alpha_true = 1.5, beta_true = 2.5, seed = 2)
)

u2_k1   <- c(rep(1, p/2), rep(0, p/2))
u2_k2   <- c(rep(0, p/2), rep(1, p/2))
U2_true <- cbind(u2_k1 / sqrt(sum(u2_k1^2)),
                 u2_k2 / sqrt(sum(u2_k2^2)))

set.seed(3)
u3_k1   <- c(rgamma(w/2, shape = 2, rate = 1), rep(0, w/2))
u3_k2   <- c(rep(0, w/2), rgamma(w/2, shape = 2, rate = 1))
U3_true <- cbind(u3_k1 / sqrt(sum(u3_k1^2)),
                 u3_k2 / sqrt(sum(u3_k2^2)))

weight   <- c(2, 1.5)
scalenpw <- sqrt(n) * sqrt(p) * sqrt(w)
sparse   <- 20

lambda_tensor <- (
  weight[1] * (U1_true[,1] %o% U2_true[,1] %o% U3_true[,1]) +
  weight[2] * (U1_true[,2] %o% U2_true[,2] %o% U3_true[,2])
) * scalenpw / sparse

set.seed(99)
X_obs <- array(rpois(n * p * w, lambda_tensor), dim = c(n, p, w))

cat(sprintf("Observed sparsity: %.1f%% zeros\n", 100 * mean(X_obs == 0)))

# =============================================================================
# Fit once (shared across all tests)
# =============================================================================

fit <- CxtEBTD(X_obs, K = 2, Xcov = X_cov,
               init = 'random_gamma', maxiter = 5,
               adj_LF_scale = FALSE, convergence_criteria = 'ELBO',
               verbose = FALSE)

# =============================================================================
# Test I — get_pip() with threshold
# =============================================================================

cat("\n===== Test I: get_pip with threshold =====\n")

pip_raw   <- get_pip(fit, mode = 'L')
pip_thresh <- get_pip(fit, mode = 'L', threshold = 0.5)
cat("Raw PIP dimensions:", dim(pip_raw), "\n")
cat("Thresholded PIP is logical:", is.logical(pip_thresh), "\n")
cat("Thresholded TRUE count:", sum(pip_thresh, na.rm = TRUE), "\n")
cat("Manual threshold matches:", all.equal(pip_thresh, pip_raw > 0.5), "\n")

# W-mode PIP
pip_W <- get_pip(fit, mode = 'W')
cat("W-mode PIP dimensions:", dim(pip_W), "\n")
cat("W-mode PIP range:", range(pip_W, na.rm = TRUE), "\n")

# F-mode PIP — should warn since ebps has no spike
cat("F-mode PIP test: ")
pip_F <- get_pip(fit, mode = 'F')
cat("F-mode PIP is NULL:", is.null(pip_F), "\n")

# =============================================================================
# Test J — get_credible_interval() all modes
# =============================================================================

cat("\n===== Test J: get_credible_interval all modes =====\n")

ci_L <- get_credible_interval(fit, mode = 'L', level = 0.95)
cat("L-mode CI: lower range", range(ci_L$lower, na.rm = TRUE),
    "upper range", range(ci_L$upper, na.rm = TRUE), "\n")
cat("Lower <= mean:", all(ci_L$lower <= ci_L$mean, na.rm = TRUE), "\n")
cat("Mean <= upper:", all(ci_L$mean <= ci_L$upper, na.rm = TRUE), "\n")

ci_W <- get_credible_interval(fit, mode = 'W', level = 0.90)
cat("W-mode 90% CI: lower range", range(ci_W$lower, na.rm = TRUE),
    "upper range", range(ci_W$upper, na.rm = TRUE), "\n")

# F-mode — should warn since ebps variance might be NA
cat("F-mode CI test: ")
ci_F <- get_credible_interval(fit, mode = 'F')
cat("F-mode CI is NULL:", is.null(ci_F), "\n")

# =============================================================================
# Test K — get_significant_patterns()
# =============================================================================

cat("\n===== Test K: get_significant_patterns =====\n")

patterns <- get_significant_patterns(fit, alpha = 0.05, mode = 'both')
cat("Number of factors:", length(patterns), "\n")
for (l in seq_along(patterns)) {
  cat(sprintf("  Factor %d: %d active times, %d active channels\n",
      patterns[[l]]$factor,
      patterns[[l]]$n_active_times,
      patterns[[l]]$n_active_channels))
}

# Strict alpha — should find fewer or no discoveries
patterns_strict <- get_significant_patterns(fit, alpha = 0.001)
cat("Strict alpha (0.001) — Factor 1 active channels:",
    patterns_strict[[1]]$n_active_channels, "\n")

# W-only mode
patterns_W <- get_significant_patterns(fit, alpha = 0.05, mode = 'W')
cat("W-only mode — active_times is NULL:", is.null(patterns_W[[1]]$active_times), "\n")

# =============================================================================
# Test L — get_posterior_quantile() on W-mode
# =============================================================================

cat("\n===== Test L: Posterior quantiles W-mode =====\n")

q_W <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = 'W')
cat("W-mode quantile names:", names(q_W), "\n")
cat("W-mode q2.5 dimensions:", dim(q_W$q2.5), "\n")
cat("W-mode q2.5 <= q97.5:", all(q_W$q2.5 <= q_W$q97.5, na.rm = TRUE), "\n")

# =============================================================================
# Test M — K=1 single component
# =============================================================================

cat("\n===== Test M: K=1 single component =====\n")

fit_K1 <- CxtEBTD(X_obs, K = 1, Xcov = X_cov,
                  init = 'random_gamma', maxiter = 5,
                  adj_LF_scale = FALSE, convergence_criteria = 'ELBO',
                  verbose = FALSE)
cat("El dimensions:", dim(fit_K1$res$ql$El), "\n")
cat("Ef dimensions:", dim(fit_K1$res$qf$Ef), "\n")
nf_K1 <- normalize_factors(fit_K1)
cat("Lambda:", round(nf_K1$lambda, 4), "\n")
recon_K1 <- reconstruct_tensor(fit_K1)
cat("Reconstruction dimensions:", dim(recon_K1), "\n")
q_K1 <- get_posterior_quantile(fit_K1, mode = 'L')
cat("Quantile dimensions:", dim(q_K1$q2.5), "\n")
cat("K=1 all passed\n")

# =============================================================================
# Test N — Rank-specific covariates with missing data
# =============================================================================

cat("\n===== Test N: Rank-specific covariates + missing data =====\n")

mask <- generate_missing_mask(X_obs, missing_rate = 0.1, seed = 7)
fit_N <- CxtEBTD_missing(mask$X_obs, K = 2,
                          Xcov = list(X_cov, NULL),
                          obs_mask = mask$obs_mask,
                          init = 'random_gamma', maxiter = 5,
                          adj_LF_scale = FALSE, verbose = FALSE)
cat("Rank 1 type:", fit_N$res$gl[[1]]$type, "\n")
cat("Rank 2 type:", fit_N$res$gl[[2]]$type, "\n")
miss_eval_N <- evaluate_missing_prediction(fit_N, mask)
cat(sprintf("Missing RMSE: %.4f\n", miss_eval_N$rmse))
cat("Rank-specific + missing passed\n")

# =============================================================================
# Test O — project_tensor dimension mismatch error
# =============================================================================

cat("\n===== Test O: project_tensor error handling =====\n")

X_wrong_p <- array(0, dim = c(5, 15, 10))  # wrong p
cat("Wrong p test: ")
tryCatch({
  project_tensor(X_wrong_p, fit)
  cat("ERROR — should have stopped\n")
}, error = function(e) cat("correctly errored:", e$message, "\n"))

X_wrong_w <- array(0, dim = c(5, 20, 8))   # wrong w
cat("Wrong w test: ")
tryCatch({
  project_tensor(X_wrong_w, fit)
  cat("ERROR — should have stopped\n")
}, error = function(e) cat("correctly errored:", e$message, "\n"))
