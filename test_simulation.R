# =============================================================================
# End-to-end simulation test for CxtEBTD
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
  lambda_i   <- as.vector(exp(X %*% gamma_true))       # n-vector of rates
  active     <- rbinom(n, 1, 1 - pi0_true)             # structural zero mask
  u          <- numeric(n)
  idx        <- which(active == 1)
  if (length(idx) > 0) {
    u[idx]   <- rgamma(length(idx), shape = alpha_true, rate = beta_true / lambda_i[idx])
  }
  u <- u / sqrt(sum(u^2))                              # normalize to L2 norm = 1
  u
}

# =============================================================================
# Simulation setup
# =============================================================================

set.seed(42)
K <- 2
n <- 100
p <- 20
w <- 10

# Covariate matrix: intercept + binary (gender) + continuous (age 20-80)
X_cov <- cbind(1, rbinom(n, 1, 0.5), runif(n, 20, 80))

# True covariate effects per component
gamma1_true <- c(0,  0.8,  0.15)
gamma2_true <- c(0, -0.5,  0.25)

# True U1 (observation mode) — covariate-dependent
U1_true <- cbind(
  gen_U1_with_covariates(n, X_cov, gamma1_true, pi0_true = 0.2,
                         alpha_true = 3,   beta_true = 2,   seed = 1),
  gen_U1_with_covariates(n, X_cov, gamma2_true, pi0_true = 0.2,
                         alpha_true = 1.5, beta_true = 2.5, seed = 2)
)

# True U2 (time mode) — block structure
u2_k1        <- c(rep(1, p/2), rep(0, p/2))
u2_k2        <- c(rep(0, p/2), rep(1, p/2))
U2_true      <- cbind(u2_k1 / sqrt(sum(u2_k1^2)),
                      u2_k2 / sqrt(sum(u2_k2^2)))

# True U3 (channel mode) — rgamma with structural zeros in complementary halves
set.seed(3)
u3_k1        <- c(rgamma(w/2, shape = 2, rate = 1), rep(0, w/2))
u3_k2        <- c(rep(0, w/2), rgamma(w/2, shape = 2, rate = 1))
U3_true      <- cbind(u3_k1 / sqrt(sum(u3_k1^2)),
                      u3_k2 / sqrt(sum(u3_k2^2)))

# Composite weights
weight      <- c(2, 1.5)
scalenpw    <- sqrt(n) * sqrt(p) * sqrt(w)
sparse      <- 20

lambda_tensor <- (
  weight[1] * (U1_true[,1] %o% U2_true[,1] %o% U3_true[,1]) +
  weight[2] * (U1_true[,2] %o% U2_true[,2] %o% U3_true[,2])
) * scalenpw / sparse

set.seed(99)
X_obs <- array(rpois(n * p * w, lambda_tensor), dim = c(n, p, w))

cat(sprintf("Observed sparsity: %.1f%% zeros\n",
            100 * mean(X_obs == 0)))

# =============================================================================
# Test A — Unsupervised
# =============================================================================

cat("\n===== Test A: Unsupervised (Xcov = NULL) =====\n")

fit_A <- CxtEBTD(X_obs, K = 2, Xcov = NULL,
                 init                 = 'random_gamma',
                 maxiter              = 30,
                 convergence_criteria = 'ELBO',
                 verbose              = TRUE)

norm_A    <- normalize_factors(fit_A)
lambda_A  <- array(0, dim = c(n, p, w))
for (k in 1:K) {
  lambda_A <- lambda_A + norm_A$lambda[k] *
    (norm_A$El[,k] %o% norm_A$Ef[,k] %o% norm_A$Ew[,k])
}
mse_A <- mean((lambda_A - lambda_tensor)^2)
cat(sprintf("Test A — Tensor MSE vs lambda_tensor: %.6f\n", mse_A))

# Per-factor RMSE via which_rank
for (k in 1:K) {
  r1 <- which_rank(U1_true[,k], norm_A$El)
  r2 <- which_rank(U2_true[,k], norm_A$Ef)
  r3 <- which_rank(U3_true[,k], norm_A$Ew)
  cat(sprintf("  Factor %d: U1 RMSE=%.4f  U2 RMSE=%.4f  U3 RMSE=%.4f\n",
              k, r1$diff, r2$diff, r3$diff))
}

# =============================================================================
# Test B — Supervised (with covariates)
# =============================================================================

cat("\n===== Test B: Supervised (Xcov supplied) =====\n")

fit_B <- CxtEBTD(X_obs, K = 2, Xcov = X_cov,
                 init                 = 'random_gamma',
                 maxiter              = 30,
                 convergence_criteria = 'ELBO',
                 verbose              = TRUE)

# Check posterior variance and PIP
cat("fit$res$ql$Varl is NULL:", is.null(fit_B$res$ql$Varl), "\n")
cat("fit$res$ql$PIPl is NULL:", is.null(fit_B$res$ql$PIPl), "\n")
if (!is.null(fit_B$res$ql$Varl)) cat("Varl range:", range(fit_B$res$ql$Varl, na.rm=TRUE), "\n")
if (!is.null(fit_B$res$ql$PIPl)) cat("PIPl range:", range(fit_B$res$ql$PIPl, na.rm=TRUE), "\n")

# Gamma estimates vs truth
cat("\nGamma estimates vs truth:\n")
cat("Component 1 — true:", gamma1_true, "\n")
cat("Component 1 — est: ", fit_B$res$gl[[1]]$gamma, "\n")
cat("Component 2 — true:", gamma2_true, "\n")
cat("Component 2 — est: ", fit_B$res$gl[[2]]$gamma, "\n")

# Tensor MSE and factor RMSE
norm_B    <- normalize_factors(fit_B)
lambda_B  <- array(0, dim = c(n, p, w))
for (k in 1:K) {
  lambda_B <- lambda_B + norm_B$lambda[k] *
    (norm_B$El[,k] %o% norm_B$Ef[,k] %o% norm_B$Ew[,k])
}
mse_B <- mean((lambda_B - lambda_tensor)^2)
cat(sprintf("\nTest B — Tensor MSE vs lambda_tensor: %.6f\n", mse_B))

for (k in 1:K) {
  r1 <- which_rank(U1_true[,k], norm_B$El)
  r2 <- which_rank(U2_true[,k], norm_B$Ef)
  r3 <- which_rank(U3_true[,k], norm_B$Ew)
  cat(sprintf("  Factor %d: U1 RMSE=%.4f  U2 RMSE=%.4f  U3 RMSE=%.4f\n",
              k, r1$diff, r2$diff, r3$diff))
}

# PIP and credible intervals for L mode
pip_B <- get_pip(fit_B, mode = 'L')
ci_B  <- get_credible_interval(fit_B, mode = 'L')
cat("CI lower range:", range(ci_B$lower, na.rm=TRUE), "\n")
cat("CI upper range:", range(ci_B$upper, na.rm=TRUE), "\n")

# =============================================================================
# Test C — Missing data
# =============================================================================

cat("\n===== Test C: Missing data =====\n")

mask <- generate_missing_mask(X_obs, missing_rate = 0.1, seed = 7)
cat(sprintf("Missing entries: %d (%.1f%% of nonzeros)\n",
            nrow(mask$missing_nonzero_indices),
            100 * nrow(mask$missing_nonzero_indices) / sum(X_obs > 0)))

fit_C <- CxtEBTD_missing(X = mask$X_obs, K = 2, obs_mask = mask$obs_mask, Xcov = X_cov,
                          init         = 'random_gamma',
                          maxiter      = 30,
                          verbose      = TRUE)

norm_C    <- normalize_factors(fit_C)
miss_eval <- evaluate_missing_prediction(fit_C, mask)
cat(sprintf("Test C — Missing-entry prediction RMSE: %.6f\n", miss_eval$rmse))
cat(sprintf("Test C — Missing-entry prediction MAE:  %.6f\n", miss_eval$mae))

# =============================================================================
# Test D — Rank-specific covariates
# =============================================================================

cat("\n===== Test D: Rank-specific covariates =====\n")

# Rank 1 gets covariates, Rank 2 is unsupervised
Xcov_list <- list(X_cov, NULL)
fit_D <- CxtEBTD(X_obs, K = 2, Xcov = Xcov_list,
                 init = 'random_gamma', maxiter = 5,
                 convergence_criteria = 'ELBO',
                 verbose = TRUE)
cat("Rank 1 gamma estimate:", fit_D$res$gl[[1]]$gamma, "\n")
cat("Rank 1 fitted_g type:", fit_D$res$gl[[1]]$type, "\n")
cat("Rank 2 fitted_g type:", fit_D$res$gl[[2]]$type, "\n")
# Rank 2 should NOT have $type == "covariate_dependent"
cat("Rank 2 is unsupervised:", !identical(fit_D$res$gl[[2]]$type, "covariate_dependent"), "\n")
norm_D <- normalize_factors(fit_D)
cat("Test D — lambda:", round(norm_D$lambda, 4), "\n")

# =============================================================================
# Test E — Posterior quantiles (using fit_B from Test B)
# =============================================================================

cat("\n===== Test E: Posterior quantiles =====\n")

q_L <- get_posterior_quantile(fit_B, probs = c(0.025, 0.5, 0.975), mode = 'L')
cat("Names:", names(q_L), "\n")
cat("Dimensions q2.5:", dim(q_L$q2.5), "\n")
cat("Dimensions q97.5:", dim(q_L$q97.5), "\n")

# Verify ordering: q2.5 <= q50 <= q97.5
cat("q2.5 <= q50 everywhere:", all(q_L$q2.5 <= q_L$q50, na.rm = TRUE), "\n")
cat("q50 <= q97.5 everywhere:", all(q_L$q50 <= q_L$q97.5, na.rm = TRUE), "\n")

# Verify elements with low PIP have quantile = 0
pip_L <- get_pip(fit_B, mode = 'L')
low_pip <- pip_L < 0.01
cat("Low PIP elements with q97.5 = 0:", sum(q_L$q97.5[low_pip] == 0, na.rm = TRUE),
    "out of", sum(low_pip, na.rm = TRUE), "\n")

# Test F-mode: should get informative error since ebps has NA shape_post
cat("F-mode quantile test: ")
tryCatch({
  get_posterior_quantile(fit_B, mode = 'F')
  cat("ERROR — should have stopped\n")
}, error = function(e) cat("correctly errored:", e$message, "\n"))

# =============================================================================
# Test F — normalize_factors (using fit_B)
# =============================================================================

cat("\n===== Test F: normalize_factors =====\n")

nf_B <- normalize_factors(fit_B)
# Verify unit norm columns
col_norms_El <- sqrt(colSums(nf_B$El^2))
col_norms_Ef <- sqrt(colSums(nf_B$Ef^2))
col_norms_Ew <- sqrt(colSums(nf_B$Ew^2))
cat("El column norms:", round(col_norms_El, 6), "\n")
cat("Ef column norms:", round(col_norms_Ef, 6), "\n")
cat("Ew column norms:", round(col_norms_Ew, 6), "\n")

# Verify lambda sorted descending
cat("Lambda sorted descending:", all(diff(nf_B$lambda) <= 0), "\n")
cat("Lambda values:", round(nf_B$lambda, 4), "\n")

# Verify reconstruction from normalized factors matches raw reconstruction
recon_raw <- reconstruct_tensor(fit_B)
recon_nf <- array(0, dim = dim(recon_raw))
for (k in 1:K) {
  lf_k <- (nf_B$El[,k] * nf_B$lambda[k]) %o% nf_B$Ef[,k] %o% nf_B$Ew[,k]
  recon_nf <- recon_nf + lf_k
}
cat("Normalized vs raw reconstruction max diff:", max(abs(recon_nf - recon_raw)), "\n")

# =============================================================================
# Test G — project_tensor (using fit_B)
# =============================================================================

cat("\n===== Test G: project_tensor =====\n")

# Project the original training tensor
proj_full <- project_tensor(X_obs, fit_B, normalize = TRUE)
cat("Projection dimensions:", dim(proj_full), "\n")

# Project without normalization
proj_raw <- project_tensor(X_obs, fit_B, normalize = FALSE)
cat("Raw projection dimensions:", dim(proj_raw), "\n")

# Single observation (matrix input)
proj_single <- project_tensor(X_obs[1,,], fit_B, normalize = TRUE)
cat("Single obs projection is vector:", is.vector(proj_single), "\n")
cat("Single obs length:", length(proj_single), "\n")
cat("Single obs matches row 1:", all.equal(proj_single, proj_full[1,]), "\n")

# =============================================================================
# Test H — reconstruct_tensor (using fit_B)
# =============================================================================

cat("\n===== Test H: reconstruct_tensor =====\n")

recon <- reconstruct_tensor(fit_B)
cat("Reconstruction dimensions:", dim(recon), "\n")
cat("Reconstruction range:", range(recon), "\n")
cat("All non-negative:", all(recon >= 0), "\n")

# MSE vs true lambda_tensor
mse_recon <- mean((recon - lambda_tensor)^2)
cat(sprintf("Reconstruction MSE vs lambda_tensor: %.6f\n", mse_recon))
