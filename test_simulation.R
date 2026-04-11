# =============================================================================
# End-to-end simulation test for CxtEBTD
# =============================================================================

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

#' L2 norm of each column of a matrix
L2_norm <- function(x) {
  apply(x, 2, function(col) sqrt(sum(col^2)))
}

#' Post-hoc normalize factor matrices from a CxtEBTD fit object.
#' Returns components sorted by descending composite weight.
post_hoc_normalize <- function(res_obj) {
  El <- res_obj$ql$El
  Ef <- res_obj$qf$Ef
  Ew <- res_obj$qw$Ew

  norm_l <- L2_norm(El)
  norm_f <- L2_norm(Ef)
  norm_w <- L2_norm(Ew)

  weights <- norm_l * norm_f * norm_w
  ord     <- order(weights, decreasing = TRUE)

  U1_normed <- sweep(El[, ord, drop = FALSE], 2, norm_l[ord], "/")
  U2_normed <- sweep(Ef[, ord, drop = FALSE], 2, norm_f[ord], "/")
  U3_normed <- sweep(Ew[, ord, drop = FALSE], 2, norm_w[ord], "/")

  list(U1_normed      = U1_normed,
       U2_normed      = U2_normed,
       U3_normed      = U3_normed,
       weights_sorted = weights[ord])
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
                 adj_LF_scale         = FALSE,
                 convergence_criteria = 'ELBO',
                 verbose              = TRUE)

norm_A    <- post_hoc_normalize(fit_A$res)
lambda_A  <- array(0, dim = c(n, p, w))
for (k in 1:K) {
  lambda_A <- lambda_A + norm_A$weights_sorted[k] *
    (norm_A$U1_normed[,k] %o% norm_A$U2_normed[,k] %o% norm_A$U3_normed[,k])
}
mse_A <- mean((lambda_A - lambda_tensor)^2)
cat(sprintf("Test A — Tensor MSE vs lambda_tensor: %.6f\n", mse_A))

# Per-factor RMSE via which_rank
for (k in 1:K) {
  r1 <- which_rank(U1_true[,k], norm_A$U1_normed)
  r2 <- which_rank(U2_true[,k], norm_A$U2_normed)
  r3 <- which_rank(U3_true[,k], norm_A$U3_normed)
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
                 adj_LF_scale         = FALSE,
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
norm_B    <- post_hoc_normalize(fit_B$res)
lambda_B  <- array(0, dim = c(n, p, w))
for (k in 1:K) {
  lambda_B <- lambda_B + norm_B$weights_sorted[k] *
    (norm_B$U1_normed[,k] %o% norm_B$U2_normed[,k] %o% norm_B$U3_normed[,k])
}
mse_B <- mean((lambda_B - lambda_tensor)^2)
cat(sprintf("\nTest B — Tensor MSE vs lambda_tensor: %.6f\n", mse_B))

for (k in 1:K) {
  r1 <- which_rank(U1_true[,k], norm_B$U1_normed)
  r2 <- which_rank(U2_true[,k], norm_B$U2_normed)
  r3 <- which_rank(U3_true[,k], norm_B$U3_normed)
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
                          adj_LF_scale = FALSE,
                          verbose      = TRUE)

miss_eval <- evaluate_missing_prediction(fit_C, mask)
cat(sprintf("Test C — Missing-entry prediction MSE: %.6f\n", miss_eval$rmse))
