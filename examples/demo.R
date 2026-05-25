# =============================================================================
# CoSparTA — end-to-end pipeline demo
# Demonstrates the full analysis workflow from raw event logs to inference.
# Data files are in the data/ directory of this package.
# =============================================================================

library(CoSparTA)
library(dplyr)
library(readr)
library(reticulate)
library(ebpm)
library(smashrgen)
library(ggplot2)
library(ggh4x)

# =============================================================================
# STEP 1: Build the tensor
# =============================================================================

# Load raw data
df             <- readRDS("data/demo_df.rds")
cov_df         <- read.csv("data/demo_covariates.csv")
channel_groups <- readRDS("data/demo_channel_groups.rds")

# Covariate matrix
Xcov_mat <- as.matrix(cov_df[, c("cov1", "cov2")])

# Reconstruct universe of labels
session_ids   <- sort(unique(df$session_id))
hour_labels   <- paste0("H", sprintf("%02d", 0:99))
channel_names <- c("Retail", "Social", "Review", "Deal", "Search")
website_names <- as.vector(
  outer(channel_names, paste0("_", sprintf("%02d", 1:10)), paste0))

# Build 1000 x 100 x 50 tensor
tensor_out <- build_tensor(
  data         = df,
  row          = "session_id",
  col          = "hour",
  slice        = "website",
  value        = "count",
  row_levels   = session_ids,
  col_levels   = hour_labels,
  slice_levels = website_names
)
X <- tensor_out$X   # 1000 x 100 x 50 integer tensor

cat("Dimensions:", dim(X), "\n")
cat("Sparsity:  ", round(mean(X == 0), 4), "\n")

# =============================================================================
# STEP 2: Fit the model
# =============================================================================


# Optional CP-APR warm-start
# Requires a Python 3.9 virtual environment with pyCP_APR and numpy:
#   python3.9 -m venv cxtebtd_env
#   source cxtebtd_env/bin/activate        # Mac/Linux
#   pip install pyCP_APR numpy
init_vals <- init_cpapr(X, K = 4, virtualenv = "cxtebtd_env")

# Supervised fit: shared covariate matrix across all factors
fit <- CoSparTA(
  X                    = X,
  K                    = 4,
  Xcov                 = Xcov_mat,
  init                 = init_vals,  # or "random_gamma" if Python unavailable
  maxiter              = 20,
  convergence_criteria = "ELBO",
  tol                  = 1e-6,
  verbose              = TRUE
)
# Key outputs:
# fit$res$ql$El   -- 1000 x 4 posterior mean loadings (U1)
# fit$res$qf$Ef   -- 100  x 4 time factors            (U2)
# fit$res$qw$Ew   -- 50   x 4 channel weights         (U3)
# fit$res$gl      -- list of K covariate coefficient vectors

# Rank-specific covariate list: both covariates, NULL, cov1 only, cov2 only
Xcov_cov1 <- Xcov_mat[, 1, drop = FALSE]
Xcov_cov2 <- Xcov_mat[, 2, drop = FALSE]
fit_rankspec <- CoSparTA(
  X                    = X,
  K                    = 4,
  Xcov                 = list(Xcov_mat, NULL, Xcov_cov1, Xcov_cov2),
  init                 = init_vals,
  maxiter              = 20,
  convergence_criteria = "ELBO",
  tol                  = 1e-6,
  verbose              = FALSE
)
cat("Rank-specific fit done\n")

# Missing data handling
mask     <- generate_missing_mask(X, missing_rate = 0.10)
fit_miss <- CoSparTA_missing(
  X                    = X,
  K                    = 4,
  Xcov                 = Xcov_mat,
  obs_mask             = mask$obs_mask,
  maxiter              = 20,
  convergence_criteria = "ELBO",
  tol                  = 1e-6,
  verbose              = FALSE
)
cat("Missing data fit done\n")

# =============================================================================
# STEP 3: Uncertainty quantification
# =============================================================================

# Exact 95% posterior credible intervals for channel mode (U3)
ci_w <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = "W")
# ci_w$q2.5  -- 50 x 4 lower bounds
# ci_w$q97.5 -- 50 x 4 upper bounds
cat("Channel mode 95% CI computed for",
    nrow(ci_w$q2.5), "websites across",
    ncol(ci_w$q2.5), "factors\n")

# Exact 95% posterior credible intervals for observation mode (U1)
ci_l <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = "L")
# ci_l$q2.5  -- 1000 x 4 lower bounds
# ci_l$q97.5 -- 1000 x 4 upper bounds
cat("Session mode 95% CI computed for",
    nrow(ci_l$q2.5), "sessions across",
    ncol(ci_l$q2.5), "factors\n")

# Delta-method CIs for covariate coefficients
ci_gamma <- get_gamma_ci(fit, method = "delta", level = 0.95)
print(ci_gamma[[1]])   # factor 1 results

# =============================================================================
# STEP 4: Covariate screening
# =============================================================================

# Fit covariate-free model
fit_unsup <- CoSparTA(
  X                    = X,
  K                    = 4,
  Xcov                 = NULL,
  init                 = init_vals,
  maxiter              = 20,
  convergence_criteria = "ELBO",
  tol                  = 1e-6,
  verbose              = FALSE
)

# Normalize and screen covariates
nf_unsup <- normalize_factors(fit_unsup)
sel <- select_covariates(
  K              = 4,
  covariate_data = as.data.frame(Xcov_mat),
  El             = nf_unsup$El
)
print(sel$selected)

# =============================================================================
# STEP 5: Post-processing
# =============================================================================

# Normalize to unit norm, sort by weight lambda
nf <- normalize_factors(fit)
cat("Lambda (weights):", round(nf$lambda, 3), "\n")

# Project new observations onto the factor space (no refit)
X_new <- X[991:1000, , ]   # held-out sessions for illustration
F_new <- project_tensor(X_new, fit)
cat("Projection dimensions:", dim(F_new), "\n")

# Reconstruct the denoised Poisson mean tensor
X_hat <- reconstruct_tensor(fit)
cat("Reconstruction dimensions:", dim(X_hat), "\n")
cat("Any negative values:      ", any(X_hat < 0), "\n")

# =============================================================================
# STEP 6: Visualization
# =============================================================================

# Faceted line plot of K time factors
plot_time_factors(
  Ef          = nf$Ef,
  time_labels = 1:100
)

# Faceted bar plot of K channel factors with category grouping
plot_channel_factors(
  Ew             = nf$Ew,
  channel_names  = website_names,
  channel_groups = channel_groups
)

cat("\nDemo complete.\n")
