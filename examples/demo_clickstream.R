# =============================================================================
# CxtEBTD — Clickstream Practical Demo
# Synthetic data generated from real comScore clickstream structure.
# Tensor: 4000 sessions x 168 hours (last 7 days) x 141 channels
# Covariates: gender, age, race, household income, education, children
# =============================================================================

library(CxtEBTD)
library(dplyr)
library(readr)
library(reticulate)
library(ebpm)
library(smashrgen)
library(ggplot2)
library(ggh4x)

# =============================================================================
# Python setup (required for CP-APR warm start)
# Requires Python 3.9+ with pyCP_APR and numpy installed:
#   python3 -m venv cxtebtd_env
#   source cxtebtd_env/bin/activate
#   pip install pyCP_APR numpy
# =============================================================================
use_virtualenv("cxtebtd_env", required = TRUE)

# =============================================================================
# STEP 1: Load data
# =============================================================================
Xtensor     <- readRDS("data/clickstream_synth_tensor.rds")
X_cov       <- read_csv("data/clickstream_synth_cov.csv", show_col_types = FALSE)
channel_info <- readRDS("data/clickstream_channel_names.rds")

cat(sprintf("Tensor dimensions: %d x %d x %d\n",
            dim(Xtensor)[1], dim(Xtensor)[2], dim(Xtensor)[3]))
cat(sprintf("Sparsity: %.4f\n", mean(Xtensor == 0)))

# =============================================================================
# STEP 2: CP-APR warm start
# =============================================================================
K       <- 10
maxiter <- 100

init_list <- init_cpapr(Xtensor, K = K, virtualenv = "cxtebtd_env",
                         method = "torch")

# =============================================================================
# STEP 3: Fit unsupervised CxtEBTD
# =============================================================================
cat("\nFitting unsupervised CxtEBTD...\n")
st <- Sys.time()
fit_ebtd <- CxtEBTD(
  X                    = Xtensor,
  K                    = K,
  Xcov                 = NULL,
  init                 = init_list,
  maxiter              = maxiter,
  convergence_criteria = "ELBO",
  tol                  = 1e-5,
  verbose              = TRUE
)
cat(sprintf("Unsupervised fit time: %s\n", format(Sys.time() - st)))

fit_ebtd_normed <- normalize_factors(fit_ebtd)

# =============================================================================
# STEP 4: Covariate selection
# =============================================================================
cat("\nSelecting covariates...\n")
result <- select_covariates(
  K              = K,
  covariate_data = X_cov,
  El             = fit_ebtd_normed$El
)

Xcov_list <- lapply(result$selected, function(idx) {
  if (length(idx) == 0) NULL else as.matrix(X_cov[, idx, drop = FALSE])
})

# =============================================================================
# STEP 5: Fit supervised CxtEBTD
# =============================================================================
cat("\nFitting supervised CxtEBTD...\n")
init_list2 <- list(fit_ebtd_normed$El, fit_ebtd_normed$Ef, fit_ebtd_normed$Ew)
st <- Sys.time()
fit_cxtebtd <- CxtEBTD(
  X                    = Xtensor,
  K                    = K,
  Xcov                 = Xcov_list,
  init                 = init_list2,
  maxiter              = maxiter,
  convergence_criteria = "ELBO",
  tol                  = 1e-5,
  verbose              = TRUE,
  ebpm.fn              = c(ebpm_point_gamma_multiplier_covariates,
                           smashrgen::ebps,
                           ebpm_point_gamma)
)
cat(sprintf("Supervised fit time: %s\n", format(Sys.time() - st)))

fit_cxtebtd_normed <- normalize_factors(fit_cxtebtd)

# =============================================================================
# STEP 6: Visualization
# =============================================================================

# --- Covariate effects (gamma) ---
all_covs <- colnames(X_cov)
df <- bind_rows(lapply(1:K, function(k) {
  gamma_full <- setNames(rep(0, length(all_covs)), all_covs)
  selected_names <- colnames(X_cov)[result$selected[[k]]]
  if (length(selected_names) > 0)
    gamma_full[selected_names] <- fit_cxtebtd_normed$gamma_list[[k]]
  data.frame(factor = paste0("F", k), covariate = all_covs, gamma = gamma_full)
})) %>%
  mutate(
    factor    = factor(factor, levels = paste0("F", 1:K)),
    covariate = factor(covariate, levels = all_covs),
    sign      = ifelse(gamma >= 0, "pos", "neg")
  )

p_gamma <- ggplot(df, aes(x = covariate, y = gamma, fill = sign)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_grid(factor ~ ., scales = "free_y") +
  scale_fill_manual(values = c("pos" = "#E57373", "neg" = "#64B5F6")) +
  scale_x_discrete(labels = c("Gender", "Age", "Race",
                               "Income", "Edu", "Children")) +
  labs(x = "Covariate", y = "Gamma") +
  theme_minimal() +
  theme(
    legend.position  = "none",
    strip.text.y     = element_text(face = "bold", size = 10),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 7),
    panel.grid.minor = element_blank()
  ) +
  force_panelsizes(rows = unit(rep(0.89, K), "cm"), TRUE)
print(p_gamma)

# --- Temporal patterns ---
p_time <- plot_time_factors(
  Ef          = fit_cxtebtd_normed$Ef,
  time_labels = 169:336,
  xlim        = c(169, 336)
)
print(p_time)

# --- Channel patterns ---
p_channel <- plot_channel_factors(
  Ew             = fit_cxtebtd_normed$Ew,
  channel_names  = channel_info$channel_names,
  channel_groups = channel_info$channel_groups
)
print(p_channel)

cat("\nDemo complete.\n")
