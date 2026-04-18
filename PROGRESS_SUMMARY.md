# CxtEBTD R Package — Progress Summary

## What is this package?
CxtEBTD (Covariate-aware Empirical Bayes Tensor Decomposition) is an R package implementing the SupEBTD algorithm for non-negative CP tensor decomposition on 3-way count tensors using Poisson likelihood, spike-and-slab priors, empirical Bayes hyperparameter estimation, and covariate-dependent factor priors. Primary applications: clickstream conversion prediction and PTSD/mental health behavioral classification. Companion software paper planned for INFORMS Journal on Computing (IJOC), Software Tools area.

## Package location
- Local: ~/Desktop/CxtEBTD/
- GitHub: https://github.com/xzhang0407/CxtEBTD (private)
- Active branch: dev (main kept clean)
- Current status: 0 errors, 0 warnings on devtools::check()

## File structure
CxtEBTD/ ├── DESCRIPTION ├── NAMESPACE ├── README.md ├── LICENSE ├── test_simulation.R ← end-to-end test script (Tests A–H, all passing) └── R/ ├── supEBTD.R ← main function CxtEBTD() ├── ebpm_covariates.R ← novel covariate-aware EBPM ├── ebpm_wrappers.R ← UQ wrappers for external EBPM functions ├── utils.R ← sparse tensor ops, normalization helpers ├── internals.R ← shared internal helpers, EM building blocks ├── inference.R ← get_pip(), get_credible_interval(), get_significant_patterns(), get_posterior_quantile() ├── postprocessing.R ← normalize_factors(), project_tensor(), reconstruct_tensor() └── missing.R ← CxtEBTD_missing(), generate_missing_mask(), evaluate_missing_prediction()

## Exported functions
| Function | Purpose |
|----------|---------|
| `CxtEBTD()` | Main decomposition — with or without covariates |
| `CxtEBTD_missing()` | Decomposition with missing data support |
| `ebpm_point_gamma_multiplier_covariates()` | Novel covariate-aware EBPM prior |
| `ebpm_point_gamma_with_uq()` | Wrapper adding var+PIP to ebpm::ebpm_point_gamma |
| `ebps_with_uq()` | Wrapper adding var to smashrgen::ebps |
| `get_pip()` | Extract posterior inclusion probabilities |
| `get_credible_interval()` | Compute credible intervals from posterior var |
| `get_significant_patterns()` | lFDR-based pattern discovery (Algorithms 1+2) |
| `get_posterior_quantile()` | Exact quantiles from spike-and-slab mixture posterior |
| `normalize_factors()` | Normalize columns to unit Frobenius norm, compute component weights λ |
| `project_tensor()` | Project new tensor data onto learned factors (Eq. 6 in ISR paper) |
| `reconstruct_tensor()` | Reconstruct denoised mean tensor from fitted factors |
| `generate_missing_mask()` | Simulate missing data for evaluation |
| `evaluate_missing_prediction()` | Evaluate imputation quality |
| `adjLF()` | Scale loadings/factors to similar norms |
| `mKL()` | Mean KL divergence |
| `poisson_to_multinom()` | Standardize Poisson factorization |
| `log_for_ebmf()` | Log-transform count matrix |
| `ebpmf_identity_smooth_control_default()` | Default smooth control parameters |

## Key architecture decisions
1. `CxtEBTD(X, K, Xcov = NULL)` — unified function, Xcov=NULL triggers unsupervised path
2. `CxtEBTD_missing()` — separate function, same interface, adds obs_mask argument
3. Default `ebpm.fn = c(ebpm_point_gamma_multiplier_covariates, ebps_with_uq, ebpm_point_gamma_with_uq)`
4. Default `init = 'random_gamma'` — Gamma(shape=100, rate=100), BPTF-style
5. Always run with `adj_LF_scale = FALSE` (default TRUE has a latent bug with Ef_smooth, partially guarded)
6. Always run with `convergence_criteria = 'ELBO'` to match dissertation pipeline
7. Rank-specific covariates: Xcov and ebpm.fn.l normalized to length-K lists early in CxtEBTD()/CxtEBTD_missing(). Per-rank dispatch with auto-fallback for unsupervised ranks.
8. Posterior Gamma shape/rate threaded through all modes via lazy-init. Enables exact quantile computation.

## Bugs fixed (as of April 2026)
1. **init='random_gamma' implemented** — replaces dead uniform/fasttopics branches. Gamma(100,100) init for L, F, W.
2. **diff_U NULL guard** — iteration tracking block now only runs when U1_true/U2_true/U3_true supplied.
3. **adj_LF_scale Ef_smooth crash** — guarded lines referencing Ef_smooth/Elogf_smooth which are never populated.
4. **Unsupervised L-mode fallback** — Xcov=NULL now correctly falls back to ebpm::ebpm_point_gamma. Previously broken.
5. **Varl/PIPl/Varf/PIPf/Varw/PIPw lazy initialization** — matrices initialized before first column assignment.
6. **Var/PIP dimension restoration** — padded back to original dims at wrap-up to match El/Ef/Ew.
7. **Rfast::rowsums namespace** — explicit qualification added.
8. **dplyr/data.table namespace** — dplyr::left_join, dplyr::group_by, dplyr::summarize, @import data.table, @importFrom dplyr %>% added.
9. **Unsupervised fallback missing in CxtEBTD_missing()** — added ebpm::ebpm_point_gamma fallback when Xcov=NULL.
10. **Var/PIP dimension restoration missing in CxtEBTD_missing()** — added zero-row padding at wrap-up.
11. **NA guard in get_posterior_quantile()** — guarded against NA-padded rows from all-zero user restoration and scalar NA shape_post from ebps F-mode.

## Test results (test_simulation.R — all 8 passing)
- Test A (unsupervised): MSE=0.000606
- Test B (supervised): MSE=0.000616, γ₁ slope ≈ 0.93
- Test C (missing): RMSE=1.178, MAE=1.029
- Test D (rank-specific covariates): rank 1 covariate_dependent, rank 2 unsupervised
- Test E (posterior quantiles): ordering invariant holds, F-mode correctly errors
- Test F (normalize_factors): unit norms, λ descending, reconstruction diff < 1e-15
- Test G (project_tensor): dims correct, single-obs consistent
- Test H (reconstruct_tensor): all non-negative, MSE=0.000616

## Known remaining issues
- adj_LF_scale=TRUE still has a latent issue: gammaF computed twice (W never scaled). Low priority since we always use FALSE.

## Not yet implemented (planned for IJOC paper)
- Bootstrap/delta method CIs for gamma coefficients
- Rank selection utility (elbow on weights, prune_rank())
- Factor stability index via bootstrap
- Contrastive trait analysis (Algorithm 3, group comparison)
- Visualization functions (trait summary plot: channel bars + time line + PIP-overlaid user loadings)
- Sparsity warning function
- Unit tests with testthat
- Vignette with worked example on real data

## Dissertation pipeline calling convention
Always called with:
- init = list(init_n, init_p, init_w)  ← from CP-APR via pyCP_APR
- adj_LF_scale = FALSE
- convergence_criteria = 'ELBO'
- U1_true, U2_true, U3_true supplied for simulation studies
- ebpm.fn = c(ebpm_point_gamma_multiplier_covariates, smashrgen::ebps, smashrgen::ebps)
