# CxtEBTD R Package — Progress Summary

## What is this package?
CxtEBTD (Covariate-aware Empirical Bayes Tensor Decomposition) is an R package implementing the SupEBTD algorithm for non-negative CP tensor decomposition on 3-way count tensors using Poisson likelihood, spike-and-slab priors, empirical Bayes hyperparameter estimation, and covariate-dependent factor priors. Primary applications: clickstream conversion prediction and PTSD/mental health behavioral classification. Companion software paper planned for INFORMS Journal on Computing (IJOC), Software Tools area.

## Package location
- Local: ~/Desktop/CxtEBTD/
- GitHub: https://github.com/xzhang0407/CxtEBTD (private)
- Active branch: dev (main kept clean)
- Current status: 0 errors, 0 warnings on devtools::check(); 22 exported functions across 10 R files + 2 C++ source files; 15 tests across 2 test scripts, all passing

## File structure
CxtEBTD/ ├── DESCRIPTION ├── NAMESPACE ├── README.md ├── LICENSE ├── test_simulation.R ← end-to-end test script (Tests A–H, all passing) ├── test_inference_edges.R ← edge-case tests (Tests I–O, all passing) └── R/ ├── supEBTD.R ← main function CxtEBTD() ├── ebpm_covariates.R ← novel covariate-aware EBPM ├── ebpm_wrappers.R ← UQ wrappers for external EBPM functions ├── utils.R ← sparse tensor ops, normalization helpers ├── internals.R ← shared internal helpers, EM building blocks ├── inference.R ← get_pip(), get_credible_interval(), get_significant_patterns(), get_posterior_quantile() ├── postprocessing.R ← normalize_factors(), project_tensor(), reconstruct_tensor(), init_cpapr(), select_covariates() ├── missing.R ← CxtEBTD_missing(), generate_missing_mask(), evaluate_missing_prediction() ├── visualization.R ← plot_time_factors(), plot_channel_factors() ├── rcpp_wrappers.R ← calc_EZ_3d_fast() thin wrapper for C++ backend ├── RcppExports.R ← auto-generated Rcpp bridge └── src/ ├── calc_EZ_3d_cpp.cpp ← C++ sparse weighted aggregation └── calc_qz_sparse_cpp.cpp ← C++ sparse softmax

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
| `plot_time_factors()` | Faceted line plot of time-mode factors |
| `plot_channel_factors()` | Faceted bar plot of channel-mode factors with optional grouping |
| `init_cpapr()` | CP-APR warm-start initialization via pyCP_APR/reticulate |
| `select_covariates()` | Two-step covariate screening: unsupervised fit → OLS regression |
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
9. Rcpp backends: `calc_EZ_3d` (sparse weighted aggregation) and `calc_qz_sparse` (softmax at non-zero entries) rewritten in C++. Eliminates dplyr overhead and avoids dense n×p×w×K array allocation. R fallbacks retained in codebase.

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
12. **Rcpp wrapper x$x vs x$v** — wrapper was passing NULL column instead of sparse count values. Fixed during validation.

## Test results (test_simulation.R — all 8 passing)
- Test A (unsupervised): MSE=0.000606
- Test B (supervised): MSE=0.000616, γ₁ slope ≈ 0.93
- Test C (missing): RMSE=1.178, MAE=1.029
- Test D (rank-specific covariates): rank 1 covariate_dependent, rank 2 unsupervised
- Test E (posterior quantiles): ordering invariant holds, F-mode correctly errors
- Test F (normalize_factors): unit norms, λ descending, reconstruction diff < 1e-15
- Test G (project_tensor): dims correct, single-obs consistent
- Test H (reconstruct_tensor): all non-negative, MSE=0.000616

## Test results (test_inference_edges.R — all passing)
- Test I (get_pip with threshold): dims correct, logical output, W-mode works, F-mode returns NULL
- Test J (get_credible_interval all modes): ordering invariant holds for L and W, F-mode returns CI
- Test K (get_significant_patterns): 10 active channels per factor, strict alpha reduces discoveries, W-only mode correct
- Test L (posterior quantiles W-mode): dims correct, ordering invariant
- Test M (K=1 single component): normalize, reconstruct, quantile all pass
- Test N (rank-specific covariates + missing data): rank 1 covariate_dependent, rank 2 unsupervised, RMSE=1.180
- Test O (project_tensor error handling): wrong-p and wrong-w correctly error

## Known remaining issues
- adj_LF_scale=TRUE still has a latent issue: gammaF computed twice (W never scaled). Low priority since we always use FALSE.

## Not yet implemented (planned for IJOC paper)
- Bootstrap/delta method CIs for gamma coefficients
- Rank selection utility (elbow on weights, prune_rank())
- Factor stability index via bootstrap
- Contrastive trait analysis (Algorithm 3, group comparison)
- User loading visualization (trait summary plot: channel bars + time line + PIP-overlaid user loadings)
- Sparsity warning function
- Pre-processing utility (build_tensor from raw event log)
- Unit tests with testthat
- Vignette with worked example on real data
- Sparse ELBO computation (calc_stm_obj Rcpp rewrite)

## Dissertation pipeline calling convention
Always called with:
- init = list(init_n, init_p, init_w)  ← from CP-APR via pyCP_APR
- adj_LF_scale = FALSE
- convergence_criteria = 'ELBO'
- U1_true, U2_true, U3_true supplied for simulation studies
- ebpm.fn = c(ebpm_point_gamma_multiplier_covariates, smashrgen::ebps, smashrgen::ebps)
