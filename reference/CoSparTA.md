# Covariate-aware Empirical Bayes Tensor Decomposition

Fits a non-negative CP tensor decomposition using empirical Bayes priors
with spike-and-slab structure and Poisson likelihood. When covariates
are supplied via `Xcov`, they enter the generative model directly
through covariate-dependent Poisson rates of the form \\\lambda_i =
\beta \cdot \exp(X_i^\top \gamma)\\, enabling simultaneous decomposition
and covariate effect estimation. Without covariates (`Xcov = NULL`), the
method reduces to unsupervised empirical Bayes Poisson tensor
decomposition.

## Usage

``` r
CoSparTA(
  X,
  K,
  Xcov = NULL,
  lib_size = NULL,
  init = "random_gamma",
  maxiter = 100,
  maxiter_init = 100,
  tol = 1e-06,
  compute_elbo_final = FALSE,
  ebpm.fn = c(ebpm_point_gamma_multiplier_covariates, ebps_with_uq,
    ebpm_point_gamma_with_uq),
  fix_L = FALSE,
  fix_F = FALSE,
  fix_W = FALSE,
  smooth_F = T,
  printevery = 10,
  verbose = TRUE,
  convergence_criteria = "factor_change",
  n_stable = 3L,
  U1_true = NULL,
  U2_true = NULL,
  U3_true = NULL
)
```

## Arguments

- X:

  A 3-dimensional non-negative integer array of dimensions `n x p x w`,
  representing observations (e.g., users), time points, and channels
  respectively.

- K:

  Integer. Number of components (CP rank) for the decomposition.

- Xcov:

  Covariate input for the observation mode. Can be: (1) `NULL` for fully
  unsupervised decomposition; (2) a numeric matrix of dimension `n x q`,
  applied identically to all K components; or (3) a list of length `K`,
  where each element is either a numeric covariate matrix (dimensions
  `n x q_k`, potentially different numbers of covariates per component)
  or `NULL` for unsupervised components. Default `NULL`.

- lib_size:

  Numeric vector of length `n` giving per-observation library sizes
  (exposure/scaling factors). Default `NULL`, which sets all sizes to 1.

- init:

  Character string specifying the initialization method. Default
  `'random_gamma'`.

- maxiter:

  Integer. Maximum number of EM iterations. Default `100`.

- maxiter_init:

  Integer. Maximum number of iterations for the initialization step.
  Default `100`.

- tol:

  Numeric. Convergence tolerance. Iterations stop when the change in the
  objective is below this threshold. Default `1e-6`.

- compute_elbo_final:

  Logical. If `TRUE`, computes the ELBO after the final iteration even
  when `convergence_criteria != 'ELBO'`. Default `FALSE`.

- ebpm.fn:

  A single function or list of three functions specifying the empirical
  Bayes prior for the L, F, and W modes respectively. Default uses
  `ebpm_point_gamma_multiplier_covariates` for L (covariate-aware),
  `ebps_with_uq` for F (smooth), and `ebpm_point_gamma_with_uq` for W
  (point-gamma), all returning posterior variance and PIP where
  applicable.

- fix_L:

  Logical. If `TRUE`, the observation-mode loadings are held fixed at
  initialization. Default `FALSE`.

- fix_F:

  Logical. If `TRUE`, the time-mode factors are held fixed. Default
  `FALSE`.

- fix_W:

  Logical. If `TRUE`, the channel-mode weights are held fixed. Default
  `FALSE`.

- smooth_F:

  Logical. If `TRUE`, applies smoothing to time-mode factors. Default
  `TRUE`.

- printevery:

  Integer. Number of iterations between progress messages when
  `verbose = TRUE`. Default `10`.

- verbose:

  Logical. If `TRUE`, prints initialization and iteration progress.
  Default `TRUE`.

- convergence_criteria:

  Character string specifying the convergence criterion:
  `'factor_change'` (maximum column-normalized change across all three
  factor matrices, default), `'mKLabs'` (mean KL divergence), or
  `'ELBO'` (evidence lower bound). `'factor_change'` is recommended: it
  avoids the O(npwK) ELBO computation and typically converges in fewer
  iterations.

- n_stable:

  Integer. Number of consecutive iterations with factor change below
  `tol` required before declaring convergence. Only used when
  `convergence_criteria = 'factor_change'`. Default `3`.

- U1_true:

  Optional matrix of true observation-mode factors, used to track
  reconstruction error during simulation studies. Default `NULL`.

- U2_true:

  Optional matrix of true time-mode factors for simulation evaluation.
  Default `NULL`.

- U3_true:

  Optional matrix of true channel-mode factors for simulation
  evaluation. Default `NULL`.

## Value

A named list with the following elements:

- elbo:

  Final ELBO value computed after the last iteration.

- obj_trace:

  Numeric vector of objective values at each iteration.

- res:

  List of variational posterior summaries. Key fields: `res$ql$El` —
  posterior mean loadings (`n x K`); `res$ql$Elogl` — posterior log-mean
  loadings (`n x K`); `res$ql$Varl` — posterior variance of loadings
  (`n x K`); `res$ql$PIPl` — posterior inclusion probabilities for
  loadings (`n x K`), giving \\P(\theta_i \neq 0 \mid x_i)\\ for each
  observation and component; `res$qf$Ef` — posterior mean time factors
  (`p x K`); `res$qf$Elogf` — posterior log-mean time factors (`p x K`);
  `res$qw$Ew` — posterior mean channel weights (`w x K`); `res$qw$Elogw`
  — posterior log-mean channel weights (`w x K`).

- diff_U:

  List of three vectors recording per-iteration reconstruction error
  relative to `U1_true`, `U2_true`, `U3_true`. Only meaningful when true
  factors are supplied.

- run_time:

  Elapsed computation time as a `difftime` object.

In addition, the following normalized and reordered fields are populated
automatically (columns scaled to unit Frobenius norm, ordered by
descending \\\lambda\\):

- res\$lambda_normed:

  Numeric vector of length `K`: component weights \\\lambda_k =
  \\l_k\\\\f_k\\\\w_k\\\\, sorted descending.

- res\$gl_normed:

  List of K gamma estimates reordered by descending \\\lambda\\.

- res\$ql\$El_normed:

  Unit-norm observation loading matrix (`n x K`), reordered by
  descending \\\lambda\\.

- res\$ql\$PIPl_normed:

  PIP matrix for L mode, reordered.

- res\$ql\$shape_post_l_normed:

  Posterior gamma shape for L, reordered.

- res\$ql\$rate_post_l_normed:

  Posterior gamma rate for L scaled by \\\\l_k\\\\, reordered, so that
  \\El\\normed\[i,k\] \sim \text{Gamma}(\text{shape},\\
  \text{rate\\normed})\\ marginally.

- res\$qf\$Ef_normed, res\$qf\$PIPf_normed,
  res\$qf\$shape_post_f_normed, res\$qf\$rate_post_f_normed:

  Same pattern for the F mode.

- res\$qw\$Ew_normed, res\$qw\$PIPw_normed,
  res\$qw\$shape_post_w_normed, res\$qw\$rate_post_w_normed:

  Same pattern for the W mode.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
X <- array(rpois(20 * 12 * 4, lambda = 1.5), dim = c(20, 12, 4))
Xcov <- matrix(rnorm(20 * 2), nrow = 20, ncol = 2)

# Supervised decomposition with covariates
fit <- CoSparTA(X, K = 3, Xcov = Xcov, maxiter = 50, verbose = FALSE)
EL <- fit$res$ql$El   # 20 x 3 loading matrix
EF <- fit$res$qf$Ef   # 12 x 3 factor matrix
EW <- fit$res$qw$Ew   # 4 x 3 weight matrix

# Unsupervised decomposition (no covariates)
fit0 <- CoSparTA(X, K = 3, Xcov = NULL, maxiter = 50, verbose = FALSE)
} # }
```
