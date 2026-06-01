# Two-Step Covariate Selection via Unsupervised Decomposition

Identifies which covariates are relevant for each tensor component by
filtering active (non-zero) observations per factor, regressing
log-transformed loadings on a candidate covariate matrix using OLS, and
selecting covariates whose coefficients are significant at level
`alpha`. Zero loadings (spike/inactive observations) are excluded before
regression to match the log-linear DGP structure \\E\[u \mid x\] \propto
\exp(x^\top \gamma)\\. The loading matrix \\E\[l\_{ik}\]\\ can be
supplied in three ways: (1) provide `El` directly, (2) provide a fitted
`fit` object to extract `El` via
[`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md),
or (3) provide `X` and `K` so the function fits an unsupervised
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
internally (original two-step behaviour).

## Usage

``` r
select_covariates(
  X = NULL,
  K,
  covariate_data,
  fit = NULL,
  El = NULL,
  alpha = 0.05,
  verbose = TRUE,
  ...
)
```

## Arguments

- X:

  A 3-dimensional non-negative integer array of dimensions `n x p x w`.
  Required when neither `fit` nor `El` is supplied.

- K:

  Integer. Number of components (CP rank). When `El` is provided
  directly, `K` must equal `ncol(El)`.

- covariate_data:

  Numeric matrix of dimension `n x q` containing candidate covariates to
  screen. Column names are used in the output if available.

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  or
  [`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).
  If supplied, `El` is extracted via `normalize_factors(fit)$El` and no
  internal `CoSparTA` call is made. Either `fit` or `El` or neither
  (with `X` and `K`) must be provided.

- El:

  Numeric matrix of dimensions `n x K` (observation loading matrix). If
  provided, used directly — no `CoSparTA` fit is run and `fit` is
  ignored.

- alpha:

  Numeric significance level for covariate selection. Default `0.05`.

- verbose:

  Logical. If `TRUE`, prints progress and selected covariates per rank.
  Default `TRUE`.

- ...:

  Additional arguments passed to
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  when running the internal unsupervised fit (e.g., `init`, `maxiter`,
  `tol`, `convergence_criteria`). Ignored when `fit` or `El` is
  supplied.

## Value

A named list with:

- selected:

  A list of length K. Each element is an integer vector of column
  indices of `covariate_data` whose coefficients are significant at
  level `alpha` for that component. Empty integer vector if no
  covariates are significant.

- summaries:

  A list of length K. Each element is the `summary.lm` object from the
  OLS regression of El\[,k\] on `covariate_data`, giving full
  coefficient estimates, standard errors, t-statistics, and p-values.

- fit_unsupervised:

  The fitted unsupervised
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  object when the internal fit was run, the supplied `fit` object when
  that was used, or `NULL` when `El` was provided directly.

- runtime_secs:

  Numeric scalar giving the total elapsed wall-clock time of the
  function call in seconds.

## Details

When running the internal unsupervised fit, `Xcov = NULL` is used and
all `...` arguments are forwarded to
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).
For each component k, rows with \\l\_{ik} = 0\\ (spike/inactive
observations) are dropped, and OLS is applied to the log-transformed
loadings of the remaining active rows: \$\$\log E\[l\_{ik}\] =
X\_{\text{cov},i} \beta_k + \epsilon\_{ik}, \quad l\_{ik} \> 0\$\$ This
matches the log-linear DGP structure of the generative model. Covariates
are selected if their two-sided p-value is below `alpha`. The intercept
is included in the regression but is never selected as a covariate (it
is excluded from the returned indices).

This is a screening procedure, not a formal statistical test. For
rigorous inference on covariate effects, fit
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
with the selected covariates and use the estimated \\\gamma\\
coefficients from the generative model.

## See also

[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md),
[`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
X <- array(rpois(100 * 20 * 10, lambda = 1.5), dim = c(100, 20, 10))
Xcov <- matrix(rnorm(100 * 5), nrow = 100)
colnames(Xcov) <- paste0("cov", 1:5)

# Style 1: internal unsupervised fit (original behaviour)
result <- select_covariates(X, K = 3, covariate_data = Xcov,
                             maxiter = 20, convergence_criteria = 'ELBO')

# Style 2: supply a fitted object
fit <- CoSparTA(X, K = 3)
result2 <- select_covariates(X, K = 3, covariate_data = Xcov, fit = fit)

# Style 3: supply El directly
nf <- normalize_factors(fit)
result3 <- select_covariates(X, K = 3, covariate_data = Xcov, El = nf$El)

# Which covariates were selected for each rank?
result$selected

# Full regression summary for rank 1
result$summaries[[1]]

# Refit with selected covariates for rank 1
Xcov_selected <- Xcov[, result$selected[[1]], drop = FALSE]
fit <- CoSparTA(X, K = 3, Xcov = Xcov_selected)
} # }
```
