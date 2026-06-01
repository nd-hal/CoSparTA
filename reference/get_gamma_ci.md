# Confidence Intervals for Covariate Coefficients

Computes confidence intervals for the covariate effect parameters
\\\gamma\\ estimated by
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).
Two methods are available:

- `"delta"`:

  Fast asymptotic intervals based on the Hessian of the marginal
  log-likelihood from the EBPM optimization. The full Hessian (covering
  \\\gamma\\ together with nuisance parameters \\\pi_0\\, \\\alpha\\,
  \\\beta\\) is inverted first, and then the \\\gamma\\ submatrix is
  extracted, implementing \\\mathrm{Var}(\hat\gamma_k) =
  \[H_k^{-1}\]\_{\gamma,\gamma}\\ per equation (21) of the paper. This
  correctly accounts for uncertainty in the nuisance parameters. These
  are **conditional** standard errors – they measure uncertainty in
  \\\gamma\\ given the current F and W estimates, and may underestimate
  the true uncertainty. Recommended for exploratory screening.

- `"bootstrap"`:

  Parametric bootstrap intervals. Generates `B` synthetic tensors from
  the fitted Poisson rates, refits
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  on each, and collects the bootstrap distribution of \\\gamma\\. This
  correctly accounts for all sources of estimation uncertainty
  (including uncertainty in F and W). Recommended for
  publication-quality inference. Computationally expensive: requires `B`
  full model refits.

## Usage

``` r
get_gamma_ci(
  fit,
  method = "delta",
  level = 0.95,
  B = 200,
  X = NULL,
  K = NULL,
  Xcov = NULL,
  normalized = TRUE,
  init_fn = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).

- method:

  Character: `"delta"` or `"bootstrap"`. Default `"delta"`.

- level:

  Numeric confidence level in `(0, 1)`. Default `0.95`.

- B:

  Integer. Number of bootstrap replicates (only used when
  `method = "bootstrap"`). Default `200`.

- X:

  The original tensor used to fit the model (required for
  `method = "bootstrap"`).

- K:

  Integer. Number of components (required for bootstrap).

- Xcov:

  Covariate matrix or list used in the original fit (required for
  bootstrap).

- normalized:

  Logical. If `TRUE` (default), reads gamma estimates from
  `fit$res$gl_normed` (reordered by descending \\\lambda\\). If `FALSE`,
  reads from `fit$res$gl` (original component order).

- init_fn:

  Optional function with signature `function(X_star, K)` that returns an
  initialization object (e.g., a list of three matrices from
  [`init_cpapr`](https://nd-hal.github.io/CoSparTA/reference/init_cpapr.md)).
  If provided, this function is called on each bootstrap replicate to
  generate a fresh initialization. If `NULL`, the `init` argument from
  `...` is reused for all replicates. Using per-replicate initialization
  (e.g., via `init_cpapr`) is recommended for publication-quality
  bootstrap inference.

- verbose:

  Logical. If `TRUE`, prints bootstrap progress. Default `TRUE`.

- ...:

  Additional arguments passed to
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  during bootstrap refitting (e.g., `maxiter`, `init`,
  `convergence_criteria`).

## Value

A list of length K (one element per component). Each element is a named
list with:

- estimate:

  Numeric vector of gamma estimates from the original fit.

- se:

  Standard errors.

- lower:

  Lower confidence bound.

- upper:

  Upper confidence bound.

- pvalue:

  Two-sided p-values (Wald test for delta method; bootstrap
  percentile-based for bootstrap, computed as
  `2 * min(prop >= 0, prop <= 0)`).

- method:

  Character string indicating which method was used.

If a component was fitted without covariates (unsupervised rank), its
entry is `NULL`.

## See also

[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md),
[`ebpm_point_gamma_multiplier_covariates`](https://nd-hal.github.io/CoSparTA/reference/ebpm_point_gamma_multiplier_covariates.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3, Xcov = Xcov, maxiter = 50)

# Fast delta method (exploratory)
ci_delta <- get_gamma_ci(fit, method = "delta")

# Parametric bootstrap (publication quality)
ci_boot <- get_gamma_ci(fit, method = "bootstrap", B = 200,
                         X = X, K = 3, Xcov = Xcov,
                         maxiter = 50, convergence_criteria = "ELBO")
} # }
```
