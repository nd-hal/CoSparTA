# Empirical Bayes Poisson Mean Estimation with Covariate-Dependent Rates

Fits a spike-and-slab point-gamma prior to a vector of counts `x` where
the Poisson rate is modulated by individual-level covariates. The
effective rate for observation \\i\\ is \\\lambda_i = \beta \cdot
\exp(X_i^\top \gamma)\\, so that the prior mean scales multiplicatively
with covariates. Prior parameters \\(\pi_0, \alpha, \beta, \gamma)\\ are
estimated by maximizing the marginal likelihood via `nlm` (with `optim`
as fallback). This function is the core building block of the
covariate-aware mode update in
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).

## Usage

``` r
ebpm_point_gamma_multiplier_covariates(
  x,
  s = 1,
  X,
  g_init = NULL,
  control = NULL
)
```

## Arguments

- x:

  Non-negative integer vector of observed counts.

- s:

  Numeric scalar or vector of length `length(x)` giving per-observation
  exposure/scaling factors. Default `1`.

- X:

  Numeric covariate matrix of dimension `n x q`, where rows correspond
  to observations in `x` and columns to covariates.

- g_init:

  Optional numeric vector of starting values in the order
  `c(pi0, alpha, beta, gamma_1, ..., gamma_q)`. If `NULL`,
  initialization uses a baseline
  [`ebpm::ebpm_point_gamma`](https://rdrr.io/pkg/ebpm/man/ebpm_point_gamma.html)
  fit with `gamma = 0`. Default `NULL`.

- control:

  Optional list of control arguments passed to `nlm`. Defaults:
  `stepmax = 1`, `gradtol = 1e-6`, `steptol = 1e-8`, `iterlim = 1000`,
  `check.analyticals = FALSE`.

## Value

A named list with:

- fitted_g:

  List of estimated prior parameters: `pi0` — spike (point-mass at zero)
  probability; `shape` — gamma shape parameter \\\alpha\\; `scale` —
  gamma scale parameter \\1/\beta\\; `gamma` — covariate coefficient
  vector of length `q`; `type` — always `"covariate_dependent"`.

- posterior:

  Data frame with one row per observation and columns: `mean` —
  posterior mean \\E\[\theta_i \mid x_i\]\\; `mean_log` — posterior
  log-mean \\E\[\log\theta_i \mid x_i\]\\ (`-Inf` for observations with
  \\x_i = 0\\); `var` — posterior variance \\\text{Var}(\theta_i \mid
  x_i)\\, decomposed as within-component variance (uncertainty in the
  gamma component) plus between-component variance (uncertainty about
  whether the spike fired); `pip` — posterior inclusion probability
  \\P(\theta_i \neq 0 \mid x_i) = 1 - \hat{\pi}\_i\\, giving the
  probability that observation \\i\\ has a truly non-zero loading;
  `shape_post` — posterior gamma shape parameter \\\alpha + x_i\\, the
  shape of the gamma kernel in the non-spike component; `rate_post` —
  posterior gamma rate parameter \\\beta / \lambda_i + s_i\\, the rate
  of the gamma kernel.

- log_likelihood:

  Maximized marginal log-likelihood.

- convergence_code:

  Optimizer convergence code. For `nlm`: 1–2 indicate convergence; 3–5
  indicate potential issues (see
  [`?nlm`](https://rdrr.io/r/stats/nlm.html)). For `optim` fallback: 0 =
  converged.

## See also

[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md),
[`ebpm_point_gamma`](https://rdrr.io/pkg/ebpm/man/ebpm_point_gamma.html)

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 200
X <- matrix(rnorm(n * 2), nrow = n)
true_gamma <- c(0.4, -0.3)
x <- rpois(n, lambda = exp(0.5 + X %*% true_gamma))

fit <- ebpm_point_gamma_multiplier_covariates(x, s = 1, X = X)
fit$fitted_g$gamma     # estimated covariate coefficients
fit$fitted_g$pi0       # estimated spike probability
head(fit$posterior)    # posterior means and log-means
} # }
```
