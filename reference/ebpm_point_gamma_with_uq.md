# Empirical Bayes Point-Gamma with Posterior Variance and PIP

A wrapper around
[`ebpm_point_gamma`](https://rdrr.io/pkg/ebpm/man/ebpm_point_gamma.html)
that extends the posterior output with posterior variance and posterior
inclusion probabilities (PIPs). All optimization is delegated to the
original function; only the posterior summary is augmented.

## Usage

``` r
ebpm_point_gamma_with_uq(x, s = 1, ...)
```

## Arguments

- x:

  Non-negative integer vector of observed counts.

- s:

  Numeric scalar or vector of length `length(x)` giving per-observation
  exposure/scaling factors. Default `1`.

- ...:

  Additional arguments passed to
  [`ebpm_point_gamma`](https://rdrr.io/pkg/ebpm/man/ebpm_point_gamma.html)
  (e.g. `g_init`, `fix_g`, `pi0`, `control`).

## Value

Same structure as
[`ebpm_point_gamma`](https://rdrr.io/pkg/ebpm/man/ebpm_point_gamma.html)
but with four additional columns in `posterior`:

- var:

  Posterior variance \\\text{Var}(\theta_i \mid x_i)\\, decomposed as
  within-component variance plus between-component variance from spike
  uncertainty.

- pip:

  Posterior inclusion probability \\P(\theta_i \neq 0 \mid x_i) = 1 -
  \hat{\pi}\_i\\.

- shape_post:

  Posterior gamma shape parameter \\a + x_i\\.

- rate_post:

  Posterior gamma rate parameter \\b + s_i\\.

## See also

[`ebpm_point_gamma`](https://rdrr.io/pkg/ebpm/man/ebpm_point_gamma.html),
[`ebpm_point_gamma_multiplier_covariates`](https://nd-hal.github.io/CoSparTA/reference/ebpm_point_gamma_multiplier_covariates.md)
