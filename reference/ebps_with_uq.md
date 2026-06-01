# Empirical Bayes Poisson Smoothing with Posterior Variance

A wrapper around [`ebps`](https://rdrr.io/pkg/smashrgen/man/ebps.html)
that augments the posterior output with posterior variance when
available. Posterior variance is computed from the log-normal formula
using the variational posterior parameters returned by `ebps`. Note that
`ebps` uses a Gaussian smoothing prior with no spike component, so
posterior inclusion probabilities (PIPs) are not applicable and returned
as `NA`.

## Usage

``` r
ebps_with_uq(x, s = NULL, ...)
```

## Arguments

- x:

  Non-negative integer vector of observed counts.

- s:

  Numeric scalar or vector of scaling factors. Default `NULL`.

- ...:

  Additional arguments passed to
  [`ebps`](https://rdrr.io/pkg/smashrgen/man/ebps.html) (e.g. `g_init`,
  `general_control`, `smooth_control`).

## Value

Same structure as [`ebps`](https://rdrr.io/pkg/smashrgen/man/ebps.html)
but with two additional columns in `posterior`:

- var:

  Posterior variance \\\text{Var}(\lambda_i \mid x_i)\\ computed via the
  log-normal formula \\\exp(2m_i + v_i)(\exp(v_i) - 1)\\ when `var_log`
  is available (i.e. `wave_trans = 'dwt'`). `NA` otherwise.

- pip:

  Not applicable for `ebps` (no spike-slab component). Always `NA`.

## See also

[`ebps`](https://rdrr.io/pkg/smashrgen/man/ebps.html),
[`ebpm_point_gamma_with_uq`](https://nd-hal.github.io/CoSparTA/reference/ebpm_point_gamma_with_uq.md)
