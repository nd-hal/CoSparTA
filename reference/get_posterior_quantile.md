# Compute Exact Posterior Quantiles from a CoSparTA fit

Computes exact quantiles of the marginal posterior distribution for each
factor element. Dispatches on the posterior family tag stored in the fit
object (`fit$res$q*$family_*`) so that user-specified `ebpm.fn`
configurations produce correct intervals regardless of prior family:

- `"point_gamma"`:

  Spike-and-slab mixture with a Gamma slab. The marginal CDF is
  \\F(\theta) = \hat\pi_i \cdot \mathbf{1}\[\theta \geq 0\] + (1 -
  \hat\pi_i) \cdot F\_{\text{Gamma}}(\theta)\\, so the quantile at
  \\\tau\\ is 0 when \\\tau \leq \hat\pi_i\\ and
  \\F\_{\text{Gamma}}^{-1}((\tau - \hat\pi_i)/(1-\hat\pi_i))\\
  otherwise.

- `"smooth_lognormal"`:

  Log-normal posterior from a smooth wavelet-based prior. Quantiles are
  computed as \\\exp(\mu\_{\log} + \sqrt{\sigma^2\_{\log}} \cdot
  \Phi^{-1}(\tau))\\, where \\\mu\_{\log}\\ and \\\sigma^2\_{\log}\\ are
  the posterior mean and variance in log space. When
  \\\sigma^2\_{\log}\\ is unavailable, a delta-method approximation
  \\\sigma^2\_{\log} \approx \text{Var}/\text{Mean}^2\\ is used.

If the family tag is absent (fit produced before this change), the
function falls back to `"point_gamma"` with a warning.

## Usage

``` r
get_posterior_quantile(
  fit,
  probs = c(0.025, 0.975),
  mode = "L",
  normalized = TRUE,
  verbose = FALSE
)
```

## Arguments

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  or
  [`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).

- probs:

  Numeric vector of probabilities in `[0, 1]` at which to evaluate the
  quantile function. Default `c(0.025, 0.975)`.

- mode:

  Character string specifying which mode to extract: `'L'` for
  observation loadings, `'F'` for time factors, `'W'` for channel
  weights. Default `'L'`.

- normalized:

  Logical. If `TRUE` (default), reads from `shape_post_*_normed`,
  `rate_post_*_normed`, and `PIP*_normed` fields (unit-norm columns,
  descending \\\lambda\\ order). If `FALSE`, reads from the
  corresponding raw fields.

- verbose:

  Logical. If `TRUE`, prints total function runtime to the console.
  Default `FALSE`.

## Value

A named list with one matrix per entry of `probs`. Each matrix has the
same dimensions as the corresponding factor matrix (`n x K` for L,
`p x K` for F, `w x K` for W). Names are formatted as
`paste0("q", round(probs * 100, 1))`, e.g. `"q2.5"` and `"q97.5"` for
the default `probs`. For modes with a `"smooth_lognormal"` family the
quantiles are strictly positive; for `"point_gamma"` modes they are
non-negative (zero where the spike mass absorbs the probability).

## See also

[`get_pip`](https://nd-hal.github.io/CoSparTA/reference/get_pip.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X, K = 3, Xcov = Xcov)

# 95% equal-tailed posterior intervals for observation loadings
q_L <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = 'L')
q_L$q2.5    # lower bound matrix (n x K)
q_L$q97.5   # upper bound matrix (n x K)

# Posterior median for channel weights
q_W <- get_posterior_quantile(fit, probs = 0.5, mode = 'W')
q_W$q50
} # }
```
