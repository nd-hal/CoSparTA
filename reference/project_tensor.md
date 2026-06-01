# Project a New Tensor onto the Learned Factor Space

Given a fitted
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
object (or raw factor matrices), projects one or more new observations
onto the learned factor space to produce a loading matrix. Each
observation is a `p x w` count slice; the function recovers the
`K`-dimensional representation without refitting.

The projection for observation \\i\\ and component \\k\\ is: \$\$F\_{ik}
= \lambda_k \\ \mathbf{f}\_k^\top X\_{\text{new}\[i,,\]} \\
\mathbf{w}\_k\$\$ when `normalize = TRUE` or raw matrices with `lambda`
are supplied, or \$\$F\_{ik} = \mathbf{f}\_k^\top X\_{\text{new}\[i,,\]}
\\ \mathbf{w}\_k\$\$ using the raw posterior means when
`normalize = FALSE`.

Computation is vectorized over observations: `X_new` is reshaped to an
`n_new x (p*w)` matrix and each Kronecker factor \\\mathbf{f}\_k \otimes
\mathbf{w}\_k\\ is a `p*w`-vector, so the full projection for component
\\k\\ reduces to a single matrix-vector multiply.

## Usage

``` r
project_tensor(
  X_new,
  fit = NULL,
  normalize = TRUE,
  Ef = NULL,
  Ew = NULL,
  lambda = NULL
)
```

## Arguments

- X_new:

  Either a 3D array of dimensions `n_new x p x w` (multiple new
  observations), or a matrix of dimensions `p x w` (a single
  observation). In the latter case the input is wrapped as
  `array(X_new, dim = c(1, p, w))` and a length-`K` vector is returned
  instead of a `1 x K` matrix.

- fit:

  A fitted object returned by
  [`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)
  or
  [`CoSparTA_missing`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA_missing.md).
  Either `fit` or all of `Ef`, `Ew`, and `lambda` must be supplied.

- normalize:

  Logical. If `TRUE` (default) and `fit` is supplied, calls
  [`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md)
  so that factor columns have unit Frobenius norm and are ordered by
  \\\lambda\\ descending. If `FALSE`, uses raw posterior means with no
  scaling or reordering. Ignored when `Ef`, `Ew`, and `lambda` are
  provided directly.

- Ef:

  Numeric matrix of dimensions `p x K` (time factor matrix, unit-norm
  columns). Provide together with `Ew` and `lambda` to bypass the `fit`
  object.

- Ew:

  Numeric matrix of dimensions `w x K` (channel factor matrix, unit-norm
  columns). Provide together with `Ef` and `lambda` to bypass the `fit`
  object.

- lambda:

  Numeric vector of length `K` (component weights). Provide together
  with `Ef` and `Ew` to bypass the `fit` object.

## Value

If `X_new` is a 3D array: an `n_new x K` numeric matrix of projected
loadings. If `X_new` is a `p x w` matrix: a numeric vector of length
`K`.

## See also

[`normalize_factors`](https://nd-hal.github.io/CoSparTA/reference/normalize_factors.md),
[`reconstruct_tensor`](https://nd-hal.github.io/CoSparTA/reference/reconstruct_tensor.md),
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- CoSparTA(X_train, K = 3)

# Project a batch of new observations (50 x p x w array)
L_new <- project_tensor(X_new, fit)          # 50 x 3 matrix

# Project a single p x w slice
l_one <- project_tensor(X_new[1, , ], fit)   # length-3 vector

# Without normalization (raw posterior means, original component order)
L_raw <- project_tensor(X_new, fit, normalize = FALSE)

# Using raw factor matrices directly (already normalized)
nf <- normalize_factors(fit)
L_new2 <- project_tensor(X_new, Ef = nf$Ef, Ew = nf$Ew, lambda = nf$lambda)
} # }
```
