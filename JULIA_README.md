# CoSparTA.jl

This repository contains a native Julia port of the CoSparTA R package for
covariate-aware, sparsity-adaptive Poisson tensor decomposition. The port keeps
the R source intact and adds the standard Julia package files (`Project.toml`,
`src/CoSparTA.jl`, and `test/runtests.jl`).

## Install and load

From the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using CoSparTA
```

Or install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/cccfran/CoSparTA.git")
using CoSparTA
```

## Quick start

```julia
using CoSparTA

sim = simulate_tensor(n=100, p=20, w=10, K=3)
model = fit(sim.X, 3; Xcov=sim.Xcov, maxiter=20,
            convergence_criteria=:factor_change)

U1 = get_loadings(model, "U1")
U2 = get_loadings(model, "U2")
U3 = get_loadings(model, "U3")
weights = get_loadings(model, "weight")
gamma = get_loadings(model, "gamma")

intervals = get_posterior_quantile(model; probs=(0.025, 0.975), mode=:W)
lower = intervals["q2.5"]
upper = intervals["q97.5"]

Xhat = reconstruct_tensor(model)
```

`fit` is the Julia equivalent of the R function `CoSparTA()`. Julia reserves
the package name `CoSparTA` for the module itself, so a same-named function
cannot be defined inside it. `fit_missing` corresponds to
`CoSparTA_missing()`; a `CoSparTA_missing` compatibility alias is also
exported.

## API mapping

| R API | Julia API |
|---|---|
| `CoSparTA(X, K, ...)` | `fit(X, K; ...)` |
| `CoSparTA_missing(X, K, ...)` | `fit_missing(X, K; ...)` |
| `$` list access | typed `CoSparTAFit` fields and `get_loadings` |
| `get_posterior_quantile(...)$q2.5` | `get_posterior_quantile(...)["q2.5"]` |
| `plot_time_factors(...)` | same name; returns a RecipesBase recipe |
| `plot_channel_factors(...)` | same name; returns a RecipesBase recipe |

The remaining exported R helpers retain their names: `build_tensor`,
`normalize_factors`, `project_tensor`, `reconstruct_tensor`, `init_cpapr`,
`select_covariates`, `match_factors`, `simulate_tensor`, `generate_missing_mask`,
`evaluate_missing_prediction`, `get_pip`, `get_significant_patterns`,
`get_gamma_ci`, and the empirical-Bayes wrappers.

## Native backend differences

The model and posterior parameterization are preserved, but two R-specific
dependencies have native Julia replacements:

- `ebpm::ebpm_point_gamma` is implemented directly by marginal-likelihood
  optimization with `Optim.jl`.
- `smashrgen::ebps` is replaced by a penalized Poisson log-rate smoother with
  a second-difference prior and Laplace posterior uncertainty. This has the
  same positive smooth-lognormal role but is not numerically identical to the
  wavelet backend.
- `init_cpapr` is a native Poisson-CP EM warm start and no longer requires
  Python or `pyCP_APR`.
- The former Rcpp sparse aggregation and responsibility kernels are ordinary
  type-stable Julia loops.

Consequently, seeded R and Julia fits are expected to agree in structure and
statistical interpretation, not bit-for-bit numerical values.

## Missing data

```julia
heldout = generate_missing_mask(sim.X; missing_rate=0.1, seed=7)
missing_model = fit_missing(heldout.X_obs, 3;
                            Xcov=sim.Xcov,
                            obs_mask=heldout.obs_mask,
                            convergence_criteria=:factor_change)
metrics = evaluate_missing_prediction(missing_model, heldout)
```

## Plotting

The plotting functions return backend-neutral RecipesBase objects. Any recipe
consumer can render them; for example:

```julia
using Plots
plot(plot_time_factors(model))
plot(plot_channel_factors(model; channel_names=string.(1:10)))
```

## Tests

```julia
using Pkg
Pkg.test()
```
