# Default control parameters for smooth factor estimation

Returns a list of default control parameters for the wavelet-based
smooth factor estimation step in
[`CoSparTA`](https://nd-hal.github.io/CoSparTA/reference/CoSparTA.md).
Pass the output as the `smooth_control` argument to override individual
settings.

## Usage

``` r
ebpmf_identity_smooth_control_default()
```

## Value

A named list with fields: `wave_trans`, `ndwt_method`, `filter.number`,
`family`, `ebnm_params`, `maxiter`, `maxiter_vga`, `make_power_of_2`,
`vga_tol`, `tol`, `warmstart`, `convergence_criteria`,
`m_init_method_for_init`.

## Examples

``` r
ctrl <- ebpmf_identity_smooth_control_default()
ctrl$filter.number <- 4  # override wavelet filter
```
