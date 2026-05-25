# Tests I-K: PIP, credible intervals, significant patterns.
# Shared objects come from helper-sim.R.

fit <- CoSparTA(X_obs, K = 2, Xcov = X_cov,
               init = 'random_gamma', maxiter = 5,
               convergence_criteria = 'factor_change', verbose = FALSE)

test_that("I: get_pip returns correct structure and respects threshold", {
  pip_raw <- get_pip(fit, mode = 'L')
  expect_equal(dim(pip_raw), c(n, K))
  expect_true(all(pip_raw >= 0 & pip_raw <= 1, na.rm = TRUE))
  pip_thresh <- get_pip(fit, mode = 'L', threshold = 0.5)
  expect_true(is.logical(pip_thresh))
  expect_equal(pip_thresh, pip_raw > 0.5)
  pip_W <- get_pip(fit, mode = 'W')
  expect_equal(dim(pip_W), c(w, K))
  pip_F <- get_pip(fit, mode = 'F')
  expect_true(is.null(pip_F))
})


test_that("K: get_significant_patterns returns one entry per factor", {
  patterns <- get_significant_patterns(fit, alpha = 0.05, mode = 'both')
  expect_equal(length(patterns), K)
  patterns_W <- get_significant_patterns(fit, alpha = 0.05, mode = 'W')
  expect_true(is.null(patterns_W[[1]]$active_times))
})
