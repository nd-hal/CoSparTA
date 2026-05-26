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

test_that("L: family tags are stored on the fit object for all three modes", {
  expect_equal(fit$res$ql$family_l, "point_gamma")
  expect_equal(fit$res$qf$family_f, "smooth_lognormal")
  expect_equal(fit$res$qw$family_w, "point_gamma")
})

test_that("M: get_posterior_quantile dispatches correctly for L mode (point_gamma)", {
  q_L <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = 'L')
  expect_equal(dim(q_L$q2.5),  c(n, K))
  expect_equal(dim(q_L$q97.5), c(n, K))
  # point_gamma: quantiles are non-negative
  expect_true(all(q_L$q2.5  >= 0, na.rm = TRUE))
  expect_true(all(q_L$q97.5 >= q_L$q2.5, na.rm = TRUE))
})

test_that("N: get_posterior_quantile dispatches correctly for W mode (point_gamma)", {
  q_W <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = 'W')
  expect_equal(dim(q_W$q2.5),  c(w, K))
  expect_equal(dim(q_W$q97.5), c(w, K))
  expect_true(all(q_W$q2.5 >= 0, na.rm = TRUE))
})

test_that("O: get_posterior_quantile dispatches correctly for F mode (smooth_lognormal)", {
  q_F <- get_posterior_quantile(fit, probs = c(0.025, 0.975), mode = 'F')
  expect_equal(dim(q_F$q2.5),  c(p, K))
  expect_equal(dim(q_F$q97.5), c(p, K))
  # smooth_lognormal: quantiles are strictly positive
  expect_true(all(q_F$q2.5  > 0, na.rm = TRUE))
  expect_true(all(q_F$q97.5 > 0, na.rm = TRUE))
  expect_true(all(q_F$q97.5 >= q_F$q2.5, na.rm = TRUE))
  # values should not be all identical (non-trivial)
  expect_true(var(as.vector(q_F$q97.5), na.rm = TRUE) > 0)
})
