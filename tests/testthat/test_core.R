# Tests A-H: core decomposition, normalization, projection, reconstruction.
# Shared objects (X_obs, X_cov, n, p, w, K, lambda_tensor, ...) come from
# helper-sim.R which testthat sources automatically.

test_that("A: unsupervised fit returns valid output", {
  fit_A <- CoSparTA(X_obs, K = 2, Xcov = NULL,
                   init = 'random_gamma', maxiter = 30,
                   convergence_criteria = 'factor_change', verbose = FALSE)
  norm_A <- normalize_factors(fit_A)
  expect_equal(dim(norm_A$El), c(n, K))
  expect_equal(dim(norm_A$Ef), c(p, K))
  expect_equal(dim(norm_A$Ew), c(w, K))
  expect_true(all(norm_A$lambda > 0))
  expect_true(all(diff(norm_A$lambda) <= 0))  # descending
})

test_that("B: supervised fit recovers gamma and returns PIP/variance", {
  fit_B <- CoSparTA(X_obs, K = 2, Xcov = X_cov,
                   init = 'random_gamma', maxiter = 30,
                   convergence_criteria = 'factor_change', verbose = FALSE)
  expect_false(is.null(fit_B$res$ql$Varl))
  expect_false(is.null(fit_B$res$ql$PIPl))
  expect_true(all(fit_B$res$ql$PIPl >= 0 & fit_B$res$ql$PIPl <= 1, na.rm = TRUE))
  expect_false(is.null(fit_B$res$gl[[1]]$gamma))
  expect_false(is.null(fit_B$res$gl[[2]]$gamma))
  expect_equal(length(fit_B$res$gl[[1]]$gamma), ncol(X_cov))
})

test_that("C: missing data fit and prediction work", {
  mask  <- generate_missing_mask(X_obs, missing_rate = 0.1, seed = 7)
  fit_C <- CoSparTA_missing(X = mask$X_obs, K = 2, obs_mask = mask$obs_mask,
                            Xcov = X_cov, init = 'random_gamma', maxiter = 30,
                            verbose = FALSE)
  eval_C <- evaluate_missing_prediction(fit_C, mask)
  expect_true(is.finite(eval_C$rmse))
  expect_true(eval_C$rmse > 0)
  expect_true(is.finite(eval_C$mae))
})

test_that("D: rank-specific covariates work correctly", {
  fit_D <- CoSparTA(X_obs, K = 2, Xcov = list(X_cov, NULL),
                   init = 'random_gamma', maxiter = 5,
                   convergence_criteria = 'factor_change', verbose = FALSE)
  expect_equal(fit_D$res$gl[[1]]$type, "covariate_dependent")
  expect_false(identical(fit_D$res$gl[[2]]$type, "covariate_dependent"))
  expect_false(is.null(fit_D$res$gl[[1]]$gamma))
})

test_that("E: posterior quantiles are ordered correctly for L and W modes", {
  fit_B <- CoSparTA(X_obs, K = 2, Xcov = X_cov,
                   init = 'random_gamma', maxiter = 30,
                   convergence_criteria = 'factor_change', verbose = FALSE)
  q_L <- get_posterior_quantile(fit_B, probs = c(0.025, 0.5, 0.975), mode = 'L')
  expect_true(all(q_L$q2.5 <= q_L$q50, na.rm = TRUE))
  expect_true(all(q_L$q50 <= q_L$q97.5, na.rm = TRUE))
  expect_equal(dim(q_L$q2.5), c(n, K))
  q_W <- get_posterior_quantile(fit_B, probs = c(0.025, 0.975), mode = 'W')
  expect_true(all(q_W$q2.5 <= q_W$q97.5, na.rm = TRUE))
  expect_equal(dim(q_W$q2.5), c(w, K))
  expect_no_error(get_posterior_quantile(fit_B, mode = 'F'))
})

test_that("F: normalize_factors returns unit-norm columns sorted by lambda", {
  fit_B <- CoSparTA(X_obs, K = 2, Xcov = X_cov,
                   init = 'random_gamma', maxiter = 30,
                   convergence_criteria = 'factor_change', verbose = FALSE)
  nf_B <- normalize_factors(fit_B)
  expect_true(all(abs(sqrt(colSums(nf_B$El^2)) - 1) < 1e-10))
  expect_true(all(abs(sqrt(colSums(nf_B$Ef^2)) - 1) < 1e-10))
  expect_true(all(abs(sqrt(colSums(nf_B$Ew^2)) - 1) < 1e-10))
  expect_true(all(diff(nf_B$lambda) <= 0))
  expect_equal(length(nf_B$lambda_order), K)
  recon_raw <- reconstruct_tensor(fit_B)
  recon_nf  <- array(0, dim = dim(recon_raw))
  for (k in 1:K) {
    recon_nf <- recon_nf + nf_B$lambda[k] *
      (nf_B$El[,k] %o% nf_B$Ef[,k] %o% nf_B$Ew[,k])
  }
  expect_lt(max(abs(recon_nf - recon_raw)), 1e-8)
})

test_that("G: project_tensor returns correct dimensions", {
  fit_B <- CoSparTA(X_obs, K = 2, Xcov = X_cov,
                   init = 'random_gamma', maxiter = 30,
                   convergence_criteria = 'factor_change', verbose = FALSE)
  proj <- project_tensor(X_obs, fit_B, normalize = TRUE)
  expect_equal(dim(proj), c(n, K))
  proj_single <- project_tensor(X_obs[1,,], fit_B, normalize = TRUE)
  expect_equal(length(proj_single), K)
  expect_equal(proj_single, proj[1,])
})

test_that("H: reconstruct_tensor returns non-negative array of correct size", {
  fit_B <- CoSparTA(X_obs, K = 2, Xcov = X_cov,
                   init = 'random_gamma', maxiter = 30,
                   convergence_criteria = 'factor_change', verbose = FALSE)
  recon <- reconstruct_tensor(fit_B)
  expect_equal(dim(recon), c(n, p, w))
  expect_true(all(recon >= 0))
})

test_that("I2: _normed fields exist, have unit-norm columns, and raw fields are unchanged", {
  fit_B <- CoSparTA(X_obs, K = 2, Xcov = X_cov,
                   init = 'random_gamma', maxiter = 30,
                   convergence_criteria = 'factor_change', verbose = FALSE)

  # raw El still present and unmodified (columns NOT unit-norm in general)
  expect_false(is.null(fit_B$res$ql$El))
  expect_equal(dim(fit_B$res$ql$El), c(n, K))

  # _normed fields present
  expect_false(is.null(fit_B$res$ql$El_normed))
  expect_false(is.null(fit_B$res$qf$Ef_normed))
  expect_false(is.null(fit_B$res$qw$Ew_normed))
  expect_false(is.null(fit_B$res$lambda_normed))
  expect_false(is.null(fit_B$res$gl_normed))

  # El_normed columns have unit Frobenius norm
  col_norms_l <- sqrt(colSums(fit_B$res$ql$El_normed^2, na.rm = TRUE))
  expect_true(all(abs(col_norms_l - 1) < 1e-10))

  # lambda_normed is descending
  expect_true(all(diff(fit_B$res$lambda_normed) <= 0))

  # raw El unchanged: lambda_normed = colnorm(El)*colnorm(Ef)*colnorm(Ew), perm applied
  nl_raw <- sqrt(colSums(fit_B$res$ql$El^2, na.rm = TRUE))
  nf_raw <- sqrt(colSums(fit_B$res$qf$Ef^2, na.rm = TRUE))
  nw_raw <- sqrt(colSums(fit_B$res$qw$Ew^2, na.rm = TRUE))
  lambda_unsorted <- nl_raw * nf_raw * nw_raw
  perm <- order(lambda_unsorted, decreasing = TRUE)
  expect_equal(fit_B$res$lambda_normed, lambda_unsorted[perm], tolerance = 1e-10)

  # shape/rate/PIP normed fields present when raw exist
  expect_false(is.null(fit_B$res$ql$PIPl_normed))
  expect_false(is.null(fit_B$res$ql$shape_post_l_normed))
  expect_false(is.null(fit_B$res$ql$rate_post_l_normed))
})

test_that("I2b: get_loadings returns correct objects", {
  fit_B <- CoSparTA(X_obs, K = 2, Xcov = X_cov,
                   init = 'random_gamma', maxiter = 30,
                   convergence_criteria = 'factor_change', verbose = FALSE)

  # U1 normalized: n x K matrix, unit-norm columns
  U1 <- get_loadings(fit_B, "U1")
  expect_true(is.matrix(U1))
  expect_equal(dim(U1), c(n, K))
  expect_true(all(abs(sqrt(colSums(U1^2, na.rm = TRUE)) - 1) < 1e-10))

  # U1 raw: matches fit$res$ql$El
  U1_raw <- get_loadings(fit_B, "U1", normalized = FALSE)
  expect_identical(U1_raw, fit_B$res$ql$El)

  # weight: numeric vector of length K, descending
  w <- get_loadings(fit_B, "weight")
  expect_true(is.numeric(w))
  expect_equal(length(w), K)
  expect_true(all(diff(w) <= 0))

  # gamma: list of length K
  g <- get_loadings(fit_B, "gamma")
  expect_true(is.list(g))
  expect_equal(length(g), K)

  # invalid what: informative error
  expect_error(get_loadings(fit_B, "invalid"),
               regexp = "should be one of")
})

test_that("I3: get_posterior_quantile uses _normed by default; normalized=FALSE uses raw", {
  fit_B <- CoSparTA(X_obs, K = 2, Xcov = X_cov,
                   init = 'random_gamma', maxiter = 30,
                   convergence_criteria = 'factor_change', verbose = FALSE)

  q_normed <- get_posterior_quantile(fit_B, probs = c(0.025, 0.975),
                                      mode = 'L', normalized = TRUE)
  q_raw    <- get_posterior_quantile(fit_B, probs = c(0.025, 0.975),
                                      mode = 'L', normalized = FALSE)

  # Both return valid quantile matrices of correct dimension
  expect_equal(dim(q_normed$q2.5), c(n, K))
  expect_equal(dim(q_raw$q2.5),    c(n, K))

  # Normalized and raw use different column orderings / scales — outputs differ
  # (unless by coincidence the factor order is identical, but values scale differently)
  nl <- sqrt(colSums(fit_B$res$ql$El^2, na.rm = TRUE))
  # If max norm > 1, normalized CIs should be smaller (divided by nl)
  if (any(nl > 1.01)) {
    expect_false(isTRUE(all.equal(q_normed$q97.5, q_raw$q97.5)))
  }
})
