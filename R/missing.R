#' Generate Missing Mask for Tensor Simulation
#'
#' @description
#' Randomly masks a proportion of non-zero entries in a tensor to simulate
#' missing data. Returns the observed tensor, observation mask, and metadata
#' needed for evaluation of imputation quality.
#'
#' @param X Non-negative integer array of dimensions \code{n x p x w}.
#' @param missing_rate Numeric proportion of non-zero entries to mask.
#'   Default \code{0.1} (10\%).
#' @param seed Optional integer random seed for reproducibility.
#'   Default \code{NULL}.
#'
#' @return A named list with:
#' \describe{
#'   \item{X_obs}{Observed tensor with masked entries set to \code{NA}.}
#'   \item{obs_mask}{Logical array of same dimensions as \code{X}, where
#'     \code{TRUE} indicates observed entries.}
#'   \item{obs_indices}{data.table of observed entry indices (V1, V2, V3).}
#'   \item{missing_nonzero_indices}{data.table of masked non-zero entry
#'     indices (V1, V2, V3).}
#'   \item{true_values}{Numeric vector of true values at masked positions.}
#'   \item{n_missing}{Number of masked entries.}
#' }
#'
#' @examples
#' \dontrun{
#' X <- array(rpois(20 * 12 * 4, lambda = 1.5), dim = c(20, 12, 4))
#' mask_info <- generate_missing_mask(X, missing_rate = 0.1, seed = 42)
#' }
#'
#' @export
generate_missing_mask <- function(X, missing_rate = 0.1, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n <- dim(X)[1]
  p <- dim(X)[2]
  w <- dim(X)[3]

  nonzero_idx <- which(X > 0, arr.ind = TRUE)
  n_nonzero <- nrow(nonzero_idx)

  n_missing <- round(n_nonzero * missing_rate)
  missing_rows <- sample(1:n_nonzero, n_missing, replace = FALSE)

  missing_nonzero_indices <- data.table(
    V1 = nonzero_idx[missing_rows, 1],
    V2 = nonzero_idx[missing_rows, 2],
    V3 = nonzero_idx[missing_rows, 3]
  )

  true_values <- X[as.matrix(missing_nonzero_indices)]

  X_obs <- X
  X_obs[as.matrix(missing_nonzero_indices)] <- NA

  obs_mask <- array(TRUE, dim = c(n, p, w))
  obs_mask[as.matrix(missing_nonzero_indices)] <- FALSE

  all_idx <- data.table(expand.grid(V1 = 1:n, V2 = 1:p, V3 = 1:w))
  missing_nonzero_indices[, missing := TRUE]
  all_idx <- merge(all_idx, missing_nonzero_indices,
                   by = c("V1", "V2", "V3"), all.x = TRUE)
  obs_indices <- all_idx[is.na(missing), .(V1, V2, V3)]

  cat(sprintf("Total entries: %d\n", n * p * w))
  cat(sprintf("Original non-zeros: %d\n", n_nonzero))
  cat(sprintf("Masked non-zeros: %d (%.1f%%)\n", n_missing, 100 * missing_rate))
  cat(sprintf("Observed entries: %d\n", nrow(obs_indices)))

  return(list(
    X_obs = X_obs,
    obs_mask = obs_mask,
    obs_indices = obs_indices,
    missing_nonzero_indices = missing_nonzero_indices[, .(V1, V2, V3)],
    true_values = true_values,
    n_missing = n_missing
  ))
}


#' Evaluate Prediction Quality on Held-out Missing Entries
#'
#' @description
#' Computes prediction metrics for held-out entries after fitting
#' \code{\link{CxtEBTD_missing}}. Reconstructs the predicted Poisson rate
#' at masked positions using the fitted factor matrices and compares against
#' the true values.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD_missing}}.
#' @param missing_info Output from \code{\link{generate_missing_mask}}.
#'
#' @return A named list with:
#' \describe{
#'   \item{rmse}{Root mean squared error.}
#'   \item{mae}{Mean absolute error.}
#'   \item{correlation}{Pearson correlation between predicted and true values.}
#'   \item{deviance}{Poisson deviance.}
#'   \item{predicted}{Numeric vector of predicted values at masked positions.}
#'   \item{true_values}{Numeric vector of true values at masked positions.}
#' }
#'
#' @examples
#' \dontrun{
#' X <- array(rpois(20 * 12 * 4, lambda = 1.5), dim = c(20, 12, 4))
#' mask_info <- generate_missing_mask(X, missing_rate = 0.1, seed = 42)
#' fit <- CxtEBTD_missing(mask_info$X_obs, K = 3,
#'                         obs_mask = mask_info$obs_mask)
#' metrics <- evaluate_missing_prediction(fit, mask_info)
#' }
#'
#' @export
evaluate_missing_prediction <- function(fit, missing_info) {

  El <- fit$res$ql$El
  Ef <- fit$res$qf$Ef
  Ew <- fit$res$qw$Ew
  K  <- ncol(El)

  missing_idx  <- as.matrix(missing_info$missing_nonzero_indices)
  n_missing    <- nrow(missing_idx)
  lambda_pred  <- numeric(n_missing)

  for (k in 1:K) {
    lambda_pred <- lambda_pred +
      El[missing_idx[, 1], k] *
      Ef[missing_idx[, 2], k] *
      Ew[missing_idx[, 3], k]
  }

  true_values <- missing_info$true_values
  rmse        <- sqrt(mean((lambda_pred - true_values)^2))
  mae         <- mean(abs(lambda_pred - true_values))
  cor_val     <- cor(lambda_pred, true_values)
  deviance    <- 2 * sum(
    true_values * log(pmax(true_values, 1e-10) / pmax(lambda_pred, 1e-10)) -
    (true_values - lambda_pred)
  )

  cat(sprintf("Missing data prediction:\n"))
  cat(sprintf("  RMSE: %.4f\n", rmse))
  cat(sprintf("  MAE: %.4f\n", mae))
  cat(sprintf("  Correlation: %.4f\n", cor_val))
  cat(sprintf("  Poisson Deviance: %.4f\n", deviance))

  return(list(
    rmse        = rmse,
    mae         = mae,
    correlation = cor_val,
    deviance    = deviance,
    predicted   = lambda_pred,
    true_values = true_values
  ))
}


#' Covariate-aware Empirical Bayes Tensor Decomposition with Missing Data
#'
#' @description
#' Extends \code{\link{CxtEBTD}} to handle tensors with missing entries.
#' Missing entries are excluded from the likelihood and the ELBO, and
#' per-observation scales are computed using only observed entries for each
#' mode. The generative model and priors are identical to
#' \code{\link{CxtEBTD}}.
#'
#' @param X A 3-dimensional non-negative integer array of dimensions
#'   \code{n x p x w}. May contain \code{NA} for missing entries.
#' @param K Integer. Number of components (CP rank).
#' @param Xcov Covariate input for the observation mode. Can be: (1) \code{NULL}
#'   for fully unsupervised decomposition; (2) a numeric matrix of dimension
#'   \code{n x q}, applied identically to all K components; or (3) a list of
#'   length \code{K}, where each element is either a numeric covariate matrix
#'   (dimensions \code{n x q_k}, potentially different numbers of covariates per
#'   component) or \code{NULL} for unsupervised components. Default \code{NULL}.
#' @param obs_mask Optional logical array of same dimensions as \code{X},
#'   where \code{TRUE} indicates observed entries. If \code{NULL}, inferred
#'   from \code{is.na(X)}. Default \code{NULL}.
#' @param lib_size Numeric vector of length \code{n}. Default \code{NULL}.
#' @param init Initialization method. Default \code{'random_gamma'}.
#' @param maxiter Maximum EM iterations. Default \code{100}.
#' @param maxiter_init Maximum initialization iterations. Default \code{100}.
#' @param tol Convergence tolerance. Default \code{1e-8}.
#' @param ebpm.fn List of three EBPM functions for L, F, W modes respectively.
#'   Default uses \code{ebpm_point_gamma_multiplier_covariates} for L,
#'   \code{ebps_with_uq} for F, and \code{ebpm_point_gamma_with_uq} for W.
#' @param fix_L Logical. Fix L at initialization. Default \code{FALSE}.
#' @param fix_F Logical. Fix F at initialization. Default \code{FALSE}.
#' @param fix_W Logical. Fix W at initialization. Default \code{FALSE}.
#' @param smooth_F Logical. Apply smoothing to F. Default \code{TRUE}.
#' @param printevery Print progress every this many iterations. Default
#'   \code{10}.
#' @param verbose Logical. Print progress. Default \code{TRUE}.
#' @param adj_LF_scale Logical. Rescale L and F each iteration. Default
#'   \code{TRUE}.
#' @param convergence_criteria Convergence criterion: \code{'ELBO'} (default)
#'   or \code{'mKLabs'}.
#' @param U1_true Optional true L matrix for simulation evaluation.
#' @param U2_true Optional true F matrix for simulation evaluation.
#' @param U3_true Optional true W matrix for simulation evaluation.
#'
#' @return Same structure as \code{\link{CxtEBTD}} with one additional field:
#' \describe{
#'   \item{obs_structure}{Internal observation structure used during fitting,
#'     required by \code{\link{evaluate_missing_prediction}}.}
#' }
#'
#' @seealso \code{\link{CxtEBTD}}, \code{\link{generate_missing_mask}},
#'   \code{\link{evaluate_missing_prediction}}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- array(rpois(20 * 12 * 4, lambda = 1.5), dim = c(20, 12, 4))
#' mask_info <- generate_missing_mask(X, missing_rate = 0.1, seed = 42)
#'
#' # With covariates
#' Xcov <- matrix(rnorm(20 * 2), nrow = 20)
#' fit <- CxtEBTD_missing(mask_info$X_obs, K = 3, Xcov = Xcov,
#'                         obs_mask = mask_info$obs_mask)
#'
#' # Without covariates
#' fit0 <- CxtEBTD_missing(mask_info$X_obs, K = 3,
#'                          obs_mask = mask_info$obs_mask)
#' }
#'
#' @export
CxtEBTD_missing <- function(X, K, Xcov = NULL,
                             obs_mask = NULL,
                             lib_size = NULL,
                             init = 'random_gamma',
                             maxiter = 100,
                             maxiter_init = 100,
                             tol = 1e-8,
                             ebpm.fn = c(ebpm_point_gamma_multiplier_covariates,
                                         ebps_with_uq,
                                         ebpm_point_gamma_with_uq),
                             fix_L = FALSE, fix_F = FALSE, fix_W = FALSE,
                             smooth_F = TRUE,
                             printevery = 10,
                             verbose = TRUE,
                             adj_LF_scale = TRUE,
                             convergence_criteria = 'ELBO',
                             U1_true = NULL, U2_true = NULL, U3_true = NULL) {

  start_time <- Sys.time()

  n_original <- dim(X)[1]
  p_original <- dim(X)[2]
  w_original <- dim(X)[3]

  # Build obs_mask from NAs if not provided
  if (is.null(obs_mask)) {
    obs_mask <- !is.na(X)
  }

  # Set missing to 0 for computation
  X[is.na(X)]   <- 0
  X[!obs_mask]  <- 0

  # Drop all-zero slices
  users_zero    <- apply(X, 1, sum) == 0
  channels_zero <- apply(X, 3, sum) == 0
  times_zero    <- apply(X, 2, sum) == 0

  X        <- X[!users_zero, , , drop = FALSE]
  X        <- X[, , !channels_zero, drop = FALSE]
  X        <- X[, !times_zero, , drop = FALSE]
  obs_mask <- obs_mask[!users_zero, , , drop = FALSE]
  obs_mask <- obs_mask[, , !channels_zero, drop = FALSE]
  obs_mask <- obs_mask[, !times_zero, , drop = FALSE]

  n        <- dim(X)[1]
  p        <- dim(X)[2]
  w        <- dim(X)[3]
  n_points <- n * p * w

  if (!is.null(Xcov)) {
    if (is.matrix(Xcov)) {
      Xcov <- rep(list(Xcov), K)
    }
    if (!is.list(Xcov) || length(Xcov) != K) {
      stop("Xcov must be NULL, a matrix, or a list of length K")
    }
    Xcov <- lapply(Xcov, function(xc) {
      if (!is.null(xc)) xc[!users_zero, , drop = FALSE] else NULL
    })
  }

  if (is.null(lib_size)) {
    lib_size <- rep(1, n)
  }

  # Sparse representation of observed non-zeros
  x        <- rbind_sparse_matrix(X, reindex = TRUE)
  non0_idx <- cbind(x$V1, x$V2, x$V3)

  # Observation structure for scale computation
  obs_idx     <- which(obs_mask, arr.ind = TRUE)
  obs_indices <- data.table(V1 = obs_idx[, 1],
                             V2 = obs_idx[, 2],
                             V3 = obs_idx[, 3])
  obs_structure <- .precompute_obs_structure(obs_indices, n, p, w)

  if (verbose) {
    cat(sprintf("Observed entries: %d / %d (%.1f%%)\n",
                nrow(obs_indices), n_points,
                100 * nrow(obs_indices) / n_points))
  }

  # EBPM function assignment
  if (length(ebpm.fn) == 1) {
    ebpm.fn.l <- ebpm.fn
    ebpm.fn.f <- ebpm.fn
    ebpm.fn.w <- ebpm.fn
  }
  if (length(ebpm.fn) == 3) {
    ebpm.fn.l <- ebpm.fn[[1]]
    ebpm.fn.f <- ebpm.fn[[2]]
    ebpm.fn.w <- ebpm.fn[[3]]
  }

  # Normalize ebpm.fn.l to a length-K list
  if (is.function(ebpm.fn.l)) {
    ebpm_fn_l_list <- rep(list(ebpm.fn.l), K)
  } else if (is.list(ebpm.fn.l) && length(ebpm.fn.l) == K) {
    ebpm_fn_l_list <- ebpm.fn.l
  } else {
    stop("ebpm.fn L-mode entry must be a single function or a list of K functions")
  }

  if (verbose) cat('Initializing loadings and factors...\n')

  res <- ebpmf_identity_init(X, K, init, maxiter_init, lib_size)

  alpha     <- res$ql$Elogl[x$V1, , drop = FALSE] +
               res$qf$Elogf[x$V2, , drop = FALSE] +
               res$qw$Elogw[x$V3, , drop = FALSE]
  exp_offset <- matrixStats::rowMaxs(alpha)
  alpha      <- alpha - outer(exp_offset, rep(1, K), FUN = '*')
  alpha      <- exp(alpha)
  alpha      <- alpha / rowSums(alpha)

  obj    <- c()
  obj[1] <- -Inf
  diff_U <- list(c(), c(), c())

  if (verbose) cat('Running iterations...\n')

  for (iter in 1:maxiter) {

    for (k in 1:K) {
      Ez <- .calc_EZ_3d_missing(x, alpha[, k], n, p, w)
      xcov_k <- if (!is.null(Xcov)) Xcov[[k]] else NULL
      fn_l_k <- ebpm_fn_l_list[[k]]
      if (is.null(xcov_k) && identical(fn_l_k, ebpm_point_gamma_multiplier_covariates)) {
        fn_l_k <- ebpm::ebpm_point_gamma
      }
      res <- .stm_update_rank1_missing(
        Ez$rs, Ez$cs, Ez$zs, k,
        fn_l_k, ebpm.fn.f, ebpm.fn.w,
        res, fix_L, fix_F, fix_W, xcov_k,
        obs_structure, n, p, w
      )
    }

    alpha      <- res$ql$Elogl[x$V1, , drop = FALSE] +
                  res$qf$Elogf[x$V2, , drop = FALSE] +
                  res$qw$Elogw[x$V3, , drop = FALSE]
    exp_offset <- matrixStats::rowMaxs(alpha)
    alpha      <- alpha - outer(exp_offset, rep(1, K), FUN = '*')
    alpha      <- exp(alpha)
    alpha      <- alpha / rowSums(alpha)

    if (convergence_criteria == 'ELBO') {
      obj[iter + 1] <- .calc_stm_obj_missing(x, n, p, w, K, res,
                                              non0_idx, obs_structure)
      if (verbose && iter %% printevery == 0) {
        cat(sprintf('Iter %d, ELBO: %f\n', iter, obj[iter + 1]))
      }
      if (is.infinite(obj[iter + 1]) && obj[iter + 1] < 0) {
        res <- res_prev
        break
      }
      if ((obj[iter + 1] - obj[iter]) / n_points < tol) break
    }

    if (convergence_criteria == 'mKLabs') {
      obj[iter + 1] <- mKL(x$v,
        (tcrossprod(res$ql$El, res$qf$Ef) * res$lib_size)[non0_idx])
      if (verbose && iter %% printevery == 0) {
        cat(sprintf('Iter %d, mKL: %f\n', iter, obj[iter + 1]))
      }
      if (abs(obj[iter + 1] - obj[iter]) <= tol) break
    }

    if (!is.null(U1_true)) {
      ret_EL <- matrix(0, nrow = n_original, ncol = K)
      ret_EF <- matrix(0, nrow = p_original, ncol = K)
      ret_EW <- matrix(0, nrow = w_original, ncol = K)
      ret_EL[!users_zero, ]    <- res$ql$El
      ret_EF[!times_zero, ]    <- res$qf$Ef
      ret_EW[!channels_zero, ] <- res$qw$Ew

      U1          <- U1_true[, 1] / norm(matrix(U1_true[, 1]), type = "F")
      U1_hat_norm <- which_rank(U1, ret_EL)
      U2          <- U2_true[, 1] / norm(matrix(U2_true[, 1]), type = "F")
      U2_hat_norm <- which_rank(U2, ret_EF)
      U3          <- U3_true[, 1] / norm(matrix(U3_true[, 1]), type = "F")
      U3_hat_norm <- which_rank(U3, ret_EW)

      diff_U[[1]] <- c(diff_U[[1]], U1_hat_norm$diff)
      diff_U[[2]] <- c(diff_U[[2]], U2_hat_norm$diff)
      diff_U[[3]] <- c(diff_U[[3]], U3_hat_norm$diff)

      if (verbose && iter %% printevery == 0) {
        cat(sprintf('  U1: %f, U2: %f, U3: %f\n',
                    U1_hat_norm$diff, U2_hat_norm$diff, U3_hat_norm$diff))
      }
    }

    res_prev <- res
  }

  if (iter == maxiter) message('Reached maximum iterations')

  elbo <- .calc_stm_obj_missing(x, n, p, w, K, res, non0_idx, obs_structure)

  ret_EL <- matrix(0, nrow = n_original, ncol = K)
  ret_EF <- matrix(0, nrow = p_original, ncol = K)
  ret_EW <- matrix(0, nrow = w_original, ncol = K)
  ret_EL[!users_zero, ]    <- res$ql$El
  ret_EF[!times_zero, ]    <- res$qf$Ef
  ret_EW[!channels_zero, ] <- res$qw$Ew

  res$ql$El <- ret_EL
  res$qf$Ef <- ret_EF
  res$qw$Ew <- ret_EW

  # Expand variance/PIP matrices back to original dimensions
  if (!is.null(res$ql$Varl)) {
    tmp <- matrix(NA_real_, n_original, K); tmp[!users_zero,] <- res$ql$Varl; res$ql$Varl <- tmp
    tmp <- matrix(NA_real_, n_original, K); tmp[!users_zero,] <- res$ql$PIPl; res$ql$PIPl <- tmp
  }
  if (!is.null(res$qf$Varf)) {
    tmp <- matrix(NA_real_, p_original, K); tmp[!times_zero,] <- res$qf$Varf; res$qf$Varf <- tmp
    tmp <- matrix(NA_real_, p_original, K); tmp[!times_zero,] <- res$qf$PIPf; res$qf$PIPf <- tmp
  }
  if (!is.null(res$qw$Varw)) {
    tmp <- matrix(NA_real_, w_original, K); tmp[!channels_zero,] <- res$qw$Varw; res$qw$Varw <- tmp
    tmp <- matrix(NA_real_, w_original, K); tmp[!channels_zero,] <- res$qw$PIPw; res$qw$PIPw <- tmp
  }
  if (!is.null(res$ql$shape_post_l)) {
    tmp <- matrix(NA_real_, n_original, K); tmp[!users_zero,] <- res$ql$shape_post_l; res$ql$shape_post_l <- tmp
    tmp <- matrix(NA_real_, n_original, K); tmp[!users_zero,] <- res$ql$rate_post_l;  res$ql$rate_post_l  <- tmp
  }
  if (!is.null(res$qf$shape_post_f)) {
    tmp <- matrix(NA_real_, p_original, K); tmp[!times_zero,] <- res$qf$shape_post_f; res$qf$shape_post_f <- tmp
    tmp <- matrix(NA_real_, p_original, K); tmp[!times_zero,] <- res$qf$rate_post_f;  res$qf$rate_post_f  <- tmp
  }
  if (!is.null(res$qw$shape_post_w)) {
    tmp <- matrix(NA_real_, w_original, K); tmp[!channels_zero,] <- res$qw$shape_post_w; res$qw$shape_post_w <- tmp
    tmp <- matrix(NA_real_, w_original, K); tmp[!channels_zero,] <- res$qw$rate_post_w;  res$qw$rate_post_w  <- tmp
  }

  return(list(
    EL           = "check res",
    EF           = "check res",
    elbo         = elbo,
    obj_trace    = obj,
    res          = res,
    diff_U       = diff_U,
    obs_structure = obs_structure,
    run_time     = difftime(Sys.time(), start_time, units = 'auto')
  ))
}


# ---- Internal helpers -------------------------------------------------------

#' @keywords internal
.precompute_obs_structure <- function(obs_indices, n, p, w) {
  user_obs_jm    <- obs_indices[, .(obs_jm = list(
                      data.table(V2 = V2, V3 = V3))), by = V1]
  setkey(user_obs_jm, V1)
  time_obs_im    <- obs_indices[, .(obs_im = list(
                      data.table(V1 = V1, V3 = V3))), by = V2]
  setkey(time_obs_im, V2)
  channel_obs_ij <- obs_indices[, .(obs_ij = list(
                      data.table(V1 = V1, V2 = V2))), by = V3]
  setkey(channel_obs_ij, V3)

  return(list(
    user_obs_jm    = user_obs_jm,
    time_obs_im    = time_obs_im,
    channel_obs_ij = channel_obs_ij,
    obs_indices    = obs_indices
  ))
}

#' @keywords internal
.compute_l_scales <- function(Ef_k, Ew_k, obs_structure, n) {
  l_scales <- numeric(n)
  for (i in 1:n) {
    user_data <- obs_structure$user_obs_jm[V1 == i]
    if (nrow(user_data) == 0 || length(user_data$obs_jm) == 0) {
      l_scales[i] <- 1e-10
    } else {
      jm_pairs    <- user_data$obs_jm[[1]]
      l_scales[i] <- sum(Ef_k[jm_pairs$V2] * Ew_k[jm_pairs$V3])
    }
  }
  return(l_scales)
}

#' @keywords internal
.compute_f_scales <- function(El_k, Ew_k, obs_structure, p) {
  f_scales <- numeric(p)
  for (j in 1:p) {
    time_data <- obs_structure$time_obs_im[V2 == j]
    if (nrow(time_data) == 0 || length(time_data$obs_im) == 0) {
      f_scales[j] <- 1e-10
    } else {
      im_pairs    <- time_data$obs_im[[1]]
      f_scales[j] <- sum(El_k[im_pairs$V1] * Ew_k[im_pairs$V3])
    }
  }
  return(f_scales)
}

#' @keywords internal
.compute_w_scales <- function(El_k, Ef_k, obs_structure, w) {
  w_scales <- numeric(w)
  for (m in 1:w) {
    channel_data <- obs_structure$channel_obs_ij[V3 == m]
    if (nrow(channel_data) == 0 || length(channel_data$obs_ij) == 0) {
      w_scales[m] <- 1e-10
    } else {
      ij_pairs    <- channel_data$obs_ij[[1]]
      w_scales[m] <- sum(El_k[ij_pairs$V1] * Ef_k[ij_pairs$V2])
    }
  }
  return(w_scales)
}

#' @keywords internal
.stm_update_rank1_missing <- function(l_seq, f_seq, w_seq, k,
                                       ebpm.fn.l, ebpm.fn.f, ebpm.fn.w,
                                       res, fix_L, fix_F, fix_W, Xcov,
                                       obs_structure, n, p, w) {
  # Update L
  if (!fix_L) {
    l_scale <- .compute_l_scales(res$qf$Ef[, k], res$qw$Ew[, k],
                                  obs_structure, n) * res$lib_size
    if (!is.null(Xcov)) {
      fit <- ebpm.fn.l(l_seq, l_scale, Xcov)
    } else {
      fit <- ebpm.fn.l(l_seq, l_scale)
    }
    res$ql$El[, k]    <- fit$posterior$mean
    res$ql$Elogl[, k] <- fit$posterior$mean_log
    if (!is.null(fit$posterior$var)) {
      if (is.null(res$ql$Varl)) res$ql$Varl <- matrix(NA_real_, nrow(res$ql$El), ncol(res$ql$El))
      if (is.null(res$ql$PIPl)) res$ql$PIPl <- matrix(NA_real_, nrow(res$ql$El), ncol(res$ql$El))
      res$ql$Varl[, k] <- fit$posterior$var
      res$ql$PIPl[, k] <- fit$posterior$pip
    }
    if (!is.null(fit$posterior$shape_post)) {
      if (is.null(res$ql$shape_post_l)) res$ql$shape_post_l <- matrix(NA_real_, nrow(res$ql$El), ncol(res$ql$El))
      if (is.null(res$ql$rate_post_l))  res$ql$rate_post_l  <- matrix(NA_real_, nrow(res$ql$El), ncol(res$ql$El))
      res$ql$shape_post_l[, k] <- fit$posterior$shape_post
      res$ql$rate_post_l[, k]  <- fit$posterior$rate_post
    }
    res$Hl[k]  <- calc_H(l_seq, l_scale, fit$log_likelihood,
                          fit$posterior$mean, fit$posterior$mean_log)
    res$gl[[k]] <- fit$fitted_g
  }

  # Update F
  if (!fix_F) {
    f_scale <- .compute_f_scales(res$ql$El[, k], res$qw$Ew[, k],
                                  obs_structure, p)
    fit <- ebpm.fn.f(f_seq, f_scale)
    res$qf$Ef[, k]    <- fit$posterior$mean
    res$qf$Elogf[, k] <- fit$posterior$mean_log
    if (!is.null(fit$posterior$var)) {
      if (is.null(res$qf$Varf)) res$qf$Varf <- matrix(NA_real_, nrow(res$qf$Ef), ncol(res$qf$Ef))
      if (is.null(res$qf$PIPf)) res$qf$PIPf <- matrix(NA_real_, nrow(res$qf$Ef), ncol(res$qf$Ef))
      res$qf$Varf[, k] <- fit$posterior$var
      res$qf$PIPf[, k] <- fit$posterior$pip
    }
    if (!is.null(fit$posterior$shape_post)) {
      if (is.null(res$qf$shape_post_f)) res$qf$shape_post_f <- matrix(NA_real_, nrow(res$qf$Ef), ncol(res$qf$Ef))
      if (is.null(res$qf$rate_post_f))  res$qf$rate_post_f  <- matrix(NA_real_, nrow(res$qf$Ef), ncol(res$qf$Ef))
      res$qf$shape_post_f[, k] <- fit$posterior$shape_post
      res$qf$rate_post_f[, k]  <- fit$posterior$rate_post
    }
    res$Hf[k]  <- calc_H(f_seq, f_scale, fit$log_likelihood,
                          fit$posterior$mean, fit$posterior$mean_log)
    res$gf[[k]] <- fit$fitted_g
  }

  # Update W
  if (!fix_W) {
    w_scale <- .compute_w_scales(res$ql$El[, k], res$qf$Ef[, k],
                                  obs_structure, w)
    fit <- ebpm.fn.w(w_seq, w_scale)
    res$qw$Ew[, k]    <- fit$posterior$mean
    res$qw$Elogw[, k] <- fit$posterior$mean_log
    if (!is.null(fit$posterior$var)) {
      if (is.null(res$qw$Varw)) res$qw$Varw <- matrix(NA_real_, nrow(res$qw$Ew), ncol(res$qw$Ew))
      if (is.null(res$qw$PIPw)) res$qw$PIPw <- matrix(NA_real_, nrow(res$qw$Ew), ncol(res$qw$Ew))
      res$qw$Varw[, k] <- fit$posterior$var
      res$qw$PIPw[, k] <- fit$posterior$pip
    }
    if (!is.null(fit$posterior$shape_post)) {
      if (is.null(res$qw$shape_post_w)) res$qw$shape_post_w <- matrix(NA_real_, nrow(res$qw$Ew), ncol(res$qw$Ew))
      if (is.null(res$qw$rate_post_w))  res$qw$rate_post_w  <- matrix(NA_real_, nrow(res$qw$Ew), ncol(res$qw$Ew))
      res$qw$shape_post_w[, k] <- fit$posterior$shape_post
      res$qw$rate_post_w[, k]  <- fit$posterior$rate_post
    }
    res$Hw[k]  <- calc_H(w_seq, w_scale, fit$log_likelihood,
                          fit$posterior$mean, fit$posterior$mean_log)
    if (!is.null(fit$fitted_g)) res$gw[[k]] <- fit$fitted_g
  }

  return(res)
}

#' @keywords internal
.calc_EZ_3d_missing <- function(x, prob, n, p, w) {
  x$prob_x <- x$v * prob

  out_x1 <- x %>% dplyr::group_by(V1) %>%
    dplyr::summarize(new = sum(prob_x), .groups = 'drop') %>% dplyr::arrange(V1)
  out_x2 <- x %>% dplyr::group_by(V2) %>%
    dplyr::summarize(new = sum(prob_x), .groups = 'drop') %>% dplyr::arrange(V2)
  out_x3 <- x %>% dplyr::group_by(V3) %>%
    dplyr::summarize(new = sum(prob_x), .groups = 'drop') %>% dplyr::arrange(V3)

  rs <- numeric(n)
  cs <- numeric(p)
  zs <- numeric(w)

  rs[out_x1$V1] <- out_x1$new
  cs[out_x2$V2] <- out_x2$new
  zs[out_x3$V3] <- out_x3$new

  return(list(rs = rs, cs = cs, zs = zs))
}

#' @keywords internal
.calc_stm_obj_missing <- function(x, n, p, w, K, res, non0_idx,
                                   obs_structure) {
  qz  <- calc_qz(n, p, w, K, res$ql, res$qf, res$qw)
  val <- 0
  sum_u <- 0

  obs_idx_mat <- as.matrix(obs_structure$obs_indices[, .(V1, V2, V3)])

  for (k in 1:K) {
    ql_k <- res$ql$El[, k]
    qf_k <- res$qf$Ef[, k]
    qw_k <- res$qw$Ew[, k]
    val  <- val + qz[, , , k] * (
      log(ql_k %o% qf_k %o% qw_k) - log(qz[, , , k])
    )
    lambda_k <- ql_k %o% qf_k %o% qw_k
    sum_u    <- sum_u + sum(lambda_k[obs_idx_mat])
  }

  E1 <- sum(x$v * val[non0_idx]) - sum_u - sum(lfactorial(x$v))
  return(E1 + sum(res$Hl) + sum(res$Hf) + sum(res$Hw))
}
