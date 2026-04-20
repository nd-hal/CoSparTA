#' Normalize and Sort Components of a CxtEBTD fit
#'
#' @description
#' Extracts the factor matrices from a fitted \code{\link{CxtEBTD}} object,
#' normalizes each column to unit Frobenius norm, computes a per-component
#' weight \eqn{\lambda_k = \|\mathbf{l}_k\| \|\mathbf{f}_k\| \|\mathbf{w}_k\|},
#' and returns components sorted by \eqn{\lambda} descending. This is the
#' canonical form for comparing decompositions and for downstream use with
#' \code{project_tensor} and \code{reconstruct_tensor}.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD}} or
#'   \code{\link{CxtEBTD_missing}}.
#'
#' @return A named list with:
#' \describe{
#'   \item{El}{Normalized and reordered observation loading matrix
#'     (\code{n x K}). Each column has unit Frobenius norm.}
#'   \item{Ef}{Normalized and reordered time factor matrix
#'     (\code{p x K}). Each column has unit Frobenius norm.}
#'   \item{Ew}{Normalized and reordered channel weight matrix
#'     (\code{w x K}). Each column has unit Frobenius norm.}
#'   \item{lambda}{Numeric vector of length \code{K} giving the component
#'     weights \eqn{\lambda_k = \|\mathbf{l}_k\| \|\mathbf{f}_k\| \|\mathbf{w}_k\|},
#'     sorted descending.}
#'   \item{order}{Integer vector of length \code{K} giving the permutation
#'     of original component indices in new (descending \eqn{\lambda}) order.}
#' }
#'
#' @seealso \code{\link{CxtEBTD}}, \code{\link{project_tensor}},
#'   \code{\link{reconstruct_tensor}}
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X, K = 3, Xcov = Xcov)
#'
#' norm_fit <- normalize_factors(fit)
#' norm_fit$lambda   # component weights, descending
#' norm_fit$order    # e.g. c(2, 1, 3) — original rank 2 is now rank 1
#' norm_fit$El       # normalized loading matrix (n x K)
#' }
#'
#' @export
normalize_factors <- function(fit) {

  El <- fit$res$ql$El
  Ef <- fit$res$qf$Ef
  Ew <- fit$res$qw$Ew

  K <- ncol(El)

  norm_l <- sqrt(colSums(El^2))
  norm_f <- sqrt(colSums(Ef^2))
  norm_w <- sqrt(colSums(Ew^2))

  # Normalize columns; leave zero-norm columns as-is
  for (k in seq_len(K)) {
    if (norm_l[k] > 0) El[, k] <- El[, k] / norm_l[k] else norm_l[k] <- 0
    if (norm_f[k] > 0) Ef[, k] <- Ef[, k] / norm_f[k] else norm_f[k] <- 0
    if (norm_w[k] > 0) Ew[, k] <- Ew[, k] / norm_w[k] else norm_w[k] <- 0
  }

  lambda <- norm_l * norm_f * norm_w

  ord <- order(lambda, decreasing = TRUE)

  list(
    El     = El[, ord, drop = FALSE],
    Ef     = Ef[, ord, drop = FALSE],
    Ew     = Ew[, ord, drop = FALSE],
    lambda = lambda[ord],
    order  = ord
  )
}


#' Project a New Tensor onto the Learned Factor Space
#'
#' @description
#' Given a fitted \code{\link{CxtEBTD}} object, projects one or more new
#' observations onto the learned factor space to produce a loading matrix.
#' Each observation is a \code{p x w} count slice; the function recovers the
#' \code{K}-dimensional representation without refitting.
#'
#' The projection for observation \eqn{i} and component \eqn{k} is:
#' \deqn{F_{ik} = \lambda_k \, \mathbf{f}_k^\top X_{\text{new}[i,,]} \, \mathbf{w}_k}
#' when \code{normalize = TRUE} (columns ordered by \eqn{\lambda} descending),
#' or
#' \deqn{F_{ik} = \mathbf{f}_k^\top X_{\text{new}[i,,]} \, \mathbf{w}_k}
#' using the raw posterior means when \code{normalize = FALSE}.
#'
#' Computation is vectorized over observations: \code{X_new} is reshaped to
#' an \code{n_new x (p*w)} matrix and each Kronecker factor
#' \eqn{\mathbf{f}_k \otimes \mathbf{w}_k} is a \code{p*w}-vector, so the
#' full projection for component \eqn{k} reduces to a single matrix-vector
#' multiply.
#'
#' @param X_new Either a 3D array of dimensions \code{n_new x p x w} (multiple
#'   new observations), or a matrix of dimensions \code{p x w} (a single
#'   observation). In the latter case the input is wrapped as
#'   \code{array(X_new, dim = c(1, p, w))} and a length-\code{K} vector is
#'   returned instead of a \code{1 x K} matrix.
#' @param fit A fitted object returned by \code{\link{CxtEBTD}} or
#'   \code{\link{CxtEBTD_missing}}.
#' @param normalize Logical. If \code{TRUE} (default), calls
#'   \code{\link{normalize_factors}} so that factor columns have unit
#'   Frobenius norm and are ordered by \eqn{\lambda} descending. If
#'   \code{FALSE}, uses raw posterior means \code{fit$res$qf$Ef} and
#'   \code{fit$res$qw$Ew} with no scaling or reordering.
#'
#' @return If \code{X_new} is a 3D array: an \code{n_new x K} numeric matrix
#'   of projected loadings. If \code{X_new} is a \code{p x w} matrix: a
#'   numeric vector of length \code{K}.
#'
#' @seealso \code{\link{normalize_factors}}, \code{\link{reconstruct_tensor}},
#'   \code{\link{CxtEBTD}}
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X_train, K = 3)
#'
#' # Project a batch of new observations (50 x p x w array)
#' L_new <- project_tensor(X_new, fit)          # 50 x 3 matrix
#'
#' # Project a single p x w slice
#' l_one <- project_tensor(X_new[1, , ], fit)   # length-3 vector
#'
#' # Without normalization (raw posterior means, original component order)
#' L_raw <- project_tensor(X_new, fit, normalize = FALSE)
#' }
#'
#' @export
project_tensor <- function(X_new, fit, normalize = TRUE) {

  # Handle single p x w matrix input
  single_obs <- is.matrix(X_new)
  if (single_obs) {
    X_new <- array(X_new, dim = c(1L, nrow(X_new), ncol(X_new)))
  }

  n_new <- dim(X_new)[1L]
  p_new <- dim(X_new)[2L]
  w_new <- dim(X_new)[3L]

  # Validate dimensions against fit
  p_fit <- nrow(fit$res$qf$Ef)
  w_fit <- nrow(fit$res$qw$Ew)

  if (p_new != p_fit) {
    stop(sprintf(
      "X_new has p = %d time points but fit was trained with p = %d.",
      p_new, p_fit
    ))
  }
  if (w_new != w_fit) {
    stop(sprintf(
      "X_new has w = %d channels but fit was trained with w = %d.",
      w_new, w_fit
    ))
  }

  # Get factors and lambda
  if (normalize) {
    nf     <- normalize_factors(fit)
    Ef     <- nf$Ef      # p x K, unit-norm columns
    Ew     <- nf$Ew      # w x K, unit-norm columns
    lambda <- nf$lambda  # length K
  } else {
    Ef     <- fit$res$qf$Ef
    Ew     <- fit$res$qw$Ew
    lambda <- NULL
  }

  K <- ncol(Ef)

  # Reshape X_new to n_new x (p*w) for vectorized projection
  X_mat <- matrix(X_new, nrow = n_new)  # n_new x (p*w)

  # For each k, Kronecker factor f_k %o% w_k is a (p*w)-vector;
  # projection for all i is X_mat %*% kron_k
  proj <- matrix(0, nrow = n_new, ncol = K)
  for (k in seq_len(K)) {
    kron_k    <- as.vector(Ef[, k] %o% Ew[, k])  # p*w vector
    proj[, k] <- X_mat %*% kron_k
  }

  if (normalize) {
    proj <- t(t(proj) * lambda)  # scale column k by lambda[k]
  }

  if (single_obs) {
    return(as.vector(proj))
  }
  proj
}


#' Reconstruct the Denoised Mean Tensor from a CxtEBTD fit
#'
#' @description
#' Reconstructs the denoised Poisson mean tensor
#' \eqn{\hat{X}[i,j,m] = \sum_{k=1}^{K} L_{ik} F_{jk} W_{mk}}
#' from the raw (unnormalized) posterior mean factor matrices stored in a
#' fitted \code{\link{CxtEBTD}} object. The result is the best rank-\code{K}
#' approximation to the observed tensor under the fitted model.
#'
#' Implementation uses a loop over \code{K} components. For each \eqn{k} the
#' outer product \eqn{\mathbf{l}_k \otimes \mathbf{f}_k \otimes \mathbf{w}_k}
#' is added to the accumulator via \code{tcrossprod} and a sweep, avoiding
#' allocation of \code{K} separate \code{n x p x w} arrays.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD}} or
#'   \code{\link{CxtEBTD_missing}}.
#'
#' @return A numeric array of dimensions \code{n x p x w} containing the
#'   reconstructed denoised mean tensor \eqn{\hat{X}}.
#'
#' @seealso \code{\link{normalize_factors}}, \code{\link{project_tensor}},
#'   \code{\link{CxtEBTD}}
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X, K = 3)
#'
#' X_hat <- reconstruct_tensor(fit)
#' dim(X_hat)          # same as dim(X)
#'
#' # Mean squared reconstruction error
#' mean((X - X_hat)^2)
#' }
#'
#' @export
reconstruct_tensor <- function(fit) {

  El <- fit$res$ql$El
  Ef <- fit$res$qf$Ef
  Ew <- fit$res$qw$Ew

  n <- nrow(El)
  p <- nrow(Ef)
  w <- nrow(Ew)
  K <- ncol(El)

  X_hat <- array(0, dim = c(n, p, w))

  for (k in seq_len(K)) {
    # outer(n×p matrix, w-vector) produces an n×p×w array per R outer() semantics
    lf_k  <- El[, k, drop = FALSE] %*% t(Ef[, k, drop = FALSE])  # n x p
    X_hat <- X_hat + outer(lf_k, Ew[, k])                         # n x p x w
  }

  X_hat
}


#' Initialize Factor Matrices via CP-APR (Poisson Tensor Factorization)
#'
#' @description
#' Computes non-negative CP-APR (Alternating Poisson Regression) factor matrices
#' to use as warm-start initialization for \code{\link{CxtEBTD}}. Calls
#' \code{pyCP_APR} via \code{reticulate}, which requires Python 3.9+ with
#' \code{pyCP_APR} and \code{numpy} installed.
#'
#' @param X A 3-dimensional non-negative integer array of dimensions
#'   \code{n x p x w}.
#' @param K Integer. Number of components (CP rank).
#' @param n_iters Integer. Maximum number of CP-APR iterations. Default \code{150}.
#' @param method Character string. Optimization backend for pyCP_APR: \code{'torch'}
#'   (default) or \code{'numpy'}.
#' @param random_state Integer. Random seed for reproducibility. Default \code{42}.
#' @param virtualenv Character string. Name of the Python virtual environment
#'   containing pyCP_APR. Default \code{'ebtd1'}.
#'
#' @return A list of three matrices \code{list(L, F, W)} with dimensions
#'   \code{n x K}, \code{p x K}, \code{w x K}, suitable for passing as the
#'   \code{init} argument to \code{\link{CxtEBTD}} or \code{\link{CxtEBTD_missing}}.
#'
#' @details
#' Requires the \code{reticulate} R package and a Python environment with
#' \code{pyCP_APR} installed. To set up:
#' \preformatted{
#' pip install pyCP_APR numpy
#' }
#' The function adds a tiny constant (1e-7) to the last tensor slice
#' \code{X[, p, w]} to prevent pyCP_APR from failing on all-zero boundary
#' slices. pyCP_APR returns factor matrices with a zero-padded first row
#' (0-indexed), which is automatically removed before returning.
#'
#' @examples
#' \dontrun{
#' X <- array(rpois(100 * 20 * 10, lambda = 1.5), dim = c(100, 20, 10))
#' init <- init_cpapr(X, K = 3)
#' fit <- CxtEBTD(X, K = 3, init = init)
#' }
#'
#' @seealso \code{\link{CxtEBTD}}
#'
#' @export
init_cpapr <- function(X, K, n_iters = 150, method = 'torch',
                       random_state = 42, virtualenv = 'ebtd1') {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for CP-APR initialization. ",
         "Install it with install.packages('reticulate').")
  }

  reticulate::use_virtualenv(virtualenv, required = TRUE)

  modules <- tryCatch({
    list(
      cpapr = reticulate::import("pyCP_APR"),
      np    = reticulate::import("numpy")
    )
  }, error = function(e) {
    stop("Failed to import pyCP_APR or numpy from virtualenv '", virtualenv, "'. ",
         "Install them with: pip install pyCP_APR numpy\n",
         "Original error: ", conditionMessage(e))
  })

  cpapr <- modules$cpapr

  # Add tiny constant to last slice to avoid all-zero boundary issues
  X[, dim(X)[2], dim(X)[3]] <- X[, dim(X)[2], dim(X)[3]] + 1e-7

  # Convert to sparse coordinate format
  X_list <- rbind_sparse_matrix(X, reindex = TRUE)
  coords  <- reticulate::r_to_py(X_list[, c("V1", "V2", "V3")])
  coords  <- coords$to_numpy()
  nnz     <- reticulate::r_to_py(X_list[, "v"])
  nnz     <- nnz[["v"]]$values

  # Fit CP-APR
  model <- cpapr$CP_APR(
    n_iters      = as.integer(n_iters),
    random_state = as.integer(random_state),
    verbose      = 0L,
    method       = method,
    return_type  = 'numpy'
  )
  Mj <- model$fit(coords = coords, values = nnz, rank = as.integer(K))

  # Extract factors — drop the zero-padded first row (0-indexed artifact)
  L     <- Mj$Factors[[1]][-1, , drop = FALSE]
  F_mat <- Mj$Factors[[2]][-1, , drop = FALSE]
  W     <- Mj$Factors[[3]][-1, , drop = FALSE]

  list(L, F_mat, W)
}


#' Two-Step Covariate Selection via Unsupervised Decomposition
#'
#' @description
#' Implements a two-step procedure for identifying which covariates are
#' relevant for each tensor component. Step 1: fit an unsupervised
#' \code{\link{CxtEBTD}} decomposition (no covariates). Step 2: for each
#' component k, regress the estimated loadings \eqn{E[l_{ik}]} on the
#' candidate covariate matrix using OLS, and select covariates whose
#' coefficients are significant at level \code{alpha}.
#'
#' @param X A 3-dimensional non-negative integer array of dimensions
#'   \code{n x p x w}.
#' @param K Integer. Number of components (CP rank).
#' @param Xcov_candidates Numeric matrix of dimension \code{n x q} containing
#'   candidate covariates to screen. Column names are used in the output if
#'   available.
#' @param alpha Numeric significance level for covariate selection. Default
#'   \code{0.05}.
#' @param verbose Logical. If \code{TRUE}, prints progress and selected
#'   covariates per rank. Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{CxtEBTD}} for the
#'   unsupervised fit (e.g., \code{init}, \code{maxiter}, \code{tol},
#'   \code{convergence_criteria}).
#'
#' @return A named list with:
#' \describe{
#'   \item{selected}{A list of length K. Each element is an integer vector of
#'     column indices of \code{Xcov_candidates} whose coefficients are
#'     significant at level \code{alpha} for that component. Empty integer
#'     vector if no covariates are significant.}
#'   \item{summaries}{A list of length K. Each element is the \code{summary.lm}
#'     object from the OLS regression of El[,k] on Xcov_candidates, giving
#'     full coefficient estimates, standard errors, t-statistics, and p-values.}
#'   \item{fit_unsupervised}{The fitted unsupervised \code{\link{CxtEBTD}}
#'     object, in case the user wants to inspect the decomposition.}
#' }
#'
#' @details
#' The unsupervised fit uses \code{Xcov = NULL} with all other arguments
#' passed via \code{...}. The OLS regression for each component k is:
#' \deqn{E[l_{ik}] = X_{\text{cov}} \beta_k + \epsilon_{ik}}
#' Covariates are selected if their two-sided p-value is below \code{alpha}.
#' The intercept is included in the regression but is never selected as a
#' covariate (it is excluded from the returned indices).
#'
#' This is a screening procedure, not a formal statistical test. For rigorous
#' inference on covariate effects, fit \code{\link{CxtEBTD}} with the selected
#' covariates and use the estimated \eqn{\gamma} coefficients from the
#' generative model.
#'
#' @examples
#' \dontrun{
#' X <- array(rpois(100 * 20 * 10, lambda = 1.5), dim = c(100, 20, 10))
#' Xcov <- matrix(rnorm(100 * 5), nrow = 100)
#' colnames(Xcov) <- paste0("cov", 1:5)
#'
#' result <- select_covariates(X, K = 3, Xcov_candidates = Xcov,
#'                              maxiter = 20, convergence_criteria = 'ELBO')
#'
#' # Which covariates were selected for each rank?
#' result$selected
#'
#' # Full regression summary for rank 1
#' result$summaries[[1]]
#'
#' # Refit with selected covariates for rank 1
#' Xcov_selected <- Xcov[, result$selected[[1]], drop = FALSE]
#' fit <- CxtEBTD(X, K = 3, Xcov = Xcov_selected)
#' }
#'
#' @seealso \code{\link{CxtEBTD}}
#' @export
select_covariates <- function(X, K, Xcov_candidates, alpha = 0.05,
                               verbose = TRUE, ...) {

  # Step 1: unsupervised decomposition
  if (verbose) cat("Step 1: fitting unsupervised CxtEBTD (K =", K, ")...\n")
  fit_unsup <- CxtEBTD(X, K = K, Xcov = NULL, ...)
  El <- fit_unsup$res$ql$El  # n x K

  cov_names <- colnames(Xcov_candidates)
  if (is.null(cov_names)) cov_names <- paste0("V", seq_len(ncol(Xcov_candidates)))

  selected_list <- vector("list", K)
  summary_list  <- vector("list", K)

  # Step 2: OLS regression of El[,k] on Xcov_candidates for each component
  for (k in seq_len(K)) {
    if (verbose) cat(sprintf("Step 2: screening covariates for rank %d...\n", k))

    lm_fit  <- lm(El[, k] ~ Xcov_candidates)
    lm_summ <- summary(lm_fit)

    # p-values for all coefficients; row 1 is the intercept — skip it
    pvals <- lm_summ$coefficients[-1, 4]

    # Selected: column indices (1-indexed into Xcov_candidates) where p < alpha
    selected_idx <- which(pvals < alpha)

    if (verbose) {
      if (length(selected_idx) == 0) {
        cat(sprintf("  Rank %d: no covariates selected at alpha = %g\n", k, alpha))
      } else {
        cat(sprintf("  Rank %d: selected covariates: %s\n", k,
                    paste(cov_names[selected_idx], collapse = ", ")))
      }
    }

    selected_list[[k]] <- selected_idx
    summary_list[[k]]  <- lm_summ
  }

  list(
    selected         = selected_list,
    summaries        = summary_list,
    fit_unsupervised = fit_unsup
  )
}


#' Match and Align Estimated Factors to Reference Factors
#'
#' @description
#' Finds the optimal permutation that aligns the columns of estimated factor
#' matrices to reference (true) factor matrices using the Tucker congruence
#' coefficient (cosine similarity) and the Hungarian algorithm. This is the
#' standard method for comparing tensor decompositions in simulation studies,
#' where the component ordering is arbitrary.
#'
#' For each pair of columns (one from reference, one from estimated), the
#' congruence coefficient is \eqn{a^\top b / (\|a\| \|b\|)}. When multiple
#' modes are provided, the joint congruence is the product of per-mode
#' congruences, so a good match must be good on all modes simultaneously.
#'
#' @param ref A list of reference factor matrices, one per mode. Each matrix
#'   has dimensions \code{d_m x K}. Typically the true factors from a
#'   simulation.
#' @param est A list of estimated factor matrices, same structure as \code{ref}.
#'   Dimensions must match: same number of modes, same \code{d_m} per mode,
#'   same K.
#' @param absolute_value Logical. If \code{TRUE} (default), uses absolute
#'   cosine similarity, ignoring sign differences. Appropriate for
#'   non-negative decompositions.
#'
#' @return A named list with:
#' \describe{
#'   \item{permutation}{Integer vector of length K giving the optimal column
#'     permutation. \code{permutation[i]} is the column index in \code{est}
#'     that best matches column \code{i} in \code{ref}.}
#'   \item{mean_congruence}{Numeric scalar: mean congruence coefficient across
#'     all K matched pairs. 1.0 = perfect recovery.}
#'   \item{per_component}{Numeric vector of length K giving the congruence
#'     for each matched pair.}
#' }
#'
#' @details
#' Requires the \code{clue} package for the Hungarian algorithm
#' (\code{clue::solve_LSAP}). If not installed, falls back to a greedy
#' matching algorithm with a warning.
#'
#' @examples
#' \dontrun{
#' # Simulation: compare estimated to true factors
#' sim <- simulate_tensor(n = 100, p = 20, w = 10, K = 3)
#' fit <- CxtEBTD(sim$X, K = 3, Xcov = sim$Xcov)
#' nf <- normalize_factors(fit)
#'
#' result <- match_factors(
#'   ref = list(sim$U1_true, sim$U2_true, sim$U3_true),
#'   est = list(nf$El, nf$Ef, nf$Ew)
#' )
#' result$mean_congruence   # overall recovery quality
#' result$permutation       # which estimated component matches which true one
#' }
#'
#' @seealso \code{\link{normalize_factors}}, \code{\link{simulate_tensor}}
#' @export
match_factors <- function(ref, est, absolute_value = TRUE) {

  # --- input validation ---
  if (!is.list(ref) || !is.list(est)) {
    stop("'ref' and 'est' must both be lists of factor matrices.")
  }
  M <- length(ref)
  if (length(est) != M) {
    stop("'ref' and 'est' must have the same number of modes.")
  }
  K <- ncol(ref[[1]])
  for (m in seq_len(M)) {
    if (!is.matrix(ref[[m]]) || !is.matrix(est[[m]])) {
      stop(sprintf("ref[[%d]] and est[[%d]] must be matrices.", m, m))
    }
    if (nrow(ref[[m]]) != nrow(est[[m]])) {
      stop(sprintf("Mode %d: ref has %d rows but est has %d rows.",
                   m, nrow(ref[[m]]), nrow(est[[m]])))
    }
    if (ncol(ref[[m]]) != K || ncol(est[[m]]) != K) {
      stop(sprintf("Mode %d: all matrices must have K = %d columns.", m, K))
    }
  }

  # --- normalize columns to unit norm, guarding against zero-norm columns ---
  normalize_cols <- function(M_mat) {
    norms <- sqrt(colSums(M_mat^2))
    norms[norms == 0] <- 1  # leave zero columns unchanged
    sweep(M_mat, 2, norms, "/")
  }

  ref_n <- lapply(ref, normalize_cols)
  est_n <- lapply(est, normalize_cols)

  # --- per-mode K x K cosine similarity matrices, then joint product ---
  S_joint <- matrix(1, nrow = K, ncol = K)
  for (m in seq_len(M)) {
    S_m <- t(ref_n[[m]]) %*% est_n[[m]]  # K x K
    if (absolute_value) S_m <- abs(S_m)
    S_joint <- S_joint * S_m
  }

  # --- optimal assignment via Hungarian algorithm (clue) or greedy fallback ---
  if (requireNamespace("clue", quietly = TRUE)) {
    perm <- as.integer(clue::solve_LSAP(S_joint, maximum = TRUE))
  } else {
    warning("Package 'clue' not available; using greedy matching (may be suboptimal). ",
            "Install with install.packages('clue') for optimal results.")
    perm    <- integer(K)
    used    <- logical(K)
    # sort rows by their maximum similarity so strongest matches go first
    row_ord <- order(apply(S_joint, 1, max), decreasing = TRUE)
    for (i in row_ord) {
      scores      <- S_joint[i, ]
      scores[used] <- -Inf
      j           <- which.max(scores)
      perm[i]     <- j
      used[j]     <- TRUE
    }
  }

  # --- per-component congruences from the matched pairs ---
  per_component <- S_joint[cbind(seq_len(K), perm)]

  list(
    permutation     = perm,
    mean_congruence = mean(per_component),
    per_component   = per_component
  )
}


#' Simulate a Count Tensor with Known Ground Truth
#'
#' @description
#' Generates a synthetic Poisson count tensor from a CP structure with known
#' factor matrices, optional covariate effects, and spike-and-slab sparsity.
#' Useful for simulation studies evaluating decomposition accuracy.
#'
#' @param n Integer. Number of observations (mode 1). Default \code{100}.
#' @param p Integer. Number of time points (mode 2). Default \code{20}.
#' @param w Integer. Number of channels (mode 3). Default \code{10}.
#' @param K Integer. Number of components. Default \code{3}.
#' @param Xcov Optional numeric covariate matrix of dimension \code{n x q}.
#'   If \code{NULL}, covariates are generated automatically: an intercept,
#'   a binary covariate, and a continuous covariate. Default \code{NULL}.
#' @param gamma_true Optional list of length K, where each element is a
#'   numeric vector of covariate coefficients for that component. If
#'   \code{NULL} and \code{Xcov} is also \code{NULL}, random coefficients
#'   are generated. If \code{Xcov} is provided, must be supplied.
#'   Set to \code{FALSE} to generate an unsupervised tensor (no covariate
#'   effects). Default \code{NULL}.
#' @param pi0 Numeric. Spike (structural zero) probability for mode-1
#'   factors. Default \code{0.2}.
#' @param alpha_true Numeric. Gamma shape parameter for mode-1 slab.
#'   Default \code{3}.
#' @param beta_true Numeric. Gamma rate parameter for mode-1 slab.
#'   Default \code{2}.
#' @param weights Numeric vector of length K giving component weights.
#'   Default \code{NULL} generates decreasing weights.
#' @param sparsity Numeric scaling factor controlling overall tensor
#'   sparsity (higher = sparser). Default \code{20}.
#' @param seed Integer random seed. Default \code{42}.
#'
#' @return A named list with:
#' \describe{
#'   \item{X}{The observed Poisson count tensor, \code{n x p x w}.}
#'   \item{lambda_true}{The true Poisson rate tensor, \code{n x p x w}.}
#'   \item{U1_true}{True mode-1 factor matrix (normalized to unit norm),
#'     \code{n x K}.}
#'   \item{U2_true}{True mode-2 factor matrix (normalized), \code{p x K}.}
#'   \item{U3_true}{True mode-3 factor matrix (normalized), \code{w x K}.}
#'   \item{weights}{Component weights, length K.}
#'   \item{Xcov}{Covariate matrix used, \code{n x q}. \code{NULL} if
#'     unsupervised.}
#'   \item{gamma_true}{List of true covariate coefficient vectors.
#'     \code{NULL} if unsupervised.}
#'   \item{sparsity_pct}{Percentage of zeros in \code{X}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Supervised simulation
#' sim <- simulate_tensor(n = 100, p = 20, w = 10, K = 3)
#' fit <- CxtEBTD(sim$X, K = 3, Xcov = sim$Xcov)
#'
#' # Unsupervised simulation
#' sim0 <- simulate_tensor(n = 100, p = 20, w = 10, K = 2, gamma_true = FALSE)
#' fit0 <- CxtEBTD(sim0$X, K = 2)
#' }
#'
#' @seealso \code{\link{CxtEBTD}}, \code{\link{match_factors}}
#' @export
simulate_tensor <- function(n = 100, p = 20, w = 10, K = 3,
                             Xcov = NULL, gamma_true = NULL,
                             pi0 = 0.2, alpha_true = 3, beta_true = 2,
                             weights = NULL, sparsity = 20, seed = 42) {

  set.seed(seed)

  unsupervised <- isFALSE(gamma_true)

  # --- covariates ---
  if (unsupervised) {
    Xcov       <- NULL
    gamma_true <- NULL
  } else {
    if (is.null(Xcov)) {
      Xcov <- cbind(1, rbinom(n, 1, 0.5), runif(n, 20, 80))
    }
    if (is.null(gamma_true)) {
      gamma_true <- lapply(seq_len(K), function(k) {
        c(0, runif(1, -1, 1), runif(1, -0.3, 0.3))
      })
    }
    if (length(gamma_true) != K) {
      stop("'gamma_true' must be a list of length K or FALSE.")
    }
    if (ncol(Xcov) != length(gamma_true[[1]])) {
      stop("Number of columns in Xcov must match length of gamma_true vectors.")
    }
  }

  # --- mode-1 factors (U1_true): covariate-dependent spike-and-slab ---
  U1_true <- matrix(0, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    if (unsupervised) {
      lambda_i <- rep(1, n)
    } else {
      lambda_i <- exp(Xcov %*% gamma_true[[k]])
    }
    active <- rbinom(n, 1, 1 - pi0)
    for (i in seq_len(n)) {
      if (active[i] == 1) {
        U1_true[i, k] <- rgamma(1, shape = alpha_true,
                                rate  = beta_true / lambda_i[i])
      }
    }
    col_norm <- sqrt(sum(U1_true[, k]^2))
    if (col_norm > 0) U1_true[, k] <- U1_true[, k] / col_norm
  }

  # --- mode-2 factors (U2_true): block structure ---
  U2_true    <- matrix(0, nrow = p, ncol = K)
  block_size <- floor(p / K)
  for (k in seq_len(K)) {
    start_idx <- (k - 1) * block_size + 1
    end_idx   <- if (k == K) p else k * block_size
    block_len <- end_idx - start_idx + 1
    U2_true[start_idx:end_idx, k] <- 1 +
      rgamma(block_len, shape = 2, rate = 2)
    col_norm <- sqrt(sum(U2_true[, k]^2))
    if (col_norm > 0) U2_true[, k] <- U2_true[, k] / col_norm
  }

  # --- mode-3 factors (U3_true): sparse group structure ---
  U3_true   <- matrix(0, nrow = w, ncol = K)
  group_size <- max(1L, floor(w / K))
  for (k in seq_len(K)) {
    start_idx <- (k - 1) * group_size + 1
    end_idx   <- if (k == K) w else k * group_size
    grp_len   <- end_idx - start_idx + 1
    U3_true[start_idx:end_idx, k] <- rgamma(grp_len, shape = 2, rate = 1)
    col_norm <- sqrt(sum(U3_true[, k]^2))
    if (col_norm > 0) U3_true[, k] <- U3_true[, k] / col_norm
  }

  # --- component weights ---
  if (is.null(weights)) {
    weights <- seq(2, 0.5, length.out = K)
  }

  # --- true rate tensor ---
  scalenpw   <- sqrt(n) * sqrt(p) * sqrt(w)
  lambda_true <- array(0, dim = c(n, p, w))
  for (k in seq_len(K)) {
    lambda_true <- lambda_true +
      weights[k] * (U1_true[, k] %o% U2_true[, k] %o% U3_true[, k])
  }
  lambda_true <- lambda_true * scalenpw / sparsity

  # --- observed tensor ---
  X <- array(rpois(n * p * w, lambda_true), dim = c(n, p, w))

  sparsity_pct <- 100 * mean(X == 0)

  list(
    X            = X,
    lambda_true  = lambda_true,
    U1_true      = U1_true,
    U2_true      = U2_true,
    U3_true      = U3_true,
    weights      = weights,
    Xcov         = Xcov,
    gamma_true   = gamma_true,
    sparsity_pct = sparsity_pct
  )
}
