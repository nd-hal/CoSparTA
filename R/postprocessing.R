#' Retrieve Factor Matrices or Weights from a CoSparTA Fit
#'
#' @description
#' Unified getter that retrieves a named output from a fitted
#' \code{\link{CoSparTA}} object by logical name, so users do not need to
#' navigate the internal slot structure (e.g.\ \code{fit$res$ql$El_normed}).
#'
#' @param fit A fitted object returned by \code{\link{CoSparTA}} or
#'   \code{\link{CoSparTA_missing}}.
#' @param what Character string specifying what to retrieve. One of:
#' \describe{
#'   \item{\code{"U1"}}{Observation loading matrix (\code{n x K}).}
#'   \item{\code{"U2"}}{Time factor matrix (\code{p x K}).}
#'   \item{\code{"U3"}}{Channel weight matrix (\code{w x K}).}
#'   \item{\code{"weight"}}{Component weight vector \eqn{\lambda} of length
#'     \code{K} (product of mode norms; sorted descending when
#'     \code{normalized = TRUE}).}
#'   \item{\code{"gamma"}}{List of length \code{K} of covariate coefficient
#'     vectors (one per component). \code{NULL} entries indicate unsupervised
#'     components.}
#' }
#' @param normalized Logical. If \code{TRUE} (default), returns the
#'   normalized and reordered version stored in the \code{_normed} fields
#'   (unit-norm columns, descending \eqn{\lambda} order). If \code{FALSE},
#'   returns the raw posterior estimates in original component order.
#'
#' @return The requested object:
#' \describe{
#'   \item{\code{"U1"}, \code{"U2"}, \code{"U3"}}{A numeric matrix.}
#'   \item{\code{"weight"}}{A numeric vector of length \code{K}.}
#'   \item{\code{"gamma"}}{A list of length \code{K}.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- CoSparTA(X, K = 3, Xcov = Xcov)
#'
#' # Normalized observation loadings (n x 3, unit-norm columns)
#' U1 <- get_loadings(fit, "U1")
#'
#' # Raw observation loadings
#' U1_raw <- get_loadings(fit, "U1", normalized = FALSE)
#'
#' # Component weights (descending)
#' w <- get_loadings(fit, "weight")
#'
#' # Covariate coefficient lists
#' gamma <- get_loadings(fit, "gamma")
#' gamma[[1]]  # gamma for component 1
#' }
#'
#' @seealso \code{\link{CoSparTA}}, \code{\link{normalize_factors}}
#' @export
get_loadings <- function(fit, what = "U1", normalized = TRUE) {
  valid <- c("U1", "U2", "U3", "weight", "gamma")
  what  <- match.arg(what, valid)

  if (normalized) {
    switch(what,
      "U1"     = fit$res$ql$El_normed,
      "U2"     = fit$res$qf$Ef_normed,
      "U3"     = fit$res$qw$Ew_normed,
      "weight" = fit$res$lambda_normed,
      "gamma"  = fit$res$gl_normed
    )
  } else {
    switch(what,
      "U1"     = fit$res$ql$El,
      "U2"     = fit$res$qf$Ef,
      "U3"     = fit$res$qw$Ew,
      "weight" = fit$res$lambda,
      "gamma"  = fit$res$gl
    )
  }
}


#' Normalize and Sort Components of a CoSparTA fit
#'
#' @description
#' Extracts the factor matrices from a fitted \code{\link{CoSparTA}} object,
#' normalizes each column to unit Frobenius norm, computes a per-component
#' weight \eqn{\lambda_k = \|\mathbf{l}_k\| \|\mathbf{f}_k\| \|\mathbf{w}_k\|},
#' and returns components sorted by \eqn{\lambda} descending. This is the
#' canonical form for comparing decompositions and for downstream use with
#' \code{project_tensor} and \code{reconstruct_tensor}.
#'
#' Since \code{\link{CoSparTA}} now produces normalized outputs directly in
#' \code{_normed} fields (e.g. \code{fit$res$ql$El_normed}), this function is
#' primarily useful for legacy code or custom normalization workflows applied
#' to fit objects that lack those fields.
#'
#' @param fit A fitted object returned by \code{\link{CoSparTA}} or
#'   \code{\link{CoSparTA_missing}}.
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
#'   \item{gamma_list}{A list of length \code{K}. Each element is the gamma
#'     coefficient vector for that component after reordering by descending
#'     lambda. \code{NULL} if no components have gamma estimates (i.e. the
#'     model was fit without covariates).}
#'   \item{lambda_order}{Integer vector of length \code{K} giving the
#'     permutation used to reorder components by descending lambda. Useful
#'     for reordering other component-indexed objects (e.g. \code{Xcov_list})
#'     to match the normalized factor ordering.}
#' }
#'
#' @seealso \code{\link{CoSparTA}}, \code{\link{project_tensor}},
#'   \code{\link{reconstruct_tensor}}
#'
#' @examples
#' \dontrun{
#' fit <- CoSparTA(X, K = 3, Xcov = Xcov)
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

  # Extract and reorder gamma coefficients if any component has them
  gl        <- fit$res$gl
  has_gamma <- !is.null(gl) &&
    any(vapply(gl, function(g) !is.null(g$gamma), logical(1)))
  if (has_gamma) {
    gamma_list         <- lapply(seq_len(K), function(k) {
      if (!is.null(gl[[k]]$gamma)) gl[[k]]$gamma else NA
    })
    gamma_list_ordered <- gamma_list[ord]
  } else {
    gamma_list_ordered <- NULL
  }

  list(
    El           = El[, ord, drop = FALSE],
    Ef           = Ef[, ord, drop = FALSE],
    Ew           = Ew[, ord, drop = FALSE],
    lambda       = lambda[ord],
    order        = ord,
    lambda_order = ord,
    gamma_list   = gamma_list_ordered
  )
}


#' Project a New Tensor onto the Learned Factor Space
#'
#' @description
#' Given a fitted \code{\link{CoSparTA}} object (or raw factor matrices),
#' projects one or more new observations onto the learned factor space to
#' produce a loading matrix. Each observation is a \code{p x w} count slice;
#' the function recovers the \code{K}-dimensional representation without
#' refitting.
#'
#' The projection for observation \eqn{i} and component \eqn{k} is:
#' \deqn{F_{ik} = \lambda_k \, \mathbf{f}_k^\top X_{\text{new}[i,,]} \, \mathbf{w}_k}
#' when \code{normalize = TRUE} or raw matrices with \code{lambda} are supplied,
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
#' @param fit A fitted object returned by \code{\link{CoSparTA}} or
#'   \code{\link{CoSparTA_missing}}. Either \code{fit} or all of \code{Ef},
#'   \code{Ew}, and \code{lambda} must be supplied.
#' @param normalize Logical. If \code{TRUE} (default) and \code{fit} is
#'   supplied, calls \code{\link{normalize_factors}} so that factor columns have
#'   unit Frobenius norm and are ordered by \eqn{\lambda} descending. If
#'   \code{FALSE}, uses raw posterior means with no scaling or reordering.
#'   Ignored when \code{Ef}, \code{Ew}, and \code{lambda} are provided directly.
#' @param Ef Numeric matrix of dimensions \code{p x K} (time factor matrix,
#'   unit-norm columns). Provide together with \code{Ew} and \code{lambda} to
#'   bypass the \code{fit} object.
#' @param Ew Numeric matrix of dimensions \code{w x K} (channel factor matrix,
#'   unit-norm columns). Provide together with \code{Ef} and \code{lambda} to
#'   bypass the \code{fit} object.
#' @param lambda Numeric vector of length \code{K} (component weights). Provide
#'   together with \code{Ef} and \code{Ew} to bypass the \code{fit} object.
#'
#' @return If \code{X_new} is a 3D array: an \code{n_new x K} numeric matrix
#'   of projected loadings. If \code{X_new} is a \code{p x w} matrix: a
#'   numeric vector of length \code{K}.
#'
#' @seealso \code{\link{normalize_factors}}, \code{\link{reconstruct_tensor}},
#'   \code{\link{CoSparTA}}
#'
#' @examples
#' \dontrun{
#' fit <- CoSparTA(X_train, K = 3)
#'
#' # Project a batch of new observations (50 x p x w array)
#' L_new <- project_tensor(X_new, fit)          # 50 x 3 matrix
#'
#' # Project a single p x w slice
#' l_one <- project_tensor(X_new[1, , ], fit)   # length-3 vector
#'
#' # Without normalization (raw posterior means, original component order)
#' L_raw <- project_tensor(X_new, fit, normalize = FALSE)
#'
#' # Using raw factor matrices directly (already normalized)
#' nf <- normalize_factors(fit)
#' L_new2 <- project_tensor(X_new, Ef = nf$Ef, Ew = nf$Ew, lambda = nf$lambda)
#' }
#'
#' @export
project_tensor <- function(X_new, fit = NULL, normalize = TRUE,
                            Ef = NULL, Ew = NULL, lambda = NULL) {

  # Handle single p x w matrix input
  single_obs <- is.matrix(X_new)
  if (single_obs) {
    X_new <- array(X_new, dim = c(1L, nrow(X_new), ncol(X_new)))
  }

  n_new <- dim(X_new)[1L]
  p_new <- dim(X_new)[2L]
  w_new <- dim(X_new)[3L]

  # Resolve factor matrices
  if (!is.null(Ef) && !is.null(Ew) && !is.null(lambda)) {
    # raw matrices supplied directly — assume already normalized
  } else if (!is.null(fit)) {
    if (normalize && !is.null(fit$res$qf$Ef_normed)) {
      Ef     <- fit$res$qf$Ef_normed
      Ew     <- fit$res$qw$Ew_normed
      lambda <- fit$res$lambda_normed
    } else if (normalize) {
      nf     <- normalize_factors(fit)
      Ef     <- nf$Ef
      Ew     <- nf$Ew
      lambda <- nf$lambda
    } else {
      Ef     <- fit$res$qf$Ef
      Ew     <- fit$res$qw$Ew
      lambda <- NULL
    }
  } else {
    stop("Either 'fit' or all of 'Ef', 'Ew', and 'lambda' must be provided.")
  }

  # Validate dimensions
  if (p_new != nrow(Ef)) {
    stop(sprintf(
      "X_new has p = %d time points but Ef has p = %d.",
      p_new, nrow(Ef)
    ))
  }
  if (w_new != nrow(Ew)) {
    stop(sprintf(
      "X_new has w = %d channels but Ew has w = %d.",
      w_new, nrow(Ew)
    ))
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

  if (!is.null(lambda)) {
    proj <- t(t(proj) * lambda)  # scale column k by lambda[k]
  }

  if (single_obs) {
    return(as.vector(proj))
  }
  proj
}


#' Reconstruct the Denoised Mean Tensor from a CoSparTA fit
#'
#' @description
#' Reconstructs the denoised Poisson mean tensor
#' \eqn{\hat{X}[i,j,m] = \sum_{k=1}^{K} L_{ik} F_{jk} W_{mk}}
#' from the raw (unnormalized) posterior mean factor matrices stored in a
#' fitted \code{\link{CoSparTA}} object. The result is the best rank-\code{K}
#' approximation to the observed tensor under the fitted model.
#'
#' Implementation uses a loop over \code{K} components. For each \eqn{k} the
#' outer product \eqn{\mathbf{l}_k \otimes \mathbf{f}_k \otimes \mathbf{w}_k}
#' is added to the accumulator via \code{tcrossprod} and a sweep, avoiding
#' allocation of \code{K} separate \code{n x p x w} arrays.
#'
#' @param fit A fitted object returned by \code{\link{CoSparTA}} or
#'   \code{\link{CoSparTA_missing}}.
#' @param normalized Logical. If \code{TRUE} (default), reads from the
#'   \code{_normed} fields (\code{El_normed}, \code{Ef_normed},
#'   \code{Ew_normed}) and scales each component by
#'   \code{lambda_normed[k]}. If \code{FALSE}, uses raw posterior means
#'   directly (current behavior).
#'
#' @return A numeric array of dimensions \code{n x p x w} containing the
#'   reconstructed denoised mean tensor \eqn{\hat{X}}.
#'
#' @seealso \code{\link{normalize_factors}}, \code{\link{project_tensor}},
#'   \code{\link{CoSparTA}}
#'
#' @examples
#' \dontrun{
#' fit <- CoSparTA(X, K = 3)
#'
#' X_hat <- reconstruct_tensor(fit)
#' dim(X_hat)          # same as dim(X)
#'
#' # Mean squared reconstruction error
#' mean((X - X_hat)^2)
#' }
#'
#' @export
reconstruct_tensor <- function(fit, normalized = TRUE) {

  if (normalized && !is.null(fit$res$ql$El_normed)) {
    El     <- fit$res$ql$El_normed
    Ef     <- fit$res$qf$Ef_normed
    Ew     <- fit$res$qw$Ew_normed
    lambda <- fit$res$lambda_normed
  } else {
    El     <- fit$res$ql$El
    Ef     <- fit$res$qf$Ef
    Ew     <- fit$res$qw$Ew
    lambda <- NULL
  }

  n <- nrow(El)
  p <- nrow(Ef)
  w <- nrow(Ew)
  K <- ncol(El)

  X_hat <- array(0, dim = c(n, p, w))

  for (k in seq_len(K)) {
    lf_k  <- El[, k, drop = FALSE] %*% t(Ef[, k, drop = FALSE])  # n x p
    term  <- outer(lf_k, Ew[, k])                                  # n x p x w
    if (!is.null(lambda)) term <- term * lambda[k]
    X_hat <- X_hat + term
  }

  X_hat
}


#' Initialize Factor Matrices via CP-APR (Poisson Tensor Factorization)
#'
#' @description
#' Computes non-negative CP-APR (Alternating Poisson Regression) factor matrices
#' to use as warm-start initialization for \code{\link{CoSparTA}}. Calls
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
#'   \code{init} argument to \code{\link{CoSparTA}} or \code{\link{CoSparTA_missing}}.
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
#' fit <- CoSparTA(X, K = 3, init = init)
#' }
#'
#' @seealso \code{\link{CoSparTA}}
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
#' Identifies which covariates are relevant for each tensor component by
#' filtering active (non-zero) observations per factor, regressing
#' log-transformed loadings on a candidate covariate matrix using OLS, and
#' selecting covariates whose coefficients are significant at level
#' \code{alpha}. Zero loadings (spike/inactive observations) are excluded
#' before regression to match the log-linear DGP structure
#' \eqn{E[u \mid x] \propto \exp(x^\top \gamma)}. The loading matrix
#' \eqn{E[l_{ik}]} can be supplied in three ways: (1) provide \code{El}
#' directly, (2) provide a fitted \code{fit} object to extract \code{El} via
#' \code{\link{normalize_factors}}, or (3) provide \code{X} and \code{K} so
#' the function fits an unsupervised \code{\link{CoSparTA}} internally
#' (original two-step behaviour).
#'
#' @param X A 3-dimensional non-negative integer array of dimensions
#'   \code{n x p x w}. Required when neither \code{fit} nor \code{El} is
#'   supplied.
#' @param K Integer. Number of components (CP rank). When \code{El} is provided
#'   directly, \code{K} must equal \code{ncol(El)}.
#' @param covariate_data Numeric matrix of dimension \code{n x q} containing
#'   candidate covariates to screen. Column names are used in the output if
#'   available.
#' @param fit A fitted object returned by \code{\link{CoSparTA}} or
#'   \code{\link{CoSparTA_missing}}. If supplied, \code{El} is extracted via
#'   \code{normalize_factors(fit)$El} and no internal \code{CoSparTA} call is
#'   made. Either \code{fit} or \code{El} or neither (with \code{X} and
#'   \code{K}) must be provided.
#' @param El Numeric matrix of dimensions \code{n x K} (observation loading
#'   matrix). If provided, used directly — no \code{CoSparTA} fit is run and
#'   \code{fit} is ignored.
#' @param alpha Numeric significance level for covariate selection. Default
#'   \code{0.05}.
#' @param verbose Logical. If \code{TRUE}, prints progress and selected
#'   covariates per rank. Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{CoSparTA}} when
#'   running the internal unsupervised fit (e.g., \code{init}, \code{maxiter},
#'   \code{tol}, \code{convergence_criteria}). Ignored when \code{fit} or
#'   \code{El} is supplied.
#'
#' @return A named list with:
#' \describe{
#'   \item{selected}{A list of length K. Each element is an integer vector of
#'     column indices of \code{covariate_data} whose coefficients are
#'     significant at level \code{alpha} for that component. Empty integer
#'     vector if no covariates are significant.}
#'   \item{summaries}{A list of length K. Each element is the \code{summary.lm}
#'     object from the OLS regression of El[,k] on \code{covariate_data}, giving
#'     full coefficient estimates, standard errors, t-statistics, and p-values.}
#'   \item{fit_unsupervised}{The fitted unsupervised \code{\link{CoSparTA}}
#'     object when the internal fit was run, the supplied \code{fit} object when
#'     that was used, or \code{NULL} when \code{El} was provided directly.}
#'   \item{runtime_secs}{Numeric scalar giving the total elapsed wall-clock time
#'     of the function call in seconds.}
#' }
#'
#' @details
#' When running the internal unsupervised fit, \code{Xcov = NULL} is used and
#' all \code{...} arguments are forwarded to \code{\link{CoSparTA}}. For each
#' component k, rows with \eqn{l_{ik} = 0} (spike/inactive observations) are
#' dropped, and OLS is applied to the log-transformed loadings of the remaining
#' active rows:
#' \deqn{\log E[l_{ik}] = X_{\text{cov},i} \beta_k + \epsilon_{ik}, \quad l_{ik} > 0}
#' This matches the log-linear DGP structure of the generative model. Covariates
#' are selected if their two-sided p-value is below \code{alpha}. The intercept
#' is included in the regression but is never selected as a covariate (it is
#' excluded from the returned indices).
#'
#' This is a screening procedure, not a formal statistical test. For rigorous
#' inference on covariate effects, fit \code{\link{CoSparTA}} with the selected
#' covariates and use the estimated \eqn{\gamma} coefficients from the
#' generative model.
#'
#' @examples
#' \dontrun{
#' X <- array(rpois(100 * 20 * 10, lambda = 1.5), dim = c(100, 20, 10))
#' Xcov <- matrix(rnorm(100 * 5), nrow = 100)
#' colnames(Xcov) <- paste0("cov", 1:5)
#'
#' # Style 1: internal unsupervised fit (original behaviour)
#' result <- select_covariates(X, K = 3, covariate_data = Xcov,
#'                              maxiter = 20, convergence_criteria = 'ELBO')
#'
#' # Style 2: supply a fitted object
#' fit <- CoSparTA(X, K = 3)
#' result2 <- select_covariates(X, K = 3, covariate_data = Xcov, fit = fit)
#'
#' # Style 3: supply El directly
#' nf <- normalize_factors(fit)
#' result3 <- select_covariates(X, K = 3, covariate_data = Xcov, El = nf$El)
#'
#' # Which covariates were selected for each rank?
#' result$selected
#'
#' # Full regression summary for rank 1
#' result$summaries[[1]]
#'
#' # Refit with selected covariates for rank 1
#' Xcov_selected <- Xcov[, result$selected[[1]], drop = FALSE]
#' fit <- CoSparTA(X, K = 3, Xcov = Xcov_selected)
#' }
#'
#' @seealso \code{\link{CoSparTA}}, \code{\link{normalize_factors}}
#' @export
select_covariates <- function(X = NULL, K, covariate_data, fit = NULL, El = NULL,
                               alpha = 0.05, verbose = TRUE, ...) {

  start_time <- Sys.time()

  # Resolve El
  if (!is.null(El)) {
    if (ncol(El) != K) {
      stop(sprintf("'El' has %d columns but K = %d.", ncol(El), K))
    }
    fit_unsup <- NULL
  } else if (!is.null(fit)) {
    if (verbose) cat("Extracting El from fit via normalize_factors()...\n")
    El        <- normalize_factors(fit)$El
    fit_unsup <- fit
  } else {
    # Step 1: unsupervised decomposition
    if (verbose) cat("Step 1: fitting unsupervised CoSparTA (K =", K, ")...\n")
    fit_unsup <- CoSparTA(X, K = K, Xcov = NULL, ...)
    El        <- fit_unsup$res$ql$El  # n x K
  }

  cov_names <- colnames(covariate_data)
  if (is.null(cov_names)) cov_names <- paste0("V", seq_len(ncol(covariate_data)))

  selected_list <- vector("list", K)
  summary_list  <- vector("list", K)

  # Step 2: OLS regression of log(El[,k]) on covariate_data for active rows
  for (k in seq_len(K)) {
    if (verbose) cat(sprintf("Step 2: screening covariates for rank %d...\n", k))

    # Filter to active (non-zero) observations; zero loadings are spike/inactive
    active_idx <- which(El[, k] > 0)
    loading_k  <- El[active_idx, k]
    cov_k      <- covariate_data[active_idx, , drop = FALSE]

    lm_fit  <- lm(log(loading_k) ~ ., data = as.data.frame(cov_k))
    lm_summ <- summary(lm_fit)

    # p-values for all coefficients; row 1 is the intercept — skip it
    pvals <- lm_summ$coefficients[-1, 4]

    # Selected: column indices (1-indexed into covariate_data) where p < alpha
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
    fit_unsupervised = fit_unsup,
    runtime_secs     = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
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
#' fit <- CoSparTA(sim$X, K = 3, Xcov = sim$Xcov)
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
#' fit <- CoSparTA(sim$X, K = 3, Xcov = sim$Xcov)
#'
#' # Unsupervised simulation
#' sim0 <- simulate_tensor(n = 100, p = 20, w = 10, K = 2, gamma_true = FALSE)
#' fit0 <- CoSparTA(sim0$X, K = 2)
#' }
#'
#' @seealso \code{\link{CoSparTA}}, \code{\link{match_factors}}
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
