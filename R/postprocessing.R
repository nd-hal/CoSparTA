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
