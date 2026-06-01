#' Covariate-aware Empirical Bayes Tensor Decomposition
#'
#' @description
#' Fits a non-negative CP tensor decomposition using empirical Bayes priors with
#' spike-and-slab structure and Poisson likelihood. When covariates are supplied
#' via \code{Xcov}, they enter the generative model directly through
#' covariate-dependent Poisson rates of the form
#' \eqn{\lambda_i = \beta \cdot \exp(X_i^\top \gamma)}, enabling simultaneous
#' decomposition and covariate effect estimation. Without covariates
#' (\code{Xcov = NULL}), the method reduces to unsupervised empirical Bayes
#' Poisson tensor decomposition.
#'
#' @param X A 3-dimensional non-negative integer array of dimensions
#'   \code{n x p x w}, representing observations (e.g., users), time points,
#'   and channels respectively.
#' @param K Integer. Number of components (CP rank) for the decomposition.
#' @param Xcov Covariate input for the observation mode. Can be: (1) \code{NULL}
#'   for fully unsupervised decomposition; (2) a numeric matrix of dimension
#'   \code{n x q}, applied identically to all K components; or (3) a list of
#'   length \code{K}, where each element is either a numeric covariate matrix
#'   (dimensions \code{n x q_k}, potentially different numbers of covariates per
#'   component) or \code{NULL} for unsupervised components. Default \code{NULL}.
#' @param lib_size Numeric vector of length \code{n} giving per-observation
#'   library sizes (exposure/scaling factors). Default \code{NULL}, which sets
#'   all sizes to 1.
#' @param init Character string specifying the initialization method.
#'   Default \code{'random_gamma'}.
#' @param maxiter Integer. Maximum number of EM iterations. Default \code{100}.
#' @param maxiter_init Integer. Maximum number of iterations for the
#'   initialization step. Default \code{100}.
#' @param tol Numeric. Convergence tolerance. Iterations stop when the change
#'   in the objective is below this threshold. Default \code{1e-6}.
#' @param compute_elbo_final Logical. If \code{TRUE}, computes the ELBO after
#'   the final iteration even when \code{convergence_criteria != 'ELBO'}.
#'   Default \code{FALSE}.
#' @param n_stable Integer. Number of consecutive iterations with factor change
#'   below \code{tol} required before declaring convergence. Only used when
#'   \code{convergence_criteria = 'factor_change'}. Default \code{3}.
#' @param ebpm.fn A single function or list of three functions specifying
#'   the empirical Bayes prior for the L, F, and W modes respectively.
#'   Default uses \code{ebpm_point_gamma_multiplier_covariates} for L
#'   (covariate-aware), \code{ebps_with_uq} for F (smooth), and
#'   \code{ebpm_point_gamma_with_uq} for W (point-gamma), all returning
#'   posterior variance and PIP where applicable.
#' @param fix_L Logical. If \code{TRUE}, the observation-mode loadings are held
#'   fixed at initialization. Default \code{FALSE}.
#' @param fix_F Logical. If \code{TRUE}, the time-mode factors are held fixed.
#'   Default \code{FALSE}.
#' @param fix_W Logical. If \code{TRUE}, the channel-mode weights are held
#'   fixed. Default \code{FALSE}.
#' @param smooth_F Logical. If \code{TRUE}, applies smoothing to time-mode
#'   factors. Default \code{TRUE}.
#' @param printevery Integer. Number of iterations between progress messages
#'   when \code{verbose = TRUE}. Default \code{10}.
#' @param verbose Logical. If \code{TRUE}, prints initialization and iteration
#'   progress. Default \code{TRUE}.
#' @param convergence_criteria Character string specifying the convergence
#'   criterion: \code{'factor_change'} (maximum column-normalized change across
#'   all three factor matrices, default), \code{'mKLabs'} (mean KL divergence),
#'   or \code{'ELBO'} (evidence lower bound). \code{'factor_change'} is
#'   recommended: it avoids the O(npwK) ELBO computation and typically converges
#'   in fewer iterations.
#' @param U1_true Optional matrix of true observation-mode factors, used to
#'   track reconstruction error during simulation studies. Default \code{NULL}.
#' @param U2_true Optional matrix of true time-mode factors for simulation
#'   evaluation. Default \code{NULL}.
#' @param U3_true Optional matrix of true channel-mode factors for simulation
#'   evaluation. Default \code{NULL}.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{elbo}{Final ELBO value computed after the last iteration.}
#'   \item{obj_trace}{Numeric vector of objective values at each iteration.}
#'   \item{res}{List of variational posterior summaries. Key fields:
#'     \code{res$ql$El} — posterior mean loadings (\code{n x K});
#'     \code{res$ql$Elogl} — posterior log-mean loadings (\code{n x K});
#'     \code{res$ql$Varl} — posterior variance of loadings (\code{n x K});
#'     \code{res$ql$PIPl} — posterior inclusion probabilities for loadings
#'     (\code{n x K}), giving \eqn{P(\theta_i \neq 0 \mid x_i)} for each
#'     observation and component;
#'     \code{res$qf$Ef} — posterior mean time factors (\code{p x K});
#'     \code{res$qf$Elogf} — posterior log-mean time factors (\code{p x K});
#'     \code{res$qw$Ew} — posterior mean channel weights (\code{w x K});
#'     \code{res$qw$Elogw} — posterior log-mean channel weights (\code{w x K}).}
#'   \item{diff_U}{List of three vectors recording per-iteration reconstruction
#'     error relative to \code{U1_true}, \code{U2_true}, \code{U3_true}.
#'     Only meaningful when true factors are supplied.}
#'   \item{run_time}{Elapsed computation time as a \code{difftime} object.}
#' }
#' In addition, the following normalized and reordered fields are populated
#' automatically (columns scaled to unit Frobenius norm, ordered by descending
#' \eqn{\lambda}):
#' \describe{
#'   \item{res$lambda_normed}{Numeric vector of length \code{K}: component
#'     weights \eqn{\lambda_k = \|l_k\|\|f_k\|\|w_k\|}, sorted descending.}
#'   \item{res$gl_normed}{List of K gamma estimates reordered by descending
#'     \eqn{\lambda}.}
#'   \item{res$ql$El_normed}{Unit-norm observation loading matrix
#'     (\code{n x K}), reordered by descending \eqn{\lambda}.}
#'   \item{res$ql$PIPl_normed}{PIP matrix for L mode, reordered.}
#'   \item{res$ql$shape_post_l_normed}{Posterior gamma shape for L, reordered.}
#'   \item{res$ql$rate_post_l_normed}{Posterior gamma rate for L scaled by
#'     \eqn{\|l_k\|}, reordered, so that \eqn{El\_normed[i,k] \sim
#'     \text{Gamma}(\text{shape},\, \text{rate\_normed})} marginally.}
#'   \item{res$qf$Ef_normed, res$qf$PIPf_normed, res$qf$shape_post_f_normed,
#'     res$qf$rate_post_f_normed}{Same pattern for the F mode.}
#'   \item{res$qw$Ew_normed, res$qw$PIPw_normed, res$qw$shape_post_w_normed,
#'     res$qw$rate_post_w_normed}{Same pattern for the W mode.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- array(rpois(20 * 12 * 4, lambda = 1.5), dim = c(20, 12, 4))
#' Xcov <- matrix(rnorm(20 * 2), nrow = 20, ncol = 2)
#'
#' # Supervised decomposition with covariates
#' fit <- CoSparTA(X, K = 3, Xcov = Xcov, maxiter = 50, verbose = FALSE)
#' EL <- fit$res$ql$El   # 20 x 3 loading matrix
#' EF <- fit$res$qf$Ef   # 12 x 3 factor matrix
#' EW <- fit$res$qw$Ew   # 4 x 3 weight matrix
#'
#' # Unsupervised decomposition (no covariates)
#' fit0 <- CoSparTA(X, K = 3, Xcov = NULL, maxiter = 50, verbose = FALSE)
#' }
#'
#' @export

CoSparTA = function(X,K,Xcov=NULL,
                                    lib_size = NULL,
                                    init = 'random_gamma',
                                    maxiter=100,
                                    maxiter_init = 100,
                                    tol=1e-6,
                                    compute_elbo_final = FALSE,
                                    #ebpm.fn=c(ebpm_point_gamma_multiplier_covariates,smashrgen::ebps,ebpm::ebpm_point_gamma),
                                    ebpm.fn=c(ebpm_point_gamma_multiplier_covariates,ebps_with_uq,ebpm_point_gamma_with_uq), # for var and pip
                                    fix_L = FALSE, fix_F = FALSE, fix_W = FALSE,
                                    smooth_F = T,
                                    #smooth_control=list(),
                                    printevery=10,
                                    verbose=TRUE,
                                    convergence_criteria = 'factor_change',
                                    n_stable = 3L,
                                    U1_true=NULL, U2_true=NULL, U3_true=NULL){

  if (convergence_criteria == 'factor_change' && tol <= 1e-8) {
    tol <- 1e-6
  }

  start_time = Sys.time()

  n_original = dim(X)[1]
  p_original = dim(X)[2]
  w_original = dim(X)[3]

  # X is n by p (time) by w (channel)

  # reduce dimension by droping 0s
  users_zero <- apply(X, 1, sum) == 0
  channels_zero <- apply(X, 3, sum) == 0
  times_zero <- apply(X, 2, sum) == 0
  X <- X[!users_zero,,]
  X <- X[,,!channels_zero]
  X <- X[,!times_zero,]

  n = dim(X)[1]
  p = dim(X)[2]
  w = dim(X)[3]
  n_points = n*p * w
  # Xcov
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

  if(is.null(lib_size)){
    lib_size = rep(1,n)
  }

  x = rbind_sparse_matrix(X,reindex=T)
  non0_idx = cbind(x$V1,x$V2,x$V3)

  if(length(ebpm.fn)==1){
    ebpm.fn.l = ebpm.fn
    ebpm.fn.f = ebpm.fn
    ebpm.fn.w = ebpm.fn
  }
  if(length(ebpm.fn)==3){
    ebpm.fn.l = ebpm.fn[[1]]
    ebpm.fn.f = ebpm.fn[[2]]
    ebpm.fn.w = ebpm.fn[[3]]
  }

  # Normalize ebpm.fn.l to a length-K list
  if (is.function(ebpm.fn.l)) {
    ebpm_fn_l_list <- rep(list(ebpm.fn.l), K)
  } else if (is.list(ebpm.fn.l) && length(ebpm.fn.l) == K) {
    ebpm_fn_l_list <- ebpm.fn.l
  } else {
    stop("ebpm.fn L-mode entry must be a single function or a list of K functions")
  }

  if(verbose){
    cat('initializing loadings and factors...')
    cat('\n')
  }

  res = ebpmf_identity_init(X,K,init,maxiter_init,lib_size)
  alpha = res$ql$Elogl[x$V1,,drop=F] + res$qf$Elogf[x$V2,,drop=F] + res$qw$Elogw[x$V3,,drop=F]
  exp_offset = matrixStats::rowMaxs(alpha)
  alpha = alpha - outer(exp_offset,rep(1,K),FUN='*')
  alpha = exp(alpha)
  alpha = alpha/Rfast::rowsums(alpha)

  obj = c()
  obj[1] = -Inf

  # list of difference in estimation and real value
  diff_U <- list(matrix(1000, nrow = maxiter, ncol=K),
                 matrix(1000, nrow = maxiter, ncol=K),
                 matrix(1000, nrow = maxiter, ncol=K))
  ret_EL_list = list()
  ret_EF_list = list()
  ret_EW_list = list()
  if(verbose){
    cat('running iterations')
    cat('\n')
  }

  El_prev <- NULL
  Ef_prev <- NULL
  Ew_prev <- NULL
  stable_count <- 0L

  for(iter in 1:maxiter){ # this is the update algo
    El_prev <- res$ql$El
    Ef_prev <- res$qf$Ef
    Ew_prev <- res$qw$Ew
    for(k in 1:K) {
      Ez = calc_EZ_3d_fast(x, alpha[,k], n, p, w)
      xcov_k <- if (!is.null(Xcov)) Xcov[[k]] else NULL
      fn_l_k <- ebpm_fn_l_list[[k]]
      if (is.null(xcov_k) && identical(fn_l_k, ebpm_point_gamma_multiplier_covariates)) {
        fn_l_k <- ebpm::ebpm_point_gamma
      }
      res = stm_update_rank1(Ez$rs, Ez$cs, Ez$zs, k, fn_l_k, ebpm.fn.f, ebpm.fn.w, res, fix_L, fix_F, fix_W, xcov_k)
    }

    # Update Z
    alpha = res$ql$Elogl[x$V1,,drop=F] + res$qf$Elogf[x$V2,,drop=F] + res$qw$Elogw[x$V3,,drop=F]
    exp_offset = matrixStats::rowMaxs(alpha)
    alpha = alpha - outer(exp_offset,rep(1,K),FUN='*')
    alpha = exp(alpha)
    alpha = alpha/Rfast::rowsums(alpha)

    if(convergence_criteria == 'factor_change'){
      norm_col <- function(M) apply(M, 2, function(x) x / sqrt(sum(x^2)))
      max_change <- max(
        max(abs(norm_col(res$ql$El) - norm_col(El_prev)), na.rm = TRUE),
        max(abs(norm_col(res$qf$Ef) - norm_col(Ef_prev)), na.rm = TRUE),
        max(abs(norm_col(res$qw$Ew) - norm_col(Ew_prev)), na.rm = TRUE)
      )
      obj[iter] <- max_change
      if (max_change < tol) {
        stable_count <- stable_count + 1L

        if (stable_count >= n_stable) break
      } else {
        stable_count <- 0L
      }
      if(verbose && iter%%printevery==0){
        cat(sprintf('At iter %d, factor_change = %e', iter, max_change))
        cat('\n')
      }
    }

    if (convergence_criteria == "recon_change") {
      G1_prev <- crossprod(El_prev)
      G2_prev <- crossprod(Ef_prev)
      G3_prev <- crossprod(Ew_prev)
      denom_sq <- sum(G1_prev * G2_prev * G3_prev)
      if (denom_sq < 1e-20) denom_sq <- 1
      max_change <- max(sapply(1:K, function(k) {
        norm_sq_curr <- sum(res$ql$El[,k]^2) * sum(res$qf$Ef[,k]^2) * sum(res$qw$Ew[,k]^2)
        norm_sq_prev <- sum(El_prev[,k]^2)   * sum(Ef_prev[,k]^2)   * sum(Ew_prev[,k]^2)
        if (norm_sq_prev < 1e-20) return(0)
        dot_product  <- sum(res$ql$El[,k] * El_prev[,k]) *
                        sum(res$qf$Ef[,k] * Ef_prev[,k]) *
                        sum(res$qw$Ew[,k] * Ew_prev[,k])
        diff_sq <- norm_sq_curr - 2 * dot_product + norm_sq_prev
        sqrt(max(diff_sq, 0)) / sqrt(denom_sq)
      }))
      obj[iter] <- max_change
      if (max_change < tol) {
        stable_count <- stable_count + 1L
        if (stable_count >= n_stable) break
      } else {
        stable_count <- 0L
      }
      if (verbose && iter %% printevery == 0) {
        cat(sprintf("At iter %d, recon_change = %e", iter, max_change))
        cat("\n")
      }
    }

    if(convergence_criteria == 'mKLabs'){
      obj[iter+1] = mKL(x$x,(tcrossprod(res$ql$El,res$qf$Ef)*res$lib_size)[non0_idx])


      if(verbose){
        if(iter%%printevery==0){
          cat(sprintf('At iter %d, mKL(X,LF) = %f',iter,obj[iter+1]))
          cat('\n')
        }
      }
      if(abs(obj[iter+1]-obj[iter])<=tol){
        break
      }
    }

    # for now I only update this; maybe revisit others later
    if(convergence_criteria=='ELBO'){
      obj[iter+1] = calc_stm_obj(x,n,p,w,K,res,non0_idx)
      if(is.infinite(obj[iter+1]) & obj[iter+1] < 0){
        res = res_prev
        break
      }
      if(verbose && iter%%printevery==0){
        cat(sprintf('At iter %d, ELBO: %f',iter,obj[iter+1]))
        cat('\n')
      }
      if((obj[iter+1]-obj[iter])/n_points<tol){ # we are maximizing elbo so not abs
        break
      }
    }

    # scale back to original dimensions; plug 0s back in
    ret_EL <- matrix(0, nrow = n_original, ncol = K)
    ret_EF <- matrix(0, nrow = p_original, ncol = K)
    ret_EW <- matrix(0, nrow = w_original, ncol = K)
    ret_EL[!users_zero,] <- res$ql$El
    ret_EF[!times_zero,] <- res$qf$Ef
    ret_EW[!channels_zero,] <- res$qw$Ew
    ret_EL_list[[iter]] = ret_EL
    ret_EF_list[[iter]] = ret_EF
    ret_EW_list[[iter]] = ret_EW
    # normalize true values
    # then calculate difference, store min diff
    if (!is.null(U1_true) && !is.null(U2_true) && !is.null(U3_true)) {
    U1 <- U1_true[,1]/norm(matrix(U1_true[,1]), type = "F")
    U1_hat_norm <- which_rank(U1, ret_EL)
    U2 <- U2_true[,1]/norm(matrix(U2_true[,1]), type = "F")
    U2_hat_norm <- which_rank(U2, ret_EF)
    U3 <- U3_true[,1]/norm(matrix(U3_true[,1]), type = "F")
    U3_hat_norm <- which_rank(U3, ret_EW)
    diff_U[[1]] <- c(diff_U[[1]], U1_hat_norm$diff)
    diff_U[[2]] <- c(diff_U[[2]], U2_hat_norm$diff)
    diff_U[[3]] <- c(diff_U[[3]], U3_hat_norm$diff)
    if(iter%%printevery==0){
      print(sprintf('At iter %d, U1: %f',iter,U1_hat_norm$diff))
      print(sprintf('At iter %d, U2: %f',iter,U2_hat_norm$diff))
      print(sprintf('At iter %d, U3: %f',iter,U3_hat_norm$diff))
    }
    } # end if (!is.null(U1_true))

    res_prev = res
  }
  ## END FOR LOOP

  if(iter==maxiter){
    message('Reached maximum iterations')
  }

  # calc elbo(approximated)
  if(verbose){
    cat('wrapping-up')
    cat('\n')
  }
  if (convergence_criteria == "ELBO") {
    elbo <- obj[iter + 1]
  } else if (compute_elbo_final) {
    elbo <- calc_stm_obj(x, n, p, w, K, res, non0_idx)
  } else {
    elbo <- NA_real_
  }

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
  res <- .add_normed_fields(res)

  fit = list(elbo=elbo,
             obj_trace=obj,
             res = res,
             diff_U = diff_U,
             run_time = difftime(Sys.time(),start_time,units='auto'))
  return(fit)
}
