#' @keywords internal
ebpmf_identity_init = function(X,
                               K,
                               init,
                               maxiter_init = 50,
                               lib_size, id1=parent.frame()$id1,
                               id2=parent.frame()$id2,id3=parent.frame()$id3,id4=parent.frame()$id4){

  n = dim(X)[1]
  p = dim(X)[2]
  w = dim(X)[3]
  # JC: p and w were switched and dropping last columns were wrong
  if(is.list(init)){
    # inital size need to match reduced X size
    L_init = init[[1]]
    F_init = init[[2]]
    W_init = init[[3]]

    L_init <- L_init[!parent.frame()$users_zero, , drop=F]
    F_init <- F_init[!parent.frame()$times_zero, , drop=F]
    W_init <- W_init[!parent.frame()$channels_zero, , drop=F]

  }else{

    if(init == 'random_gamma'){
      L_init = matrix(rgamma(n*K, shape=100, rate=100), nrow=n, ncol=K)
      F_init = matrix(rgamma(p*K, shape=100, rate=100), nrow=p, ncol=K)
      W_init = matrix(rgamma(w*K, shape=100, rate=100), nrow=w, ncol=K)
    }
  }
  # adjust scale of L and F, mainly for stability.
  # here we need to change; for now maybe ignore?
  #ratio = poisson_to_libsize(F_init,L_init,lib_size)
  #L_init = ratio$L
  #F_init = ratio$FF

  gl = list()
  gf = list()
  gw = list()

  ql = list(El = L_init, Elogl = log(L_init+1e-18))
  qf = list(Ef = F_init, Elogf = log(F_init+1e-18))
  qw = list(Ew = W_init, Elogw = log(W_init+1e-18))


  return(list(ql=ql,
              qf=qf,
              qw=qw,
              gl=gl,
              gf=gf,
              gw=gw,
              Hl = rep(0,K),
              Hf = rep(0,K),
              Hw = rep(0,K),
              lib_size=lib_size))

}

#' @keywords internal
calc_qz = function(n,p,w,K,ql,qf,qw){ # 28

  qz = array(dim = c(n,p,w,K))
  for(k in 1:K){
    ql_k <- ql$Elogl[,k]
    qf_k <- qf$Elogf[,k]
    qw_k <- qw$Elogw[,k]
    qz_k <- outer(outer(ql_k, qf_k, "+"), qw_k, "+")
    qz[,,,k] <- qz_k
  }

  qz = exp(qz)
  qz_sum <- apply(qz, 1:3, sum)
  dim(qz_sum) <- c(n, p, w)
  for(k in 1:K) {
    qz[,,,k] <- qz[,,,k] / qz_sum
  }

  qz = pmax(qz,1e-16) # parallel maxima

  return(qz) # softmax is for posterior mean
}


#' @keywords internal
stm_update_rank1 = function(l_seq, f_seq, w_seq, k, ebpm.fn.l, ebpm.fn.f, ebpm.fn.w,
                             res, fix_L, fix_F, fix_W, Xcov = NULL){

  # update l
  if(!fix_L) {
    resrev = res
    l_scale = sum(resrev$qf$Ef[,k])*sum(resrev$qw$Ew[,k])*resrev$lib_size #s
    if (!is.null(Xcov)) {
      fit = ebpm.fn.l(l_seq, l_scale, Xcov)
    } else {
      fit = ebpm.fn.l(l_seq, l_scale)
    }

    #fit = ebpm.fn.l(l_seq,l_scale,Xcov)
    res$ql$El[,k] = fit$posterior$mean
    res$ql$Elogl[,k] = fit$posterior$mean_log
    if (!is.null(fit$posterior$var)) {
      if (is.null(res$ql$Varl)) res$ql$Varl <- matrix(NA_real_, nrow(res$ql$El), ncol(res$ql$El))
      if (is.null(res$ql$PIPl)) res$ql$PIPl <- matrix(NA_real_, nrow(res$ql$El), ncol(res$ql$El))
      res$ql$Varl[,k] = fit$posterior$var
      res$ql$PIPl[,k] = fit$posterior$pip
    }
    if (!is.null(fit$posterior$shape_post)) {
      if (is.null(res$ql$shape_post_l)) res$ql$shape_post_l <- matrix(NA_real_, nrow(res$ql$El), ncol(res$ql$El))
      if (is.null(res$ql$rate_post_l))  res$ql$rate_post_l  <- matrix(NA_real_, nrow(res$ql$El), ncol(res$ql$El))
      res$ql$shape_post_l[, k] <- fit$posterior$shape_post
      res$ql$rate_post_l[, k]  <- fit$posterior$rate_post
    }
    res$Hl[k] = calc_H(l_seq,l_scale,fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)
    res$gl[[k]] = fit$fitted_g
  }

  # update F
  if(!fix_F){
    f_scale = sum(res$ql$El[,k])*sum(res$qw$Ew[,k])

    fit = ebpm.fn.f(f_seq,f_scale)
    res$qf$Ef[,k] = fit$posterior$mean
    res$qf$Elogf[,k] = fit$posterior$mean_log
    if (!is.null(fit$posterior$var)) {
      if (is.null(res$qf$Varf)) res$qf$Varf <- matrix(NA_real_, nrow(res$qf$Ef), ncol(res$qf$Ef))
      if (is.null(res$qf$PIPf)) res$qf$PIPf <- matrix(NA_real_, nrow(res$qf$Ef), ncol(res$qf$Ef))
      res$qf$Varf[,k] = fit$posterior$var
      res$qf$PIPf[,k] = fit$posterior$pip
    }
    if (!is.null(fit$posterior$shape_post)) {
      if (is.null(res$qf$shape_post_f)) res$qf$shape_post_f <- matrix(NA_real_, nrow(res$qf$Ef), ncol(res$qf$Ef))
      if (is.null(res$qf$rate_post_f))  res$qf$rate_post_f  <- matrix(NA_real_, nrow(res$qf$Ef), ncol(res$qf$Ef))
      res$qf$shape_post_f[, k] <- fit$posterior$shape_post
      res$qf$rate_post_f[, k]  <- fit$posterior$rate_post
    }
    res$Hf[k] = calc_H(f_seq,f_scale,fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)
    res$gf[[k]] = fit$fitted_g
  }

  # Update W
  if(!fix_W){
    w_scale = sum(res$ql$El[,k])*sum(res$qf$Ef[,k])
    fit = ebpm.fn.w(w_seq,w_scale)
    res$qw$Ew[,k] = fit$posterior$mean
    res$qw$Elogw[,k] = fit$posterior$mean_log
    if (!is.null(fit$posterior$var)) {
      if (is.null(res$qw$Varw)) res$qw$Varw <- matrix(NA_real_, nrow(res$qw$Ew), ncol(res$qw$Ew))
      if (is.null(res$qw$PIPw)) res$qw$PIPw <- matrix(NA_real_, nrow(res$qw$Ew), ncol(res$qw$Ew))
      res$qw$Varw[,k] = fit$posterior$var
      res$qw$PIPw[,k] = fit$posterior$pip
    }
    if (!is.null(fit$posterior$shape_post)) {
      if (is.null(res$qw$shape_post_w)) res$qw$shape_post_w <- matrix(NA_real_, nrow(res$qw$Ew), ncol(res$qw$Ew))
      if (is.null(res$qw$rate_post_w))  res$qw$rate_post_w  <- matrix(NA_real_, nrow(res$qw$Ew), ncol(res$qw$Ew))
      res$qw$shape_post_w[, k] <- fit$posterior$shape_post
      res$qw$rate_post_w[, k]  <- fit$posterior$rate_post
    }
    res$Hw[k] = calc_H(w_seq,w_scale,fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)
    if(!is.null(fit$fitted_g)) res$gw[[k]] <- fit$fitted_g
  }

  return(res)

}


#' @keywords internal
calc_approx_elbo_F = function(x,alpha,K,ebpm.fn.f,res,ebps_control){

  for(k in 1:K){
    Ez = calc_EZ(x, alpha[,k])
    fit = ebpm.fn.f(Ez$cs,sum(res$lib_size*res$ql$El[,k]),
                    g_init = list(sigma2 = res$gf$sigma2[k]),
                    q_init = list(m=res$qf$Elogf[,k],smooth = res$qf$Elogf_smooth[,k]),
                    general_control = list(maxiter=1,
                                           maxiter_vga = ebps_control$maxiter_vga,
                                           make_power_of_2=ebps_control$make_power_of_2,
                                           vga_tol=ebps_control$vga_tol,
                                           tol = ebps_control$tol),
                    smooth_control = list(wave_trans='dwt',
                                          filter.number = ebps_control$filter.number,
                                          family = ebps_control$family,
                                          ebnm_params=ebps_control$ebnm_params,
                                          warmstart=ebps_control$warmstart))
    res$Hf[k] = calc_H(Ez$cs,sum(res$lib_size*res$ql$El[,k]),fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)
  }

  res

}


#' Default control parameters for smooth factor estimation
#'
#' @description
#' Returns a list of default control parameters for the wavelet-based smooth
#' factor estimation step in \code{\link{CxtEBTD}}. Pass the output as the
#' \code{smooth_control} argument to override individual settings.
#'
#' @return A named list with fields: \code{wave_trans}, \code{ndwt_method},
#'   \code{filter.number}, \code{family}, \code{ebnm_params},
#'   \code{maxiter}, \code{maxiter_vga}, \code{make_power_of_2},
#'   \code{vga_tol}, \code{tol}, \code{warmstart},
#'   \code{convergence_criteria}, \code{m_init_method_for_init}.
#'
#' @examples
#' ctrl <- ebpmf_identity_smooth_control_default()
#' ctrl$filter.number <- 4  # override wavelet filter
#'
#' @export
ebpmf_identity_smooth_control_default = function(){
  list(wave_trans='ndwt',
       ndwt_method = "ti.thresh",
       filter.number = 1,
       family = 'DaubExPhase',
       ebnm_params=list(),
       maxiter=1,
       maxiter_vga = 10,
       make_power_of_2='extend',
       vga_tol=1e-3,
       tol = 1e-2,
       warmstart=TRUE,
       convergence_criteria = 'nugabs',
       m_init_method_for_init = 'vga')
}

#' @keywords internal
calc_stm_obj = function(x,n,p,w,K,res,non0_idx){
  qz = calc_qz(n,p,w, K,res$ql,res$qf,res$qw) # update; qz is n*p*w*K

  val = 0
  sum_u <- 0
  for(k in 1:K){
    ql_k <- res$ql$El[,k]
    qf_k <- res$qf$Ef[,k]
    qw_k <- res$qw$Ew[,k]
    val = val + qz[,,,k]*( log(ql_k %o% qf_k %o% qw_k) - log(qz[,,,k]))
    sum_u <- sum_u + sum(ql_k %o% qf_k %o% qw_k)
  }

  E1 = sum(x$v*val[non0_idx]) - sum_u - sum(lfactorial(x$v))
  # sum(x$v*val[non0_idx]): Esumsum Z_ijk*log(uuu) - Esumsum Z*log(pi)
  # sum_u: Esumsum (u_il u_jl u_kl)
  # sum(lfactorial(x$v)): Esum log(T_ijk!)

  return(E1+sum(res$Hl)+sum(res$Hf)+sum(res$Hw))
}
