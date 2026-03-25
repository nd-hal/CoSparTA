#print("source JC ebptd function")
CxtEBTD = function(X,K,Xcov,
                                    lib_size = NULL,
                                    init = 'fasttopics',
                                    maxiter=100,
                                    maxiter_init = 100,
                                    tol=1e-8,
                                    ebpm.fn=c(ebpm::ebpm_point_gamma,smashrgen::ebps), # need add another function
                                    fix_L = FALSE, fix_F = FALSE, fix_W = FALSE,
                                    smooth_F = T,
                                    #smooth_control=list(),
                                    printevery=10,
                                    verbose=TRUE,
                                    adj_LF_scale = TRUE,
                                    convergence_criteria = 'mKLabs',
                                    U1_true=NULL, U2_true=NULL, U3_true=NULL){

  # remove first/last hours / channels that are all 0, and are at the start or end of the matrices
  start_time = Sys.time()

  n_original = dim(X)[1]
  p_original = dim(X)[2]
  w_original = dim(X)[3]

  #browser()
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
  ######## Xcov

  Xcov <- Xcov[!users_zero,]
  if(is.null(lib_size)){
    lib_size = rep(1,n) # what is this?
  }

  #X = Matrix(X,sparse = T)
  #x = summary(X)
  #non0_idx = cbind(x$i,x$j) # i is first dim, j is second dim, x is value
  # replace above with tensor functions
  x = rbind_sparse_matrix(X,reindex=T)
  non0_idx = cbind(x$V1,x$V2,x$V3)

  #smooth_control = modifyList(ebpmf_identity_smooth_control_default(),smooth_control,keep.null = TRUE)
  if(length(ebpm.fn)==1){
    ebpm.fn.l = ebpm.fn
    ebpm.fn.f = ebpm.fn
    ebpm.fn.w = ebpm.fn
  }
  if(length(ebpm.fn)==3){
    ebpm.fn.l = ebpm.fn[[1]]
    ebpm.fn.f = ebpm.fn[[2]]
    ebpm.fn.w = ebpm.fn[[3]]
    #print(ebpm.fn.w)
  }

  if(verbose){
    cat('initializing loadings and factors...')
    cat('\n')
  }

  res = ebpmf_identity_init(X,K,init,maxiter_init,lib_size)
  # what are we doing here? is alpha a scale?
  alpha = res$ql$Elogl[x$V1,,drop=F] + res$qf$Elogf[x$V2,,drop=F] + res$qw$Elogw[x$V3,,drop=F]
  exp_offset = matrixStats::rowMaxs(alpha)
  alpha = alpha - outer(exp_offset,rep(1,K),FUN='*')
  alpha = exp(alpha)
  alpha = alpha/rowsums(alpha)


  obj = c()
  obj[1] = -Inf

  # list of difference in estimation and real value
  #diff_U <- list(c(), c(), c())
  # record diff_U for each K
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

  # ######################################
  for(iter in 1:maxiter){ # this is the update algo
    #print(iter)
    for(k in 1:K) {
      Ez = calc_EZ_3d(x, alpha[,k],n,p,w)
      res = stm_update_rank1(Ez$rs,Ez$cs,Ez$zs, k,ebpm.fn.l,ebpm.fn.f,ebpm.fn.w, res,fix_L,fix_F,fix_W,Xcov) # here we need to update
    }

    # Update Z
    # EZ = Calc_EZ(X,K,EZ,res$ql,res$qf)
    alpha = res$ql$Elogl[x$V1,,drop=F] + res$qf$Elogf[x$V2,,drop=F] + res$qw$Elogw[x$V3,,drop=F]
    exp_offset = matrixStats::rowMaxs(alpha)
    alpha = alpha - outer(exp_offset,rep(1,K),FUN='*')
    alpha = exp(alpha)
    alpha = alpha/rowsums(alpha)

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
      #browser()
      obj[iter+1] = calc_stm_obj(x,n,p,w,K,res,non0_idx)
      #print(obj)
      print(obj[iter+1])
      if(is.infinite(obj[iter+1]) & obj[iter+1] < 0){
        res = res_prev
        break
      }
      #print(obj[iter+1]-obj[iter])
      if(verbose){
        if(iter%%printevery==0){
          print(sprintf('At iter %d, ELBO: %f',iter,obj[iter+1]))
        }
      }
      # print("diff:")
      # print((obj[iter+1]-obj[iter]))
      if((obj[iter+1]-obj[iter])/n_points<tol){ # we are maximizing elbo so not abs
        #if(abs((obj[iter+1]-obj[iter])/n_points)<tol){ # why not abs?
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
    U1 <- U1_true[,1]/norm(matrix(U1_true[,1]), type = "F")
    U1_hat_norm <- which_rank(U1, ret_EL)
    #
    U2 <- U2_true[,1]/norm(matrix(U2_true[,1]), type = "F")
    U2_hat_norm <- which_rank(U2, ret_EF)
    #
    U3 <- U3_true[,1]/norm(matrix(U3_true[,1]), type = "F")
    U3_hat_norm <- which_rank(U3, ret_EW)

    # U1_hat_norm <- which_rank_modified(U1_true, ret_EL)
    # U2_hat_norm <- which_rank_modified(U2_true, ret_EF)
    # U3_hat_norm <- which_rank_modified(U3_true, ret_EW)
    # for (k in 1:K) {
    #   diff_U[[1]][iter, k] <- U1_hat_norm[[k]]$diff
    #   diff_U[[2]][iter, k] <- U2_hat_norm[[k]]$diff
    #   diff_U[[3]][iter, k] <- U3_hat_norm[[k]]$diff
    # }
    diff_U[[1]] <- c(diff_U[[1]], U1_hat_norm$diff)
    diff_U[[2]] <- c(diff_U[[2]], U2_hat_norm$diff)
    diff_U[[3]] <- c(diff_U[[3]], U3_hat_norm$diff)

    # plot(U2_hat_norm$u_hat_norm, col = "red")
    # points(U2)

    if(iter%%printevery==0){
      print(sprintf('At iter %d, U1: %f',iter,U1_hat_norm$diff))
      print(sprintf('At iter %d, U2: %f',iter,U2_hat_norm$diff))
      print(sprintf('At iter %d, U3: %f',iter,U3_hat_norm$diff))
    }

    if(adj_LF_scale){ # we don't use it now; but likely will use it when K>1; this is like the lambda
      gammaL = colSums(res$ql$El)
      gammaF = colSums(res$qf$Ef)
      gammaF = colSums(res$qf$Ef)
      adjScale = sqrt(gammaL*gammaF)
      sl = adjScale/gammaL
      sf = adjScale/gammaF
      res$ql$El = t(t(res$ql$El) * sl)
      res$ql$Elogl = res$ql$Elogl + outer(rep(1,n),log(sl))
      res$qf$Ef = t(t(res$qf$Ef) * sf)
      res$qf$Ef_smooth = t(t(res$qf$Ef_smooth) * sf)
      res$qf$Elogf = res$qf$Elogf + outer(rep(1,p),log(sf))
      res$qf$Elogf_smooth = res$qf$Elogf_smooth + outer(rep(1,p),log(sf))
    }
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
  elbo = calc_stm_obj(x,n,p,w,K,res,non0_idx)

  ### We want to add this part later
  #ldf = poisson_to_multinom(res$qf$Ef,res$ql$El)
  #EL = ldf$L
  #EF = ldf$FF
  EF_smooth = NULL

  # add: now we finished looping, we find for each U1/U2/U3, on each K, at which iter the diff is smallest
  min_indices <- list()
  # Find the index of the minimum value for each column in each matrix of diff_U
  # for (i in 1:length(diff_U)) {
  #   min_indices[[i]] <- find_min_index(diff_U[[i]])
  # }
  # best_EL <- extract_best_columns(ret_EL_list, min_indices[[1]])
  # best_EF <- extract_best_columns(ret_EF_list, min_indices[[2]])
  # best_EW <- extract_best_columns(ret_EW_list, min_indices[[3]])
  #
  res$ql$El <- ret_EL
  res$qf$Ef <- ret_EF
  res$qw$Ew <- ret_EW
  #res$ql$El <- best_EL
  # res$qf$Ef <- best_EF
  # res$qw$Ew <- best_EW

  fit = list(EL ="check res", # EL = EL,
             EF ="check res", # EF = EF,
             EF_smooth = EF_smooth,
             elbo=elbo,
             #d=ldf$s,
             obj_trace=obj,
             res = res,
             diff_U = diff_U,
             run_time = difftime(Sys.time(),start_time,units='auto'))
  return(fit)
}
