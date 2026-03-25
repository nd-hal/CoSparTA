#'@title Make Init L and F on similar scale
#'@param L N by K loading matrix
#'@param FF p by K factor matrix
#'@export
adjLF = function(L,FF){
  gammaL = colSums(L)
  gammaF = colSums(FF)
  adjScale = sqrt(gammaL*gammaF)
  L = t(t(L) * (adjScale/gammaL))
  FF = t(t(FF) * (adjScale/gammaF))
  return(list(L_init = L, F_init=FF))
}

#'@title Log transformation of scRNA-seq count matrix for EBMF
#'@export
log_for_ebmf = function(Y){
  log(1+median(rowSums(Y))/0.5*Y/rowSums(Y))
}


#'@title calculate mean KL divergence of 2 nonnegative matrices
#'@export
#'

mKL = function(A,B){
  D = A*log(A/B)-A+B
  mean(as.matrix(D),na.rm=T)
}


#'@title standard loadings and factors from Poisson matrix Factorization
#'@param L: n by k matrix
#'@param FF: k by p matrix
#'@export
#'

poisson_to_multinom <- function (FF, L) {
  L <- t(t(L) * colSums(FF))
  s <- rowSums(L)
  L <- L / s
  FF <- scale.cols(FF)
  return(list(FF = FF,L = L,s = s))
}

#poisson_to_libsize(F_init,L_init,lib_size)
poisson_to_libsize <- function (FF, L, lib_size) {
  res = poisson_to_multinom(FF,L)
  size = res$s/lib_size
  multinom_to_poisson(res$FF,res$L,size)
}

multinom_to_poisson <- function (FF, L,size) {
  L = L * size
  res = adjLF(L,FF)
  return(list(FF = res$F_init,L = res$L_init))
}

scale.cols <- function (A)
  apply(A,2,function (x) x/sum(x))


calc_EZ = function(x, prob){ # this is what's used in ebpmf_identity
  Ez = sparseMatrix(i = x$i, j = x$j, x = x$x * prob)
  return(list(rs = Matrix::rowSums(Ez), cs = Matrix::colSums(Ez)))
}

softmax3d=function(x){ # updated
  #x = x - array(apply(x,c(1,2),max),dim=dim(x))
  x = exp(x)
  p=as.vector(x)/as.vector(rowSums(x,dims=3)) # from dims=2
  p = pmax(p,1e-10)
  dim(p) <- dim(x)
  return(p)
}

rowMax = function(X){
  do.call(pmax.int, c(na.rm = TRUE, as.data.frame(X)))
}


# to replace Matrix & summary command for tensor
rbind_sparse_matrix <- function(X, reindex = F) {
  # X is input matrix
  convs <- NULL
  row_cnt <- 0

  conv <- as.simple_sparse_array(X)
  #V1: index; ie, U1/L
  #V2: channel; U2/F
  #V3: time; U3/W
  row_cnt <- 0
  ret <- data.table(conv$i)
  ret$v <- conv$v
  if(reindex) {
    unique_i <- unique(ret$V1)
    unique_i <- unique_i[order(unique_i)]
    indx_tbl <- data.frame(new = 1:length(unique_i), V1 = unique_i)
    ret <- left_join(ret, indx_tbl, by = "V1")
    ret <- ret[, .(V1 = new, V2, V3, v)]
    # if(i > 1) ret[, V1 := V1 + max(convs$V1)]
  } else {
    if(i > 1) row_cnt <- row_cnt + max(conv[[i-1]]$i[,1])
    ret[, V1 := V1 + row_cnt]
  }
  ret <- replace(ret,is.na(ret),0)
  ret # out: number * 4; V1,V2,V3,i
}

calc_H = function(x,s,loglik,pm,pmlog){ # this calculates entropy
  if(is.null(loglik)){
    H = 0
  }else{
    H = loglik - sum(x*log(s)+x*pmlog-pm*s-lfactorial(x))
  }
  H
}

calc_EZ_3d <- function(x,prob,n,p,w){ #29

  x$prob_x = x$v*prob

  out_x1= x %>% group_by(V1) %>% summarize(new = sum(prob_x))
  out_x2= x %>% group_by(V2) %>% summarize(new = sum(prob_x))
  out_x3= x %>% group_by(V3) %>% summarize(new = sum(prob_x))

  return(list(rs=out_x1$new,
              cs=out_x2$new,
              zs=out_x3$new))
}


# This is for a part calculating ELBO; replace the tcrossprod only part
cr3d <- function(ew,res=parent.frame()$res){ # check
  sum(ew*tcrossprod(res$lib_size*res$ql$El,res$qf$Ef))
}


scale_to_norm1 <- function(df) {
  # Function to calculate the norm of a vector
  norm <- function(x) sqrt(sum(x^2))

  # Apply the scaling to each column
  scaled_df <- apply(df, 2, function(x) x / norm(x))

  return((scaled_df))
}


which_rank <- function(U, Uhat) {
  U_norm <- U/norm(matrix(U), type = "F")
  Uhat_norm <- scale_to_norm1(Uhat)
  Uhat_norm = as.matrix(Uhat_norm)
  diffs <- apply(Uhat_norm, 2, function(u_hat) {
    norm(as.matrix(U_norm - u_hat), "F")
  })
  r = which.min(diffs)

  list(r = r,
       u_hat = Uhat[, r],
       u_hat_norm = Uhat_norm[, r],
       diff = diffs[r])
}
normalize_columns <- function(matrix) {
  t(apply(matrix, 2, function(col) col / norm(matrix(col), type = "F")))
}

which_rank_modified <- function(U, Uhat) {
  U_norm <- scale_to_norm1(U)
  Uhat_norm <- scale_to_norm1(Uhat)
  available_columns <- 1:ncol(U) # this is 1:K

  results <- vector("list", ncol(U))

  for (i in 1:ncol(U)) {
    diffs <- sapply(available_columns, function(j) {
      norm(matrix(U_norm[, i] - Uhat_norm[, j]), type = "F")
    })

    best_match <- which.min(diffs)
    matched_column <- available_columns[best_match]

    results[[i]] <- list(
      u_index = i,
      r = best_match,
      u_hat = Uhat[, matched_column],
      u_hat_norm = Uhat_norm[, matched_column],
      diff = diffs[best_match]
    )

    available_columns <- setdiff(available_columns, matched_column)
  }

  return(results) # a list of length K
}

find_min_index <- function(matrix) {
  apply(matrix, 2, which.min)
}

extract_best_columns <- function(result_list, min_indices) {
  K <- length(min_indices)
  extracted_columns <- matrix(0, nrow = nrow(result_list[[1]]), ncol = K)
  for (k in 1:K) {
    iter_index <- min_indices[k]
    extracted_columns[, k] <- result_list[[iter_index]][, k]
  }
  return(extracted_columns)
}
