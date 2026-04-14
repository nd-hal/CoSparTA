#' Rescale loadings and factors to similar column norms
#'
#' @description
#' Adjusts the column scales of L and F so that their geometric mean is
#' preserved. Used for numerical stability when K > 1.
#'
#' @param L Numeric matrix of dimensions \code{n x K} (observation loadings).
#' @param FF Numeric matrix of dimensions \code{p x K} (time factors).
#'
#' @return A named list with:
#' \describe{
#'   \item{L_init}{Rescaled loading matrix, same dimensions as \code{L}.}
#'   \item{F_init}{Rescaled factor matrix, same dimensions as \code{FF}.}
#' }
#' @export

adjLF = function(L,FF){
  gammaL = colSums(L)
  gammaF = colSums(FF)
  adjScale = sqrt(gammaL*gammaF)
  L = t(t(L) * (adjScale/gammaL))
  FF = t(t(FF) * (adjScale/gammaF))
  return(list(L_init = L, F_init=FF))
}

#' Log-transformation of count matrix for empirical Bayes matrix factorization
#'
#' @description
#' Applies a library-size normalized log transformation suitable as
#' preprocessing input for EBMF methods.
#'
#' @param Y Non-negative integer count matrix (\code{n x p}).
#'
#' @return A numeric matrix of the same dimensions as \code{Y} with
#'   log-transformed, library-size normalized values.
#' @export
log_for_ebmf = function(Y){
  log(1+median(rowSums(Y))/0.5*Y/rowSums(Y))
}


#' Mean KL divergence between two non-negative matrices
#'
#' @description
#' Computes the mean generalized KL divergence
#' \eqn{D(A \| B) = A \log(A/B) - A + B}, averaged over all elements.
#' Used as a convergence criterion in \code{\link{CxtEBTD}}.
#'
#' @param A Non-negative numeric matrix or vector (observed).
#' @param B Non-negative numeric matrix or vector (fitted), same dimensions as \code{A}.
#'
#' @return A single numeric value giving the mean KL divergence.
#' @export

mKL = function(A,B){
  D = A*log(A/B)-A+B
  mean(as.matrix(D),na.rm=T)
}


#' Standardize Poisson matrix factorization to multinomial parameterization
#'
#' @description
#' Converts a Poisson factorization \eqn{X \approx L F^\top} into a
#' standardized form where columns of F sum to 1 and L absorbs the scale.
#'
#' @param FF Numeric matrix of dimensions \code{p x K} (factors).
#' @param L Numeric matrix of dimensions \code{n x K} (loadings).
#'
#' @return A named list with:
#' \describe{
#'   \item{FF}{Column-normalized factor matrix (\code{p x K}).}
#'   \item{L}{Rescaled loading matrix (\code{n x K}).}
#'   \item{s}{Numeric vector of length \code{n} giving row scales.}
#' }
#' @export

poisson_to_multinom <- function (FF, L) {
  L <- t(t(L) * colSums(FF))
  s <- rowSums(L)
  L <- L / s
  FF <- scale.cols(FF)
  return(list(FF = FF,L = L,s = s))
}

#' @keywords internal
poisson_to_libsize <- function (FF, L, lib_size) {
  res = poisson_to_multinom(FF,L)
  size = res$s/lib_size
  multinom_to_poisson(res$FF,res$L,size)
}

#' @keywords internal
multinom_to_poisson <- function (FF, L,size) {
  L = L * size
  res = adjLF(L,FF)
  return(list(FF = res$F_init,L = res$L_init))
}

#' @keywords internal
scale.cols <- function (A)
  apply(A,2,function (x) x/sum(x))

#' @keywords internal
calc_EZ = function(x, prob){ # this is what's used in ebpmf_identity
  Ez = sparseMatrix(i = x$i, j = x$j, x = x$x * prob)
  return(list(rs = Matrix::rowSums(Ez), cs = Matrix::colSums(Ez)))
}

#' @keywords internal
softmax3d=function(x){ # updated
  #x = x - array(apply(x,c(1,2),max),dim=dim(x))
  x = exp(x)
  p=as.vector(x)/as.vector(rowSums(x,dims=3)) # from dims=2
  p = pmax(p,1e-10)
  dim(p) <- dim(x)
  return(p)
}

#' @keywords internal
rowMax = function(X){
  do.call(pmax.int, c(na.rm = TRUE, as.data.frame(X)))
}


#' @keywords internal
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
    ret <- data.table::as.data.table(dplyr::left_join(ret, indx_tbl, by = "V1"))
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

#' @keywords internal
calc_EZ_3d <- function(x,prob,n,p,w){ #29

  x$prob_x = x$v*prob

  out_x1= x %>% dplyr::group_by(V1) %>% dplyr::summarize(new = sum(prob_x))
  out_x2= x %>% dplyr::group_by(V2) %>% dplyr::summarize(new = sum(prob_x))
  out_x3= x %>% dplyr::group_by(V3) %>% dplyr::summarize(new = sum(prob_x))

  return(list(rs=out_x1$new,
              cs=out_x2$new,
              zs=out_x3$new))
}



#' @keywords internal
cr3d <- function(ew,res=parent.frame()$res){ # check
  sum(ew*tcrossprod(res$lib_size*res$ql$El,res$qf$Ef))
}

#' @keywords internal
scale_to_norm1 <- function(df) {
  # Function to calculate the norm of a vector
  norm <- function(x) sqrt(sum(x^2))

  # Apply the scaling to each column
  scaled_df <- apply(df, 2, function(x) x / norm(x))

  return((scaled_df))
}

#' @keywords internal
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

#' @keywords internal
normalize_columns <- function(matrix) {
  t(apply(matrix, 2, function(col) col / norm(matrix(col), type = "F")))
}

#' @keywords internal
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

#' @keywords internal
find_min_index <- function(matrix) {
  apply(matrix, 2, which.min)
}

#' @keywords internal
extract_best_columns <- function(result_list, min_indices) {
  K <- length(min_indices)
  extracted_columns <- matrix(0, nrow = nrow(result_list[[1]]), ncol = K)
  for (k in 1:K) {
    iter_index <- min_indices[k]
    extracted_columns[, k] <- result_list[[iter_index]][, k]
  }
  return(extracted_columns)
}
