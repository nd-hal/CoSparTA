// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

//' Sparse Softmax (Responsibility) Computation
//'
//' @description
//' Computes the softmax responsibility matrix alpha (nnz x K) at nonzero
//' positions of a sparse 3-way tensor, using the log-sum-exp trick for
//' numerical stability. Replaces the five-line R block:
//' \preformatted{
//'   alpha = Elogl[V1,] + Elogf[V2,] + Elogw[V3,]
//'   exp_offset = rowMaxs(alpha)
//'   alpha = alpha - outer(exp_offset, rep(1,K))
//'   alpha = exp(alpha)
//'   alpha = alpha / rowsums(alpha)
//' }
//'
//' @param V1 Integer vector of user (L-mode) indices (1-based), length nnz.
//' @param V2 Integer vector of time (F-mode) indices (1-based), length nnz.
//' @param V3 Integer vector of channel (W-mode) indices (1-based), length nnz.
//' @param Elogl Numeric matrix (n x K): expected log loadings for L-mode.
//' @param Elogf Numeric matrix (p x K): expected log loadings for F-mode.
//' @param Elogw Numeric matrix (w x K): expected log loadings for W-mode.
//'
//' @return A numeric matrix of size (nnz x K): the row-normalized softmax
//'   responsibilities at each nonzero position.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix calc_qz_sparse_cpp(IntegerVector V1, IntegerVector V2,
                                  IntegerVector V3,
                                  NumericMatrix Elogl,
                                  NumericMatrix Elogf,
                                  NumericMatrix Elogw) {
  int nnz = V1.size();
  int K   = Elogl.ncol();

  NumericMatrix alpha(nnz, K);

  // Step 1: fill log-responsibility matrix
  for (int i = 0; i < nnz; i++) {
    int li = V1[i] - 1;
    int fi = V2[i] - 1;
    int wi = V3[i] - 1;
    for (int k = 0; k < K; k++) {
      alpha(i, k) = Elogl(li, k) + Elogf(fi, k) + Elogw(wi, k);
    }
  }

  // Step 2: log-sum-exp softmax row by row
  for (int i = 0; i < nnz; i++) {
    double row_max = alpha(i, 0);
    for (int k = 1; k < K; k++) {
      if (alpha(i, k) > row_max) row_max = alpha(i, k);
    }
    double row_sum = 0.0;
    for (int k = 0; k < K; k++) {
      alpha(i, k) = std::exp(alpha(i, k) - row_max);
      row_sum += alpha(i, k);
    }
    for (int k = 0; k < K; k++) {
      alpha(i, k) /= row_sum;
    }
  }

  return alpha;
}
