// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

//' Sparse Weighted Aggregation for Tensor EM
//'
//' @description
//' For a sparse 3-way tensor stored in coordinate format, computes the
//' rank-k sufficient statistics (row sums, column sums, slice sums) needed
//' for the EM update of factor matrices L, F, W. Replaces the dplyr
//' group-by approach in \code{calc_EZ_3d()}.
//'
//' @param V1 Integer vector of user (L-mode) indices (1-based).
//' @param V2 Integer vector of time (F-mode) indices (1-based).
//' @param V3 Integer vector of channel (W-mode) indices (1-based).
//' @param x_vals Numeric vector of observed counts at each nonzero position.
//' @param alpha_k Numeric vector of length \code{nnz}: the rank-k responsibility
//'   weights (alpha[, k]).
//' @param n Integer. Number of users (L-mode dimension).
//' @param p Integer. Number of time points (F-mode dimension).
//' @param w Integer. Number of channels (W-mode dimension).
//'
//' @return A named list with three numeric vectors:
//'   \item{rs}{Length-\code{n} vector: row sums (L-mode aggregation).}
//'   \item{cs}{Length-\code{p} vector: column sums (F-mode aggregation).}
//'   \item{zs}{Length-\code{w} vector: slice sums (W-mode aggregation).}
//'
//' @export
// [[Rcpp::export]]
List calc_EZ_3d_cpp(IntegerVector V1, IntegerVector V2, IntegerVector V3,
                    NumericVector x_vals, NumericVector alpha_k,
                    int n, int p, int w) {
  int nnz = V1.size();

  NumericVector rs(n, 0.0);
  NumericVector cs(p, 0.0);
  NumericVector zs(w, 0.0);

  for (int i = 0; i < nnz; i++) {
    double contrib = x_vals[i] * alpha_k[i];
    rs[V1[i] - 1] += contrib;
    cs[V2[i] - 1] += contrib;
    zs[V3[i] - 1] += contrib;
  }

  return List::create(
    Named("rs") = rs,
    Named("cs") = cs,
    Named("zs") = zs
  );
}
