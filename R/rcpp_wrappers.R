#' Fast Sparse EZ Aggregation Wrapper
#'
#' @description
#' Thin R wrapper around \code{\link{calc_EZ_3d_cpp}} that accepts the same
#' sparse coordinate-format data.table \code{x} used by \code{calc_EZ_3d()},
#' making the Rcpp replacement a drop-in at the call site in
#' \code{CoSparTA()} and \code{CoSparTA_missing()}.
#'
#' @param x A data.table with integer columns \code{V1}, \code{V2}, \code{V3}
#'   (1-based indices) and a numeric column \code{x} (observed counts).
#' @param alpha_k Numeric vector of length \code{nrow(x)}: rank-k
#'   responsibilities.
#' @param n Integer. L-mode (user) dimension.
#' @param p Integer. F-mode (time) dimension.
#' @param w Integer. W-mode (channel) dimension.
#'
#' @return Named list with \code{rs} (length n), \code{cs} (length p),
#'   \code{zs} (length w).
#'
#' @keywords internal
calc_EZ_3d_fast <- function(x, alpha_k, n, p, w) {
  calc_EZ_3d_cpp(
    V1      = as.integer(x$V1),
    V2      = as.integer(x$V2),
    V3      = as.integer(x$V3),
    x_vals  = as.numeric(x$v),
    alpha_k = alpha_k,
    n = n, p = p, w = w
  )
}
