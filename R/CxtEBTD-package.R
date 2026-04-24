#' @docType package
#' @name CxtEBTD-package
#' @import data.table
#' @importFrom dplyr %>% group_by summarize arrange left_join
#' @importFrom matrixStats rowMaxs
#' @importFrom methods new
#' @importFrom Rcpp sourceCpp
#' @importFrom stats ave cor dnbinom lm median nlm optim plogis pnorm
#'   qgamma qnorm quantile rbinom rgamma rpois runif sd
#' @importFrom utils globalVariables
#' @importFrom slam as.simple_sparse_array
#' @importFrom Matrix sparseMatrix
#' @useDynLib CxtEBTD, .registration = TRUE
#' @keywords internal
NULL

# Suppress R CMD check NOTEs for data.table and ggplot2 non-standard evaluation
utils::globalVariables(c(
  ".", "V1", "V2", "V3", "v", "i", "new", "prob_x",
  "r", "s", "count",
  "time", "value", "rank", "within_idx", "channel_group", "channel_name"
))
