#' @docType package
#' @name CxtEBTD-package
#' @import data.table
#' @importFrom dplyr %>% group_by summarize arrange left_join
#' @importFrom matrixStats rowMaxs
#' @importFrom methods new
#' @importFrom Rcpp sourceCpp
#' @importFrom stats cor dnbinom median nlm optim plogis qnorm rgamma
#' @importFrom slam as.simple_sparse_array
#' @importFrom Matrix sparseMatrix
#' @useDynLib CxtEBTD, .registration = TRUE
#' @keywords internal
NULL

# Suppress R CMD check NOTEs for data.table non-standard evaluation variables
utils::globalVariables(c(".", "V1", "V2", "V3", "v", "i", "new", "prob_x"))
