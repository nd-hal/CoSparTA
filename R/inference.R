#' Extract Posterior Inclusion Probabilities from a CxtEBTD fit
#'
#' @description
#' Extracts the posterior inclusion probability (PIP) matrices from a fitted
#' \code{\link{CxtEBTD}} object. PIPs give the probability that each factor
#' element is truly non-zero vs. noise, derived from the spike-and-slab prior.
#' A PIP close to 1 indicates strong evidence for a non-zero loading; a PIP
#' close to 0 indicates the loading is likely noise.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD}}.
#' @param mode Character string specifying which mode to extract:
#'   \code{'L'} for observation loadings, \code{'F'} for time factors,
#'   \code{'W'} for channel weights. Default \code{'L'}.
#' @param threshold Numeric value in \code{[0, 1]}. If supplied, returns a
#'   logical matrix indicating which elements exceed the threshold. If
#'   \code{NULL}, returns the raw PIP matrix. Default \code{NULL}.
#'
#' @return If \code{threshold = NULL}, a numeric matrix of PIPs with the same
#'   dimensions as the corresponding factor matrix (\code{n x K} for L,
#'   \code{p x K} for F, \code{w x K} for W). If \code{threshold} is
#'   supplied, a logical matrix of the same dimensions where \code{TRUE}
#'   indicates PIP > threshold. Returns \code{NULL} with a warning if PIPs
#'   are not available for the requested mode (e.g. when \code{ebps_with_uq}
#'   is used for F, which has no spike component).
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X, K = 3, Xcov = Xcov)
#'
#' # Raw PIP matrix for observation loadings
#' pip_L <- get_pip(fit, mode = 'L')
#'
#' # Which loadings exceed 0.9 PIP threshold
#' sig_L <- get_pip(fit, mode = 'L', threshold = 0.9)
#' }
#'
#' @export
get_pip <- function(fit, mode = 'L', threshold = NULL) {

  pip_mat <- switch(mode,
    'L' = fit$res$ql$PIPl,
    'F' = fit$res$qf$PIPf,
    'W' = fit$res$qw$PIPw,
    stop("mode must be one of 'L', 'F', or 'W'")
  )

  if (is.null(pip_mat) || all(is.na(pip_mat))) {
    warning(paste("PIPs not available for mode", mode,
                  "-- the prior used for this mode has no spike component."))
    return(NULL)
  }

  if (!is.null(threshold)) {
    if (threshold < 0 || threshold > 1) {
      stop("threshold must be between 0 and 1")
    }
    return(pip_mat > threshold)
  }

  return(pip_mat)
}


#' Compute Credible Intervals from a CxtEBTD fit
#'
#' @description
#' Computes approximate credible intervals for factor elements using the
#' posterior mean and variance stored in a fitted \code{\link{CxtEBTD}}
#' object. Intervals are computed under a normal approximation to the
#' posterior, which is reasonable for elements with non-negligible PIP.
#' For elements with low PIP (likely zero), intervals should be interpreted
#' cautiously.
#'
#' @param fit A fitted object returned by \code{\link{CxtEBTD}}.
#' @param mode Character string specifying which mode to extract:
#'   \code{'L'} for observation loadings, \code{'F'} for time factors,
#'   \code{'W'} for channel weights. Default \code{'L'}.
#' @param level Numeric credible level in \code{(0, 1)}. Default \code{0.95}
#'   for 95\% credible intervals.
#'
#' @return A list with three matrices of the same dimensions as the
#'   corresponding factor matrix:
#' \describe{
#'   \item{mean}{Posterior mean matrix.}
#'   \item{lower}{Lower bound of the credible interval.}
#'   \item{upper}{Upper bound of the credible interval.}
#' }
#' Returns \code{NULL} with a warning if posterior variance is not available
#' for the requested mode.
#'
#' @examples
#' \dontrun{
#' fit <- CxtEBTD(X, K = 3, Xcov = Xcov)
#'
#' # 95% credible intervals for observation loadings
#' ci_L <- get_credible_interval(fit, mode = 'L', level = 0.95)
#' ci_L$mean    # posterior means
#' ci_L$lower   # lower bounds
#' ci_L$upper   # upper bounds
#'
#' # 90% credible intervals for channel weights
#' ci_W <- get_credible_interval(fit, mode = 'W', level = 0.90)
#' }
#'
#' @export
get_credible_interval <- function(fit, mode = 'L', level = 0.95) {

  if (level <= 0 || level >= 1) {
    stop("level must be between 0 and 1")
  }

  mean_mat <- switch(mode,
    'L' = fit$res$ql$El,
    'F' = fit$res$qf$Ef,
    'W' = fit$res$qw$Ew,
    stop("mode must be one of 'L', 'F', or 'W'")
  )

  var_mat <- switch(mode,
    'L' = fit$res$ql$Varl,
    'F' = fit$res$qf$Varf,
    'W' = fit$res$qw$Varw
  )

  if (is.null(var_mat) || all(is.na(var_mat))) {
    warning(paste("Posterior variance not available for mode", mode,
                  "-- credible intervals cannot be computed."))
    return(NULL)
  }

  z <- qnorm((1 + level) / 2)
  sd_mat <- sqrt(var_mat)

  return(list(
    mean  = mean_mat,
    lower = pmax(0, mean_mat - z * sd_mat),  # non-negative since loadings >= 0
    upper = mean_mat + z * sd_mat
  ))
}
