
#' Find robust multivariate outliers
#' 
#' Implements the multivariate covariance estimation and outlier rejection 
#' approach advocated by Leys et al. (2018; \code{doi.org/10.1016/j.jesp.2017.09.011}). 
#' Used by \code{\link{RMW}} to manage outliers through robust multivariate Winsorization.
#'
#' @param dat Data frame of numeric (or coerced to numeric) variables
#' @inheritParams RMW
#'
#' @return A list containing:\describe{
#'    \item{\code{excluded}}{Logical vector indicating whether each index is an outlier.}
#'    \item{\code{robCov}}{Robust estimate of covariance; output of \code{\link[MASS]{cov.mcd}}.}
#'    \item{\code{mDist}}{Numeric vector; robust Mahalanobis distance of each index from the multivariate center.}
#'    \item{\code{cutoff}}{Numeric scalar; given the \code{reject_alpha}, this is the Mahalanobis distance beyond which an index is considered to be an outlier.}
#' }
#'
#' @export
#'
#' @examples
#' iris_noOutliers <- RobMulOutliers(iris[,1:4],.9,.01)
RobMulOutliers <- function(dat, dat_proportion, reject_alpha) {
  require(MASS)
  covar_fit <- cov.mcd(dat, quantile.used = nrow(dat) * 
                         dat_proportion, nsamp = "best")
  mhmcd <- mahalanobis(dat, covar_fit$center, covar_fit$cov)
  cutoff <- (qchisq(p = 1 - reject_alpha, df = ncol(dat)))
  excluded <- rep(F, length(mhmcd))
  excluded[which(mhmcd > cutoff)] <- T
  names(excluded) <- rownames(dat)
  return(list(excluded=excluded, robCov=covar_fit, mDist=mhmcd,cutoff=cutoff))
}