#' Robust Univariate Winsorization
#'
#' Winsorization involves the setting of some threshold for probable outliers 
#' (e.g., top 2 percent of cases or cases more than 4 MAD from the median) and replacing them
#' with values that correspond to that threshold (e.g., replacing a value 6 MAD from the median
#' with one 4 MAD from the median). This function implements a robust version of this
#' approach, following an asymmetric adaptation of the median-MAD approach advocated by Leys et al. 
#' (\code{doi.org/10.1016/j.jesp.2013.03.013}).
#' 
#' With \code{method=='simple'}: Given a numeric vector, find the median, the median deviation of the lower 50% of cases,
#' and the median deviation of the upper 50% of cases, and scale each asymmetric median
#' deviation to be equivalent to SD scale (the default of \code{\link[stats]{mad}}). 
#' Using some rejection alpha (i.e., tail density), identify
#' the number of Gaussian standard deviations \code{qThresh} corresponding to that alpha. 
#' For all values above \code{median + qThresh*upperMAD}, replace them 
#' with \code{median + qThresh*upperMAD}. For all values 
#' below \code{median - qThresh*lowerMAD}, replace them 
#' with \code{median - qThresh*lowerMAD}.
#' 
#' With \code{method=='distr'}: \emph{Not yet implemented}
#'
#' @param x Numeric vector
#' @param reject_alpha The Gaussian alpha threshold for identifying outliers (a higher alpha leads to more shifted outliers). Defaults to an alpha corresponding to the density in both tails of a Gaussian cutoff of 3 standard deviations.
#'
#' @returns Vector the same length as \code{x}, with the outliers of \code{x} Winsorized.
#'
#' @export
#'
#' @examples
#' sepalWidthWins <- RUW(iris$Sepal.Width, reject_alpha = .01)
#' 1-mean(sepalWidthWins == iris$Sepal.Width) # less than 2 percent of values Winsorized
#' 
RUW <- function(x
                ,reject_alpha = pnorm(-3)*2
                ,method='simple'
                ){
  
  if(length(unique(na.omit(x)))<4){warning('There are too few unique values for RUW') ; return(x)}
  
  qThresh <- abs(qnorm(reject_alpha/2)) # divide tail density by 2 for the two tails
  
  if(method =='simple'){
  xMed <- median(x,na.rm = T)
  
  lowMAD <- (median(x,na.rm = T) - median(x[x < median(x,na.rm = T)]) ) * 1/qnorm(3/4)
  highMAD <- (median(x[x > median(x,na.rm=T)]) - median(x,na.rm = T)) * 1/qnorm(3/4)
  
  upper <- xMed + qThresh*highMAD
  lower <- xMed - qThresh*lowMAD
  xPrime <- x
  xPrime[!is.na(xPrime) & xPrime < lower] <- lower
  xPrime[!is.na(xPrime) & xPrime > upper] <- upper
  return(xPrime)
  }
  if(method=='distr'){
    
    # First, split the distribution into top 50% and bottom 50%
    # Second, take the quantiles between 25% and 50% and separately
    # the quantiles between 50% and 75%, of the distribution, and 
    # find the Gaussian mean and SD that correspond to these
    # quantiles (i.e., "left half Gaussian" and "right half Gaussian").
    # Third, identify values beyond the qThresh*SD, and use their 
    # empirical quantile to determine the replacement value on the qThresh*SD scale.
    
    # NOTE: this procedure might get topsy-turvy if 1/nObs (which limits the
    # empirical quantile) is smaller than reject_alpha. For example, if
    # there are 200 observations and a reject_alpha of .05, then we actual
    # expect 10 values beyond the cutoff. 
    # We need to make sure that ranks are always preserved.
    # Actually small nObs [and scaling edge correction] is a problem.
    # Because, for example, empirical quantile of .8 could have a Gaussian
    # quantile of .99, and be shifted to Gaussian quantile of .8
    # Whereas empirical quantile of .79 might have had a Gaussian quantile
    # of .9, and thereby end up being larger than the first [shifted] value.
    # It's not really clear how Boudt managed this, right?
    
  }
  stop('"method" must be "simple" or "distr"')
}