
#' Robust Multivariate Winsorization
#' 
#' Winsorization involves the setting of some threshold for probable outliers 
#' (e.g., top 2 percent of cases or cases more than 4 MAD from the median) and replacing them
#' with values that correspond to that threshold (e.g., replacing a value 6 MAD from the median
#' with one 4 MAD from the median). This function implements a multivariate extension of this
#' approach, following the multivariate covariance estimation and outlier rejection 
#' approach advocated by Leys et al. (\code{doi.org/10.1016/j.jesp.2017.09.011}). That is, this function
#' uses robust Mahalanobis distance multivariate outlier identification. Given that the 
#' Mahalanobis distance from the multivariate distribution's center is directionless \emph{per se}, this function
#' takes each identified outlier and iteratively shifts it by reducing its distance from the center in a 
#' "straight line" in the multivariate space, thereby minimizing warping of the multivariate
#' distribution. Once all identified outliers' Mahalanobis distances from the center have 
#' been shifted and reduced to under the original outlier cutoff, the resulting data frame is returned.
#' In summary, this means that 
#' \strong{multivariate outlier cases will be shifted directly toward the data's center.}
#' 
#' A very small (e.g., <.5) \code{dat_proportion} or a very large (e.g., >.1) 
#' \code{reject_alpha} may lead to unstable or biased behavior. Generally, 
#' it is recommended to keep \code{dat_proportion > .667} and \code{reject_alpha < .05}, 
#' but systematic testing outside these ranges may provide useful insights. 
#' The default is \code{reject_alpha} corresponding to the density in both tails of a Gaussian 
#' cutoff of 3 standard deviations.
#' 
#' Note that this method's influence on multivariate Gaussian distributions' bivariate product-moment
#' correlations has been shown to be minimal and unbiased (e.g., in testing of data
#' simulated from true multivariate Gaussian distributions, the bulk of simulated distributions
#' have correlations change between -.005 and +.005). However, consequences for inherently skewed
#' distributions (e.g., lognormal), multimodal distributions, or for multivariate statistics (e.g., factor analysis) have
#' not been systematically investigated.
#' 
#' See also: Dixon and Tukey (\code{doi.org/10.2307/1266226}) ; 
#' Yale and Forsythe (\code{doi.org/10.1080/00401706.1976.10489449})
#'
#' @param dat A data frame
#' @param dat_proportion The proportion of data used for robust covariance estimation. Inversely related to number of shifted outliers.
#' @param reject_alpha The chi-squared alpha threshold for identifying outliers (a higher alpha leads to more shifted outliers).
#' @param returnOrigDat Logical. TRUE returns original data as an attribute of transformed data. FALSE does not and therefore avoids duplicating objects and using memory unnecessarily.
#'
#' @return A data frame corresponding to \code{dat} with outliers shifted toward multivariate center, with the attributes:\describe{
#'    \item{\code{shiftedCases}}{Logical vector indicating whether each index is [was] an outlier.}
#'    \item{\code{robCov}}{Robust estimate of covariance; output of \code{\link[MASS]{cov.mcd}}.}
#'    \item{\code{dat_proportion}}{Numeric scalar \code{dat_proportion}.}
#'    \item{\code{reject_alpha}}{Numeric scalar \code{reject_alpha}.}
#'    \item{\code{origDat}}{\emph{optional} \code{dat}.}
#' }
#'
#' @export
#'
#' @examples
#' ##-##-##-##-##-##-##-##-##
#' ### Simple example with the iris dataset ###
#' iris_wins <- RMW(iris)
#' mean(attr(iris_wins,'shiftedCases')) # fewer than 5 percent of cases Winsorized
#' cor(cbind(iris[,1:4],iris_wins[,1:4])) # correlation with original data >.99
#' 
#' ##-##-##-##-##-##-##-##-##
#' ### visualize changes to a simulated bivariate Gaussian: ###
#'
#' set.seed(11)
#'
#' # generate data:
#' d <- data.frame(x = rnorm(100)) 
#' d$y = scale(rnorm(100,sd=2) + d$x)
#' 
#' # shift outliers:
#' dFixed <- RMW(d,reject_alpha = .01)
#' 
#' # plot:
#' 
#' d$color <- dFixed$color <- 'black'
#'
#' d[attr(dFixed,'shiftedCases'),'color'] <- 
#'   dFixed[attr(dFixed,'shiftedCases'),'color'] <- 
#'   rainbow(sum(attr(dFixed,'shiftedCases')),start = .1)
#' 
#' plot(d$x,d$y,col=d$color
#'      ,sub = paste0('new pearson: '
#'                    ,signif(cor(dFixed$x,dFixed$y,method='pearson'),3)
#'                    ,'; original pearson: '
#'                    ,signif(cor(d$x,d$y,method='pearson'),3))
#'      ,xlab = '',ylab='')
#' points(dFixed$x, dFixed$y, col = dFixed$color,pch = 19)
#' points(attr(dFixed , 'robCov')$center['x']
#'        ,attr(dFixed , 'robCov')$center['y']
#'        ,pch = 9, col='red') # visualize multivariate center
#' rm(.Random.seed, envir = globalenv())
#' 
#' ##-##-##-##-##-##-##-##-##-##-##-##
#' 
RMW <- function(dat
                , dat_proportion = .8
                , reject_alpha = pnorm(-3)*2
                , returnOrigDat = T){
  
  library(MASS)
  
  # to do:
  # 
  # - use an automatic (e.g., dat_proportion = "check" ) ?
  
  isNonNumeric <- 
    ( as.numeric(sapply(dat,class) =='character') +
        as.numeric(sapply(dat,class) =='factor') +
        as.numeric(apply(dat
                         ,2
                         ,function(x){length(unique(x))<4}
        ))
    )>0
  
  if(any(isNonNumeric)){
    message('Only numeric variables with more than 3 unique values are included in this outlier rejection')
    if(sum(!isNonNumeric)<2){stop('not enough numeric variables with at least 3 unique values')}
    dat_nonNum <- data.frame(dat[,isNonNumeric])
    colnames(dat_nonNum) <- colnames(dat)[isNonNumeric]
    dat <- dat[,!isNonNumeric]
  }
  
  outDat <- dat
  
  ## not currently used "automatic proportion"
  if(is.numeric(dat_proportion)){
    outlList <- RobMulOutliers(dat,dat_proportion,reject_alpha)
  }else{
    prop05 <- optimize(function(curProp){(.05 - mean(
      RobMulOutliers(dat=dat,curProp,reject_alpha=reject_alpha)$excluded
    ))^2}
    ,interval = c(.5,(nrow(dat)-1)/nrow(dat)))
    cat(prop05$minimum)
    outlList <- RobMulOutliers(dat,prop05$minimum,reject_alpha)
  }
  ### ^
  
  ## Winsorize:
  for(curExcl in which(outlList$excluded)){
    
    optimVals <- optimise(function(propShift
                                   ,outlList = outlList
                                   ,dat=dat
                                   ,curExcl=curExcl
    ){
      rowDat <- (dat[curExcl,] - outlList$robCov$center)*propShift + outlList$robCov$center
      
      abs(mahalanobis(rowDat,outlList$robCov$center,cov = outlList$robCov$cov) - 
            outlList$cutoff )
    }
    ,interval = c(.25,.999)
    ,outlList = outlList
    ,dat=dat
    ,curExcl=curExcl
    )
    
    outDat[curExcl,] <- (dat[curExcl,] - outlList$robCov$center)*optimVals$minimum + outlList$robCov$center
  }
  
  if(any(isNonNumeric)){
    dat <- cbind(dat,dat_nonNum)
    outDat <- cbind(outDat,dat_nonNum)
  }
  
  if(returnOrigDat){
    attr(outDat,'origDat') <- dat
  }
  
  attr(outDat,'robCov') <- outlList$robCov
  attr(outDat,'shiftedCases') <- outlList$excluded
  attr(outDat,'reject_alpha') <- reject_alpha
  attr(outDat,'dat_proportion') <- dat_proportion
  
  return(outDat)
  
  if(F){ # for testing
    
    par(mfrow = c(2,2)) ; for(curPlot in 1:4){
      
      d <- data.frame(
        x = rnorm(100)
      ) ; d$y = scale(rnorm(100,sd=2) + d$x)
      # d[1:2,] <- data.frame(x=c(-2,2.57),y=c(2,2.75))
      
      dFixed <- RMW(d)
      
      d$color <- dFixed$color <- 'black'
      
      d[attr(dFixed,'shiftedCases'),'color'] <- 
        dFixed[attr(dFixed,'shiftedCases'),'color'] <- 
        rainbow(sum(attr(dFixed,'shiftedCases')),start = .1)
      
      plot(d$x,d$y,col=d$color
           ,sub = paste0('new pearson: '
                         ,signif(cor(dFixed$x,dFixed$y,method='pearson'),3)
                         ,'; original pearson: '
                         ,signif(cor(d$x,d$y,method='pearson'),3))
           ,xlab = '',ylab='')
      points(dFixed$x, dFixed$y, col = dFixed$color,pch = 19)
      points(attr(dFixed , 'robCov')$center['x']
             ,attr(dFixed , 'robCov')$center['y']
             ,pch = 9, col='red')
      
    } ; par(mfrow = c(1,1))
    
    deltaPearson <- replicate(200,{
      d <- data.frame(
        x = rnorm(100)
      ) ; d$y = scale(rnorm(100,sd=2) + d$x)
      
      dFixed <- RMW(d)
      
      cor(dFixed$x,dFixed$y,method='pearson') -
        cor(d$x,d$y,method='pearson')
    })
    hist(deltaPearson,breaks = 25)
    
    iris1 <- RMW(iris)
    
  }
  
}
