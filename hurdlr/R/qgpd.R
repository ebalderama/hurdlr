#' Quantile function of GPD
#' 
#' @description \code{PE} is used to calculate the likelihood of a user-defined
#' 'extreme' value count observation in a double-hurdle regression model.
#' 
#' @param p vector of zero-count probabilities.
#' 
#' @param q vector of 'extreme' count probabilities.
#' 
#' @param log logical operator. If \code{TRUE}, probabilities \code{p} 
#' and \code{q} are given as \code{log(p)}, \code{log(q)}.
#' 
#' @details  
#' 
#' @return 
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Dr. Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @seealso \code{\link{hurdle}} 
#' 

#Quantile function of GPD
qgpd <- function(p, mu = 0, sigma = 1, xi = 1, lower.tail = T){
  options(warn = -1)
  if(lower.tail == F){p <- 1 - p}
  if(xi != 0){inv.cdf <- sigma*((1-p)^(-xi) - 1)/xi + mu}
    else{inv.cdf <- -sigma*log(1-p) + mu}
    if(xi < 0 & inv.cdf > mu - sigma/xi){inv.cdf <- 1}
  return(inv.cdf)
}