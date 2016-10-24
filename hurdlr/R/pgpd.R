#' Cumulative distribution function of GPD
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
pgpd <- function(q, mu = 0, sigma = 1, xi = 1, lower.tail = T){
  options(warn = -1)
  if(xi != 0) {cdf <- 1 - (1 + xi*(q - mu)/sigma)^(-1/xi)}
  else{cdf <- 1 - exp(-(q - mu)/sigma)}
  cdf[q < mu] <- 0
  if(xi < 0){cdf[q > mu - sigma/xi] <- 1}
    if(lower.tail){return(cdf)}
    else{return(1 - cdf)}
}
