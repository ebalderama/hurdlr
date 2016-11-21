#________________________________________________
#Documentation


#' Zero-inflated Negative Binomial Data Likelihood
#' 
#' @description Data likelihood fuction for zero-inflated negative binomial 
#' model regression using \code{\link{zero_nb}}.
#' 
#' @param y numeric response vector.
#' 
#' @param z vector of binary operators. \code{z == 0} for observations 
#' considered belnging to the negative binomial distribution, \code{z == 1} 
#' for observations considered to be 'extra' zeros.
#' 
#' @param mu current value for the negative binomial likelihood mu parameter.
#' 
#' @param size size parameter for negative binomial distribution.
#' 
#' @param p vector of 'extra' zero-count probabilities.
#' 
#' 
#' @return The log-likelihood of the zero-inflated negative binomial fit 
#' for the current iteration of the MCMC algorithm.
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @seealso \code{\link{zero_nb}}
#' 
#' 

#________________________________________________
#Source code

#Define zero-inflated negative binomial data loglikelihood at z = 1, z = 0
loglik_zinb <- function(y, z, mu, size, p) {
  
  ll <- log(1 - p) + dnbinom(y, mu = mu, size = size, log = T)
  ll[z==1] <- log(p[z==1])
  return(ll)
}