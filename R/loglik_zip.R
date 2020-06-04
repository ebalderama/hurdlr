#________________________________________________
#Documentation


#' Zero-inflated Poisson Data Likelihood
#' 
#' @description Data likelihood fuction for zero-inflated Poisson model 
#' regression using \code{\link{zero_poisson}}.
#' 
#' @param y numeric response vector.
#' 
#' @param z vector of binary operators. \code{z == 0} for observations 
#' considered belnging to the negative binomial distribution, \code{z == 1} 
#' for observations considered to be 'extra' zeros.
#' 
#' @param lam current value for the Poisson likelihood lambda parameter.
#' 
#' @param p vector of 'extra' zero-count probabilities.
#' 
#' 
#' @return The log-likelihood of the zero-inflated Poisson fit for the
#' current iteration of the MCMC algorithm.
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @seealso \code{\link{zero_poisson}}
#' 
#' 

#________________________________________________
#Source code

#Define zero-inflated Poisson data loglikelihood at z = 1, z = 0 
loglik_zip <- function(y, z, lam, p) {
  
  ll <- log(1 - p) + dpois(y, lam, log = T)
  ll[z==1] <- log(p[z==1])
  return(ll)
}