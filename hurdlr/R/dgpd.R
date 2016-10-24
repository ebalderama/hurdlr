#________________________________________________
#Documentation


#' Probability distribution function of GPD
#' 
#' @description \code{zero_nb} is used to fit zero-inflated 
#' negative binomial regression models to count data via Bayesian inference.
#' 
#' @param y numeric response vector.
#' 
#' @param x numeric predictor matrix.
#' 
#' @param a shape parameter for gamma prior distributions.
#' 
#' @param b rate parameter for gamma prior distributions.
#' 
#' @param mu.start initial value for mu parameter.
#' 
#' @param beta.prior.mean mu parameter for normal prior distributions.
#' 
#' @param beta.prior.sd standard deviation for normal prior distributions.
#' 
#' @param iters number of iterations for the Markov chain to run.
#' 
#' @param burn numeric burn-in length.
#' 
#' @param nthin numeric thinning rate.
#' 
#' @param plots logical operator. \code{TRUE} to output plots.
#' 
#' @param progress.bar logical operator. \code{TRUE} to print progress bar.
#' 
#' @details 
#' 
#' @return 
#' 
#' @author 
#' Dr. Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @example 
#' 

dgpd <- function(x, mu = 0, sigma = 1, xi = 1, log = F){
  options(warn = -1)
  log.den <- (log(sigma) - (xi + 1)*log(sigma + xi*(x - mu)))/xi
  log.den[x < mu] <- -Inf
  if(xi < 0){log.den[x > (mu - sigma/xi)] <- -Inf}
    if(log){return(log.den)}
    else{return(exp(log.den))}
}

