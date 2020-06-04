#________________________________________________
#Documentation


#' Control Parameters for Hurdle Model Count Data Regression
#' 
#' @description Various parameters for fitting control of hurdle 
#' model regression using \code{\link{hurdle}}.
#' 
#' @param a shape parameter for gamma prior distributions.
#' 
#' @param b rate parameter for gamma prior distributions.
#' 
#' @param size size parameter for negative binomial likelihood distributions.
#' 
#' @param beta.prior.mean mu parameter for normal prior distributions.
#' 
#' @param beta.prior.sd standard deviation for normal prior distributions.
#' 
#' @param beta.tune Markov-chain tuning for regression coefficient estimation.
#' 
#' @param pars.tune Markov chain tuning for parameter estimation of 'extreme' 
#' observations distribution.
#' 
#' @param lam.start initial value for the poisson likelihood lambda parameter.
#' 
#' @param mu.start initial value for the negative binomial or log normal 
#' likelihood mu parameter.
#' 
#' @param sigma.start initial value for the generalized pareto likelihood
#' sigma parameter.
#' 
#' @param xi.start initial value for the generalized pareto likelihood 
#' xi parameter.
#' 
#' @return A list of all input values.
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @seealso \code{\link{hurdle}}
#' 
#' 

#________________________________________________
#Source code

#Control parameters for hurdle() model-building
hurdle_control <- function(a = 1, b = 1, size = 1,
                           beta.prior.mean = 0, beta.prior.sd = 1000,
                           beta.tune = 1, pars.tune = 0.2,
                           lam.start = 1, mu.start = 1,
                           sigma.start = 1, xi.start = 1){

  output <- list(a = a, b = b, size = size,
                 beta.prior.mean = beta.prior.mean, beta.prior.sd = beta.prior.sd,
                 pars.tune = pars.tune, beta.tune = beta.tune,
                 lam.start = lam.start, mu.start = mu.start,
                 sigma.start = sigma.start, xi.start = xi.start)

  return(output)
}
