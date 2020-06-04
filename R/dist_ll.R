#________________________________________________
#Documentation


#' Distributional Likelihood for Hurdle Model Count Data Regression
#' 
#' @description \code{dist_ll} is the data likelihood fuction for hurdle model 
#' regression using \code{\link{hurdle}}.
#' 
#' @param y numeric response vector.
#' 
#' @param hurd numeric threshold for 'extreme' observations of two-hurdle 
#' models. \code{Inf} for one-hurdle models.
#' 
#' @param lam current value for the poisson likelihood lambda parameter.
#' 
#' @param size size parameter for negative binomial likelihood distributions.
#' 
#' @param mu current value for the negative binomial or log normal likelihood 
#' mu parameter.
#' 
#' @param xi current value for the generalized pareto likelihood xi parameter.
#' 
#' @param sigma current value for the generalized pareto likelihood sigma 
#' parameter.
#' 
#' @param dist character specification of response distribution.
#' 
#' @param log logical operator. if \code{TRUE}, probabilities p are given 
#' as log(p).
#' 
#' @param g.x logical operator. \code{TRUE} if operating within the third 
#' component of the likelihood function (the likelihood of 'extreme' observations).
#' 
#' @details Currently, Poisson, Negative Binomial, log-Normal, 
#' and Generalized Pareto distributions are available.
#' 
#' @return The log-likelihood of the zero-inflated Poisson fit for the
#' current iteration of the MCMC algorithm.
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @seealso \code{\link{hurdle}}
#' @import stats
#' 


#________________________________________________
#Source code

#Hurdle model data distribution likelihood
dist_ll <- function(y, hurd = Inf, lam = NULL, size = 1, mu = NULL, xi = NULL, sigma = NULL,
                    dist = c("poisson", "nb", "lognormal", "gpd"),
                    g.x = F, log = T){

  trunc1 <- ifelse(g.x, hurd - 1, 0)
  trunc2 <- ifelse(g.x, Inf, hurd - 1)

  f <- match.arg(dist)
  ll <- switch(f,
               poisson = dpois(y, lam, log = T) - log(ppois(trunc2, lam) - ppois(trunc1, lam)),
               nb = dnbinom(y, size = size, mu = mu, log = T) - log(pnbinom(trunc2, size = size, mu = mu) - pnbinom(trunc1, size = size, mu = mu)),
               lognormal = mlnorm(y, meanlog = mu, sdlog = sqrt(mu+mu^2/size), log = T) - log(plnorm(trunc2 + 0.5, meanlog = mu, sdlog = sqrt(mu+mu^2/size)) - plnorm(trunc1 + 0.5, meanlog = mu, sdlog = sqrt(mu+mu^2/size))),
               gpd = mgpd(y, trunc1 + 1, sigma, xi, log = T) - log(pgpd(trunc2 + 0.5, trunc1 + 1, sigma, xi)))

  if(log){return(ll)}
  else{return(exp(ll))}
}

