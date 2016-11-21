#________________________________________________
#Documentation


#' MCMC Third-Component Parameter Update Function for Hurdle Model Count Data
#' Regression
#' 
#' @description MCMC algorithm for updating the third-component likelihood 
#' parameters in hurdle model regression using \code{\link{hurdle}}.
#' 
#' @param y numeric response vector of observations within the bounds of the 
#' third component of the likelihood function, \eqn{y[y \ge hurd]}.
#' 
#' @param hurd numeric threshold for 'extreme' observations of two-hurdle models.
#' 
#' @param dist character specification of response distribution for the third 
#' component of the likelihood function.
#' 
#' @param like.part numeric vector of the current third-component likelihood values.
#' 
#' @param a shape parameter for gamma prior distributions.
#' 
#' @param b rate parameter for gamma prior distributions.
#' 
#' @param size size parameter for negative binomial likelihood distributions.
#' 
#' @param lam current value for the poisson likelihood lambda parameter.
#' 
#' @param mu current value for the negative binomial or log normal likelihood 
#' mu parameter.
#' 
#' @param xi current value for the generalized pareto likelihood xi parameter.
#' 
#' @param sigma current value for the generalized pareto likelihood sigma 
#' parameter.
#' 
#' @param lam.acc,mu.acc,xi.acc,sigma.acc current MCMC values for third-component 
#' parameter acceptance rates.
#' 
#' @param lam.tune,mu.tune,xi.tune,sigma.tune current MCMC tuning values for 
#' each third-component parameter. 
#' 
#' @param g.x logical operator. \code{TRUE} if operating within the third component 
#' of the likelihood function (the likelihood of 'extreme' observations).
#' 
#' 
#' @return A list of MCMC-updated likelihood estimator(s) for the third-component
#' parameter(s) and each parameter's MCMC acceptance ratio.  
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @seealso 
#' \code{\link{hurdle}} \cr
#' \code{\link{dist_ll}}
#' 
#' 

#________________________________________________
#Source code

update_pars <- function(y, hurd, dist, like.part,
                        a, b, size, lam, mu, xi, sigma,
                        lam.acc, mu.acc, xi.acc, sigma.acc,
                        lam.tune, mu.tune, xi.tune, sigma.tune,
                        g.x = F){

  like.cur <- like.part

  if(dist == "poisson"){

    #Update Poisson lambda
    lam.new <- exp(rnorm(1, log(lam), lam.tune))
    like.new <- dist_ll(y, hurd, lam = lam.new, dist = dist, g.x = g.x, log = T)

    lam.ratio <- dgamma(lam.new, a, b, log = T) -
      dgamma(lam, a, b, log = T) +
      sum(like.new) - sum(like.cur)

    if(is.finite(lam.ratio))if(log(runif(1)) < lam.ratio){
      like.cur <- like.new
      lam <- lam.new
      lam.acc <- lam.acc + 1
    }

    output <- list(lam = lam,
                   like.cur = like.new,
                   lam.acc = lam.acc)
  }

  if(dist == "nb" | dist == "lognormal"){

    #Update NB or lognormal mu
    mu.new <- exp(rnorm(1, log(mu), mu.tune))
    like.new <- dist_ll(y, hurd, mu = mu.new, dist = dist, g.x = g.x, log = T)

    mu.ratio <- dgamma(mu.new, a, b, log = T) -
      dgamma(mu, a, b, log = T) +
      sum(like.new) - sum(like.cur)

    if(is.finite(mu.ratio))if(log(runif(1)) < mu.ratio){
      like.cur <- like.new
      mu <- mu.new
      mu.acc <- mu.acc + 1
    }

    output <- list(mu = mu,
                   like.cur = like.new,
                   mu.acc = mu.acc)
  }

  if(dist == "gpd"){

    #Update GPD xi
    xi.new <- rnorm(1, xi, xi.tune)
    like.new <- dist_ll(y, hurd, xi = xi.new, sigma = sigma,
                        dist = dist, g.x = g.x, log = T)

    xi.ratio <- dnorm(xi.new, 0, 1, log = T) - dnorm(xi, 0, 1, log = T) +
      sum(like.new) - sum(like.cur)

    if(log(runif(1)) < xi.ratio){
      like.cur <- like.new
      xi <- xi.new
      xi.acc <- xi.acc + 1
    }

    #Update GPD sigma
    sigma.new <- exp(rnorm(1, log(sigma), sigma.tune))
    like.new <- dist_ll(y, hurd, xi = xi, sigma = sigma.new,
                        dist = dist, g.x = g.x, log = T)

    sigma.ratio <- dgamma(sigma.new, a, b, log = T) -
      dgamma(sigma, a, b, log = T) +
      sum(like.new) - sum(like.cur)

    if(log(runif(1)) < sigma.ratio){
      like.cur <- like.new
      sigma <- sigma.new
      sigma.acc <- sigma.acc + 1
    }

    output <- list(xi = xi,
                   sigma = sigma,
                   xi.acc = xi.acc,
                   sigma.acc = sigma.acc)
  }

  return(output)
}



