#________________________________________________
#Documentation


#' @name GenPareto
#' 
#' @aliases dgpd
#' @aliases mgpd
#' @aliases pgpd
#' @aliases qgpd
#' @aliases rgpd
#' 
#' @title The Generalized Pareto Distribution
#' 
#' @description Density, distribution function, quantile function and random 
#' generation for the Generalized Pareto distribution with parameters
#' \code{mu}, \code{sigma}, and \code{xi}.
#' 
#' @param x,q vector of quantiles.
#' 
#' @param p numeric predictor matrix.
#' 
#' @param n number of random values to return.
#' 
#' @param mu location parameter.
#' 
#' @param sigma (non-negative) scale parameter.
#' 
#' @param xi shape parameter.
#' 
#' @param log logical; if \code{TRUE}, probabilities p are given 
#' as log(p).
#' 
#' @param lower.tail logical; if \code{TRUE}, probabilities are 
#' \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' 
#' @details 
#' The generalized pareto distribution has density
#' \deqn{f(x) = \frac{\sigma^{\frac{1}{\xi}}}{(\sigma + \xi(x-\mu))^{\frac{1}{\xi}+1}}}{%
#' f(x) = sigma^(1/xi) / (sigma + xi*(x-mu))^(1/xi + 1)}
#' 
#' @return \code{dgpd} gives the continuous density, \code{pgpd} gives the distribution 
#' function, \code{qgpd} gives the quantile function, and \code{rgpd} 
#' generates random deviates. 
#'
#' \code{mgpd} gives a probability mass function for a discretized version of GPD.
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @examples
#' dexp(1,rate=.5) #Exp(rate) equivalent to gpd with mu=0 AND xi=0, and sigma=1/rate.
#' dgpd(1,mu=0,sigma=2,xi=0) #cannot take xi=0.
#' dgpd(1,mu=0,sigma=2,xi=0.0000001) #but can get close.
#' 
#' ##"mass" function of GPD
#' mgpd(8) == pgpd(8.5) - pgpd(7.5)

#________________________________________________
#dgpd()

#' @rdname GenPareto
#' @export

#Probability distribution function of GPD
dgpd <- function(x, mu = 0, sigma = 1, xi = 1, log = F){
  options(warn = -1)
  log.den <- (log(sigma) - (xi + 1)*log(sigma + xi*(x - mu)))/xi
  log.den[x < mu] <- -Inf
  if(xi < 0){log.den[x > (mu - sigma/xi)] <- -Inf}
    if(log){return(log.den)}
    else{return(exp(log.den))}
}

#________________________________________________
#mgpd()

#' @rdname GenPareto
#' @export

#Probability distribution function of discrete GPD
mgpd <- function(x, mu = 0, sigma = 1, xi = 1, log = F){
  pmf <- pgpd(x + 0.5, mu, sigma, xi) - pgpd(x - 0.5, mu, sigma, xi)
  if(log){return(log(pmf))}
  else{return(pmf)}
}

#________________________________________________
#pgpd()

#' @rdname GenPareto
#' @export

#Cumulative distribution function of GPD
pgpd <- function(q, mu = 0, sigma = 1, xi = 1, lower.tail = T){
  options(warn = -1)
  if(xi != 0) {cdf <- 1 - (1 + xi*(q - mu)/sigma)^(-1/xi)}
  else{cdf <- 1 - exp(-(q - mu)/sigma)}
  cdf[q < mu] <- 0
  if(xi < 0){cdf[q > mu - sigma/xi] <- 1}
  if(lower.tail){return(cdf)}
  else{return(1 - cdf)}
}

#________________________________________________
#qgpd()

#' @rdname GenPareto
#' @export

#Quantile function of GPD
qgpd <- function(p, mu = 0, sigma = 1, xi = 1, lower.tail = T){
  options(warn = -1)
  if(lower.tail == F){p <- 1 - p}
  if(xi != 0){inv.cdf <- sigma*((1-p)^(-xi) - 1)/xi + mu}
  else{inv.cdf <- -sigma*log(1-p) + mu}
  if(xi < 0 & inv.cdf > mu - sigma/xi){inv.cdf <- 1}
  return(inv.cdf)
}

#________________________________________________
#rgpd()

#' @rdname GenPareto
#' @export

#GPD random number generator
rgpd <- function(n, mu = 0, sigma = 1, xi = 1){
  options(warn = -1)
  mu + sigma*(runif(n)^(-xi) - 1)/xi
}
