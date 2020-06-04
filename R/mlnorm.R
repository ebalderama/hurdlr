#________________________________________________
#Documentation


#' Density Function for Discrete Log Normal Distribution
#' 
#' @description Density function of the discrete log normal distribution
#' whose logarithm has mean equal to \code{meanlog} and standard deviation 
#' equal to \code{sdlog}.
#' 
#' @param x vector of quantiles.
#' 
#' @param meanlog mean of the distribution on the log scale.
#' 
#' @param sdlog standard deviation of the distribution on the log scale.
#' 
#' @param log logical; if \code{TRUE}, probabilities p are given 
#' as log(p).
#' 
#' 
#' @return Discrete log-normal distributional density.
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' 

#________________________________________________
#mlnorm()

#Probability distribution function of discrete log normal
mlnorm <- function(x, meanlog = 0, sdlog = 1, log = T){
  pmf <- plnorm(x + 0.5, meanlog, sdlog) - plnorm(x - 0.5, meanlog, sdlog)
  if(log){return(log(pmf))}
    else{return(pmf)}
}
