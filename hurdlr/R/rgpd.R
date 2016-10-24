#' GPD random number generator
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
rgpd <- function(n, mu = 0, sigma = 1, xi = 1){
  options(warn = -1)
  mu + sigma*(runif(n)^(-xi) - 1)/xi
}
