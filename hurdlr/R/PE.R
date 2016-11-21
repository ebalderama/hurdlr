#________________________________________________
#Documentation


#' Extreme Count Probability Likelihood
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
#' 
#' @return A vector of probabilities.
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @seealso \code{\link{hurdle}} 
#' 

#________________________________________________
#Source code

#Extreme count probability likelihood
PE <- function(p, q, log = T){
  if(log){log(1-p) + log(q)}
  else{(1-p)*q}
}