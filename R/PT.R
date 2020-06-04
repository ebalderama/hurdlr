#________________________________________________
#Documentation


#' Typical Count Probability Likelihood
#' 
#' @description \code{PT} is used to calculate the likelihood of a user-defined
#' 'typical' value count observation in a double-hurdle regression model.
#' 
#' @param p vector of zero-count probabilities.
#' 
#' @param q vector of 'typical' count probabilities.
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

#Typical count probability likelihood
PT <- function(p, q, log = T){
  if(log){log(1-p) + log(1-q)}
  else{(1-p)*(1-q)}
}