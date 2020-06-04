#________________________________________________
#Documentation


#' Zero Count Probability Likelihood
#' 
#' @description \code{PZ} is used to calculate the likelihood of a 
#' zero-value count observation in a single or double-hurdle regression model.
#' 
#' @param p vector of zero-count probabilities.
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

#Zero count probability likelihood
PZ <- function(p, log = T){
  if(log){log(p)}
  else{p}
}
