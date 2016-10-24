#' Probability distribution function of discrete log normal
#' 
#' @description \code{hurdle} is used to fit single or 
#' double-hurdle regression models to count data via Bayesian inference.
#' 
#' @param y numeric response vector.
#' 
#' @param x numeric predictor matrix.
#' 
#' @param hurd numeric threshold for 'extreme' observations of two-hurdle models. 
#' \code{Inf} for one-hurdle models.
#' 
#' @param dist character specification of response distribution.
#' 
#' @param dist.2 character specification of response distribution for 
#' 'extreme' observations of two-hurdle models.
#' 
#' @param control list of parameters for controlling the fitting process, 
#' specified by \code{\link{hurdle_control}}.
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
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Dr. Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @example 
#' 
mlnorm <- function(y, meanlog = 0, sdlog = 1, log = T){
  pmf <- plnorm(y + 0.5, meanlog, sdlog) - plnorm( y - 0.5, meanlog, sdlog)
  if(log){return(log(pmf))}
    else{return(pmf)}
}
