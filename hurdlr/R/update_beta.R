#________________________________________________
#Documentation


#' MCMC Second-Component Parameter Update Function for Hurdle Model Count Data
#' Regression
#' 
#' @description MCMC algorithm for updating the second-component likelihood 
#' parameters in hurdle model regression using \code{\link{hurdle}}.
#' 
#' @param y numeric response vector of observations within the bounds of the 
#' second component of the likelihood function, \eqn{y[0 < y \& y < hurd]}
#' 
#' @param x optional numeric predictor matrix for response observations within 
#' the bounds of the second component of the likelihood function, 
#' \eqn{y[0 < y \& y < hurd]}.
#' 
#' @param hurd numeric threshold for 'extreme' observations of two-hurdle models.
#' 
#' @param dist character specification of response distribution for the third 
#' component of the likelihood function.
#' 
#' @param like.part numeric vector of the current third-component likelihood values.
#' 
#' @param beta.prior.mean mu parameter for normal prior distributions.
#' 
#' @param beta.prior.sd standard deviation for normal prior distributions.
#' 
#' @param beta numeric matrix of current regression coefficient parameter values.
#' 
#' @param XB \eqn{x*beta[,1]} product matrix for response observations within the bounds 
#' of the second component of the likelihood function, \eqn{y[0 < y \& y < hurd]}.
#' 
#' @param beta.acc numeric matrix of current MCMC acceptance rates for 
#' regression coefficient parameters.
#' 
#' @param beta.tune numeric matrix of current MCMC tuning values for regression 
#' coefficient estimation.
#' 
#' @param g.x logical operator. \code{TRUE} if operating within the third component 
#' of the likelihood function (the likelihood of 'extreme' observations).
#' 
#' 
#' @return A list of MCMC-updated regression coefficients for the estimation of 
#' the second-component likelihood parameter as well as each coefficient's MCMC
#' acceptance ratio.
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

update_beta <- function(y, x, hurd, dist, like.part,
                        beta.prior.mean, beta.prior.sd,
                        beta, XB, beta.acc, beta.tune,
                        g.x = F){

  beta1 <- beta[,1]
  beta1.acc <- beta.acc[,1]
  beta1.tune <- beta.tune[,1]
  x <- as.matrix(x)

  like.cur <- like.part

  if(dist == "poisson"){

    for(j in 1:nrow(beta)){

      #Update Poisson lambda
      beta.new <- rnorm(1, beta1[j], beta1.tune[j])
      XB.new <- XB + x[,j]*(beta.new - beta1[j])
      lam.new <- exp(XB.new)
      like.new <- dist_ll(y, hurd, lam = lam.new, dist = dist, g.x = g.x, log = T)

      lam.ratio <- dnorm(beta.new, beta.prior.mean, beta.prior.sd, log = T) -
        dnorm(beta1[j], beta.prior.mean, beta.prior.sd, log = T) +
        sum(like.new) - sum(like.cur)

      if(is.finite(lam.ratio))if(log(runif(1)) < lam.ratio){
        like.cur <- like.new
        beta1[j] <- beta.new
        beta1.acc[j] <- beta1.acc[j] + 1
      }
    }
  }

  if(dist == "nb" | dist == "lognormal"){

    for(j in 1:nrow(beta)){

      #Update NB or lognormal mu
      beta.new <- rnorm(1, beta1[j], beta1.tune[j])
      XB.new <- XB + x[,j]*(beta.new - beta1[j])
      mu.new <- exp(XB.new)
      like.new <- dist_ll(y, hurd, mu = mu.new, dist = dist, g.x = g.x, log = T)

      mu.ratio <- dnorm(beta.new, beta.prior.mean, beta.prior.sd, log = T) -
        dnorm(beta1[j], beta.prior.mean, beta.prior.sd, log = T) +
        sum(like.new) - sum(like.cur)

      if(is.finite(mu.ratio))if(log(runif(1)) < mu.ratio){
        like.cur <- like.new
        beta1[j] <- beta.new
        beta1.acc[j] <- beta1.acc[j] + 1
      }
    }
  }

### Comment out GPD from the typical distribution
#
#  if(dist == "gpd"){
#
#    for(j in 1:nrow(beta)){
#
#      #Update GPD sigma
#      beta.new <- rnorm(1, beta1[j], beta1.tune[j])
#      XB.new <- XB + x[,j]*(beta.new - beta1[j])
#      sigma.new <- exp(XB.new)
#      like.new <- dist_ll(y, hurd, xi = xi, sigma = sigma.new, dist = dist, g.x = F, log = T)
#
#      sigma.ratio <- dnorm(beta.new, beta.prior.mean, beta.prior.sd, log = T) -
#        dnorm(beta1[j], beta.prior.mean, beta.prior.sd, log = T) +
#        sum(like.new) - sum(like.cur)
#
#      if(is.finite(sigma.ratio))if(log(runif(1)) < sigma.ratio){
#        like.cur <- like.new
#        beta1[j] <- beta.new
#        beta1.acc[j] <- beta1.acc[j] + 1
#      }
#    }
#  }

  output <- list(beta1 = beta1,
                 beta1.acc = beta1.acc)

  return(output)
}





