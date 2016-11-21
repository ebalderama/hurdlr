#Zero-Inflated Negative Binomial model

#Zero-Infladed Negative Binomial likelihood:

#              {   p,                   zi = 1
#    Yi|zi =   {
#              {   1 - p * f(x|mu),     zi = 0

# logit(p) = g(p) = BX = B0 + B1*x1 + B2*x2 ...
# f(x|mu, size) ~ NB(mu, size)
# zi ~ Binomial(1, p)
# mu ~ Gamma(a, b)


#________________________________________________
#Documentation


#' Zero-Inflated Negative Binomial Regression Model
#' 
#' @description \code{zero_nb} is used to fit zero-inflated 
#' negative binomial regression models to count data via Bayesian inference.
#' 
#' @param y numeric response vector.
#' 
#' @param x numeric predictor matrix.
#' 
#' @param size size parameter for negative binomial likelihood distributions.
#' 
#' @param a shape parameter for gamma prior distributions.
#' 
#' @param b rate parameter for gamma prior distributions.
#' 
#' @param mu.start initial value for mu parameter.
#' 
#' @param beta.prior.mean mu parameter for normal prior distributions.
#' 
#' @param beta.prior.sd standard deviation for normal prior distributions.
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
#' @details Fits a zero-inflated negative binomial (ZINB) model.
#' 
#' @return \code{zero_nb} returns a list which includes the items
#' \describe{
#'    \item{mu}{numeric vector; posterior distribution of mu parameter}
#'    \item{beta}{numeric matrix; posterior distributions of regression coefficients}
#'    \item{p}{numeric vector; posterior distribution of parameter 'p', the 
#'    probability of a given zero observation belonging to the model's zero component}
#'    \item{ll}{numeric vector; posterior log-likelihood}    
#' }
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' 

#________________________________________________
#Source code

zero_nb <- function(y, x, size, a = 1, b = 1, mu.start = 1, 
                    beta.prior.mean = 0, beta.prior.sd = 1,
                    iters = 1000, burn = 500, nthin = 1,
                    plots = T, progress.bar = T) {
  
  
  #Initial values
  pb <- txtProgressBar(min = 0, max = iters, style = 3)
  x <- cbind(1, x)
  x.col <- ncol(x)
  z <- rep(1, length(y))
  B <- rep(beta.prior.mean, length.out = x.col)
  XB <- x%*%B
  p <- 1/(1 + exp(-XB))
  mu <- mu.start
  ll <- loglik_zinb(y = y, z = z, mu = mu, size = size, p = p)
  
  #Tuning
  beta.tune <- rep(0.02, x.col)
  mu.tune <- 0.02
  beta.acc <- rep(0, x.col)
  mu.acc <- 0
  att <- 0
  
  #Create matrix of MCMC values for all parameters B, mu
  keep.beta <- matrix(0, iters, length(B))
  keep.mu <- rep(mu, iters)
  keep.ll <- rep(sum(ll), iters)
  keep.p <- rep(mean(p), iters)
  
  #MCMC start
  for(i in 2:iters) {for(thin in 1:nthin){
    
    att <- att + 1
    
    #Update p (Update Betas)
    
    for(j in 1:x.col){
      
      beta.current <- B[j]
      beta.new <- rnorm(1, beta.current, beta.tune[j])
      XB.new <- XB + x[,j]*(beta.new - beta.current)
      p.new <- 1/(1 + exp(-XB.new))
      ll.new <- loglik_zinb(y = y, z = z, mu = mu, size = size, p = p.new)
      
      p.ratio <- sum(ll.new) + sum(dnorm(beta.new, beta.prior.mean, beta.prior.sd, log = T)) -
        sum(ll) - sum(dnorm(beta.current, beta.prior.mean, beta.prior.sd, log = T))
      #if(is.finite(p.ratio))
      if(log(runif(1)) < p.ratio){
        
        B[j] <- beta.new
        XB <- XB.new
        p <- p.new
        ll <- ll.new
        beta.acc[j] <- beta.acc[j] + 1
      }
    }
    
    #Update z
    
    P1 <- ifelse(y==0, 1, 0)*p
    P0 <- dnbinom(y, mu = mu, size = size)*(1 - p)
    z  <- rbinom(length(y), 1, P1/(P1 + P0))
    ll <- loglik_zinb(y = y, z = z, mu = mu, size = size, p = p)
    
    #Update mu
    
    mu.current <- mu
    mu.new <- exp(rnorm(1, log(mu.current), mu.tune))
    
    mu.ratio <- dgamma(mu.new, a, b, log = T) + sum(dpois(y[z==0], mu.new, log = T)) -
      dgamma(mu.current, a, b, log = T) - sum(dpois(y[z==0], mu.current, log = T))
    
    if(log(runif(1)) < mu.ratio) {
      mu <- mu.new
      mu.acc <- mu.acc + 1
      ll <- loglik_zinb(y = y, z = z, mu = mu.new, size = size, p = p)
    }
    
  }#End thinning
    
    
    keep.beta[i,] <- B
    keep.mu[i] <- mu
    keep.ll[i] <- sum(ll)
    keep.p[i] <- mean(p)
    
    
    #Update tuning parameters
    
    if(i < 0.75*burn & att > 50){
      beta.tune <- ifelse(beta.acc/att < 0.20, 0.8*beta.tune, beta.tune)
      beta.tune <- ifelse(beta.acc/att > 0.50, 1.2*beta.tune, beta.tune)
      mu.tune <- ifelse(mu.acc/att < 0.20, 0.8*mu.tune, mu.tune)
      mu.tune <- ifelse(mu.acc/att > 0.50, 1.2*mu.tune, mu.tune)
      
      beta.acc <- rep(0, x.col)
      mu.acc <- att <- 0
    }
    
    #Plots
    
    if(plots==T){
      if(i > burn & i%%100==0){
        
        par(mfrow = c(1,3))
        plot(keep.mu[burn:i], type = "s",
             ylab = "", main = bquote(Mean~mu))
        abline(h = mean(keep.mu[burn:i]), col = 2)
        plot(keep.p[burn:i],type="s",
             ylab = "", main = bquote(beta[0]~Zero~Prob))
        abline(h = mean(keep.p[burn:i]), col = 2)
        plot(keep.ll[burn:i], type = "s",
             ylab = "", main = bquote(Log~Likelihood))
        abline(h = mean(keep.ll[burn:i]), col = 2)
      }
    }
    
    if(progress.bar){setTxtProgressBar(pb, i)}
  }
  
  return(list(mu = keep.mu[(burn + 1):iters],
              beta = keep.beta[(burn + 1):iters,],
              p = keep.p[(burn + 1):iters,],
              ll = keep.ll[(burn + 1):iters,]
  ))
  
  
}
