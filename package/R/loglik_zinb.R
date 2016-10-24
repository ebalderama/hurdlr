#Define zero-inflated negative binomial data loglikelihood at z = 1, z = 0
loglik_zinb <- function(y, z, mu, size, p) {
  
  ll <- log(1 - p) + dnbinom(y, mu = mu, size = size, log = T)
  ll[z==1] <- log(p[z==1])
  return(ll)
}