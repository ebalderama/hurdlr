#Define zero-inflated poisson data loglikelihood at z = 1, z = 0
loglik_zip <- function(y, z, lam, p) {
  
  ll <- log(1 - p) + dpois(y, lam, log = T)
  ll[z==1] <- log(p[z==1])
  return(ll)
}