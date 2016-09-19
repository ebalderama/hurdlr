
#Probability distribution function of GPD
dgpd <- function(x, mu = 0, sigma = 1, xi = 1, log = F){
  options(warn = -1)
  log.den <- (log(sigma) - (xi + 1)*log(sigma + xi*(x - mu)))/xi
  log.den[x < mu] <- -Inf
  if(xi < 0){log.den[x > (mu - sigma/xi)] <- -Inf}
    if(log){return(log.den)}
    else{return(exp(log.den))}
}

