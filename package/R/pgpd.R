
#Cumulative distribution function of GPD
pgpd <- function(q, mu = 0, sigma = 1, xi = 1, lower.tail = T){
  options(warn = -1)
  if(xi != 0) {cdf <- 1 - (1 + xi*(q - mu)/sigma)^(-1/xi)}
  else{cdf <- 1 - exp(-(q - mu)/sigma)}
  cdf[q < mu] <- 0
  if(xi < 0){cdf[q > mu - sigma/xi] <- 1}
    if(lower.tail){return(cdf)}
    else{return(1 - cdf)}
}
