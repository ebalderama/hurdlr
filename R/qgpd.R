
#Quantile function of GPD
qgpd <- function(p, mu = 0, sigma = 1, xi = 1, lower.tail = T){
  options(warn = -1)
  if(lower.tail == F){p <- 1 - p}
  if(xi != 0){inv.cdf <- sigma*((1-p)^(-xi) - 1)/xi + mu}
    else{inv.cdf <- -sigma*log(1-p) + mu}
    if(xi < 0 & inv.cdf > mu - sigma/xi){inv.cdf <- 1}
  return(inv.cdf)
}