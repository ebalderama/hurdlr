
#Probability distribution function of discrete GPD
mgpd <- function(x, mu = 0, sigma = 1, xi = 1, log = F){
  pmf <- pgpd(x + 0.5, mu, sigma, xi) - pgpd(x - 0.5, mu, sigma, xi)
  if(log){return(log(pmf))}
    else{return(pmf)}
}
