
#GPD random number generator
rgpd <- function(n, mu = 0, sigma = 1, xi = 1){
  options(warn = -1)
  mu + sigma*(runif(n)^(-xi) - 1)/xi
}
