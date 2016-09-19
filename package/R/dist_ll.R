
#Data distribution likelihood
dist_ll <- function(y, hurd = Inf, lam = NULL, size = 1, mu = NULL, xi = NULL, sigma = NULL,
                    dist = c("poisson", "nb", "lognormal", "gpd", "none"),
                    g.x = F, log = T){

  trunc1 <- ifelse(g.x, hurd - 1, 0)
  trunc2 <- ifelse(g.x, Inf, hurd - 1)

  f <- match.arg(dist)
  ll <- switch(f,
               poisson = dpois(y, lam, log = T) - log(ppois(trunc2, lam) - ppois(trunc1, lam)),
               nb = dnbinom(y, size = size, mu = mu, log = T) - log(pnbinom(trunc2, size = size, mu = mu) - pnbinom(trunc1, size = size, mu = mu)),
               lognormal = mlnorm(y, meanlog = mu, sdlog = sqrt(mu+mu^2/n), log = T) - log(plnorm(trunc2 + 0.5, meanlog = mu, sdlog = sqrt(mu+mu^2/n)) - plnorm(trunc1 + 0.5, meanlog = mu, sdlog = sqrt(mu+mu^2/n))),
               gpd = mgpd(y, trunc1 + 1, sigma, xi, log = T) - log(pgpd(trunc2 + 0.5, trunc1 + 1, sigma, xi)),
               none = 0)

  if(log){return(ll)}
  else{return(exp(ll))}
}

