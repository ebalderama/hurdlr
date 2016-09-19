hurdle_control <- function(a = 1, b = 1, size = 1,
                           beta.prior.mean = 0, beta.prior.sd = 1000,
                           beta.tune = 1, pars.tune = 0.2,
                           lam.start = 1, mu.start = 1,
                           sigma.start = 1, xi.start = 1){

  output <- list(a = a, b = b, size = size,
                 beta.prior.mean = beta.prior.mean, beta.prior.sd = beta.prior.sd,
                 pars.tune = pars.tune, beta.tune = beta.tune,
                 lam.start = lam.start, mu.start = mu.start,
                 sigma.start = sigma.start, xi.start = xi.start)

  return(output)
}
