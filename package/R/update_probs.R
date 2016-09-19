
update_probs <- function(y, x, hurd, p, q,
                         beta.prior.mean, beta.prior.sd,
                         pZ, pT, pE,
                         beta, XB2, XB3, beta.acc, beta.tune){

  pZ.cur <- pZ
  pT.cur <- pT
  pE.cur <- pE

  #Update p
  beta2 <- beta[,2]
  beta2.acc <- beta.acc[,2]
  beta2.tune <- beta.tune[,2]



  for(j in 1:nrow(beta)){

    beta.new <- rnorm(1, beta2[j], beta2.tune[j])
    XB.new <- XB2 + x[,j]*(beta.new - beta2[j])
    p.new <- 1/(1 + exp(-XB.new))
    pZ.new <- PZ(p.new[y == 0], log = T)
    pT.new <- PT(p.new[0 < y & y <hurd], q[0 < y & y <hurd], log = T)
    pE.new <- PE(p.new[y >=hurd], q[y >=hurd], log = T)

    p.ratio <- dnorm(beta.new, beta.prior.mean, beta.prior.sd, log = T) -
      dnorm(beta2[j], beta.prior.mean, beta.prior.sd, log = T) +
      sum(pZ.new, pT.new, pE.new) - sum(pZ.cur, pT.cur, pE.cur)

    if(is.finite(p.ratio))if(log(runif(1)) < p.ratio){
      beta2[j] <- beta.new
      beta2.acc[j] <- beta2.acc[j] + 1
      pZ.cur <- pZ.new
      pT.cur <- pT.new
      pE.cur <- pE.new
    }
  }

  #Update q
  beta3 <- beta[,3]
  beta3.acc <- beta.acc[,3]
  beta3.tune <- beta.tune[,3]

  if(hurd != Inf){

    for(j in 1:nrow(beta)){

      beta.new <- rnorm(1, beta3[j], beta3.tune[j])
      XB.new <- XB3 + x[,j]*(beta.new - beta3[j])
      q.new <- 1/(1 + exp(-XB.new))
      pT.new <- PT(p[0 < y & y <hurd], q.new[0 < y & y <hurd], log = T)
      pE.new <- PE(p[y >=hurd], q.new[y >=hurd], log = T)

      q.ratio <- dnorm(beta.new, beta.prior.mean, beta.prior.sd, log = T) -
        dnorm(beta3[j], beta.prior.mean, beta.prior.sd, log = T) +
        sum(pT.new, pE.new) - sum(pT.cur, pE.cur)

      if(is.finite(q.ratio))if(log(runif(1)) < q.ratio){
        beta3[j] <- beta.new
        beta3.acc[j] <- beta3.acc[j] + 1
        pT.cur <- pT.new
        pE.cur <- pE.new
      }
    }
  }

  output <- list(beta2 = beta2,
                 beta3 = beta3,
                 beta2.acc = beta2.acc,
                 beta3.acc = beta3.acc
  )

  return(output)
}


