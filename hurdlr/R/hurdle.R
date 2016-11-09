#Truncated distribution function

#One-hurdle likelihood:

#                 {   p,                            yi = 0
#    Yi|theta =   {
#                 {   (1-p)*f(x|theta),             yi > 0

#Two-hurdle liklihood:

#                 {   p,                            yi = 0
#                 {
#    Yi|theta =   {   (1-q)*(1-p)*f(x|theta),       0 < yi < psi
#                 {
#                 {   q*(1-p)*g(x|theta),           yi >= psi


#________________________________________________
#Documentation


#' Hurdle Model Count Data Regression
#' 
#' @description \code{hurdle} is used to fit single or 
#' double-hurdle regression models to count data via Bayesian inference.
#' 
#' @param y numeric response vector.
#' 
#' @param x numeric predictor matrix.
#' 
#' @param hurd numeric threshold for 'extreme' observations of two-hurdle models. 
#' \code{Inf} for one-hurdle models.
#' 
#' @param dist character specification of response distribution.
#' 
#' @param dist.2 character specification of response distribution for 
#' 'extreme' observations of two-hurdle models.
#' 
#' @param control list of parameters for controlling the fitting process, 
#' specified by \code{\link{hurdle_control}}.
#' 
#' @param iters number of iterations for the Markov chain to run.
#' 
#' @param burn numeric burn-in length.
#' 
#' @param nthin numeric thinning rate.
#' 
#' @param plots logical operator. \code{TRUE} to output plots.
#' 
#' @param progress.bar logical operator. \code{TRUE} to print progress bar.
#' 
#' @details 
#' 
#' @return \code{Hurdle} returns a list which includes the items
#' \describe{
#'    \item{pD}{measure of model dimensionality \eq{p_D} where 
#'    \eq{p_D = \bar{D} - D(\bar{\theta}}) is the \eq{"mean posterior deviance - 
#'    deviance of posterior means"}}
#'    \item{DIC}{Deviance Information Criterion where \eq{DIC = \bar{D} - p_D}}
#'    \item{PPO}{Posterior Predictive Ordinate (PPO) measure of fit}
#'    \item{CPO}{Conditional Predictive Ordinate (CPO) measure of fit}
#'    \item{pars.means}{posterior mean(s) of third-component parameter(s) if 
#'    \code{hurd != Inf}}  
#'    \item{ll.means}{posterior means of the log-likelihood distributions of 
#'    all model components}
#'    \item{beta.means}{posterior means regression coefficients}
#'    \item{dev}{posterior deviation where \eq{D = -2LogL}}
#'    \item{beta}{posterior distributions of regression coefficients}
#'    \item{pars}{posterior distribution(s) of third-component parameter(s) if 
#'    \code{hurd != Inf}}  
#'  }
#' 
#' @author 
#' Taylor Trippe <\email{ttrippe@@luc.edu}> \cr
#' Earvin Balderama <\email{ebalderama@@luc.edu}>
#' 
#' @examples
#' #Generate some data:
#' p=0.5; q=0.25; lam=3;
#' mu=10; sigma=7; xi=0.75;
#' n=200
#'
#' set.seed(2016)
#' y <- rbinom(n,1,p)
#' nz <- sum(1-y)
#' extremes <- rbinom(sum(y),1,q)
#' ne <- sum(extremes)
#' nt <- n-nz-ne
#' yt <- sample(mu-1,nt,replace=T,prob=dpois(1:(mu-1),3)/(ppois(mu-1,lam)-ppois(0,lam)))
#' yz <- round(rgpd(nz,mu,sigma,xi))
#' y[y==1] <- c(yt,yz)
#' 
#' m1 <- hurdle(y)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 

#________________________________________________
#Source code

hurdle <- function(y, x = NULL, hurd = Inf,
                   dist = c("poisson", "nb", "lognormal", "gpd"),
                   dist.2 = c("none", "gpd", "poisson", "lognormal", "nb"),
                   control = hurdle_control(),
                   iters = 1000, burn = 500, nthin = 1,
                   plots = T, progress.bar = T){

  #Initial values
  N <- length(y)
  dist <- match.arg(dist)
  dist.2 <- match.arg(dist.2)
  hurd <- ifelse(dist.2 == "none", Inf, hurd)
  pb <- txtProgressBar(min = 0, max = iters, style = 3)
  x <- cbind(rep(1, N), x)
  x <- as.matrix(x)
  x.col <- ncol(x)
  regressions <- 3
  out.fx <- out.gx <- out.probs <- NULL

  beta <- matrix(0, x.col, regressions)
  XB1 <- c(x%*%beta[,1])
  XB2 <- c(x%*%beta[,2])
  XB3 <- c(x%*%beta[,3])

  lam <- mu <- sigma <- exp(XB1)
  xi <- control$xi.start
  p <- 1/(1 + exp(-XB2))
  q <- 1/(1 + exp(-XB3))

  #Priors
  a <- control$a
  b <- control$b
  size <- control$size
  beta.prior.mean <- control$beta.prior.mean
  beta.prior.sd <- control$beta.prior.sd

  #Tuning
  beta.tune <- matrix(rep(c(control$pars.tune, rep(control$probs.tune, 2)), each = x.col), x.col, regressions)
  beta.acc <- matrix(0, x.col, regressions)
  att <- 0

  #Keep vectors
  keep.beta <- array(0, c(iters, x.col, regressions))
  keep.pars <- NULL

  #Diagnostics
  keep.dev <- track.ppo <- track.icpo <- 0

  #If using a two-hurdle model
  if(hurd != Inf){

    #Initial values
    lam2 <- control$lam.start
    mu2 <- control$mu.start
    xi2 <- control$xi.start
    sigma2 <- control$sigma.start

    #Tuning
    lam2.tune <- mu2.tune <- xi2.tune <- sigma2.tune <- control$pars.tune
    lam2.acc <- mu2.acc <- xi2.acc <- sigma2.acc <- 0

    #Keep vectors
    keep.pars <- matrix(0, iters, 4)
  }


  #Likelihood
  pZ <- PZ(p[y == 0], log = T)
  pT <- PT(p[0 < y & y < hurd], q[0 < y & y < hurd], log = T)
  LT <- dist_ll(y = y[0 < y & y < hurd], hurd = hurd, size = size,
                lam = lam[0 < y & y < hurd], mu = mu[0 < y & y < hurd],
                xi = xi, sigma = sigma[0 < y & y < hurd], dist = dist, g.x = F, log = T)
  pE <- PE(p[y >= hurd], q[y >= hurd], log = T)
  LE <- dist_ll(y = y[y >= hurd], hurd = hurd, size = size,
                lam = control$lam.start, mu = control$mu.start,
                xi = xi, sigma = control$sigma.start, dist = dist.2, g.x = T, log = T)
  ll <- NULL
  ll[y == 0] <- pZ
  ll[0 < y & y < hurd] <- pT + LT
  ll[y >= hurd] <- pE + LE

  #MCMC start
  for(i in 2:iters){
    for(thin in 1:nthin){

      att <- att + 1

      #Update f(x) parameters
      out.fx <- update_beta(y = y[0 < y & y < hurd], x[which(0 < y & y < hurd),],
                            hurd, dist = dist, like.part = LT,
                            beta.prior.mean, beta.prior.sd,
                            beta, XB1[0 < y & y < hurd], beta.acc, beta.tune,
                            g.x = F)

      if(dist == "poisson"){
        beta[,1] <- out.fx$beta1; beta.acc[,1] <- out.fx$beta1.acc
        XB1 <- x%*%beta[,1]
        lam <- exp(XB1)
      }

      if(dist == "nb" | dist == "lognormal"){
        beta[,1] <- out.fx$beta1; beta.acc[,1] <- out.fx$beta1.acc
        XB1 <- x%*%beta[,1]
        mu <- exp(XB1)
      }

      if(dist == "gpd"){
        beta[,1] <- out.fx$beta1; beta.acc[,1] <- out.fx$beta1.acc
        XB1 <- x%*%beta[,1]
        sigma <- exp(XB1)
      }

      LT <- dist_ll(y = y[0 < y & y < hurd], hurd = hurd, size = size,
                    lam = lam[0 < y & y < hurd], mu = mu[0 < y & y < hurd], xi = xi,
                    sigma = sigma[0 < y & y < hurd], dist = dist, g.x = F, log = T)


      if(hurd != Inf){

        #Update g(x) parameters
        out.gx <- update_pars(y = y[y >= hurd], hurd,
                              dist = dist.2, like.part = LE,
                              a, b, size, lam2, mu2, xi2, sigma2,
                              lam2.acc, mu2.acc, xi2.acc, sigma2.acc,
                              lam2.tune, mu2.tune, xi2.tune, sigma2.tune,
                              g.x = T)
        if(dist.2 == "poisson"){
          lam2 <- out.gx$lam; lam2.acc <- out.gx$lam.acc
        }

        if(dist.2 == "nb" | dist.2 == "lognormal"){
          mu2 <- out.gx$mu; mu2.acc <- out.gx$mu.acc
        }

        if(dist.2 == "gpd"){
          xi2 <- out.gx$xi; xi2.acc <- out.gx$xi.acc
          sigma2 <- out.gx$sigma; sigma2.acc <- out.gx$sigma.acc
        }

        LE <- dist_ll(y = y[y >= hurd], hurd = hurd, size = size,
                      lam = lam2, mu = mu2, xi = xi2,
                      sigma = sigma2, dist = dist.2, g.x = T, log = T)
      }

      out.probs <- update_probs(y, x, hurd, p, q,
                                beta.prior.mean, beta.prior.sd,
                                pZ, pT, pE,
                                beta, XB2, XB3, beta.acc, beta.tune)

      beta[,2] <- out.probs$beta2; beta.acc[,2] <- out.probs$beta2.acc
      XB2 <- x%*%beta[,2]
      p <- 1/(1 + exp(-XB2))

      beta[,3] <- out.probs$beta3; beta.acc[,3] <- out.probs$beta3.acc
      XB3 <- x%*%beta[,3]
      q <- 1/(1 + exp(-XB3))

      #Update likelihood
      pZ <- PZ(p[y == 0], log = T)
      pT <- PT(p[0 < y & y < hurd], q[0 < y & y < hurd], log = T)
      pE <- PE(p[y >= hurd], q[y >= hurd], log = T)

      ll[y == 0] <- pZ
      ll[0 < y & y < hurd] <- pT + LT
      ll[y >= hurd] <- pE + LE


    }#End thinning

    keep.dev[i] <- -2*sum(ll)
    keep.beta[i,,] <- beta

    if(hurd != Inf){
      keep.pars[i,] <- c(lam2, mu2, sigma2, xi2)
    }

    if(i > burn){
      track.ppo <- track.ppo + exp(ll)
      track.icpo <- track.icpo + 1/exp(ll)
    }

    #Update tuning parameters
    if(i < 0.75*burn & att > 50){

      for(k in 1:regressions){
        beta.tune[,k] <- ifelse(beta.acc[,k]/att < 0.20, 0.8*beta.tune[,k], beta.tune[,k])
        beta.tune[,k] <- ifelse(beta.acc[,k]/att > 0.50, 1.2*beta.tune[,k], beta.tune[,k])
      }

      beta.acc <- matrix(0, x.col, regressions)

      if(hurd != Inf){
        lam2.tune <- ifelse(lam2.acc/att < 0.20, 0.8*lam2.tune, lam2.tune)
        lam2.tune <- ifelse(lam2.acc/att > 0.50, 1.2*lam2.tune, lam2.tune)
        mu2.tune <- ifelse(mu2.acc/att < 0.20, 0.8*mu2.tune, mu2.tune)
        mu2.tune <- ifelse(mu2.acc/att > 0.50, 1.2*mu2.tune, mu2.tune)
        xi2.tune <- ifelse(xi2.acc/att < 0.20, 0.8*xi2.tune, xi2.tune)
        xi2.tune <- ifelse(xi2.acc/att > 0.50, 1.2*xi2.tune, xi2.tune)
        sigma2.tune <- ifelse(sigma2.acc/att < 0.20, 0.8*sigma2.tune, sigma2.tune)
        sigma2.tune <- ifelse(sigma2.acc/att > 0.50, 1.2*sigma2.tune, sigma2.tune)

        lam2.acc <- mu2.acc <- xi2.acc <- sigma2.acc <- 0
      }
      att <- 0
    }

    #Plots
    if(plots == T){

      if(i > burn & i%%200 == 0){

        num.plots <- ifelse(hurd != Inf, par(mfrow = c(2,3)), par(mfrow = c(1,3))); num.plots
        plot(keep.beta[burn:i,1,1], type = "l",
             ylab = "", main = bquote(beta[0]~Typical~Mean))
        abline(h = mean(keep.beta[burn:i,1,1]), col = 2)
        plot(keep.beta[burn:i,1,2], type = "l",
             ylab = "", main = bquote(beta[0]~Zero~Prob))
        abline(h = mean(keep.beta[burn:i,1,2]), col = 2)

        if(hurd != Inf){

          plot(keep.beta[burn:i,1,3], type = "l",
               ylab = "", main = bquote(beta[0]~Extreme~Prob))
          abline(h = mean(keep.beta[burn:i,1,3]), col = 2)

          if(dist.2 == "poisson"){
            plot(keep.pars[burn:i,1], type = "l",
                 ylab = "", main = bquote(Extreme~Mean~lambda))
            abline(h = mean(keep.pars[burn:i,1]), col = 2)
          }

          if(dist.2 == "nb" | dist.2 == "lognormal"){
            plot(keep.pars[burn:i,2], type = "l",
                 ylab = "", main = bquote(Extreme~Mean~mu))
            abline(h = mean(keep.pars[burn:i,2]), col = 2)
          }

          if(dist.2 == "gpd"){
            plot(keep.pars[burn:i,4], type = "l",
                 ylab = "", main = bquote(xi[gpd]))
            abline(h = mean(keep.pars[burn:i,4]), col = 2)
            plot(keep.pars[burn:i,3], type = "l",
                 ylab = "", main = bquote(sigma[gpd]))
            abline(h = mean(keep.pars[burn:i,3]), col = 2)
          }
        }

        plot(burn:i, keep.dev[burn:i], type="l",
             ylab = "", xlab = "Index", main = bquote(D==-2*sum(loglike)))
        abline(h = mean(keep.dev[burn:i]), col = 2)
      }
    }

    if(progress.bar){setTxtProgressBar(pb, i)}

  }#End loop


  #Diagnostics (plots, pD, DIC)
  beta.bar <- apply(array(keep.beta[(burn):iters,,], c(iters-burn, x.col, regressions)), 2:3, mean)

  XB1 <- c(x%*%beta.bar[,1])
  XB2 <- c(x%*%beta.bar[,2])
  XB3 <- c(x%*%beta.bar[,3])

  typ.mean.bar <- pz.bar <- pe.bar <- NULL
  typ.mean.bar <- exp(XB1)
  pz.bar <- 1/(1 + exp(-XB2))
  pe.bar <- 1/(1 + exp(-XB3))

  if(hurd != Inf){
    ext.lam.bar <- mean(keep.pars[burn:i,1])
    ext.mu.bar <- mean(keep.pars[burn:i,2])
    ext.sigma.bar <- mean(keep.pars[burn:i,3])
    ext.xi.bar <- mean(keep.pars[burn:i,4])
  }

  PZ.bar <- PT.bar <- PE.bar <- LT.bar <- LE.bar <- NULL
  PZ.bar <- PZ(pz.bar[y == 0], log = T)
  PT.bar <- PT(pz.bar[0 < y & y < hurd], pe.bar[0 < y & y < hurd], log = T)
  PE.bar <- PE(pz.bar[y >= hurd], pe.bar[y >= hurd], log = T)
  LT.bar <- dist_ll(y[0 < y & y < hurd], hurd = hurd,
                    lam = typ.mean.bar, mu = typ.mean.bar,
                    dist = dist, g.x = F, log = T)

  ll.means <- c(PZ.bar, PT.bar, LT.bar)

  if(hurd != Inf){
    LE.bar <- dist_ll(y[y > hurd], hurd = hurd,
                      lam = ext.lam.bar, mu = ext.mu.bar,
                      sigma = ext.sigma.bar, xi = ext.xi.bar,
                      dist = dist.2, g.x = T, log = T)

    ll.means <- c(PZ.bar, PT.bar, PE.bar, LT.bar, LE.bar)
  }

  #pD (effective # of parameters) = mean deviance - deviance at posterior means
  pD <- mean(keep.dev[(burn + 1):iters]) - (-2*sum(ll.means))

  #DIC = mean deviance + pD
  DIC <- mean(keep.dev[(burn + 1):iters]) + pD


  ll.means <- list(PZ.bar = PZ.bar, PT.bar = PT.bar, LT.bar = LT.bar)
  pars.means <- list(typ.mean.bar = typ.mean.bar, pz.bar = pz.bar)
  beta.means <- list(typ.mean = beta.bar[,1], pz = beta.bar[,2])

  if(hurd != Inf){

    if(dist.2 == "poisson"){
      pars.means <- list(typ.mean.bar = typ.mean.bar, pz.bar = pz.bar,
                         pe.bar = pe.bar, ext.lam.bar = ext.lam.bar)
    }

    if(dist.2 == "nb" | dist.2 == "lognormal"){
      pars.means <- list(typ.mean.bar = typ.mean.bar, pz.bar = pz.bar,
                         pe.bar = pe.bar, ext.mu.bar = ext.mu.bar)
    }

    if(dist.2 == "gpd"){
      pars.means <- list(typ.mean.bar = typ.mean.bar, pz.bar = pz.bar,
                         pe.bar = pe.bar, ext.sigma.bar = ext.sigma.bar, ext.xi.bar = ext.xi.bar)
    }

    ll.means <- list(PZ.bar = PZ.bar, PT.bar = PT.bar, PE.bar = PE.bar,
                     LT.bar = LT.bar, LE.bar = LE.bar)
    beta.means <- list(typ.mean = beta.bar[,1], pz = beta.bar[,2], pe = beta.bar[,3])
  
    #Output for two-hurdle model
    output <- list(pD = pD,
                   DIC = DIC,
                   PPO = track.ppo/(iters - burn),
                   CPO = (iters - burn)/(track.icpo),
                   pars.means = pars.means,
                   ll.means = ll.means,
                   beta.means = beta.means,
                   dev = keep.dev[(burn + 1):iters],
                   beta = keep.beta[(burn + 1):iters,,],
                   pars = keep.pars[(burn + 1):iters,]
      )
  }

  if(hurd == Inf){
    #Keep betas
    keep.beta <- keep.beta[,,-3]
    
    #Output for one-hurdle model
    output <- list(pD = pD,
                   DIC = DIC,
                   PPO = track.ppo/(iters - burn),
                   CPO = (iters - burn)/(track.icpo),
                   ll.means = ll.means,
                   beta.means = beta.means,
                   dev = keep.dev[(burn + 1):iters],
                   beta = keep.beta[(burn + 1):iters,,]
        )
    }
  return(output)
}
