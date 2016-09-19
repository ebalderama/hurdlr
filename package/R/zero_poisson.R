#Zero-Inflated Poisson model

#Zero-Inflated Poisson likelihood:

#              {   p,                   zi = 1
#    Yi|zi =   {
#              {   1 - p * f(x|lam),    zi = 0 

# logit(p) = g(p) = BX = B0 + B1*x1 + B2*x2 ...
# f(x|lam) ~ Poisson(lam)
# zi ~ Binomial(1, p)
# lambda ~ Gamma(a, b)



#Y data: zero-inflated poisson
rzip <- function(n, lambda, p = 1) {
	y  <- rbinom(n, 1, 1 - p)
	y[y==1] <- rpois(sum(y==1), lambda)
	return(y)
}


#Define data loglikelihood at z = 1, z = 0
loglike <- function(y, z, lam, p) {
  
	ll <- log(1 - p) + dpois(y, lam, log = T)
	ll[z==1] <- log(p[z==1])
	return(ll)
}

#MCMC start
zero_poisson <- function(y, x, a = 1, b = 1,
                lam.start = 1, lam.sd = 0.05,
                beta.prior.mean = 0, beta.prior.sd = 1,
                iters = 1000, burn = 500, nthin = 1,
                plots = T, progress.bar = T) {
  
    
  #Initial values
  pb <- txtProgressBar(min = 0, max = iters, style = 3)
  x <- cbind(1, x)
  x.col <- ncol(x)
  z <- rep(1, length(y))
  B <- rep(beta.prior.mean, length.out = x.col)
  XB <- x%*%B
  p <- 1/(1 + exp(-XB))
  lam <- lam.start
  ll <- loglike(y = y, z = z, lam = lam, p = p)
  
  #Tuning
  beta.tune <- rep(0.02, x.col)
  lam.tune <- 0.02
  beta.acc <- rep(0, x.col)
  lam.acc <- 0
  att <- 0

  #Create matrix of MCMC values for all parameters B, lambda  
  keep.beta <- matrix(0, iters, length(B))
  keep.lam <- rep(lam, iters)
  keep.ll <- rep(sum(ll), iters)
  keep.p <- rep(mean(p), iters)

  
for(i in 2:iters) {for(thin in 1:nthin){
	
  att <- att + 1
  
  #Update p (Update Betas)
  
  for(j in 1:x.col){
    
  	beta.current <- B[j]
  	beta.new <- rnorm(1, beta.current, beta.tune[j])
  	XB.new <- XB + x[,j]*(beta.new - beta.current)
  	p.new <- 1/(1 + exp(-XB.new))
  	ll.new <- loglike(y = y, z = z, lam = lam, p = p.new)
  	
  	p.ratio <- sum(ll.new) + sum(dnorm(beta.new, beta.prior.mean, beta.prior.sd, log = T)) - 
  	           sum(ll) - sum(dnorm(beta.current, beta.prior.mean, beta.prior.sd, log = T))
  	#if(is.finite(p.ratio))
  	if(log(runif(1)) < p.ratio){
  		
  	  B[j] <- beta.new
  		XB <- XB.new
  		p <- p.new
  		ll <- ll.new
  		beta.acc[j] <- beta.acc[j] + 1
  	}
  }
  
  #Update z
  
    P1 <- ifelse(y==0, 1, 0)*p
    P0 <- dpois(y, lam)*(1 - p)
    z  <- rbinom(length(y), 1, P1/(P1 + P0))
    ll <- loglike(y = y, z = z, lam = lam, p = p)
    
  #Update lambda
    
    lam.current <- lam
    lam.new <- exp(rnorm(1, log(lam.current), lam.tune))

    lam.ratio <- dgamma(lam.new, a, b, log = T) + sum(dpois(y[z==0], lam.new, log = T)) - 
                      dgamma(lam.current, a, b, log = T) - sum(dpois(y[z==0], lam.current, log = T))
    
    if(log(runif(1)) < lam.ratio) {
      lam <- lam.new
      lam.acc <- lam.acc + 1
      ll <- loglike(y = y, z = z, lam = lam.new, p = p)
    }
  
  }#End thinning
  
  
  keep.beta[i,] <- B
  keep.lam[i] <- lam
  keep.ll[i] <- sum(ll)
  keep.p[i] <- mean(p)
  
  
  #Update tuning parameters
  
  if(i < 0.75*burn & att > 50){
    beta.tune <- ifelse(beta.acc/att < 0.20, 0.8*beta.tune, beta.tune)
    beta.tune <- ifelse(beta.acc/att > 0.50, 1.2*beta.tune, beta.tune)
    lam.tune <- ifelse(lam.acc/att < 0.20, 0.8*lam.tune, lam.tune)
    lam.tune <- ifelse(lam.acc/att > 0.50, 1.2*lam.tune, lam.tune)
    
    beta.acc <- rep(0, x.col)
    lam.acc <- att <- 0
  }
  
  
  #Plots
  if(plots==T){
  if(i%%100==0){
  	
  	par(mfrow = c(2,4))
  	plot(keep.lam[1:i], type = "s")
  	  abline(h = mean(keep.lam[1:i]), col = 2)
  	plot(y, p, col = z+1, ylim = 0:1)
  	plot(keep.p[49:i],type="s")
  	  abline(h = mean(keep.p[49:i]), col = 2)
  	plot(keep.ll[49:i], type = "s")
  	  abline(h = mean(keep.ll[49:i]), col = 2)
  	plot(keep.beta[1:i, 1], type = "s")
  	  abline(h = mean(keep.beta[1:i, 1]), col = 2)
  	plot(keep.beta[1:i, 2], type = "s")
  	  abline(h = mean(keep.beta[1:i, 2]), col = 2)
  	plot(keep.beta[1:i, 3], type = "s")
  	  abline(h = mean(keep.beta[1:i, 3]), col = 2)
  	plot(keep.beta[1:i, 4], type = "s")
  	  abline(h = mean(keep.beta[1:i, 4]), col = 2)
  }
  }
  
  if(progress.bar){setTxtProgressBar(pb, i)}
}


  return(list(lam = keep.lam[(burn + 1):iters],
              beta = keep.beta[(burn + 1):iters,],
              ll = sum(ll)
              ))


}


#Example

if(F){

#Generate X and Y data
y <- rzip(100, lambda = 3, p = 0.20)
sum(y==0)

#X data: random whatever
n <- length(y)
x1 <- rexp(n)
x2 <- runif(n, 0, 1)
x3 <- rnorm(n, 0, 0.1)
x4 <- x1 + x2 + rnorm(n, 0, 0.2)
x5 <- 2*x1 - 2*x2 + runif(n, -0.5, 0.5)
x <- cbind(x1, x2, x3, x4, x5)

g <- zero_poisson(y, x = x, iters = 1000, burn = 200, nthin = 5, beta.prior.sd = 1)


} 
