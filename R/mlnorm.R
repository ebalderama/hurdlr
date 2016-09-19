
#Probability distribution function of discrete log normal
mlnorm <- function(y, meanlog = 0, sdlog = 1, log = T){
  pmf <- plnorm(y + 0.5, meanlog, sdlog) - plnorm( y - 0.5, meanlog, sdlog)
  if(log){return(log(pmf))}
    else{return(pmf)}
}
