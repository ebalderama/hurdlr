
#Typical count probability likelihood
PT <- function(p, q, log = T){
  if(log){log(1-p) + log(1-q)}
  else{(1-p)*(1-q)}
}