
#Extreme count probability likelihood
PE <- function(p, q, log = T){
  if(log){log(1-p) + log(q)}
  else{(1-p)*q}
}