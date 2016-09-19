
#Zero count probability likelihood
PZ <- function(p, log = T){
  if(log){log(p)}
  else{p}
}