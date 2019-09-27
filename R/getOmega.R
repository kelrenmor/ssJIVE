getOmega <- function(p){
  # Use smoothing Omega (from li2015supervised)
  Q <- matrix(0,nrow=p,ncol=p)
  R <- matrix(0,nrow=p,ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){
        Q[i,j] <- -2
        R[i,j] <- 2/3
      }
      if(abs(i-j)==1){
        Q[i,j] <- 1
        R[i,j] <- 1/6
      }
    }
  }
  Q <- Q[,2:(p-1)]
  R <- R[1:(p-2),1:(p-2)]
  Omega <- Q %*% solve(R) %*% t(Q) # pxp
  return(Omega)
}
