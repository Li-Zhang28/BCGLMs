#' Simulation for continuous response and compositional predictors
#' @description Simulation of data with dependent compositional predictors and continuous response.
#' The complex dependence structure can be inferred from paper.
#' @param n sample size
#' @param p dimension of compositional predictors
#' @return
#' @export
#' @import phyloseq
#' @import MASS
#' @references \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.9946}
sim_c <- function(n,p) {
  gammatrue=rep(0,p)
  true_index=seq(18,40,by=2)
  gammatrue[true_index]=1
  b1=c(2.0800,-1.3900,2.1200,1.3100,-1.8600,-1.3400)
  b2=c(-1.4100,-1.1500,2.5100,-1.9500,1.9300,-0.8500)
  b=rep(0,p)
  index1=seq(1,length(true_index),by=2)
  index2=seq(2,length(true_index),by=2)
  b[true_index[index1]]=b1
  b[true_index[index2]]=b2
  sigmaX=1
  #sigma=1/snr*mean(abs(c(b1,b2)))
  xcor=matrix(0,p,p)
  for (i in seq(18,38,by=2)){
    for (j in seq(i+2,40,by=2)){
      xcor[i,j]=0.75-0.015*abs(i-j)
    }
  }
  Xcor=xcor+t(xcor)+diag(p)
  theta=rep(0,p)
  theta[true_index]=log(0.5*p)
  beta=gammatrue*b
  s=diag(p)*sigmaX
  X=mvrnorm(n,theta,s%*%Xcor%*%s)
  x=exp(X)
  x1 = x/rowSums(x)
  x2=log(x1)
  epsilon=rnorm(n,0,1.6)
  y=x2%*%beta+epsilon

  return(list(x=x1,y=y))
}
