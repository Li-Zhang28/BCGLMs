#' Simulation for continuous response and mixed effects compositional predictors
#' @description Simulation of data with dependent compositional predictors and continuous response.
#' The simulation setting and dependence structure can be inferred from paper.
#' @param n sample size
#' @param p dimension of compositional predictors
#' @param snr proportion of small effects to the total number of predictors p
#' @return
#' @export
#' @import phyloseq
#' @import MASS
#' @references \url{https://link.springer.com/article/10.1186/s12859-025-06114-3}

sim_cmm <-
  function(n,p,snr) {
    gammatrue=rep(1,p)
    np=snr*p
    true_index=seq(18,28,by=2)
    b1=c(2.0800,-1.1600,-2.1200)
    b2=c(1.5000,-0.8600,0.56)
    b=rep(0,p)
    index1=seq(1,length(true_index),by=2)
    index2=seq(2,length(true_index),by=2)
    b[true_index[index1]]=b1
    b[true_index[index2]]=b2
    b[(p-np):(p-1)]=rnorm(np,0,0.2)
    b[p]=-sum(b[1:(p-1)])
    sigmaX=1
    xcor=matrix(0,p,p)
    for (i in seq(18,26,by=2)){
      for (j in seq(i+2,28,by=2)){
        xcor[i,j]=0.75-0.015*abs(i-j)
      }
    }

    for (i in (p-np):(p-1)){
      for (j in (p-np+1):(p)){
        xcor[i,j]=0.25-0.00015*abs(i-j)
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

