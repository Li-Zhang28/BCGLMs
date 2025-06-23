#' Simulation for survival response and compositional predictors
#' @description Simulation of data with dependent compositional predictors and survival response.
#' The simulation setting and dependence structure can be inferred from paper.
#' @param n sample size
#' @param p dimension of compositional predictors
#' @return
#' @export
#' @import phyloseq
#' @import MASS
#' @import BhGLM
#' @import survival
#' @references \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.9946}
sim_s <- function(n,p) {

  gammatrue=rep(0,p)
  true_index=seq(18,28,by=2)
  gammatrue[true_index]=1
  b1=c(2.0800,-1.3900,2.1200)
  b2=c(1.3100,-1.8600,-2.26)
  b=rep(0,p)
  index1=seq(1,length(true_index),by=2)
  index2=seq(2,length(true_index),by=2)
  b[true_index[index1]]=b1
  b[true_index[index2]]=b2
  sigmaX=1
  #sigma=1/snr*mean(abs(c(b1,b2)))
  xcor=matrix(0,p,p)
  for (i in seq(18,26,by=2)){
    for (j in seq(i+2,28,by=2)){
      xcor[i,j]=0.75-0.015*abs(i-j)
    }
  }
  Xcor=xcor+t(xcor)+diag(p)
  s=diag(p)
  theta=rep(0,p)
  theta[true_index]=log(0.5*p)
  beta=gammatrue*b

  #########
  rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
                       method = c("svd", "eigen", "chol"), pre0.9_9994 = FALSE)
  {
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE)) {
      stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma))
      stop("mean and sigma have non-conforming size")
    method <- match.arg(method)
    R <- if (method == "eigen") {
      ev <- eigen(sigma, symmetric = TRUE)
      if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
        warning("sigma is numerically not positive definite")
      }
      t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
    }
    else if (method == "svd") {
      s. <- svd(sigma)
      if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
        warning("sigma is numerically not positive definite")
      }
      t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
    }
    else if (method == "chol") {
      R <- chol(sigma, pivot = TRUE)
      R[, order(attr(R, "pivot"))]
    }
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
    retval <- sweep(retval, 2, mean, "+")
    nm <- names(mean)
    if (is.null(nm) && !is.null(colnames(sigma))) nm <- colnames(sigma)
    colnames(retval) <- nm
    if (n == 1) drop(retval)
    else retval
  }


  sim.x <- function(n, m, corr = 0.6, v = rep(1, m), p = 0.5, genotype = NULL,
                    method = c("svd", "chol", "eigen"), joint = TRUE, verbose = FALSE)
  {
    start.time <- Sys.time()
    vars <- paste("x", 1:m, sep = "")
    V=s%*%Xcor%*%s
    rownames(V) <- colnames(V) <- unique(vars)

    method <- method[1]
    x <- matrix(0, n, m)
    colnames(x) <- unique(vars)
    x <- rmvnorm(n = n, mean =theta, sigma = V,method = method)

    stop.time <- Sys.time()
    minutes <- round(difftime(stop.time, start.time, units = "min"), 3)
    if (verbose) cat("simulation time:", minutes, "minutes \n")

    return(data.frame(x))
  }


  k = sim.x(n=n, m=p)

  colMeans(k)
  cor.test(k$x18,k$x20)
  sd(k$x1)
  x1=exp(k)
  x2 = x1/rowSums(x1)
  rowSums(x2)
  x3=log(x2)

  yy=sim.y(x=x3, mu=0, coefs=beta, sigma=1.6, theta=2)
  yy$coefs
  y = yy$y.surv;
  d=table(y[,2]);d[1]/sum(d)##censoring proportion

  dat=data.frame(y,x3)


  return(list(x=x2,y=y))
}
