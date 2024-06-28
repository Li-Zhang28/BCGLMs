#' Bayesian Compositional Ordinal Regression
#' @param x abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param y outcome (ordinal)
#' @param plot if TRUE, shows the plots (default = TRUE)
#' @import phyloseq
#' @import brms
#' @import BhGLM
#' @importFrom stats gaussian quantile
#' @return
#' @export
#' @author Li Zhang
#' @references \url{https://journals.sagepub.com/doi/abs/10.1177/09622802241247730}
bco <- function(x,y,plot=TRUE) {

  dist = c("bray", "jaccard", "jsd", "unifrac", "wunifrac", "dpcoa")
  m = dist[1]
  otu=otu_table(x,taxa_are_rows = F)
  dis.taxa = distance(otu, method=m, type="taxa")
  dis.taxa = as.matrix(dis.taxa); dim(dis.taxa)
  ###similarity matrix
  simi.mat <- function(dis.mat)
  {
    n <- ncol(dis.mat)
    D2 <- dis.mat^2
    I <- diag(1, n)
    II <- array(1, c(n,n))
    K <- -0.5 * (I - II/n) %*% D2 %*% (I - II/n)
    # to ensure a positive semi-definite matrix
    # method 1
    ev <- eigen(K)
    v <- ev$vectors
    e <- abs(ev$values)
    K <- v %*% diag(e) %*% t(v)
    K
  }

  K.taxa = simi.mat(dis.taxa); dim(K.taxa)

  taxa = as.numeric(K.taxa)
  node1 = node2 = w = NULL
  for (i in 1:(ncol(K.taxa)-1))
    for (j in (i+1):ncol(K.taxa))
    {
      if(abs(K.taxa[i,j]) > quantile(abs(taxa))[4]) ###0.18
      {
        node1 = c(node1,i);
        node2 = c(node2,j);
        w = c(w, K.taxa[i,j])
      }
    }
  w = sqrt(abs(w)); length(w)

  log_x=log(x)

  fm = bf(y ~ ., center=T)

  dat=data.frame(y,log_x)
  ### In real data analysis, the number of identified taxa in univariate generalized linear models can serve as a reasonable
  ### prior guess for 𝑚0

  bp4 = set_prior("horseshoe(df=3, df_global=3)", class="b")
  bp4= bp4 + set_prior("target += normal_lpdf(mean(b) | 0, 0.001)", check=F) # sensitive to the variance of mean(b)
  bp4 = bp4 + set_prior("target += -0.5*dot_self(w .* (log(hs_local[node1])-log(hs_local[node2])))", check=F)

  ln = length(node1)
  stanvars =stanvar(x=ln, name="ln", scode="int ln;", block="data") +
    stanvar(x=node1, name="node1", scode="int node1[ln];", block="data") +
    stanvar(x=node2, name="node2", scode="int node2[ln];", block="data") +
    stanvar(x=w, name="w", scode="vector[ln] w;", block="data")

  f4= brm(fm, data=dat,family = cumulative("logit"),prior=bp4, stanvars=stanvars,init_r=0.1, control = list(adapt_delta = 0.99,max_treedepth= 20),
          chains=2, iter=2000)
}
