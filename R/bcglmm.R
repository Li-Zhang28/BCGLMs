#' Bayesian Compositional GLMMs
#' @description Bayesian compositional GLMMs methods for microbiome data with continuous or binary outcome.
#' @param x abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param y outcome (binary or continuous)
#' @param family gaussian or bernoulli
#' @param df_local the degree of freedom in half-t prior for local scale in horseshoe prior (default=1)
#' @param df_global the degree of freedom in half-t prior for global scale in horseshoe prior (default=1)
#' @param dist method to calculate dissimilarity matrix, default is Bray-curtis
#' @import phyloseq
#' @import brms
#' @import BhGLM
#' @importFrom stats gaussian quantile
#' @return
#' @export
#' @examples dat=sim_cmm(n=400,p=100,snr=0.2)
#'           otu=otu_table(dat$x,taxa_are_rows = F)
#'           fit=bcglmm(x=dat$x,y=dat$y,otu=otu,family=gaussian,df_local=3,df_global=3,dist="bray")
#'           summary(fit)
#'           fixef(fit)
#'           mcmc_plot(fit,variable = "^b_X", regex = TRUE)
#'           plot(fit,variable=c("b_X18","b_X20","b_X22","b_X24","b_X26","b_X28"),nvariables = 6)
#'
#' @author Li Zhang
#' @references \url{https://link.springer.com/article/10.1186/s12859-025-06114-3}
bcglmm <- function(x,y,otu,family=gaussian,df_local=1,df_global=1,dist="bray") {
  # --- Input handling ---
  # If x is already a phyloseq object, use directly
  if (inherits(otu, "phyloseq")) {
    phy <- otu
  }
  # If x is an otu_table, wrap it in a minimal phyloseq object
  else if (inherits(otu, "otu_table")) {
    phy <- phyloseq(otu)
  }
  else {
    stop("Input must be a phyloseq or otu_table object.")
  }

  dis.taxa = phyloseq::distance(phy, method=dist, type="taxa")
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
    e=e+0.00000002
    K <- v %*% diag(e) %*% t(v)
    K
  }

  K.taxa = simi.mat(dis.taxa); dim(K.taxa)

  taxa = as.numeric(K.taxa)
  node1 = node2 = w = NULL
  for (i in 1:(ncol(K.taxa)-1))
    for (j in (i+1):ncol(K.taxa))
    {
      if(abs(K.taxa[i,j]) > 0.18)#quantile(abs(taxa))[3])
      {
        node1 = c(node1,i);
        node2 = c(node2,j);
        w = c(w, K.taxa[i,j])
      }
    }
  w = sqrt(abs(w)); length(w)

  ####
  dis.sample = distance(phy, method=dist, type="sample")
  dis.sample = as.matrix(dis.sample); dim(dis.sample)

  # construct similarity matrix from distance matrix
  # see Zhao et al. (2015)
  K.sample = simi.mat(dis.sample); dim(K.sample)

  sample = as.factor(1:length(y))
  A = K.sample
  rownames(A) = sample
  colnames(A) = sample

  log_x=log(x)

  dat = data.frame(y, log_x, sample)

  fm = bf(y ~ . - sample + (1|gr(sample, cov=A)))

  ### In real data analysis, the number of identified taxa in univariate generalized linear models can serve as a reasonable
  ### prior guess for 𝑚0


  bp4 = set_prior(
    paste0("horseshoe(df=",df_local,", df_global=", df_global, ")"),
    class = "b"
  )
    bp4= bp4 + set_prior("target += normal_lpdf(mean(b) | 0, 0.001)", check=F) # sensitive to the variance of mean(b)
    bp4 = bp4 + set_prior("target += -0.5*dot_self(w .* (log(hs_local[node1])-log(hs_local[node2])))", check=F)

    ln = length(node1)
    stanvars =stanvar(x=ln, name="ln", scode="int ln;", block="data") +
      stanvar(x=node1, name="node1", scode="int node1[ln];", block="data") +
      stanvar(x=node2, name="node2", scode="int node2[ln];", block="data") +
      stanvar(x=w, name="w", scode="vector[ln] w;", block="data")

    f5=brm(fm, data=dat, data2=list(A=A), family=family, prior=bp4, stanvars=stanvars, control = list(adapt_delta = 0.999,max_treedepth= 20), chains=4, iter=4000)


  return(f5)
}



