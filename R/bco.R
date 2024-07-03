#' Bayesian Compositional Ordinal Regression
#' @param x abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param y outcome (ordinal)
#' @import phyloseq
#' @import brms
#' @import BhGLM
#' @importFrom stats gaussian quantile
#' @return
#' @export
#' @examples dat=sim_o(200,100)
#'           fit=bco(dat$x,dat$y)
#'           summary(fit)
#'           fixef(fit)
#'           mcmc_plot(fit)
#'
#' @author Li Zhang
#' @references \url{https://journals.sagepub.com/doi/abs/10.1177/09622802241247730}

bco <- function(x,y,similarity=NULL) {

  log_x=log(x)

  fm = bf(y ~ ., center=T)

  dat=data.frame(y,log_x)
  ### In real data analysis, the number of identified taxa in univariate generalized linear models can serve as a reasonable
  ### prior guess for 𝑚0
  if(is.null(similarity)){

    bp4 = set_prior("horseshoe(df=3, df_global=3)", class="b")
    bp4= bp4 + set_prior("target += normal_lpdf(mean(b) | 0, 0.001)", check=F) # sensitive to the variance of mean(b)


    f4= brm(fm, data=dat,family = cumulative("logit"),prior=bp4, init_r=0.1, control = list(adapt_delta = 0.99,max_treedepth= 20),
            chains=2, iter=2000)

  }

  if(!is.null(similarity)){
    node1=similarity$node1
    node2=similarity$node2
    w=similarity$w

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
  return(f4)
}
