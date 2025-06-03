#' Bayesian Compositional Ordinal Regression
#' @description Bayesian compositional method for microbiome data with ordinal outcome.
#' @param x abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param y outcome (ordinal)
#' @param df_local the degree of freedom in half-t prior for local scale in horseshoe prior (default=1)
#' @param df_global the degree of freedom in half-t prior for global scale in horseshoe prior (default=1)
#' @param similarity measures the relatedness among taxa (default=NULL)
#' @import phyloseq
#' @import brms
#' @import BhGLM
#' @importFrom stats gaussian quantile
#' @return
#' @export
#' @examples dat=sim_o(n=400,p=100,q=c(0.2,0.5))
#'           fit=bco(x=dat$x,y=dat$y,df_local=1,df_global=1)
#'           summary(fit)
#'           fixef(fit)
#'           mcmc_plot(fit,variable = "^b_X", regex = TRUE)
#'           plot(fit,variable=c("b_X18","b_X20","b_X22","b_X24","b_X26","b_X28"),nvariables = 6)
#'
#' @author Li Zhang
#' @references \url{https://journals.sagepub.com/doi/abs/10.1177/09622802241247730}

bco <- function(x,y,df_local=1,df_global=1,similarity=NULL) {

  log_x=log(x)

  fm = bf(y ~ ., center=T)

  dat=data.frame(y,log_x)
  ### In real data analysis, the number of identified taxa in univariate generalized linear models can serve as a reasonable
  ### prior guess for 𝑚0
  if(is.null(similarity)){

    bp4 = set_prior(
      paste0("horseshoe(df=",df_local,", df_global=", df_global, ")"),
      class = "b"
    )
    bp4= bp4 + set_prior("target += normal_lpdf(mean(b) | 0, 0.001)", check=F) # sensitive to the variance of mean(b)


    f4= brm(fm, data=dat,family = cumulative("logit"),prior=bp4, init_r=0.1, control = list(adapt_delta = 0.99,max_treedepth= 20),
            chains=4, iter=4000)

  }

  if(!is.null(similarity)){
    node1=similarity$node1
    node2=similarity$node2
    w=similarity$w

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

  f4= brm(fm, data=dat,family = cumulative("logit"),prior=bp4, stanvars=stanvars,init_r=0.1, control = list(adapt_delta = 0.99,max_treedepth= 20),
          chains=4, iter=4000)
  }
  return(f4)
}
