#' Bayesian Compositional survival analysis
#' @description Bayesian compositional survival analysis for microbiome data.
#' @param x abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param y survival outcome
#' @param df_local the degree of freedom in half-t prior for local scale in horseshoe prior (default=1)
#' @param df_global the degree of freedom in half-t prior for global scale in horseshoe prior (default=1)
#' @param similarity measures the relatedness among taxa (default=NULL)
#' @import phyloseq
#' @import brms
#' @import BhGLM
#' @importFrom stats gaussian quantile
#' @return
#' @export
#' @examples dat=sim_s(n=400,p=100)
#'           sim=similarity(dat$x)
#'           fit=bccoxph(x=dat$x,y=dat$y,df_local=1,df_global=1,similarity=sim)
#'           summary(fit)
#'           fixef(fit)
#'           mcmc_plot(fit,variable = "^b_X", regex = TRUE)
#'           plot(fit,variable=c("b_X18","b_X20","b_X22","b_X24","b_X26","b_X28"),nvariables = 6)
#'
#' @author Zhenying Ding

bccoxph <- function(x,y,df_local=1,df_global=1,similarity=NULL) {

  log_x=log(x)

  fm = bf(y[,1]| cens(1-y[,2]) ~  ., center=T)

  dat=data.frame(y,log_x)
  ### In real data analysis, the number of identified taxa in univariate generalized linear models can serve as a reasonable
  ### prior guess for 𝑚0
  if(is.null(similarity)){
    bp4 = set_prior(
      paste0("horseshoe(df=",df_local,", df_global=", df_global, ")"),
      class = "b"
    )
    bp4= bp4 + set_prior("target += normal_lpdf(mean(b) | 0, 0.001)", check=F) # sensitive to the variance of mean(b)


    f4= brm(fm, data=dat,family = cox(),prior=bp4, control = list(adapt_delta = 0.95,max_treedepth= 18),
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

    f4= brm(fm, data=dat,family = cox(),prior=bp4, stanvars=stanvars,control = list(adapt_delta = 0.95,max_treedepth= 18),
            chains=4, iter=4000)

  }
  return(f4)
}



